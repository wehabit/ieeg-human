"""Build a scalp-only sidecar for simultaneous scalp-versus-iEEG comparisons.

This is intentionally *not* a replacement for ``cache_lc_series.py``.  It pins and
reuses the exact v8 HUP cache interval (``night_s``, ``hours``, and ``sf``), streams
only explicitly allow-listed scalp channels from that same portal snapshot, and
stores reversible band-power support plus pre-threshold slow-oscillation candidates.
It never redetects ECG beats, chooses another night, or restages sleep.

The downstream comparison fixes the non-EEG inputs, but its EEG measurement arms
differ jointly in modality, location, reference, and contact aggregation.  It is
therefore not a pure sensor-modality intervention or an exact Lecci/Naji
replication.  The online scalp reference is undocumented, and HUP staging/ECG
processing also differ from those papers.

Examples
--------
Preview the frozen channel plan without contacting the portal::

    .venv/bin/python analysis/cache_paired_scalp.py --show-plan

Acquire the nine datasets with active C3/C03 scalp recordings::

    .venv/bin/python analysis/cache_paired_scalp.py --force

Acquire only the two datasets that also contain F3::

    .venv/bin/python analysis/cache_paired_scalp.py --subjects 160,187 --force
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import traceback
from typing import Iterable

import numpy as np
from scipy import signal

from cache_lc_series import (
    FILTER_EDGE_S,
    SIGMA_FIXED,
    SO_BAND_NAJI,
    SWA_BAND_L,
    _write_multichannel_power,
    detect_so_candidates,
    empty_channel_activity_extrema,
    finalize_channel_activity_qc,
    prepare_continuous_signal,
    update_channel_activity_extrema,
)
from cohort_stages_3ABD import CHUNK_S, band_sos
from infraslow_rr_sigma_coherence import notch, pull_continuous, sess
from pipeline_version import (
    atomic_json_dump,
    atomic_savez,
    file_sha256,
    git_is_dirty,
    git_revision,
    npz_scalar_text,
    runtime_versions,
    utc_now,
)
from results_3A_tutorial_style import ied_clean_mask


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_IEEG_CACHE = os.path.join(ROOT, "data", "derived", "lc_infraslow")
DEFAULT_OUTPUT = os.path.join(ROOT, "data", "derived", "paired_scalp")
LEGACY_IEEG_SCHEMA = "2026-07-neutral-per-contact-gap-aware-source-pin-v8"
SCALP_CACHE_SCHEMA = "2026-07-paired-scalp-sidecar-v1"
PIPELINE = "cache_paired_scalp"
HISTORICAL_POWER_SUPPORT = 0.5

# This is an explicit, audit-derived acquisition plan.  C3/C03 is the Lecci-aligned
# primary scalp sensor.  F3 is an incomplete unilateral Naji sensitivity (the paper
# derived a cardiac curve per referenced F3/A2 and F4/A1 electrode, then averaged
# the electrode-specific HR-maximum/RR-minimum times).  Fz was not used by Naji.
#
# HUP182 is retained in the acquisition plan because it has an active C3 scalp
# channel, even though its pinned iEEG 3A endpoint may remain unavailable.  That
# distinction is made downstream rather than selected here from scalp outcomes.
SCALP_CHANNEL_PLAN = {
    "HUP160_phaseII": {"c3": "C3", "f3": "F3", "fz": "Fz"},
    "HUP182_phaseII": {"c3": "C3"},
    "HUP185_phaseII": {"c3": "C3", "fz": "Fz"},
    "HUP187_phaseII": {"c3": "C3", "f3": "F3", "fz": "Fz"},
    "HUP191_phaseII": {"c3": "C3", "fz": "Fz"},
    "HUP199_phaseII": {"c3": "C3", "fz": "Fz"},
    "HUP205_phaseII": {"c3": "C3", "fz": "Fz"},
    "HUP211_phaseII": {"c3": "C3", "fz": "Fz"},
    "HUP212_phaseII": {"c3": "C03", "fz": "Fz"},
}

_DEPENDENCY_FILES = (
    "analysis/cache_paired_scalp.py",
    "analysis/cache_lc_series.py",
    "analysis/cohort_stages_3ABD.py",
    "analysis/infraslow_rr_sigma_coherence.py",
    "analysis/pipeline_version.py",
    "analysis/results_3A_tutorial_style.py",
    "env/requirements.txt",
    "env/PYTHON_VERSION",
)


def _hash_files(relative_paths: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for relative in sorted(relative_paths):
        path = os.path.join(ROOT, *relative.split("/"))
        if not os.path.isfile(path):
            raise FileNotFoundError(f"paired-scalp dependency is missing: {path}")
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        with open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        digest.update(b"\0")
    return digest.hexdigest()


def cache_dependency_sha256() -> str:
    """Digest the exact source/config bytes that can change a scalp sidecar."""
    return _hash_files(_DEPENDENCY_FILES)


def _normalize_scalp_label(label: str) -> str:
    """Normalize case and harmless zero padding, without fuzzy channel matching."""
    value = str(label).strip().upper()
    head = value.rstrip("0123456789")
    tail = value[len(head):]
    if tail:
        tail = str(int(tail))
    return head + tail


def _resolve_channels(labels, role_plan):
    """Resolve every explicit plan label uniquely against portal labels."""
    labels = [str(value) for value in labels]
    normalized = {}
    for label in labels:
        normalized.setdefault(_normalize_scalp_label(label), []).append(label)
    resolved = {}
    for role, planned_label in role_plan.items():
        matches = normalized.get(_normalize_scalp_label(planned_label), [])
        if len(matches) != 1:
            raise RuntimeError(
                f"planned scalp channel {planned_label!r} ({role}) resolves to "
                f"{matches!r}; exactly one portal label is required")
        resolved[role] = matches[0]
    if len(set(resolved.values())) != len(resolved):
        raise RuntimeError("two paired-scalp roles resolved to the same portal channel")
    return resolved


def _source_node_map(ds):
    nodes = {}
    for node in getattr(ds, "ts_array", []):
        label = node.findtext("channelLabel")
        if label is not None:
            nodes[str(label)] = node
    return nodes


def _channel_identity(ds, label, nodes=None):
    detail = ds.get_time_series_details(label)
    nodes = _source_node_map(ds) if nodes is None else nodes
    node = nodes.get(label)
    return {
        "revision_id": str(getattr(detail, "portal_id", "") or ""),
        "data_check": None if node is None else node.findtext("dataCheck"),
        "start_time_us": int(detail.start_time),
        "end_time_us": int(detail.end_time),
        "duration_us": float(detail.duration),
        "number_of_samples": int(detail.number_of_samples),
        "sample_rate_hz": float(detail.sample_rate),
    }


def _same_series_geometry(actual, reference):
    exact_keys = (
        "start_time_us",
        "end_time_us",
        "number_of_samples",
    )
    if any(actual[key] != reference[key] for key in exact_keys):
        return False
    return bool(
        np.isclose(actual["duration_us"], reference["duration_us"], rtol=0, atol=1)
        and np.isclose(
            actual["sample_rate_hz"], reference["sample_rate_hz"], rtol=0, atol=1e-9)
    )


def validate_pinned_ieeg_cache(subject, cache_dir=DEFAULT_IEEG_CACHE):
    """Validate the exact legacy cache bytes using its terminal run manifest.

    The check deliberately does not compare the v8 cache-builder digest with the
    current source tree.  That would make an immutable historical cache unreadable
    after a bug fix.  Instead, the terminal v8 manifest pins the exact NPZ bytes,
    and this sidecar records both hashes as immutable inputs.
    """
    manifest_path = os.path.join(cache_dir, "RUN_MANIFEST.json")
    cache_path = os.path.join(cache_dir, f"{subject}.npz")
    if not os.path.isfile(manifest_path) or not os.path.isfile(cache_path):
        raise FileNotFoundError(f"missing pinned iEEG cache input for {subject}")
    with open(manifest_path) as handle:
        manifest = json.load(handle)
    if manifest.get("run_state") != "complete":
        raise RuntimeError(f"{manifest_path} is not a completed run")
    if manifest.get("pipeline") != "cache_lc_series":
        raise RuntimeError(f"{manifest_path} is not a cache_lc_series run")
    if manifest.get("cache_schema_version") != LEGACY_IEEG_SCHEMA:
        raise RuntimeError(
            f"paired comparison requires legacy schema {LEGACY_IEEG_SCHEMA!r}; "
            f"manifest has {manifest.get('cache_schema_version')!r}")
    if subject not in set(manifest.get("completed", [])):
        raise RuntimeError(f"{subject} is not completed in the pinned cache manifest")
    expected = manifest.get("result_files_sha256", {}).get(subject)
    actual = file_sha256(cache_path)
    if not expected or actual != expected:
        raise RuntimeError(f"{subject} cache bytes differ from the terminal v8 manifest")
    with np.load(cache_path, allow_pickle=False) as cache:
        if npz_scalar_text(cache, "subject") != subject:
            raise RuntimeError(f"{cache_path} embeds a different participant")
        if npz_scalar_text(cache, "status") != "ok":
            raise RuntimeError(f"{cache_path} is not an OK cache")
        if npz_scalar_text(cache, "cache_schema_version") != LEGACY_IEEG_SCHEMA:
            raise RuntimeError(f"{cache_path} is not a v8 neutral cache")
        failed = json.loads(npz_scalar_text(cache, "failed_chunks_json", "[]"))
        if failed:
            raise RuntimeError(f"{cache_path} contains failed acquisition chunks")
        required = ("night_s", "hours", "sf", "source_identity_json")
        missing = [key for key in required if key not in cache.files]
        if missing:
            raise RuntimeError(f"{cache_path} lacks frozen fields: {missing}")
        pinned = {
            "night_s": float(np.asarray(cache["night_s"]).item()),
            "hours": float(np.asarray(cache["hours"]).item()),
            "sf": float(np.asarray(cache["sf"]).item()),
            "source_identity": json.loads(
                npz_scalar_text(cache, "source_identity_json")),
        }
    return {
        "path": cache_path,
        "sha256": actual,
        "manifest_path": manifest_path,
        "manifest_sha256": file_sha256(manifest_path),
        "manifest_run_id": manifest.get("run_id"),
        **pinned,
    }


def _portal_source_identity(ds, resolved, pinned):
    old_identity = pinned["source_identity"]
    expected_snapshot = str(old_identity.get("snapshot_id", "") or "")
    actual_snapshot = str(getattr(ds, "snap_id", "") or "")
    if actual_snapshot != expected_snapshot:
        raise RuntimeError(
            f"portal snapshot changed: expected {expected_snapshot!r}, "
            f"got {actual_snapshot!r}")
    old_channels = old_identity.get("channels", {})
    reference_label = next(
        (value for value in old_identity.get("cortical_channels", [])
         if value in old_channels),
        None,
    )
    if reference_label is None:
        raise RuntimeError("pinned source identity has no iEEG reference channel")
    labels = list(ds.get_channel_labels())
    if reference_label not in labels:
        raise RuntimeError(
            f"pinned iEEG reference channel {reference_label!r} is absent")
    nodes = _source_node_map(ds)
    reference = _channel_identity(ds, reference_label, nodes)
    expected_reference = old_channels[reference_label]
    for key in (
        "revision_id",
        "data_check",
        "start_time_us",
        "end_time_us",
        "number_of_samples",
    ):
        if reference.get(key) != expected_reference.get(key):
            raise RuntimeError(
                f"pinned iEEG reference {reference_label} changed at {key}")
    for key, tolerance in (("duration_us", 1), ("sample_rate_hz", 1e-9)):
        if not np.isclose(
                reference[key], expected_reference[key], rtol=0, atol=tolerance):
            raise RuntimeError(
                f"pinned iEEG reference {reference_label} changed at {key}")
    if not np.isclose(reference["sample_rate_hz"], pinned["sf"], rtol=0, atol=1e-9):
        raise RuntimeError("portal reference sample rate differs from frozen cache sf")

    scalp = {}
    for role, label in resolved.items():
        identity = _channel_identity(ds, label, nodes)
        if not _same_series_geometry(identity, reference):
            raise RuntimeError(
                f"scalp channel {label!r} does not share the pinned iEEG time base")
        scalp[label] = identity
    return {
        "dataset_name": str(getattr(ds, "name", "") or ""),
        "snapshot_id": actual_snapshot,
        "reference_ieeg_channel": reference_label,
        "reference_ieeg_identity": reference,
        "scalp_role_to_channel": dict(resolved),
        "scalp_channels": scalp,
        "geometry_requirement": (
            "same start/end/sample count/duration/sample rate as the pinned "
            "simultaneously recorded iEEG reference channel"),
    }


def build_subject(subject, *, cache_dir=DEFAULT_IEEG_CACHE,
                  output_dir=DEFAULT_OUTPUT, force=False):
    if subject not in SCALP_CHANNEL_PLAN:
        raise ValueError(f"{subject} is not in the explicit paired-scalp plan")
    pinned = validate_pinned_ieeg_cache(subject, cache_dir)
    output_path = os.path.join(output_dir, f"{subject}.npz")
    dependency_digest = cache_dependency_sha256()
    if os.path.exists(output_path) and not force:
        with np.load(output_path, allow_pickle=False) as existing:
            reusable = (
                npz_scalar_text(existing, "status") == "ok"
                and npz_scalar_text(
                    existing, "cache_schema_version") == SCALP_CACHE_SCHEMA
                and npz_scalar_text(
                    existing, "ieeg_cache_sha256") == pinned["sha256"]
                and npz_scalar_text(
                    existing, "cache_dependency_sha256") == dependency_digest
            )
        if reusable:
            print(f"[{subject}] validated existing scalp sidecar", flush=True)
            return output_path, "reused"
        raise RuntimeError(
            f"{output_path} is stale or unvalidated; rerun with --force")

    session = sess()
    ds = None
    try:
        ds = session.open_dataset(subject)
        labels = list(ds.get_channel_labels())
        resolved = _resolve_channels(labels, SCALP_CHANNEL_PLAN[subject])
        source_identity = _portal_source_identity(ds, resolved, pinned)
        channels = list(resolved.values())
        channel_indices = [labels.index(channel) for channel in channels]
        sf = float(pinned["sf"])
        night_s = float(pinned["night_s"])
        hours = float(pinned["hours"])
        total_s = int(round(hours * 3600.0))
        if total_s <= 0 or not np.isclose(total_s, hours * 3600.0):
            raise RuntimeError("frozen cache duration is not a positive whole second")

        n_channels = len(channels)
        sigma_values = np.full((n_channels, total_s), np.nan)
        swa_values = np.full((n_channels, total_s), np.nan)
        sigma_numerator = np.full((n_channels, total_s), np.nan)
        swa_numerator = np.full((n_channels, total_s), np.nan)
        sigma_denominator = np.zeros((n_channels, total_s), dtype=np.uint32)
        swa_denominator = np.zeros((n_channels, total_s), dtype=np.uint32)
        activity = empty_channel_activity_extrema(n_channels)
        candidates = {channel: [] for channel in channels}
        failed_chunks = []

        sigma_sos = band_sos(SIGMA_FIXED, sf)
        swa_sos = band_sos(SWA_BAND_L, sf, 3)
        so_sos = signal.butter(
            3, list(SO_BAND_NAJI), btype="band", fs=sf, output="sos")
        print(
            f"[{subject}] {sf:g} Hz | {channels} | frozen "
            f"{night_s / 3600:.3f}-{(night_s + total_s) / 3600:.3f} h",
            flush=True,
        )

        t = 0.0
        while t < total_s:
            duration = min(float(CHUNK_S), total_s - t)
            pull_start = max(0.0, t - FILTER_EDGE_S)
            pull_stop = min(float(total_s), t + duration + FILTER_EDGE_S)
            try:
                data = pull_continuous(
                    ds,
                    channel_indices,
                    night_s + pull_start,
                    pull_stop - pull_start,
                )
                if data.ndim != 2 or data.shape[1] != n_channels:
                    raise RuntimeError(
                        f"portal returned unexpected block shape {data.shape}")
            except Exception as exc:
                failed_chunks.append({
                    "start_s": float(t),
                    "duration_s": float(duration),
                    "error": f"{type(exc).__name__}: {exc}",
                })
                raise

            core_start = int(round((t - pull_start) * sf))
            requested_core = int(round(duration * sf))
            core_count = min(requested_core, max(0, len(data) - core_start))
            if core_count != requested_core:
                raise RuntimeError(
                    f"portal returned {core_count} core samples; expected "
                    f"{requested_core} at t={t:g}")
            core_stop = core_start + core_count
            update_channel_activity_extrema(
                activity, data[core_start:core_stop].T)

            prepared = [
                prepare_continuous_signal(data[:, index], sf)
                for index in range(n_channels)
            ]
            filtered = np.asarray([
                notch(signal.detrend(filled), sf)
                for filled, _ in prepared
            ])
            measured = np.asarray([mask for _, mask in prepared])
            for index, channel in enumerate(channels):
                waveform = filtered[index]
                if not np.isfinite(waveform).all() or np.std(waveform) < 1e-9:
                    continue
                so = signal.sosfiltfilt(so_sos, waveform)
                clean = (
                    ied_clean_mask(waveform, sf, pad_s=5.0)
                    & measured[index]
                )
                found = detect_so_candidates(so, sf)
                for trough, down, up, peak_to_peak in found:
                    trough = int(trough)
                    if (
                        core_start <= trough < core_stop
                        and clean[trough]
                    ):
                        candidates[channel].append((
                            trough / sf + pull_start,
                            down,
                            up,
                            peak_to_peak,
                        ))

            offset = int(round(t))
            n_seconds = int(core_count // int(round(sf)))
            _write_multichannel_power(
                data,
                sigma_sos,
                sf,
                n_seconds,
                total_s,
                offset,
                sigma_values,
                core_start_sample=core_start,
                numerator_dest=sigma_numerator,
                clean_sample_count_dest=sigma_denominator,
            )
            _write_multichannel_power(
                data,
                swa_sos,
                sf,
                n_seconds,
                total_s,
                offset,
                swa_values,
                core_start_sample=core_start,
                numerator_dest=swa_numerator,
                clean_sample_count_dest=swa_denominator,
            )
            t += duration

        activity_qc = finalize_channel_activity_qc(activity)
        if failed_chunks:
            raise RuntimeError(
                f"{subject} has {len(failed_chunks)} failed acquisition chunks")
        payload = {
            "status": "ok",
            "cache_schema_version": SCALP_CACHE_SCHEMA,
            "cache_dependency_sha256": dependency_digest,
            "generated_at_utc": utc_now(),
            "code_revision": git_revision(ROOT),
            "code_dirty": git_is_dirty(ROOT),
            "runtime_versions_json": json.dumps(
                runtime_versions(), sort_keys=True),
            "subject": subject,
            "source_dataset": subject,
            "source_kind": "iEEG.org API",
            "source_identity_json": json.dumps(
                source_identity, sort_keys=True),
            "channel_roles_json": json.dumps(resolved, sort_keys=True),
            "scalp_chans": np.asarray(channels, dtype="<U16"),
            "sf": sf,
            "night_s": night_s,
            "hours": hours,
            "ieeg_cache_sha256": pinned["sha256"],
            "ieeg_cache_manifest_sha256": pinned["manifest_sha256"],
            "ieeg_cache_manifest_run_id": str(
                pinned["manifest_run_id"] or ""),
            "ieeg_cache_schema_version": LEGACY_IEEG_SCHEMA,
            "failed_chunks_json": json.dumps(
                failed_chunks, sort_keys=True),
            "filter_edge_seconds": float(FILTER_EDGE_S),
            "chunk_seconds": float(CHUNK_S),
            "sigma_band_hz": np.asarray(SIGMA_FIXED, float),
            "swa_band_hz": np.asarray(SWA_BAND_L, float),
            "so_band_hz": np.asarray(SO_BAND_NAJI, float),
            "sigma_fixed_by_channel": sigma_values,
            "swa_by_channel": swa_values,
            "sigma_fixed_power_numerator_by_channel": sigma_numerator,
            "sigma_fixed_clean_sample_count_by_channel": sigma_denominator,
            "swa_power_numerator_by_channel": swa_numerator,
            "swa_clean_sample_count_by_channel": swa_denominator,
            "power_samples_per_second": int(round(sf)),
            "power_historical_minimum_clean_fraction_per_second": (
                HISTORICAL_POWER_SUPPORT),
            "power_support_semantics": (
                "numerator=sum of IED-clean squared-envelope samples in each "
                "nominal 1-s bin; denominator=count of admitted samples"),
            "scalp_signal_nonflat_mask": activity_qc["nonflat_mask"],
            "scalp_signal_raw_minimum": activity_qc["minimum"],
            "scalp_signal_raw_maximum": activity_qc["maximum"],
            "scalp_signal_raw_dynamic_range": activity_qc["dynamic_range"],
            "scalp_signal_numerical_flat_tolerance": activity_qc[
                "numerical_flat_tolerance"],
            "scalp_signal_finite_sample_count": activity_qc["finite_count"],
            "scalp_signal_activity_qc_method": (
                "raw full-interval range must exceed 64 float64 eps times "
                "max(1, maximum absolute raw value)"),
            "ecg_reused_not_redetected": True,
            "staging_reused_not_recomputed": True,
        }
        for channel in channels:
            values = np.asarray(
                sorted(candidates[channel], key=lambda value: value[0]),
                float,
            ).reshape(-1, 4)
            payload[f"so_candidate_t_{channel}"] = values[:, 0]
            payload[f"so_candidate_down_{channel}"] = values[:, 1]
            payload[f"so_candidate_up_{channel}"] = values[:, 2]
            payload[f"so_candidate_p2p_{channel}"] = values[:, 3]
        atomic_savez(output_path, **payload)
        print(
            f"[{subject}] saved {sum(len(v) for v in candidates.values())} "
            f"pre-threshold SO candidates -> {output_path}",
            flush=True,
        )
        return output_path, "completed"
    finally:
        try:
            if ds is not None and hasattr(session, "close_dataset"):
                session.close_dataset(ds)
        finally:
            if hasattr(session, "close"):
                session.close()


def _parse_subjects(value):
    if value is None or not str(value).strip():
        return list(SCALP_CHANNEL_PLAN)
    result = []
    for item in str(value).split(","):
        token = item.strip()
        if not token:
            continue
        subject = (
            token if token.startswith("HUP")
            else f"HUP{int(token)}_phaseII"
        )
        if subject.startswith("HUP") and not subject.endswith("_phaseII"):
            subject += "_phaseII"
        if subject not in SCALP_CHANNEL_PLAN:
            raise ValueError(f"{subject} is not in the explicit scalp plan")
        result.append(subject)
    if len(result) != len(set(result)):
        raise ValueError("subjects contain duplicates")
    return result


def _write_manifest(output_dir, *, requested, completed, reused, failed,
                    cache_dir):
    os.makedirs(output_dir, exist_ok=True)
    result_hashes = {
        subject: file_sha256(os.path.join(output_dir, f"{subject}.npz"))
        for subject in completed + reused
    }
    base_manifest = os.path.join(cache_dir, "RUN_MANIFEST.json")
    payload = {
        "schema_version": SCALP_CACHE_SCHEMA,
        "pipeline": PIPELINE,
        "run_state": "complete" if not failed else "failed",
        "generated_at_utc": utc_now(),
        "code_revision": git_revision(ROOT),
        "code_dirty": git_is_dirty(ROOT),
        "runtime_versions": runtime_versions(),
        "cache_dependency_sha256": cache_dependency_sha256(),
        "requested": requested,
        "completed": completed,
        "reused": reused,
        "failed": failed,
        "channel_plan": {
            subject: SCALP_CHANNEL_PLAN[subject] for subject in requested
        },
        "config": {
            "ieeg_cache_directory_relative": os.path.relpath(
                os.path.abspath(cache_dir), ROOT).replace(os.sep, "/"),
            "ieeg_cache_manifest_sha256": file_sha256(base_manifest),
            "legacy_ieeg_schema": LEGACY_IEEG_SCHEMA,
            "scalp_cache_schema": SCALP_CACHE_SCHEMA,
            "filter_edge_seconds": FILTER_EDGE_S,
            "chunk_seconds": CHUNK_S,
            "ecg_reused_not_redetected": True,
            "staging_reused_not_recomputed": True,
        },
        "result_files_sha256": result_hashes,
    }
    atomic_json_dump(payload, os.path.join(output_dir, "RUN_MANIFEST.json"))
    return payload


def main():
    parser = argparse.ArgumentParser(
        description="Cache simultaneous scalp channels on frozen v8 HUP intervals")
    parser.add_argument("--subjects", help="comma-separated HUP numbers/names")
    parser.add_argument("--ieeg-cache", default=DEFAULT_IEEG_CACHE)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--show-plan", action="store_true",
        help="print the frozen acquisition plan and exit without portal access")
    args = parser.parse_args()
    subjects = _parse_subjects(args.subjects)
    if args.show_plan:
        print(json.dumps(
            {subject: SCALP_CHANNEL_PLAN[subject] for subject in subjects},
            indent=2,
            sort_keys=True,
        ))
        return

    completed, reused, failed = [], [], []
    for subject in subjects:
        try:
            _, status = build_subject(
                subject,
                cache_dir=os.path.abspath(args.ieeg_cache),
                output_dir=os.path.abspath(args.output_dir),
                force=args.force,
            )
            (reused if status == "reused" else completed).append(subject)
        except Exception as exc:
            traceback.print_exc()
            failed.append({
                "subject": subject,
                "error": f"{type(exc).__name__}: {exc}",
            })
    _write_manifest(
        os.path.abspath(args.output_dir),
        requested=subjects,
        completed=completed,
        reused=reused,
        failed=failed,
        cache_dir=os.path.abspath(args.ieeg_cache),
    )
    if failed:
        print(
            f"{len(failed)} paired-scalp acquisition(s) failed; manifest is "
            "marked failed",
            file=sys.stderr,
        )
        raise SystemExit(1)
    print(
        f"{len(completed)} built, {len(reused)} reused -> "
        f"{os.path.abspath(args.output_dir)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
