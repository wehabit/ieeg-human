"""Build a scalp-only sidecar for simultaneous scalp-versus-iEEG comparisons.

This is intentionally *not* a replacement for ``cache_lc_series.py``.  It pins and
reuses the exact HUP cache interval (``night_s``, ``hours``, and ``sf``), streams
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

Acquire the ten datasets with unique, geometry-matched C3/C03 recordings::

    .venv/bin/python analysis/cache_paired_scalp.py --force

Acquire two example datasets that also contain F3::

    .venv/bin/python analysis/cache_paired_scalp.py --subjects 160,187 --force
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import traceback
import uuid
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
    require_complete_acquisition,
    update_channel_activity_extrema,
)
from hup_portal import (
    portal_core_sample_geometry,
    pull_continuous_exact,
    same_series_geometry,
)
from staging_helpers import CHUNK_S, band_sos
from infraslow_rr_sigma_coherence import notch, sess
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    atomic_json_dump,
    atomic_savez,
    file_sha256,
    git_is_dirty,
    git_revision,
    npz_scalar_text,
    require_integer_sample_rate,
    runtime_versions,
    utc_now,
    validate_completed_cache_failures,
    validate_full_interval_acquisition,
    validate_terminal_run_manifest,
)
from paired_scalp_sidecar_validation import (
    TERMINAL_SIDECAR_MANIFEST_FIELDS,
    SidecarValidationContract,
    _is_lower_hex,
    _required_array,
    _terminal_sidecar_partitions,
    _valid_manifest_run_id,
    _valid_utc_timestamp,
    _validate_activity,
    _validate_power_support as _validate_power_support_impl,
    _validate_so_candidates,
    dependency_sha256_at_revision as _dependency_sha256_at_revision_impl,
    recorded_commit_exists as _recorded_commit_exists_impl,
    validate_reusable_sidecar_run as _validate_reusable_sidecar_run_impl,
    validate_scalp_sidecar_payload as _validate_scalp_sidecar_payload_impl,
    validate_terminal_sidecar_manifest as _validate_terminal_manifest_impl,
)
from signal_qc import ied_clean_mask


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_IEEG_CACHE = os.path.join(ROOT, "data", "derived", "lc_infraslow")
DEFAULT_OUTPUT = os.path.join(ROOT, "data", "derived", "paired_scalp")
IEEG_CACHE_SCHEMA = CACHE_SCHEMA_VERSION
SCALP_CACHE_SCHEMA = "2026-07-paired-scalp-sidecar-v3"
PIPELINE = "cache_paired_scalp"
HISTORICAL_POWER_SUPPORT = 0.5
# Compatibility alias for the inventory auditor; implementation is centralized.
_same_series_geometry = same_series_geometry

# This is an explicit, audit-derived acquisition plan. C3/C03 is the
# Lecci-motivated scalp location used by the adapted 3A comparison. The current
# production 3B sensitivity remains limited to F3/Fz roles; HUP138's F4 is
# acquired for a future bilateral analysis and is not yet admitted to the
# estimator or group claims. The paper
# derived a cardiac curve per referenced F3/A2 and F4/A1 electrode, then averaged
# the electrode-specific HR-maximum/RR-minimum times).  Fz was not used by Naji.
#
# HUP138 and HUP182 are retained because channel acquisition is intentionally
# broader than the downstream paired intersection. HUP138 became spectrum-
# eligible after the v9 rebuild; HUP182 has C3 but no frozen iEEG spectrum.
# Endpoint eligibility and full-interval activity are decided only after these
# sidecars exist, rather than being hard-coded from an earlier result release.
SCALP_CHANNEL_PLAN = {
    "HUP138_phaseII": {
        "c3": "C3",
        "f3": "F3",
        "f4": "F4",
        "fz": "Fz",
    },
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

# The validator consumes existing bytes but cannot change bytes emitted by this
# producer, so paired_scalp_sidecar_validation.py is intentionally not included.
_DEPENDENCY_FILES = (
    "analysis/cache_paired_scalp.py",
    "analysis/cache_lc_series.py",
    "analysis/hup_portal.py",
    "analysis/infraslow_rr_sigma_coherence.py",
    "analysis/pipeline_version.py",
    "analysis/signal_qc.py",
    "analysis/staging_helpers.py",
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


def _read_pinned_ieeg_manifest(cache_dir):
    """Read and validate one immutable terminal HUP cache manifest."""
    manifest_path = os.path.join(cache_dir, "RUN_MANIFEST.json")
    if not os.path.isfile(manifest_path):
        raise FileNotFoundError(
            f"missing pinned iEEG cache manifest: {manifest_path}")
    with open(manifest_path, "rb") as handle:
        manifest_bytes = handle.read()
    try:
        manifest = json.loads(manifest_bytes)
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"{manifest_path} is not valid JSON") from exc
    partition = validate_terminal_run_manifest(
        manifest, source=manifest_path)
    config = manifest.get("config")
    producer_digest = (
        config.get("cache_code_sha256") if isinstance(config, dict) else None)
    if (
        partition["run_state"] != "complete"
        or manifest.get("pipeline") != "cache_lc_series"
        or manifest.get("analysis_version") != ANALYSIS_VERSION
        or manifest.get("cache_schema_version") != IEEG_CACHE_SCHEMA
        or not isinstance(producer_digest, str)
        or len(producer_digest) != 64
        or any(
            character not in "0123456789abcdef"
            for character in producer_digest
        )
        or config.get("cache_schema_version") != IEEG_CACHE_SCHEMA
    ):
        raise RuntimeError(
            f"{manifest_path} is not an exact completed current HUP cache run")
    return {
        "path": manifest_path,
        "sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "manifest": manifest,
        "partition": partition,
        "producer_digest": producer_digest,
    }


def _validate_pinned_ieeg_subject(subject, cache_dir, frozen_manifest):
    """Validate one cache against already frozen terminal-manifest bytes."""
    partition = frozen_manifest["partition"]
    if subject not in set(partition["completed"]):
        raise RuntimeError(f"{subject} is not completed in the pinned cache manifest")
    cache_path = os.path.join(cache_dir, f"{subject}.npz")
    if not os.path.isfile(cache_path):
        raise FileNotFoundError(f"missing pinned iEEG cache input for {subject}")
    expected = partition["result_files_sha256"][subject]
    actual = file_sha256(cache_path)
    if actual != expected:
        raise RuntimeError(f"{subject} cache bytes differ from the terminal manifest")
    with np.load(cache_path, allow_pickle=False) as cache:
        if npz_scalar_text(cache, "subject") != subject:
            raise RuntimeError(f"{cache_path} embeds a different participant")
        if npz_scalar_text(cache, "status") != "ok":
            raise RuntimeError(f"{cache_path} is not an OK cache")
        if npz_scalar_text(cache, "cache_schema_version") != IEEG_CACHE_SCHEMA:
            raise RuntimeError(f"{cache_path} is not a v9 neutral cache")
        if (
            npz_scalar_text(cache, "cache_code_sha256")
            != frozen_manifest["producer_digest"]
        ):
            raise RuntimeError(
                f"{cache_path} producer digest differs from its terminal manifest")
        validate_completed_cache_failures(
            cache, source=cache_path, require_ecg=True)
        required = (
            "night_s",
            "hours",
            "sf",
            "source_identity_json",
            "acquisition_sample_counts_json",
        )
        missing = [key for key in required if key not in cache.files]
        if missing:
            raise RuntimeError(f"{cache_path} lacks frozen fields: {missing}")
        acquisition = validate_full_interval_acquisition(
            cache,
            source=cache_path,
            core_purpose="analysis_core",
            subrequest_purpose="analysis_subrequest",
            core_chunk_s=CHUNK_S,
            filter_edge_s=FILTER_EDGE_S,
        )
        try:
            source_identity = json.loads(
                npz_scalar_text(cache, "source_identity_json"))
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            raise RuntimeError(
                f"{cache_path} has malformed source identity") from exc
        if not isinstance(source_identity, dict):
            raise RuntimeError(
                f"{cache_path} source identity must be a JSON object")
        pinned = {
            "night_s": float(np.asarray(cache["night_s"]).item()),
            "hours": float(np.asarray(cache["hours"]).item()),
            "sf": float(acquisition["sample_rate_hz"]),
            "source_identity": source_identity,
        }
    return {
        "path": cache_path,
        "sha256": actual,
        "manifest_path": frozen_manifest["path"],
        "manifest_sha256": frozen_manifest["sha256"],
        "manifest_run_id": frozen_manifest["manifest"].get("run_id"),
        "cache_code_sha256": frozen_manifest["producer_digest"],
        "validated_analysis_core_count": acquisition["core_count"],
        **pinned,
    }


def freeze_pinned_ieeg_inputs(subjects, cache_dir=DEFAULT_IEEG_CACHE):
    """Freeze one base manifest and every requested HUP cache identity."""
    cache_dir = os.path.abspath(cache_dir)
    subjects = list(subjects)
    if not subjects or len(subjects) != len(set(subjects)):
        raise RuntimeError(
            "pinned iEEG subjects must be a nonempty list without duplicates")
    frozen_manifest = _read_pinned_ieeg_manifest(cache_dir)
    pinned = {
        subject: _validate_pinned_ieeg_subject(
            subject, cache_dir, frozen_manifest)
        for subject in subjects
    }
    return {
        "cache_dir": cache_dir,
        "manifest_path": frozen_manifest["path"],
        "manifest_sha256": frozen_manifest["sha256"],
        "manifest_run_id": frozen_manifest["manifest"].get("run_id"),
        "subjects": pinned,
    }


def assert_pinned_ieeg_inputs_unchanged(frozen, subjects=None):
    """Fail if a frozen base manifest or any selected cache changed in place."""
    selected = (
        list(frozen["subjects"])
        if subjects is None
        else list(subjects)
    )
    if file_sha256(frozen["manifest_path"]) != frozen["manifest_sha256"]:
        raise RuntimeError(
            "pinned iEEG terminal manifest changed during paired-scalp work")
    for subject in selected:
        pinned = frozen["subjects"].get(subject)
        if pinned is None:
            raise RuntimeError(f"{subject} is absent from frozen iEEG inputs")
        if file_sha256(pinned["path"]) != pinned["sha256"]:
            raise RuntimeError(
                f"{subject} pinned iEEG cache changed during paired-scalp work")
    return True


def assert_pinned_ieeg_subject_unchanged(pinned):
    """Check one frozen subject for direct ``build_subject`` callers."""
    if file_sha256(pinned["manifest_path"]) != pinned["manifest_sha256"]:
        raise RuntimeError(
            "pinned iEEG terminal manifest changed during paired-scalp work")
    if file_sha256(pinned["path"]) != pinned["sha256"]:
        raise RuntimeError(
            "pinned iEEG cache changed during paired-scalp work")
    return True


def validate_pinned_ieeg_cache(subject, cache_dir=DEFAULT_IEEG_CACHE):
    """Validate one exact current cache without requiring today's source tree."""
    frozen = freeze_pinned_ieeg_inputs([subject], cache_dir)
    return frozen["subjects"][subject]


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
        if not same_series_geometry(identity, reference):
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


def _manifest_config(cache_dir, *, pinned_manifest_sha256=None):
    """Return the complete current sidecar-run configuration contract."""
    base_manifest = os.path.join(os.path.abspath(cache_dir), "RUN_MANIFEST.json")
    if not os.path.isfile(base_manifest):
        raise FileNotFoundError(
            f"paired-scalp base manifest is missing: {base_manifest}")
    actual_manifest_sha256 = file_sha256(base_manifest)
    if (
        pinned_manifest_sha256 is not None
        and actual_manifest_sha256 != pinned_manifest_sha256
    ):
        raise RuntimeError(
            "pinned iEEG terminal manifest changed before paired-scalp "
            "configuration was frozen")
    return {
        "ieeg_cache_directory_relative": os.path.relpath(
            os.path.abspath(cache_dir), ROOT).replace(os.sep, "/"),
        "ieeg_cache_manifest_sha256": actual_manifest_sha256,
        "ieeg_cache_schema": IEEG_CACHE_SCHEMA,
        "scalp_cache_schema": SCALP_CACHE_SCHEMA,
        "filter_edge_seconds": FILTER_EDGE_S,
        "chunk_seconds": CHUNK_S,
        "ecg_reused_not_redetected": True,
        "staging_reused_not_recomputed": True,
    }


def _recorded_commit_exists(value):
    """Compatibility wrapper for callers that inspect recorded revisions."""
    return _recorded_commit_exists_impl(value, root=ROOT)


def _dependency_sha256_at_revision(revision):
    """Hash dependency bytes stored by one recorded Git commit."""
    return _dependency_sha256_at_revision_impl(
        revision,
        root=ROOT,
        dependency_files=_DEPENDENCY_FILES,
        commit_exists=_recorded_commit_exists,
    )


def _validation_contract():
    """Bind validator behavior to producer-owned byte-affecting constants."""
    return SidecarValidationContract(
        analysis_version=ANALYSIS_VERSION,
        ieeg_cache_schema=IEEG_CACHE_SCHEMA,
        scalp_cache_schema=SCALP_CACHE_SCHEMA,
        pipeline=PIPELINE,
        channel_plan=SCALP_CHANNEL_PLAN,
        filter_edge_s=FILTER_EDGE_S,
        chunk_s=CHUNK_S,
        historical_power_support=HISTORICAL_POWER_SUPPORT,
        sigma_band_hz=tuple(SIGMA_FIXED),
        swa_band_hz=tuple(SWA_BAND_L),
        so_band_hz=tuple(SO_BAND_NAJI),
    )


def validate_terminal_sidecar_manifest(
        manifest, *, source, output_dir, cache_dir,
        expected_requested=None):
    """Prove one exact, complete current paired-scalp terminal manifest."""
    return _validate_terminal_manifest_impl(
        manifest,
        source=source,
        output_dir=output_dir,
        cache_dir=cache_dir,
        expected_requested=expected_requested,
        contract=_validation_contract(),
        dependency_sha256=cache_dependency_sha256,
        manifest_config=_manifest_config,
        current_runtime_versions=runtime_versions,
        commit_exists=_recorded_commit_exists,
        dependency_sha256_for_revision=_dependency_sha256_at_revision,
    )


def _validate_power_support(cache, *, path, n_channels, total_s, sf):
    """Compatibility wrapper around the contract-bound NPZ validator."""
    return _validate_power_support_impl(
        cache,
        path=path,
        n_channels=n_channels,
        total_s=total_s,
        sf=sf,
        historical_power_support=HISTORICAL_POWER_SUPPORT,
    )


def validate_scalp_sidecar_payload(
        cache, *, path, subject, pinned, dependency_digest):
    """Validate one exact result-producing scalp sidecar payload."""
    return _validate_scalp_sidecar_payload_impl(
        cache,
        path=path,
        subject=subject,
        pinned=pinned,
        dependency_digest=dependency_digest,
        contract=_validation_contract(),
        normalize_scalp_label=_normalize_scalp_label,
    )


def _validate_reused_sidecar_payload(
        cache, *, path, subject, pinned, dependency_digest):
    """Backward-compatible reuse wrapper around the shared exact validator."""
    return validate_scalp_sidecar_payload(
        cache,
        path=path,
        subject=subject,
        pinned=pinned,
        dependency_digest=dependency_digest,
    )


def validate_reusable_sidecar_run(
        output_dir, requested, *, cache_dir=DEFAULT_IEEG_CACHE,
        frozen_inputs=None):
    """Prove a complete prior sidecar run before reusing any existing bytes."""
    return _validate_reusable_sidecar_run_impl(
        output_dir,
        requested,
        cache_dir=cache_dir,
        frozen_inputs=frozen_inputs,
        freeze_pinned_inputs=freeze_pinned_ieeg_inputs,
        assert_pinned_inputs_unchanged=assert_pinned_ieeg_inputs_unchanged,
        validate_terminal_manifest=validate_terminal_sidecar_manifest,
        validate_sidecar_payload=_validate_reused_sidecar_payload,
    )


def build_subject(subject, *, cache_dir=DEFAULT_IEEG_CACHE,
                  output_dir=DEFAULT_OUTPUT, force=False,
                  validated_reuse=None, pinned=None,
                  dependency_digest=None):
    if subject not in SCALP_CHANNEL_PLAN:
        raise ValueError(f"{subject} is not in the explicit paired-scalp plan")
    if pinned is None:
        pinned = validate_pinned_ieeg_cache(subject, cache_dir)
    elif pinned.get("path") != os.path.join(
            os.path.abspath(cache_dir), f"{subject}.npz"):
        raise RuntimeError(
            f"{subject} frozen iEEG input does not match the cache directory")
    assert_pinned_ieeg_subject_unchanged(pinned)
    output_path = os.path.join(output_dir, f"{subject}.npz")
    if dependency_digest is None:
        dependency_digest = cache_dependency_sha256()
    elif dependency_digest != cache_dependency_sha256():
        raise RuntimeError(
            "paired-scalp source changed after its dependency digest was frozen")
    if os.path.exists(output_path) and not force:
        if validated_reuse is None:
            validated_reuse = validate_reusable_sidecar_run(
                output_dir, [subject], cache_dir=cache_dir)
        proven = validated_reuse.get("subjects", {}).get(subject)
        if (
            proven is None
            or proven.get("path") != output_path
            or proven.get("sha256") != file_sha256(output_path)
        ):
            raise RuntimeError(
                f"{output_path} was not proven by the exact prior terminal "
                "manifest; rerun with --force")
        print(f"[{subject}] validated existing scalp sidecar", flush=True)
        return output_path, "reused"

    session = sess()
    ds = None
    try:
        ds = session.open_dataset(subject)
        labels = list(ds.get_channel_labels())
        resolved = _resolve_channels(labels, SCALP_CHANNEL_PLAN[subject])
        source_identity = _portal_source_identity(ds, resolved, pinned)
        channels = list(resolved.values())
        channel_indices = [labels.index(channel) for channel in channels]
        sf = float(require_integer_sample_rate(
            pinned["sf"], source=f"{subject} pinned iEEG cache"))
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
        acquisition_sample_counts = []

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
                data = pull_continuous_exact(
                    ds,
                    channel_indices,
                    night_s + pull_start,
                    pull_stop - pull_start,
                    sf,
                    records=acquisition_sample_counts,
                    purpose="paired_scalp_subrequest",
                )
                if data.ndim != 2 or data.shape[1] != n_channels:
                    raise RuntimeError(
                        f"portal returned unexpected block shape {data.shape}")
                core_start, requested_core = portal_core_sample_geometry(
                    t, pull_start, duration, sf)
                returned_core = max(
                    0, min(len(data), core_start + requested_core) - core_start)
                acquisition_sample_counts.append({
                    "purpose": "paired_scalp_core",
                    "analysis_start_s": float(t),
                    "analysis_duration_s": float(duration),
                    "requested_sample_count": requested_core,
                    "returned_sample_count": int(returned_core),
                    "status": (
                        "ok" if returned_core == requested_core
                        else "sample_count_mismatch"
                    ),
                })
                if returned_core != requested_core:
                    raise RuntimeError(
                        f"portal returned {returned_core} core samples; "
                        f"expected {requested_core} at t={t:g}")
            except Exception as exc:
                failure = {
                    "start_s": float(t),
                    "duration_s": float(duration),
                    "error": f"{type(exc).__name__}: {exc}",
                }
                failed_chunks.append(failure)
                require_complete_acquisition(failed_chunks)

            core_count = requested_core
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
        require_complete_acquisition(failed_chunks)
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
            "ieeg_cache_schema_version": IEEG_CACHE_SCHEMA,
            "failed_chunks_json": json.dumps(
                failed_chunks, sort_keys=True),
            "acquisition_sample_counts_json": json.dumps(
                acquisition_sample_counts, sort_keys=True),
            "acquisition_sample_count_semantics": (
                "every bounded portal subrequest and non-overlapping scalp "
                "analysis core records exact requested/returned counts; any "
                "mismatch is fatal"),
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
        assert_pinned_ieeg_subject_unchanged(pinned)
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


def _manifest_base(
        *, requested, config, run_id, run_state, dependency_digest=None):
    """Build shared paired-scalp run identity without terminal claims."""
    requested = list(requested)
    dependency_digest = (
        cache_dependency_sha256()
        if dependency_digest is None
        else str(dependency_digest)
    )
    return {
        "run_id": str(run_id),
        "schema_version": SCALP_CACHE_SCHEMA,
        "pipeline": PIPELINE,
        "run_state": run_state,
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": IEEG_CACHE_SCHEMA,
        "generated_at_utc": utc_now(),
        "code_revision": git_revision(ROOT),
        "code_dirty": git_is_dirty(ROOT),
        "runtime_versions": runtime_versions(),
        "cache_dependency_sha256": dependency_digest,
        "requested": requested,
        "channel_plan": {
            subject: SCALP_CHANNEL_PLAN[subject] for subject in requested
        },
        "config": dict(config),
    }


def _write_in_progress_manifest(
        output_dir, *, requested, config, run_id, dependency_digest=None):
    """Invalidate any prior terminal sidecar manifest before file writes."""
    os.makedirs(output_dir, exist_ok=True)
    payload = _manifest_base(
        requested=requested,
        config=config,
        run_id=run_id,
        run_state="in_progress",
        dependency_digest=dependency_digest,
    )
    payload.update({
        "completed": [],
        "reused": [],
        "failed": [],
        "result_files_sha256": {},
    })
    atomic_json_dump(payload, os.path.join(output_dir, "RUN_MANIFEST.json"))
    return payload


def _write_manifest(output_dir, *, requested, completed, reused, failed,
                    config, run_id, dependency_digest=None):
    os.makedirs(output_dir, exist_ok=True)
    requested = list(requested)
    completed = list(completed)
    reused = list(reused)
    failed = list(failed)
    failed_subjects = [
        value.get("subject") if isinstance(value, dict) else None
        for value in failed
    ]
    partitions = (set(completed), set(reused), set(failed_subjects))
    if (
        len(requested) != len(set(requested))
        or len(completed) != len(set(completed))
        or len(reused) != len(set(reused))
        or any(value is None for value in failed_subjects)
        or len(failed_subjects) != len(set(failed_subjects))
        or any(
            partitions[left] & partitions[right]
            for left, right in ((0, 1), (0, 2), (1, 2))
        )
        or set.union(*partitions) != set(requested)
    ):
        raise RuntimeError(
            "paired-scalp terminal subject partitions are invalid")
    result_hashes = {
        subject: file_sha256(os.path.join(output_dir, f"{subject}.npz"))
        for subject in completed + reused
    }
    payload = _manifest_base(
        requested=requested,
        config=config,
        run_id=run_id,
        run_state="complete" if not failed else "failed",
        dependency_digest=dependency_digest,
    )
    payload.update({
        "completed": completed,
        "reused": reused,
        "failed": failed,
        "result_files_sha256": result_hashes,
    })
    atomic_json_dump(payload, os.path.join(output_dir, "RUN_MANIFEST.json"))
    return payload


def main():
    parser = argparse.ArgumentParser(
        description="Cache simultaneous scalp channels on frozen HUP intervals")
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

    output_dir = os.path.abspath(args.output_dir)
    cache_dir = os.path.abspath(args.ieeg_cache)
    dependency_digest = cache_dependency_sha256()
    frozen_inputs = freeze_pinned_ieeg_inputs(subjects, cache_dir)
    config = _manifest_config(
        cache_dir,
        pinned_manifest_sha256=frozen_inputs["manifest_sha256"],
    )
    existing = [
        subject for subject in subjects
        if os.path.exists(os.path.join(output_dir, f"{subject}.npz"))
    ]
    validated_reuse = None
    if existing and not args.force:
        validated_reuse = validate_reusable_sidecar_run(
            output_dir,
            subjects,
            cache_dir=cache_dir,
            frozen_inputs=frozen_inputs,
        )

    run_id = str(uuid.uuid4())
    _write_in_progress_manifest(
        output_dir,
        requested=subjects,
        config=config,
        run_id=run_id,
        dependency_digest=dependency_digest,
    )
    completed, reused, failed = [], [], []
    for subject in subjects:
        try:
            _, status = build_subject(
                subject,
                cache_dir=cache_dir,
                output_dir=output_dir,
                force=args.force,
                validated_reuse=validated_reuse,
                pinned=frozen_inputs["subjects"][subject],
                dependency_digest=dependency_digest,
            )
            (reused if status == "reused" else completed).append(subject)
        except Exception as exc:
            traceback.print_exc()
            failed.append({
                "subject": subject,
                "error": f"{type(exc).__name__}: {exc}",
            })
    assert_pinned_ieeg_inputs_unchanged(frozen_inputs, subjects)
    if cache_dependency_sha256() != dependency_digest:
        raise RuntimeError(
            "paired-scalp source changed during acquisition; terminal "
            "publication is refused")
    _write_manifest(
        output_dir,
        requested=subjects,
        completed=completed,
        reused=reused,
        failed=failed,
        config=config,
        run_id=run_id,
        dependency_digest=dependency_digest,
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
        f"{output_dir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
