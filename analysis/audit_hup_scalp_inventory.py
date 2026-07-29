"""Freeze a complete scalp-channel inventory for the pinned HUP cohort.

The paired scalp analysis must not infer that an unrequested channel was absent
from a sidecar containing only positively selected channels.  This audit starts
from the complete frozen HUP cache-manifest universe, opens each exact pinned
portal snapshot, and records the full ordered label list plus geometry for every
relevant 10-20 scalp/reference candidate.

This is a metadata audit: it does not stream EEG or inspect endpoint values.
Numerical activity for result-producing channels is established separately by
the byte-pinned full-interval scalp sidecars.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import traceback

import numpy as np

from artifact_contracts import SCALP_INVENTORY_SCHEMA
from cache_paired_scalp import (
    DEFAULT_IEEG_CACHE,
    DEFAULT_OUTPUT as DEFAULT_SCALP_CACHE,
    IEEG_CACHE_SCHEMA,
    SCALP_CACHE_SCHEMA,
    _channel_identity,
    _normalize_scalp_label,
    _same_series_geometry,
    _source_node_map,
    cache_dependency_sha256,
)
from infraslow_rr_sigma_coherence import sess
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    atomic_json_dump,
    file_sha256,
    git_is_dirty,
    git_revision,
    runtime_versions,
    source_tree_sha256,
    utc_now,
)
from qc_profiles import load_qc_profile, qc_profile_sha256


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_QC_GRID = os.path.join(
    ROOT, "outputs", "qc_grid_public", "locked",
    "overlap11_endpoint_local__hup.json")
DEFAULT_SOURCE_PIN = os.path.join(ROOT, "analysis", "hup_ieeg_source_pin.json")
DEFAULT_OUTPUT = os.path.join(ROOT, "outputs", "paired_scalp_inventory")
PROFILE_ID = "overlap11_endpoint_local"
SCHEMA = SCALP_INVENTORY_SCHEMA
PIPELINE = "audit_hup_scalp_inventory"

# Exact normalized labels only.  This does not treat shaft names such as LC3 as
# C3, and normalization removes harmless numerical zero padding (C03 -> C3).
ROLE_TARGETS = {
    "c3": "C3",
    "c4": "C4",
    "f3": "F3",
    "f4": "F4",
    "fz": "FZ",
    "a1": "A1",
    "a2": "A2",
    "m1": "M1",
    "m2": "M2",
}


def _load_json(path):
    with open(path) as handle:
        return json.load(handle)


def _role_matches(labels):
    normalized = {}
    for label in labels:
        normalized.setdefault(_normalize_scalp_label(label), []).append(label)
    return {
        role: list(normalized.get(target, []))
        for role, target in ROLE_TARGETS.items()
    }


def _frozen_profile_subjects(payload):
    matches = [
        profile for profile in payload.get("profiles", [])
        if profile.get("qc_profile_id") == PROFILE_ID
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"expected exactly one frozen QC-grid profile {PROFILE_ID!r}")
    subjects = matches[0].get("subjects", [])
    if len(subjects) != len({
            value.get("subject") for value in subjects}):
        raise RuntimeError("frozen QC-grid profile contains duplicate subjects")
    return {
        value["subject"]: value
        for value in subjects
    }, matches[0]


def _verify_pinned_snapshot(ds, subject, source_pin):
    expected = source_pin["datasets"][subject]
    actual_snapshot = str(getattr(ds, "snap_id", "") or "")
    if actual_snapshot != str(expected["snapshot_id"]):
        raise RuntimeError(
            f"{subject} snapshot mismatch: expected "
            f"{expected['snapshot_id']!r}, got {actual_snapshot!r}")
    labels = [str(value) for value in ds.get_channel_labels()]
    nodes = _source_node_map(ds)
    for label, expected_compact in expected["channels"].items():
        if label not in labels:
            raise RuntimeError(
                f"{subject} pinned channel {label!r} is absent")
        identity = _channel_identity(ds, label, nodes)
        actual_compact = [
            identity["revision_id"],
            identity["data_check"],
        ]
        if actual_compact != expected_compact:
            raise RuntimeError(
                f"{subject} pinned channel {label!r} revision/dataCheck "
                "changed")
    reference = next(
        iter(expected.get("cortical_channels", [])),
        expected.get("ekg"),
    )
    if reference is None:
        raise RuntimeError(f"{subject} source pin has no reference channel")
    return labels, nodes, reference, _channel_identity(ds, reference, nodes)


def _load_sidecar_manifest(sidecar_dir, ieeg_manifest_sha256):
    path = os.path.join(sidecar_dir, "RUN_MANIFEST.json")
    if not os.path.isfile(path):
        return None, path
    manifest = _load_json(path)
    if (
        manifest.get("pipeline") != "cache_paired_scalp"
        or manifest.get("run_state") != "complete"
        or manifest.get("schema_version") != SCALP_CACHE_SCHEMA
    ):
        raise RuntimeError("paired scalp sidecar manifest is not terminal/current")
    if manifest.get("cache_dependency_sha256") != cache_dependency_sha256():
        raise RuntimeError("paired scalp sidecar manifest has a stale source digest")
    if (
        manifest.get("config", {}).get("ieeg_cache_manifest_sha256")
        != ieeg_manifest_sha256
    ):
        raise RuntimeError(
            "paired scalp sidecars do not reference the frozen iEEG manifest")
    return manifest, path


def _sidecar_evidence(subject, sidecar_dir, manifest):
    if manifest is None:
        return None
    expected = manifest.get("result_files_sha256", {}).get(subject)
    if not expected:
        return None
    path = os.path.join(sidecar_dir, f"{subject}.npz")
    if not os.path.isfile(path) or file_sha256(path) != expected:
        raise RuntimeError(f"{subject} sidecar bytes do not match its manifest")
    with np.load(path, allow_pickle=False) as cache:
        roles = json.loads(str(np.asarray(
            cache["channel_roles_json"]).item()))
        channels = [str(value) for value in cache["scalp_chans"]]
        nonflat = np.asarray(cache["scalp_signal_nonflat_mask"], bool)
        finite = np.asarray(cache["scalp_signal_finite_sample_count"], int)
        dynamic = np.asarray(cache["scalp_signal_raw_dynamic_range"], float)
        if not (
            len(channels) == len(nonflat) == len(finite) == len(dynamic)
        ):
            raise RuntimeError(f"{subject} sidecar activity arrays do not align")
        by_channel = {
            channel: {
                "numerically_nonflat_full_interval": bool(nonflat[index]),
                "finite_sample_count": int(finite[index]),
                "raw_dynamic_range": float(dynamic[index]),
            }
            for index, channel in enumerate(channels)
        }
        return {
            "path_relative": os.path.relpath(path, ROOT).replace(os.sep, "/"),
            "sha256": expected,
            "roles": roles,
            "channels": by_channel,
        }


def _classification(record):
    if record.get("query_status") != "ok":
        return "inventory_error"
    c3 = record["roles"]["c3"]
    if len(c3["matches"]) == 0:
        return "no_c3_or_c03_label"
    if len(c3["matches"]) != 1:
        return "ambiguous_c3_or_c03_labels"
    if not c3["matches"][0]["same_geometry_as_pinned_reference"]:
        return "c3_geometry_mismatch"
    if not record["frozen_3a"]["spectrum_available"]:
        return "c3_present_but_frozen_3a_spectrum_unavailable"
    sidecar = record.get("sidecar")
    if sidecar is None:
        return "requires_full_interval_activity_sidecar"
    label = c3["matches"][0]["label"]
    activity = sidecar["channels"].get(label)
    if activity is None:
        return "c3_not_streamed_in_sidecar"
    if not activity["numerically_nonflat_full_interval"]:
        return "c3_numerically_flat"
    return "paired_3a_eligible"


def build_inventory(
    *,
    cache_dir=DEFAULT_IEEG_CACHE,
    sidecar_dir=DEFAULT_SCALP_CACHE,
    qc_grid_path=DEFAULT_QC_GRID,
    source_pin_path=DEFAULT_SOURCE_PIN,
):
    cache_manifest_path = os.path.join(cache_dir, "RUN_MANIFEST.json")
    cache_manifest = _load_json(cache_manifest_path)
    if (
        cache_manifest.get("run_state") != "complete"
        or cache_manifest.get("pipeline") != "cache_lc_series"
        or cache_manifest.get("cache_schema_version") != IEEG_CACHE_SCHEMA
    ):
        raise RuntimeError(
            "HUP cache manifest is not a terminal current-schema "
            "cache_lc_series run")
    requested = [str(value) for value in cache_manifest.get("requested", [])]
    if len(requested) != len(set(requested)) or not requested:
        raise RuntimeError("frozen HUP cache requested universe is invalid")

    source_pin = _load_json(source_pin_path)
    if set(source_pin.get("datasets", {})) != set(requested):
        raise RuntimeError(
            "source-pin subjects do not exactly match the frozen HUP universe")
    qc_grid = _load_json(qc_grid_path)
    qc_subjects, frozen_profile = _frozen_profile_subjects(qc_grid)
    completed = [str(value) for value in cache_manifest.get("completed", [])]
    if set(qc_subjects) != set(completed):
        raise RuntimeError(
            "frozen QC-grid subjects do not exactly match completed cache "
            "subjects")
    skipped = {
        str(value.get("subject"))
        for value in cache_manifest.get("skipped", [])
    }
    if set(requested) != set(completed) | skipped:
        raise RuntimeError(
            "completed and explicitly skipped cache subjects do not partition "
            "the frozen requested universe")
    cache_manifest_sha256 = file_sha256(cache_manifest_path)
    if qc_grid.get("cache_manifest_sha256") != cache_manifest_sha256:
        raise RuntimeError(
            "frozen QC grid was not computed from this exact cache manifest")
    current_profile_sha256 = qc_profile_sha256(
        load_qc_profile(PROFILE_ID))
    if frozen_profile.get("qc_profile_sha256") != current_profile_sha256:
        raise RuntimeError(
            "frozen QC-grid profile differs from the current locked profile")
    sidecar_manifest, sidecar_manifest_path = _load_sidecar_manifest(
        sidecar_dir, cache_manifest_sha256)

    session = sess()
    records = []
    try:
        for subject in requested:
            print(f"[{subject}] scalp-label inventory", flush=True)
            record = {
                "subject": subject,
                "query_status": "error",
                "roles": {},
                "frozen_3a": {
                    "record_available": subject in qc_subjects,
                    "spectrum_available": bool(
                        qc_subjects.get(subject, {})
                        .get("result_3a", {})
                        .get("endpoint_availability", {})
                        .get("spectrum", False)
                    ),
                    "coherence_available": bool(
                        qc_subjects.get(subject, {})
                        .get("result_3a", {})
                        .get("endpoint_availability", {})
                        .get("fixed_0p02_coherence", False)
                    ),
                    "cross_correlation_available": bool(
                        qc_subjects.get(subject, {})
                        .get("result_3a", {})
                        .get("endpoint_availability", {})
                        .get("cross_correlation", False)
                    ),
                },
            }
            ds = None
            try:
                ds = session.open_dataset(subject)
                labels, nodes, reference_label, reference = (
                    _verify_pinned_snapshot(
                        ds, subject, source_pin))
                matches = _role_matches(labels)
                role_records = {}
                for role, role_labels in matches.items():
                    role_records[role] = {
                        "target_normalized_label": ROLE_TARGETS[role],
                        "matches": [
                            {
                                "label": label,
                                "identity": _channel_identity(
                                    ds, label, nodes),
                                "same_geometry_as_pinned_reference": (
                                    _same_series_geometry(
                                        _channel_identity(ds, label, nodes),
                                        reference,
                                    )
                                ),
                            }
                            for label in role_labels
                        ],
                    }
                record.update({
                    "query_status": "ok",
                    "dataset_name": str(getattr(ds, "name", "") or ""),
                    "snapshot_id": str(getattr(ds, "snap_id", "") or ""),
                    "ordered_channel_labels": labels,
                    "ordered_channel_labels_sha256": (
                        hashlib.sha256(
                            json.dumps(
                                labels,
                                ensure_ascii=False,
                                separators=(",", ":"),
                            ).encode("utf-8")
                        ).hexdigest()
                    ),
                    "pinned_reference_channel": reference_label,
                    "pinned_reference_identity": reference,
                    "roles": role_records,
                    "sidecar": _sidecar_evidence(
                        subject, sidecar_dir, sidecar_manifest),
                })
            except Exception as exc:
                record["error"] = f"{type(exc).__name__}: {exc}"
                record["traceback"] = traceback.format_exc()
            finally:
                if ds is not None and hasattr(session, "close_dataset"):
                    session.close_dataset(ds)
            record["selection_classification"] = _classification(record)
            records.append(record)
    finally:
        if hasattr(session, "close"):
            session.close()

    classes = {}
    for record in records:
        classes.setdefault(record["selection_classification"], []).append(
            record["subject"])
    errors = [
        record["subject"] for record in records
        if record["query_status"] != "ok"
    ]
    metadata = {
        "schema_version": SCHEMA,
        "pipeline": PIPELINE,
        "run_state": "complete" if not errors else "failed",
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "generated_at_utc": utc_now(),
        "code_revision": git_revision(ROOT),
        "code_dirty": git_is_dirty(ROOT),
        "source_tree_sha256": source_tree_sha256(ROOT),
        "audit_script_sha256": file_sha256(os.path.abspath(__file__)),
        "runtime_versions": runtime_versions(),
        "selection_was_blind_to_scalp_endpoint_values": True,
        "selection_inputs": (
            "exact frozen cohort membership; pinned portal snapshot and "
            "channel metadata; exact C3/C03 label and geometry; frozen 3A "
            "spectrum availability; prespecified numerical activity QC from "
            "full-interval sidecars"),
        "reference_warning": (
            "A channel label does not identify its online reference. Presence "
            "of A1/A2/M1/M2 also does not prove a reconstructable Naji montage."),
        "requested_subjects": requested,
        "n_requested": len(requested),
        "frozen_cache_manifest_path_relative": os.path.relpath(
            cache_manifest_path, ROOT).replace(os.sep, "/"),
        "frozen_cache_manifest_sha256": cache_manifest_sha256,
        "source_pin_path_relative": os.path.relpath(
            source_pin_path, ROOT).replace(os.sep, "/"),
        "source_pin_sha256": file_sha256(source_pin_path),
        "qc_grid_path_relative": os.path.relpath(
            qc_grid_path, ROOT).replace(os.sep, "/"),
        "qc_grid_sha256": file_sha256(qc_grid_path),
        "qc_profile_id": PROFILE_ID,
        "qc_profile_sha256": frozen_profile["qc_profile_sha256"],
        "current_locked_qc_profile_sha256": current_profile_sha256,
        "sidecar_manifest_path_relative": (
            None if sidecar_manifest is None
            else os.path.relpath(
                sidecar_manifest_path, ROOT).replace(os.sep, "/")
        ),
        "sidecar_manifest_sha256": (
            None if sidecar_manifest is None
            else file_sha256(sidecar_manifest_path)
        ),
        "classification_counts": {
            key: len(value) for key, value in sorted(classes.items())
        },
        "classification_subjects": {
            key: value for key, value in sorted(classes.items())
        },
        "query_error_subjects": errors,
    }
    return {**metadata, "subjects": records}


def main():
    parser = argparse.ArgumentParser(
        description="Audit relevant scalp labels across all pinned HUP snapshots")
    parser.add_argument("--ieeg-cache", default=DEFAULT_IEEG_CACHE)
    parser.add_argument("--scalp-cache", default=DEFAULT_SCALP_CACHE)
    parser.add_argument("--qc-grid", default=DEFAULT_QC_GRID)
    parser.add_argument("--source-pin", default=DEFAULT_SOURCE_PIN)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    payload = build_inventory(
        cache_dir=os.path.abspath(args.ieeg_cache),
        sidecar_dir=os.path.abspath(args.scalp_cache),
        qc_grid_path=os.path.abspath(args.qc_grid),
        source_pin_path=os.path.abspath(args.source_pin),
    )
    inventory_path = os.path.join(output_dir, "hup_scalp_channel_inventory.json")
    atomic_json_dump(payload, inventory_path)
    manifest = {
        key: payload[key]
        for key in (
            "schema_version",
            "pipeline",
            "run_state",
            "analysis_version",
            "cache_schema_version",
            "generated_at_utc",
            "code_revision",
            "code_dirty",
            "source_tree_sha256",
            "audit_script_sha256",
            "runtime_versions",
            "requested_subjects",
            "n_requested",
            "frozen_cache_manifest_sha256",
            "source_pin_sha256",
            "qc_grid_sha256",
            "qc_profile_id",
            "qc_profile_sha256",
            "current_locked_qc_profile_sha256",
            "sidecar_manifest_sha256",
            "classification_counts",
            "classification_subjects",
            "query_error_subjects",
        )
    }
    manifest["inventory_file"] = os.path.basename(inventory_path)
    manifest["inventory_file_sha256"] = file_sha256(inventory_path)
    atomic_json_dump(
        manifest, os.path.join(output_dir, "RUN_MANIFEST.json"))
    print(
        json.dumps(payload["classification_counts"], sort_keys=True),
        flush=True,
    )
    if payload["run_state"] != "complete":
        sys.exit(1)


if __name__ == "__main__":
    main()
