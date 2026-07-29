"""Fast integrity checks for the simultaneous scalp–iEEG sidecar and results."""
from __future__ import annotations

import argparse
import json
import os
import tempfile

import numpy as np

from compact_qc_grid_artifacts import validate_public_artifacts
from paired_artifact_validation import validate_checked_in_paired_evidence
from paired_reporting import _validate_csv_artifact
from cache_paired_scalp import (
    HISTORICAL_POWER_SUPPORT,
    IEEG_CACHE_SCHEMA,
    SCALP_CACHE_SCHEMA,
    SCALP_CHANNEL_PLAN,
    FILTER_EDGE_S,
    SIGMA_FIXED,
    SO_BAND_NAJI,
    SWA_BAND_L,
    _manifest_config,
    _normalize_scalp_label,
    _resolve_channels,
    cache_dependency_sha256,
)
from paired_scalp_ieeg_comparison import (
    STAGES_3B,
    _group_summary,
    _normalized_csv_rows,
    _endpoint_local_3a_materializations,
    _parse_subjects,
    _role_pair_csv_rows,
    _validate_scalp_inventory,
    _validate_scalp_cache,
)
from pipeline_version import (
    ANALYSIS_VERSION,
    atomic_savez,
    file_sha256,
    git_revision,
    runtime_versions,
    source_tree_sha256,
    utc_now,
)
from staging_helpers import CHUNK_S


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--artifact-validation",
        choices=("offline", "publication"),
        default="publication",
        help=(
            "offline validates every checked-in artifact and skips only absent "
            "data/derived manifests; publication also requires the live private "
            "cache lineage"
        ),
    )
    return parser.parse_args()


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


ARGS = parse_args()
PUBLIC_QC_EVIDENCE = validate_public_artifacts()
check(
    "all compact QC summaries and the full locked-profile snapshot validate",
    len(PUBLIC_QC_EVIDENCE) == 9,
)
PUBLIC_EVIDENCE = validate_checked_in_paired_evidence(ROOT)
check(
    "checked-in paired inventory and result evidence is complete and byte-pinned",
    bool(PUBLIC_EVIDENCE["subject_results"]["subjects"]),
)
if ARGS.artifact_validation == "offline":
    print(
        "  INFO  non-publication CI skipped absent ignored lineage: "
        + ", ".join(PUBLIC_EVIDENCE["skipped_derived_lineage"])
    )
else:
    check(
        "publication validation has every private derived manifest",
        not PUBLIC_EVIDENCE["skipped_derived_lineage"],
    )


check(
    "scalp labels normalize case and harmless zero padding only",
    _normalize_scalp_label(" c03 ") == "C3"
    and _normalize_scalp_label("Fz") == "FZ",
)
check(
    "the explicit role plan resolves a unique portal label",
    _resolve_channels(
        ["C03", "F3", "Fz"], {"c3": "C3", "f3": "F3", "fz": "FZ"}
    )
    == {"c3": "C03", "f3": "F3", "fz": "Fz"},
)
check(
    "no-argument comparison defers subject selection to the audited inventory",
    _parse_subjects(None) is None,
)
try:
    _resolve_channels(["C3", "C03"], {"c3": "C3"})
    ambiguous_rejected = False
except RuntimeError:
    ambiguous_rejected = True
check("ambiguous zero-padded scalp labels fail closed", ambiguous_rejected)


n_seconds = 180
base = np.linspace(1.0, 2.0, n_seconds)
hr = np.linspace(55.0, 65.0, n_seconds)
ieeg_sigma = base.copy()
ieeg_swa = (base * 2).copy()
scalp_sigma = (base * 3).copy()
scalp_swa = (base * 4).copy()
ieeg_sigma[10:20] = np.nan
ieeg_swa[30:40] = np.nan
scalp_sigma[50:60] = np.nan
scalp_swa[70:80] = np.nan
hr[90:100] = np.nan
materialized_stub = {
    "sigma_parietal": ieeg_sigma,
    "swa_parietal": ieeg_swa,
    "hr_1": hr,
    "stage_lab": np.asarray(["N2"] * (n_seconds // 30)),
}
scalp_stub = {
    **materialized_stub,
    "sigma_parietal": scalp_sigma,
    "swa_parietal": scalp_swa,
}
nonzero_hr_profile = {"hr": {"minimum_coverage": 0.8}}
endpoint_materializations, support = _endpoint_local_3a_materializations(
    materialized_stub, scalp_stub, nonzero_hr_profile)
left_eeg = endpoint_materializations["ieeg"]["eeg"]
right_eeg = endpoint_materializations["scalp"]["eeg"]
left_cardiac = endpoint_materializations["ieeg"]["cardiac"]
right_cardiac = endpoint_materializations["scalp"]["cardiac"]
left_eeg_mask = (
    np.isfinite(left_eeg["sigma_parietal"])
    & np.isfinite(left_eeg["swa_parietal"])
)
right_eeg_mask = (
    np.isfinite(right_eeg["sigma_parietal"])
    & np.isfinite(right_eeg["swa_parietal"])
)
left_cardiac_mask = (
    np.isfinite(left_cardiac["sigma_parietal"])
    & np.isfinite(left_cardiac["swa_parietal"])
    & np.isfinite(left_cardiac["hr_1"])
)
right_cardiac_mask = (
    np.isfinite(right_cardiac["sigma_parietal"])
    & np.isfinite(right_cardiac["swa_parietal"])
    & np.isfinite(right_cardiac["hr_1"])
)
check(
    "paired 3A uses exact endpoint-local EEG and cardiac masks",
    np.array_equal(left_eeg_mask, right_eeg_mask)
    and np.array_equal(left_cardiac_mask, right_cardiac_mask)
    and int(left_eeg_mask.sum()) == 140
    and int(left_cardiac_mask.sum()) == 130
    and support["eeg_common"]["n_seconds"] == 140
    and support["cardiac_common"]["n_seconds"] == 130
    and support["eeg_common"]["support_mask_sha256"]
    != support["cardiac_common"]["support_mask_sha256"]
    and support["cardiac_common_is_subset_of_eeg_common"]
    and support["cardiac_common"]["minimum_hr_coverage"] == 0.8
    and not support["cardiac_common"]["passes_hr_coverage_profile"]
    and not left_cardiac["hr_meets_profile"]
    and not right_cardiac["hr_meets_profile"],
)

inventory_path = os.path.join(
    ROOT, "outputs", "paired_scalp_inventory",
    "hup_scalp_channel_inventory.json")
if ARGS.artifact_validation == "publication":
    inventory = _validate_scalp_inventory(inventory_path)
    audited_subjects = inventory["paired_3a_eligible_subjects"]
    hup138_record = next(
        record
        for record in PUBLIC_EVIDENCE["inventory"]["subjects"]
        if record["subject"] == "HUP138_phaseII"
    )
    check(
        "all-cohort inventory proves a fully activity-audited 3A intersection",
        inventory["n_frozen_hup_participants"] == 25
        and bool(audited_subjects)
        and PUBLIC_EVIDENCE["result_manifest"]["requested_subjects"]
        == audited_subjects
        and inventory["participants_with_unique_f3_and_f4_labels"]
        == ["HUP138_phaseII"]
        and hup138_record["selection_classification"]
        in {"paired_3a_eligible", "c3_numerically_flat"},
    )


with tempfile.TemporaryDirectory() as directory:
    subject = "HUP160_phaseII"
    base_dir = os.path.join(directory, "base")
    os.makedirs(base_dir)
    base_manifest_path = os.path.join(base_dir, "RUN_MANIFEST.json")
    with open(base_manifest_path, "w") as handle:
        json.dump({"synthetic": True}, handle)
    cache_path = os.path.join(directory, f"{subject}.npz")
    current_dependency = cache_dependency_sha256()
    sf = 100
    duration_s = 1
    roles = SCALP_CHANNEL_PLAN[subject]
    channels = list(roles.values())
    n_channels = len(channels)
    geometry = {
        "revision_id": "synthetic-revision",
        "data_check": "synthetic-data-check",
        "start_time_us": 0,
        "end_time_us": 1_000_000,
        "duration_us": 1_000_000.0,
        "number_of_samples": sf,
        "sample_rate_hz": float(sf),
    }
    pinned = {
        "sha256": "frozen-ieeg-cache-sha256",
        "manifest_path": base_manifest_path,
        "manifest_sha256": file_sha256(base_manifest_path),
        "manifest_run_id": "frozen-ieeg-run",
        "night_s": 123.0,
        "hours": duration_s / 3600.0,
        "sf": float(sf),
        "source_identity": {
            "snapshot_id": "snapshot-1",
            "channels": {"A1": geometry},
        },
    }
    source_identity = {
        "dataset_name": subject,
        "snapshot_id": "snapshot-1",
        "reference_ieeg_channel": "A1",
        "reference_ieeg_identity": geometry,
        "scalp_role_to_channel": roles,
        "scalp_channels": {channel: geometry for channel in channels},
    }
    denominator = np.full((n_channels, duration_s), sf, dtype=np.int64)
    numerator = np.full((n_channels, duration_s), 2.0)
    power = numerator / denominator
    candidates = {
        f"so_candidate_{field}_{channel}": np.asarray([], dtype=float)
        for channel in channels
        for field in ("t", "down", "up", "p2p")
    }
    atomic_savez(
        cache_path,
        status="ok",
        cache_schema_version=SCALP_CACHE_SCHEMA,
        cache_dependency_sha256=current_dependency,
        subject=subject,
        source_dataset=subject,
        source_kind="iEEG.org API",
        ieeg_cache_sha256=pinned["sha256"],
        ieeg_cache_manifest_sha256=pinned["manifest_sha256"],
        ieeg_cache_manifest_run_id=pinned["manifest_run_id"],
        ieeg_cache_schema_version=IEEG_CACHE_SCHEMA,
        night_s=pinned["night_s"],
        hours=pinned["hours"],
        sf=pinned["sf"],
        failed_chunks_json="[]",
        acquisition_sample_counts_json=json.dumps([
            {
                "purpose": "paired_scalp_subrequest",
                "request_start_s": pinned["night_s"],
                "request_duration_s": duration_s,
                "requested_sample_count": sf,
                "returned_sample_count": sf,
                "requested_channel_count": n_channels,
                "returned_channel_count": n_channels,
                "status": "ok",
            },
            {
                "purpose": "paired_scalp_core",
                "analysis_start_s": 0.0,
                "analysis_duration_s": duration_s,
                "requested_sample_count": sf,
                "returned_sample_count": sf,
                "status": "ok",
            },
        ]),
        source_identity_json=json.dumps(source_identity),
        channel_roles_json=json.dumps(roles),
        scalp_chans=np.asarray(channels),
        filter_edge_seconds=float(FILTER_EDGE_S),
        chunk_seconds=float(CHUNK_S),
        sigma_band_hz=np.asarray(SIGMA_FIXED, float),
        swa_band_hz=np.asarray(SWA_BAND_L, float),
        so_band_hz=np.asarray(SO_BAND_NAJI, float),
        sigma_fixed_by_channel=power,
        swa_by_channel=power,
        sigma_fixed_power_numerator_by_channel=numerator,
        sigma_fixed_clean_sample_count_by_channel=denominator,
        swa_power_numerator_by_channel=numerator,
        swa_clean_sample_count_by_channel=denominator,
        power_samples_per_second=sf,
        power_historical_minimum_clean_fraction_per_second=(
            HISTORICAL_POWER_SUPPORT),
        scalp_signal_nonflat_mask=np.ones(n_channels, dtype=bool),
        scalp_signal_raw_minimum=np.zeros(n_channels),
        scalp_signal_raw_maximum=np.ones(n_channels),
        scalp_signal_raw_dynamic_range=np.ones(n_channels),
        scalp_signal_numerical_flat_tolerance=np.full(
            n_channels, 64 * np.finfo(float).eps),
        scalp_signal_finite_sample_count=np.full(
            n_channels, sf, dtype=np.int64),
        ecg_reused_not_redetected=True,
        staging_reused_not_recomputed=True,
        **candidates,
    )
    manifest_path = os.path.join(directory, "RUN_MANIFEST.json")
    manifest = {
        "run_id": "11111111-1111-4111-8111-111111111111",
        "pipeline": "cache_paired_scalp",
        "run_state": "complete",
        "schema_version": SCALP_CACHE_SCHEMA,
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": IEEG_CACHE_SCHEMA,
        "generated_at_utc": utc_now(),
        "code_revision": git_revision(ROOT),
        "code_dirty": False,
        "runtime_versions": runtime_versions(),
        "cache_dependency_sha256": current_dependency,
        "requested": [subject],
        "completed": [subject],
        "reused": [],
        "failed": [],
        "channel_plan": {subject: roles},
        "config": _manifest_config(base_dir),
        "result_files_sha256": {
            subject: file_sha256(cache_path),
        },
    }
    with open(manifest_path, "w") as handle:
        json.dump(manifest, handle)
    validated = _validate_scalp_cache(subject, directory, pinned)
    check(
        "a byte-pinned current scalp sidecar validates",
        validated["sha256"] == file_sha256(cache_path)
        and validated["cache_dependency_sha256"] == current_dependency,
    )
    manifest["reused"] = [subject]
    with open(manifest_path, "w") as handle:
        json.dump(manifest, handle)
    try:
        _validate_scalp_cache(subject, directory, pinned)
        overlapping_partition_rejected = False
    except RuntimeError:
        overlapping_partition_rejected = True
    check(
        "live comparison rejects overlapping sidecar manifest partitions",
        overlapping_partition_rejected,
    )
    manifest["reused"] = []
    with open(manifest_path, "w") as handle:
        json.dump(manifest, handle)

    with np.load(cache_path, allow_pickle=False) as cache:
        malformed = {key: np.asarray(cache[key]) for key in cache.files}
    malformed.pop("so_candidate_t_C3")
    atomic_savez(cache_path, **malformed)
    manifest["result_files_sha256"][subject] = file_sha256(cache_path)
    with open(manifest_path, "w") as handle:
        json.dump(manifest, handle)
    try:
        _validate_scalp_cache(subject, directory, pinned)
        malformed_payload_rejected = False
    except RuntimeError:
        malformed_payload_rejected = True
    check(
        "live comparison rejects a rehashed incomplete sidecar payload",
        malformed_payload_rejected,
    )

    atomic_savez(cache_path, **{
        **malformed,
        "so_candidate_t_C3": np.asarray([], dtype=float),
    })
    manifest["result_files_sha256"][subject] = file_sha256(cache_path)
    manifest["cache_dependency_sha256"] = "stale-source-digest"
    with open(manifest_path, "w") as handle:
        json.dump(manifest, handle)
    try:
        _validate_scalp_cache(subject, directory, pinned)
        stale_rejected = False
    except RuntimeError:
        stale_rejected = True
    check("a stale scalp-cache source digest fails closed", stale_rejected)


result_dir = os.path.join(ROOT, "outputs", "paired_scalp_ieeg")
result_manifest = PUBLIC_EVIDENCE["result_manifest"]
paired_results = PUBLIC_EVIDENCE["subject_results"]

role_pair_rows = _role_pair_csv_rows(paired_results["subjects"])
normalized_rows = _normalized_csv_rows(
    paired_results["subjects"], role_pair_rows)


actual_normalized_rows = _validate_csv_artifact(
    os.path.join(result_dir, "paired_metrics.csv"),
    normalized_rows,
)
actual_role_rows = _validate_csv_artifact(
    os.path.join(result_dir, "role_pair_metrics.csv"),
    role_pair_rows,
)
check(
    "checked-in normalized CSV exactly matches subject_results-derived rows",
    bool(actual_normalized_rows),
)
check(
    "checked-in role-pair CSV exactly matches subject_results-derived rows",
    bool(actual_role_rows),
)

normalized_keys = [
    (
        row["subject"],
        row["question"],
        row["stage"],
        row["modality"],
    )
    for row in normalized_rows
]
check(
    "normalized CSV has no duplicate subject/question/stage/modality",
    len(normalized_keys) == len(set(normalized_keys))
    and len(actual_normalized_rows) == len({
        (
            row["subject"],
            row["question"],
            row["stage"],
            row["modality"],
        )
        for row in actual_normalized_rows
    }),
)
normalized_ieeg_3b = [
    row for row in normalized_rows
    if row["question"] == "3B" and row["modality"] == "iEEG"
]
check(
    "normalized CSV has exactly one iEEG 3B row per subject-stage",
    len(normalized_ieeg_3b)
    == len(paired_results["subjects"]) * len(STAGES_3B)
    and all(
        row["comparison_role"] == "subject_stage"
        for row in normalized_ieeg_3b
    ),
)
available_scalp_role_keys = {
    (row["subject"], row["stage"], row["comparison_role"])
    for row in role_pair_rows
    if (
        row["question"] == "3B"
        and row["modality"] != "iEEG"
        and row["sensor_available"]
    )
}
ieeg_role_keys = {
    (row["subject"], row["stage"], row["comparison_role"])
    for row in role_pair_rows
    if row["question"] == "3B" and row["modality"] == "iEEG"
}
check(
    "role-pair CSV retains one iEEG mate for every available scalp role",
    available_scalp_role_keys == ieeg_role_keys,
)
check(
    "CSV normalization leaves authoritative group summary unchanged",
    _group_summary(paired_results["subjects"])
    == paired_results["group_summary"],
)

hashes_match = all(
    file_sha256(os.path.join(result_dir, name)) == expected
    for name, expected in result_manifest["result_files_sha256"].items()
)
check("paired result files match their terminal manifest", hashes_match)
check(
    "paired results identify the exact current analysis source tree",
    result_manifest.get("source_tree_sha256") == source_tree_sha256(ROOT),
)
common_support_ok = True
for subject_result in paired_results["subjects"]:
    shared = subject_result["shared_inputs"]
    geometry = shared["paired_3a_geometry_verification"]
    endpoint_support = shared["paired_3a_endpoint_support"]
    spectrum_geometry = geometry["spectrum_peak_negative_control"]
    cardiac_geometry = geometry["coherence_cross_correlation"]
    left = subject_result["ieeg"]["result_3a"]
    right = subject_result["scalp"]["result_3a"]
    common_support_ok &= (
        geometry["status"] == "exact_endpoint_local_geometry"
        and spectrum_geometry["status"] == "exact_shared_eeg_geometry"
        and cardiac_geometry["status"] == "exact_shared_cardiac_geometry"
        and endpoint_support["eeg_common"]["n_seconds"] > 0
        and endpoint_support["cardiac_common"]["n_seconds"]
        <= endpoint_support["eeg_common"]["n_seconds"]
        and endpoint_support[
            "cardiac_common_is_subset_of_eeg_common"]
        and bool(endpoint_support["eeg_common"]["support_mask_sha256"])
        and bool(endpoint_support["cardiac_common"]["support_mask_sha256"])
        and left["n_bouts"] == right["n_bouts"]
        and left["bout_seconds"] == right["bout_seconds"]
    )
    left_coherence = left.get("coherence")
    right_coherence = right.get("coherence")
    common_support_ok &= (
        (left_coherence is None) == (right_coherence is None)
    )
    if left_coherence is not None:
        common_support_ok &= (
            left_coherence["K"] == right_coherence["K"]
            and left_coherence["n_valid"] == right_coherence["n_valid"]
            and np.isclose(
                left_coherence["analytic_threshold"],
                right_coherence["analytic_threshold"],
                rtol=0,
                atol=1e-12,
            )
        )
    left_xcorr = left.get("cross_correlation")
    right_xcorr = right.get("cross_correlation")
    common_support_ok &= (
        (left_xcorr is None) == (right_xcorr is None)
    )
    if left_xcorr is not None:
        common_support_ok &= (
            left_xcorr["n_intervals"] == right_xcorr["n_intervals"]
            and cardiac_geometry[
                "cross_correlation_retained_window_count"]
            == left_xcorr["n_intervals"]
            and bool(cardiac_geometry[
                "cross_correlation_retained_window_starts_sha256"])
        )
check(
    "every saved primary 3A pair has exact endpoint-local support geometry",
    common_support_ok,
)

print("ALL PAIRED-SCALP CHECKS PASSED")
