"""Fast integrity checks for the simultaneous scalp–iEEG sidecar and results."""
from __future__ import annotations

import json
import os
import tempfile

import numpy as np

from cache_paired_scalp import (
    SCALP_CACHE_SCHEMA,
    _normalize_scalp_label,
    _resolve_channels,
    cache_dependency_sha256,
)
from paired_scalp_ieeg_comparison import (
    STAGES_3B,
    _group_summary,
    _normalized_csv_rows,
    _role_pair_csv_rows,
    _shared_3a_materializations,
    _validate_scalp_inventory,
    _validate_scalp_cache,
)
from pipeline_version import atomic_savez, file_sha256, source_tree_sha256


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


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
left_shared, right_shared, support = _shared_3a_materializations(
    materialized_stub,
    scalp_stub,
    {"hr": {"minimum_coverage": 0.0}},
)
left_mask = (
    np.isfinite(left_shared["sigma_parietal"])
    & np.isfinite(left_shared["swa_parietal"])
    & np.isfinite(left_shared["hr_1"])
)
right_mask = (
    np.isfinite(right_shared["sigma_parietal"])
    & np.isfinite(right_shared["swa_parietal"])
    & np.isfinite(right_shared["hr_1"])
)
check(
    "disjoint missingness is intersected into one exact paired 3A mask",
    np.array_equal(left_mask, right_mask)
    and int(left_mask.sum()) == n_seconds - 50
    and support["n_shared_seconds"] == int(left_mask.sum()),
)

inventory_path = os.path.join(
    ROOT, "outputs", "paired_scalp_inventory",
    "hup_scalp_channel_inventory.json")
if os.path.isfile(inventory_path):
    audited_subjects = [
        "HUP160_phaseII",
        "HUP185_phaseII",
        "HUP187_phaseII",
        "HUP191_phaseII",
        "HUP199_phaseII",
        "HUP205_phaseII",
        "HUP211_phaseII",
        "HUP212_phaseII",
    ]
    inventory = _validate_scalp_inventory(
        inventory_path, audited_subjects)
    check(
        "all-cohort inventory proves the exact eight-subject 3A intersection",
        inventory["n_frozen_hup_participants"] == 25
        and inventory["paired_3a_eligible_subjects"] == audited_subjects
        and inventory["participants_with_unique_f3_and_f4_labels"]
        == ["HUP138_phaseII"],
    )


with tempfile.TemporaryDirectory() as directory:
    subject = "HUPTEST_phaseII"
    cache_path = os.path.join(directory, f"{subject}.npz")
    current_dependency = cache_dependency_sha256()
    pinned = {
        "sha256": "frozen-ieeg-cache-sha256",
        "manifest_sha256": "frozen-ieeg-manifest-sha256",
        "night_s": 123.0,
        "hours": 7.0,
        "sf": 512.0,
    }
    atomic_savez(
        cache_path,
        status="ok",
        cache_schema_version=SCALP_CACHE_SCHEMA,
        cache_dependency_sha256=current_dependency,
        subject=subject,
        ieeg_cache_sha256=pinned["sha256"],
        ieeg_cache_manifest_sha256=pinned["manifest_sha256"],
        night_s=pinned["night_s"],
        hours=pinned["hours"],
        sf=pinned["sf"],
        failed_chunks_json="[]",
        ecg_reused_not_redetected=True,
        staging_reused_not_recomputed=True,
    )
    manifest_path = os.path.join(directory, "RUN_MANIFEST.json")
    manifest = {
        "pipeline": "cache_paired_scalp",
        "run_state": "complete",
        "schema_version": SCALP_CACHE_SCHEMA,
        "cache_dependency_sha256": current_dependency,
        "completed": [subject],
        "reused": [],
        "config": {
            "ieeg_cache_manifest_sha256": pinned["manifest_sha256"],
        },
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
result_manifest_path = os.path.join(result_dir, "RUN_MANIFEST.json")
if os.path.isfile(result_manifest_path):
    with open(result_manifest_path) as handle:
        result_manifest = json.load(handle)
    with open(os.path.join(result_dir, "subject_results.json")) as handle:
        paired_results = json.load(handle)

    role_pair_rows = _role_pair_csv_rows(paired_results["subjects"])
    normalized_rows = _normalized_csv_rows(
        paired_results["subjects"], role_pair_rows)
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
        len(normalized_keys) == len(set(normalized_keys)),
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
        left = subject_result["ieeg"]["result_3a"]
        right = subject_result["scalp"]["result_3a"]
        common_support_ok &= (
            geometry["status"] == "exact_shared_geometry"
            and shared["paired_3a_common_support"]["n_shared_seconds"] > 0
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
                and geometry[
                    "cross_correlation_retained_window_count"]
                == left_xcorr["n_intervals"]
                and bool(geometry[
                    "cross_correlation_retained_window_starts_sha256"])
            )
    check(
        "every saved primary 3A pair has exact shared support geometry",
        common_support_ok,
    )

print("ALL PAIRED-SCALP CHECKS PASSED")
