"""Synthetic regressions for compact QC-grid publication integrity."""
from __future__ import annotations

import json
import os
import tempfile

from compact_qc_grid_artifacts import (
    LOCKED_RELATIVE_PATH,
    MANIFEST_NAME,
    build_public_artifacts,
    validate_public_artifacts,
)
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    atomic_json_dump,
    file_sha256,
    runtime_versions,
    source_tree_sha256,
)
from qc_profiles import expand_qc_grid, manifest_qc_config


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


def unavailable_subject():
    return {
        "subject": "SYNTHETIC",
        "support_diagnostics": {},
        "result_3a": {
            "endpoint_availability": {
                "spectrum": False,
                "fixed_0p02_coherence": False,
                "cross_correlation": False,
            },
        },
        "result_3b": {
            "stages": {
                stage: {"available_under_profile": False}
                for stage in ("N2", "N3", "NREM")
            },
        },
        "result_3d": {
            "descriptive_effect_available": False,
            "support_passes_profile": None,
        },
    }


def empty_summary(minimum_participants=5):
    return {
        "n_subjects": 1,
        "minimum_participants": minimum_participants,
        "n_3a_spectrum": 0,
        "n_3a_coherence": 0,
        "n_3a_cross_correlation": 0,
        "n_3b_N2": 0,
        "n_3b_N3": 0,
        "n_3b_pooled_NREM_exploratory": 0,
        "n_3d_descriptive_effect": 0,
        "n_3d_support_evaluable": 0,
        "meets_minimum_participants": {
            "3a_spectrum": False,
            "3a_coherence": False,
            "3a_cross_correlation": False,
            "3b_N2": False,
            "3b_N3": False,
            "3b_pooled_NREM_exploratory": False,
            "3d_descriptive_effect": False,
            "3d_support_evaluable": False,
        },
    }


def write_full_grids(input_root):
    for grid_id in (
        "coverage_oat_v1",
        "staging_window_support_v1",
        "auxiliary_window_support_v1",
        "event_count_oat_v1",
    ):
        for cohort in ("hup", "respect"):
            profiles = []
            for profile in expand_qc_grid(grid_id):
                profiles.append({
                    **manifest_qc_config(profile),
                    "subjects": [unavailable_subject()],
                    "summary": empty_summary(
                        profile["cohort"]["minimum_participants"]),
                })
            payload = {
                "analysis_version": ANALYSIS_VERSION,
                "cache_schema_version": CACHE_SCHEMA_VERSION,
                "source_tree_sha256": source_tree_sha256(ROOT),
                "grid_id": grid_id,
                "grid_role": "synthetic artifact-contract test",
                "cache_pipeline": (
                    "cache_lc_series" if cohort == "hup"
                    else "stage_ds003848"),
                "cache_directory_relative": f"data/derived/{cohort}",
                "cache_manifest_sha256": (
                    "a" * 64 if cohort == "hup" else "b" * 64),
                "cache_manifest_run_id": f"synthetic-{cohort}",
                "cache_files_sha256": {"SYNTHETIC": "c" * 64},
                "cache_runtime_versions": runtime_versions(),
                "runtime_versions": runtime_versions(),
                "generated_at_utc": "2026-07-01T00:00:00+00:00",
                "requested_subjects": ["SYNTHETIC"],
                "analysed_subjects": ["SYNTHETIC"],
                "explicitly_skipped_subjects": [],
                "staging_calibration_provenance": {"synthetic": True},
                "staging_calibration_application_to_this_cache": (
                    "synthetic"),
                "limitations": [],
                "profiles": profiles,
            }
            atomic_json_dump(
                payload,
                os.path.join(
                    input_root, grid_id, f"{cohort}_qc_grid.json"),
            )


def rehash_public_file(output_root, relative):
    manifest_path = os.path.join(output_root, MANIFEST_NAME)
    with open(manifest_path, encoding="utf-8") as handle:
        manifest = json.load(handle)
    manifest["result_files_sha256"][relative] = file_sha256(
        os.path.join(output_root, *relative.split("/")))
    atomic_json_dump(manifest, manifest_path)


with tempfile.TemporaryDirectory() as directory:
    input_root = os.path.join(directory, "physical-input")
    output_root = os.path.join(directory, "staged-public")
    write_full_grids(input_root)
    built = build_public_artifacts(
        input_root=input_root,
        output_root=output_root,
    )
    check(
        "custom output root owns its locked snapshot and terminal manifest",
        len(built) == 9
        and os.path.isfile(os.path.join(output_root, LOCKED_RELATIVE_PATH))
        and os.path.isfile(os.path.join(output_root, MANIFEST_NAME))
        and len(validate_public_artifacts(output_root)) == 9,
    )
    summary_relative = "coverage_oat_v1/hup_qc_grid_summary.json"
    summary_path = os.path.join(output_root, *summary_relative.split("/"))
    with open(summary_path, encoding="utf-8") as handle:
        summary = json.load(handle)
    check(
        "public provenance records a stable logical full-grid path",
        summary["full_grid_path_relative"]
        == "outputs/qc_grid/coverage_oat_v1/hup_qc_grid.json",
    )

    summary["profiles"][0]["qc_profile"]["power"][
        "minimum_contacts"] = 999
    atomic_json_dump(summary, summary_path)
    rehash_public_file(output_root, summary_relative)
    try:
        validate_public_artifacts(output_root)
    except RuntimeError:
        profile_tamper_rejected = True
    else:
        profile_tamper_rejected = False
    check(
        "manifest-rehashed QC-profile tampering still fails semantic validation",
        profile_tamper_rejected,
    )

    build_public_artifacts(input_root=input_root, output_root=output_root)
    with open(summary_path, encoding="utf-8") as handle:
        summary = json.load(handle)
    summary["profiles"] = summary["profiles"][:-1]
    atomic_json_dump(summary, summary_path)
    rehash_public_file(output_root, summary_relative)
    try:
        validate_public_artifacts(output_root)
    except RuntimeError:
        dropped_profile_rejected = True
    else:
        dropped_profile_rejected = False
    check(
        "manifest-rehashed profile deletion fails the exact grid contract",
        dropped_profile_rejected,
    )

    build_public_artifacts(input_root=input_root, output_root=output_root)
    manifest_path = os.path.join(output_root, MANIFEST_NAME)
    with open(manifest_path, encoding="utf-8") as handle:
        manifest = json.load(handle)
    manifest["run_state"] = "in_progress"
    atomic_json_dump(manifest, manifest_path)
    try:
        validate_public_artifacts(output_root)
    except RuntimeError:
        interrupted_rejected = True
    else:
        interrupted_rejected = False
    check(
        "an interrupted public-artifact set cannot validate as terminal",
        interrupted_rejected,
    )

    write_full_grids(input_root)
    injected_source = os.path.join(
        input_root, "coverage_oat_v1", "hup_qc_grid.json")
    with open(injected_source, encoding="utf-8") as handle:
        injected = json.load(handle)
    injected["unexpected_top_level_raw_vector"] = [1, 2, 3]
    atomic_json_dump(injected, injected_source)
    try:
        build_public_artifacts(
            input_root=input_root,
            output_root=os.path.join(directory, "top-level-injection"),
        )
    except RuntimeError:
        top_level_injection_rejected = True
    else:
        top_level_injection_rejected = False
    check(
        "unexpected full-grid fields cannot leak into compact public artifacts",
        top_level_injection_rejected,
    )

    write_full_grids(input_root)
    with open(injected_source, encoding="utf-8") as handle:
        injected = json.load(handle)
    injected["profiles"][0]["unexpected_profile_numeric_vector"] = [4, 5, 6]
    atomic_json_dump(injected, injected_source)
    try:
        build_public_artifacts(
            input_root=input_root,
            output_root=os.path.join(directory, "profile-injection"),
        )
    except RuntimeError:
        profile_injection_rejected = True
    else:
        profile_injection_rejected = False
    check(
        "unexpected profile fields cannot leak into compact public artifacts",
        profile_injection_rejected,
    )

    write_full_grids(input_root)
    with open(injected_source, encoding="utf-8") as handle:
        injected = json.load(handle)
    injected["profiles"][0]["subjects"][0][
        "unexpected_subject_numeric_vector"
    ] = [7, 8, 9]
    atomic_json_dump(injected, injected_source)
    try:
        build_public_artifacts(
            input_root=input_root,
            output_root=os.path.join(directory, "subject-injection"),
        )
    except RuntimeError:
        subject_injection_rejected = True
    else:
        subject_injection_rejected = False
    check(
        "unexpected subject fields cannot leak into the locked snapshot",
        subject_injection_rejected,
    )

    write_full_grids(input_root)
    build_public_artifacts(input_root=input_root, output_root=output_root)
    with open(summary_path, encoding="utf-8") as handle:
        summary = json.load(handle)
    summary["profiles"][0]["summary"]["minimum_participants"] = 999
    summary["profiles"][0]["summary"]["meets_minimum_participants"] = {
        "fabricated": True,
    }
    atomic_json_dump(summary, summary_path)
    rehash_public_file(output_root, summary_relative)
    try:
        validate_public_artifacts(output_root)
    except RuntimeError:
        fabricated_reporting_gate_rejected = True
    else:
        fabricated_reporting_gate_rejected = False
    check(
        "manifest-rehashed reporting-gate fabrication fails semantic validation",
        fabricated_reporting_gate_rejected,
    )

print("ALL PUBLIC-ARTIFACT CHECKS PASSED")
