"""Regression checks for staging-calibration cache lineage."""
from __future__ import annotations

import os
import tempfile

import numpy as np

from calibrate_staging_windows import (
    _manifest_subjects,
    calibration_source_files_sha256,
    calibration_records,
    update_staging_calibration_pin,
)
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    atomic_json_dump,
    atomic_savez,
    cache_code_sha256,
    file_sha256,
    write_run_manifest,
)


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


check(
    "calibration records exact purpose-specific source hashes",
    set(calibration_source_files_sha256()) == {
        "analysis/calibrate_staging_windows.py",
        "analysis/pipeline_version.py",
    }
    and all(
        len(value) == 64
        for value in calibration_source_files_sha256().values()
    ),
)


def write_cache(path, subject, *, embedded_subject=None, extra=False):
    shape = (1, 2, 14)
    payload = dict(
        subject=embedded_subject or subject,
        status="ok",
        cache_schema_version=CACHE_SCHEMA_VERSION,
        cache_code_sha256=cache_code_sha256(ROOT),
        cortical_chans=np.asarray(["A1"]),
        ep_valid_welch_window_mask_by_contact=np.ones(shape, np.uint8),
        ep_window_swa_power_by_contact=np.ones(shape, float),
        ep_window_total_power_by_contact=np.full(shape, 2.0),
    )
    if extra:
        payload["tampered"] = np.asarray(1)
    atomic_savez(path, **payload)


def write_manifest(cache_dir, subject):
    write_run_manifest(
        cache_dir,
        pipeline="stage_ds003848",
        requested=[subject],
        completed=[subject],
        skipped=[],
        failed=[],
        config={
            "cache_schema_version": CACHE_SCHEMA_VERSION,
            "cache_code_sha256": cache_code_sha256(ROOT),
        },
        run_id="synthetic-calibration-cache",
        run_state="complete",
        result_files_sha256={
            subject: file_sha256(
                os.path.join(cache_dir, f"{subject}.npz")),
        },
    )


with tempfile.TemporaryDirectory() as cache_dir:
    subject = "sub-TEST0001"
    path = os.path.join(cache_dir, f"{subject}.npz")
    write_cache(path, subject)
    write_manifest(cache_dir, subject)
    subjects, manifest = _manifest_subjects(cache_dir)
    records, hashes = calibration_records(
        cache_dir, subjects, manifest["result_files_sha256"])
    check(
        "calibration accepts an exact manifest-pinned neutral cache",
        subjects == [subject]
        and len(records) > 0
        and hashes == manifest["result_files_sha256"],
    )

    write_cache(path, subject, extra=True)
    try:
        calibration_records(
            cache_dir, subjects, manifest["result_files_sha256"])
    except RuntimeError as exc:
        tamper_rejected = "differ" in str(exc)
    else:
        tamper_rejected = False
    check(
        "calibration rejects cache bytes changed after the terminal manifest",
        tamper_rejected,
    )

    write_cache(path, subject, embedded_subject="sub-WRONG")
    write_manifest(cache_dir, subject)
    try:
        _manifest_subjects(cache_dir)
    except RuntimeError as exc:
        identity_rejected = "identity" in str(exc)
    else:
        identity_rejected = False
    check(
        "calibration rejects a manifest-hashed cache with the wrong embedded subject",
        identity_rejected,
    )


with tempfile.TemporaryDirectory() as synthetic_root:
    artifact_path = os.path.join(
        synthetic_root, "outputs", "calibration.json")
    profile_path = os.path.join(
        synthetic_root, "analysis", "profiles.json")
    artifact = {
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "cache_manifest_sha256": "a" * 64,
        "cache_pipeline": "stage_ds003848",
        "cache_run_id": "synthetic-current-cache",
        "recommended_minimum_valid_windows": 11,
        "support_results": [{
            "minimum_valid_windows": 11,
            "exact_support": {"meets_calibration_targets": True},
        }],
    }
    profile_set = {
        "schema_version": "qc-profile-set-v1",
        "method_config": {
            "staging_calibration": {
                "analysis_version": "old",
                "artifact_relative_path": "old.json",
                "artifact_sha256": "0" * 64,
                "cache_manifest_sha256": "1" * 64,
                "cache_pipeline": "stage_ds003848",
                "cache_run_id": "old",
                "cache_schema_version": "old",
                "calibration_scope": "synthetic outcome-blind calibration",
                "hup_transport_status": "not calibrated for HUP",
                "not_a_paper_requirement": True,
                "recommended_minimum_valid_windows": 11,
            },
        },
    }
    atomic_json_dump(artifact, artifact_path)
    atomic_json_dump(profile_set, profile_path)
    expected_artifact_sha = file_sha256(artifact_path)
    pin = update_staging_calibration_pin(
        artifact_path,
        profile_path,
        root=synthetic_root,
    )
    first_profile_sha = file_sha256(profile_path)
    second_pin = update_staging_calibration_pin(
        artifact_path,
        profile_path,
        root=synthetic_root,
    )
    check(
        "calibration pin update is exact, contextual, and deterministic",
        pin == second_pin
        and file_sha256(profile_path) == first_profile_sha
        and pin["artifact_relative_path"] == "outputs/calibration.json"
        and pin["artifact_sha256"] == expected_artifact_sha
        and pin["cache_run_id"] == "synthetic-current-cache"
        and pin["calibration_scope"]
        == "synthetic outcome-blind calibration",
    )

    changed = dict(artifact)
    changed["recommended_minimum_valid_windows"] = 12
    changed["support_results"] = [{
        "minimum_valid_windows": 12,
        "exact_support": {"meets_calibration_targets": True},
    }]
    atomic_json_dump(changed, artifact_path)
    try:
        update_staging_calibration_pin(
            artifact_path,
            profile_path,
            root=synthetic_root,
        )
    except RuntimeError as exc:
        changed_recommendation_rejected = "recommendation changed" in str(exc)
    else:
        changed_recommendation_rejected = False
    check(
        "pinning cannot silently change the locked calibration recommendation",
        changed_recommendation_rejected,
    )
