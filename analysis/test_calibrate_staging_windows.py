"""Regression checks for staging-calibration cache lineage."""
from __future__ import annotations

import os
import tempfile

import numpy as np

from calibrate_staging_windows import (
    _manifest_subjects,
    calibration_source_files_sha256,
    calibration_records,
)
from pipeline_version import (
    CACHE_SCHEMA_VERSION,
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
