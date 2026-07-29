"""Finalize explicitly retried subjects from one failed HUP cache batch.

A failed terminal manifest records the exact pre-run bytes for each failed
subject.  Retry acquisition is performed with the ordinary producer without
overwriting that manifest (for example by calling ``cache_lc_series.run``).
This finalizer then:

* preserves the exact failed manifest in ``recovery_manifests/``;
* proves every retried file is new or differs from its recorded pre-run bytes;
* requires a later generation time, exact acquisition counts, and empty fatal
  acquisition/ECG arrays;
* proves every non-retried result is byte-identical to the failed manifest; and
* writes one complete replacement manifest atomically, with recovery provenance
  included in that single write.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import tempfile

import numpy as np

from cache_lc_series import FILTER_EDGE_S
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    cache_code_sha256,
    file_sha256,
    npz_scalar_text,
    runtime_versions,
    validate_completed_cache_failures,
    validate_full_interval_acquisition,
    validate_terminal_run_manifest,
    write_run_manifest,
)
from staging_helpers import CHUNK_S


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECOVERY_DIRECTORY = "recovery_manifests"


def _parse_utc(value, *, field):
    try:
        parsed = dt.datetime.fromisoformat(str(value))
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"{field} is not a valid ISO timestamp") from exc
    if parsed.tzinfo is None:
        raise RuntimeError(f"{field} must include a timezone")
    return parsed.astimezone(dt.timezone.utc)


def _atomic_write_bytes(path, content):
    """Create or verify an immutable byte-for-byte manifest archive."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path):
        with open(path, "rb") as handle:
            if handle.read() != content:
                raise RuntimeError(
                    f"recovery archive collision at {path}")
        return
    fd, temporary = tempfile.mkstemp(
        prefix=".tmp-", suffix=".json", dir=os.path.dirname(path))
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def _archive_failed_manifest(cache_dir, manifest_path, prior_sha256):
    with open(manifest_path, "rb") as handle:
        original_bytes = handle.read()
    archive_name = f"{prior_sha256}.failed.RUN_MANIFEST.json"
    archive_path = os.path.join(
        cache_dir, RECOVERY_DIRECTORY, archive_name)
    _atomic_write_bytes(archive_path, original_bytes)
    if file_sha256(archive_path) != prior_sha256:
        raise RuntimeError("failed-run manifest archive is not byte-identical")
    return {
        "relative_path": os.path.relpath(
            archive_path, cache_dir).replace(os.sep, "/"),
        "sha256": prior_sha256,
    }


def _subject_partition(prior, retried_subjects):
    requested = list(prior.get("requested", []))
    completed = list(prior.get("completed", []))
    skipped = list(prior.get("skipped", []))
    failed = list(prior.get("failed", []))
    skipped_subjects = [
        value.get("subject") if isinstance(value, dict) else None
        for value in skipped
    ]
    failed_subjects = [
        value.get("subject") if isinstance(value, dict) else None
        for value in failed
    ]
    if (
        len(requested) != len(set(requested))
        or len(completed) != len(set(completed))
        or len(skipped_subjects) != len(set(skipped_subjects))
        or len(failed_subjects) != len(set(failed_subjects))
        or any(subject is None for subject in skipped_subjects + failed_subjects)
        or set(completed) & set(skipped_subjects)
        or set(completed) & set(failed_subjects)
        or set(skipped_subjects) & set(failed_subjects)
        or set(completed) | set(skipped_subjects) | set(failed_subjects)
        != set(requested)
        or set(failed_subjects) != set(retried_subjects)
    ):
        raise RuntimeError(
            "prior manifest has an invalid or unexpected subject partition")
    return requested, completed, skipped, failed


def _validate_identity(cache, *, path, subject, expected_status, digest):
    if (
        npz_scalar_text(cache, "subject") != subject
        or npz_scalar_text(cache, "status") != expected_status
        or npz_scalar_text(
            cache, "cache_schema_version") != CACHE_SCHEMA_VERSION
        or npz_scalar_text(cache, "cache_code_sha256") != digest
    ):
        raise RuntimeError(
            f"{path} embedded identity/status/schema/source digest is invalid")


def _failed_retry_baselines(failed):
    baselines = {}
    for record in failed:
        baseline = record.get("pre_run_output")
        if not isinstance(baseline, dict):
            raise RuntimeError(
                "failed record lacks pre-run output identity; a retry cannot "
                "be proven from this batch")
        existed = baseline.get("existed")
        sha256 = baseline.get("sha256")
        if (
            not isinstance(existed, bool)
            or (existed and (
                not isinstance(sha256, str) or len(sha256) != 64))
            or (not existed and sha256 is not None)
        ):
            raise RuntimeError("failed record has malformed pre-run output identity")
        baselines[record["subject"]] = {
            "existed": existed,
            "sha256": sha256,
        }
    return baselines


def finalize(cache_dir, retried_subjects):
    cache_dir = os.path.abspath(cache_dir)
    manifest_path = os.path.join(cache_dir, "RUN_MANIFEST.json")
    prior_manifest_sha256 = file_sha256(manifest_path)
    with open(manifest_path) as handle:
        prior = json.load(handle)

    retried_subjects = list(retried_subjects)
    if (
        not retried_subjects
        or len(retried_subjects) != len(set(retried_subjects))
    ):
        raise RuntimeError(
            "retried subjects must be a nonempty list without duplicates")
    current_digest = cache_code_sha256(ROOT)
    partition = validate_terminal_run_manifest(
        prior, source=manifest_path)
    if (
        partition["run_state"] != "failed"
        or prior.get("pipeline") != "cache_lc_series"
        or prior.get("analysis_version") != ANALYSIS_VERSION
        or prior.get("cache_schema_version") != CACHE_SCHEMA_VERSION
        or prior.get("runtime_versions") != runtime_versions()
        or prior.get("config", {}).get("cache_code_sha256") != current_digest
    ):
        raise RuntimeError(
            "prior manifest is not an exact current failed cache batch")

    requested, prior_completed, skipped, failed = _subject_partition(
        prior, retried_subjects)
    baselines = _failed_retry_baselines(failed)
    prior_hashes = prior.get("result_files_sha256")
    non_retried = prior_completed + [value["subject"] for value in skipped]
    if (
        not isinstance(prior_hashes, dict)
        or set(prior_hashes) != set(non_retried)
    ):
        raise RuntimeError(
            "prior failed manifest does not pin every non-retried result")

    prior_generated_at = _parse_utc(
        prior.get("generated_at_utc"), field="prior generated_at_utc")
    hashes = {}
    retry_proof = {}
    skipped_subjects = {value["subject"] for value in skipped}
    for subject in requested:
        path = os.path.join(cache_dir, f"{subject}.npz")
        if not os.path.isfile(path):
            raise RuntimeError(f"requested cache is missing: {path}")
        current_sha256 = file_sha256(path)
        hashes[subject] = current_sha256
        if subject not in retried_subjects:
            if current_sha256 != prior_hashes[subject]:
                raise RuntimeError(
                    f"non-retried cache changed after the failed batch: {subject}")
            with np.load(path, allow_pickle=False) as cache:
                expected_status = (
                    "skip" if subject in skipped_subjects else "ok")
                _validate_identity(
                    cache,
                    path=path,
                    subject=subject,
                    expected_status=expected_status,
                    digest=current_digest,
                )
                if expected_status == "ok":
                    validate_completed_cache_failures(
                        cache, source=path, require_ecg=True)
            continue

        baseline = baselines[subject]
        if baseline["existed"] and current_sha256 == baseline["sha256"]:
            raise RuntimeError(
                f"retried cache bytes did not change for {subject}")
        with np.load(path, allow_pickle=False) as cache:
            _validate_identity(
                cache,
                path=path,
                subject=subject,
                expected_status="ok",
                digest=current_digest,
            )
            validate_completed_cache_failures(
                cache, source=path, require_ecg=True)
            generated_at_text = npz_scalar_text(cache, "generated_at_utc")
            generated_at = _parse_utc(
                generated_at_text,
                field=f"{subject} generated_at_utc",
            )
            if generated_at <= prior_generated_at:
                raise RuntimeError(
                    f"{subject} was not generated after the failed batch")
            acquisition = validate_full_interval_acquisition(
                cache,
                source=path,
                core_purpose="analysis_core",
                subrequest_purpose="analysis_subrequest",
                core_chunk_s=CHUNK_S,
                filter_edge_s=FILTER_EDGE_S,
            )
        retry_proof[subject] = {
            "pre_run_output_existed": baseline["existed"],
            "pre_run_output_sha256": baseline["sha256"],
            "retried_output_sha256": current_sha256,
            "retried_output_generated_at_utc": generated_at_text,
            "validated_acquisition_count_records": len(
                acquisition["records"]),
            "validated_analysis_core_count": acquisition["core_count"],
            "validated_subrequest_count": acquisition["subrequest_count"],
            "validated_analysis_duration_s": acquisition["duration_s"],
        }

    archive = _archive_failed_manifest(
        cache_dir, manifest_path, prior_manifest_sha256)
    completed = prior_completed + retried_subjects
    recovery_provenance = {
        "prior_batch_run_id": prior["run_id"],
        "prior_batch_manifest_sha256": prior_manifest_sha256,
        "prior_batch_manifest_archive": archive,
        "prior_failed_records": failed,
        "retried_subjects": retried_subjects,
        "retry_proof_by_subject": retry_proof,
        "non_retried_result_files_sha256": {
            subject: prior_hashes[subject] for subject in non_retried
        },
        "retry_method": (
            "current cache producer rebuilt only the failed subjects; this "
            "finalizer proved changed/new bytes, later generation timestamps, "
            "exact portal sample counts, and empty fatal failure arrays"),
    }
    return write_run_manifest(
        cache_dir,
        pipeline="cache_lc_series",
        requested=requested,
        completed=completed,
        skipped=skipped,
        failed=[],
        config=prior["config"],
        run_state="complete",
        result_files_sha256=hashes,
        extra_fields={"recovery_provenance": recovery_provenance},
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument(
        "--retried-subject", action="append", required=True)
    args = parser.parse_args()
    manifest = finalize(args.cache_dir, args.retried_subject)
    print(
        f"finalized {args.cache_dir}: {len(manifest['completed'])} completed, "
        f"{len(manifest['skipped'])} skipped, 0 failed",
        flush=True,
    )


if __name__ == "__main__":
    main()
