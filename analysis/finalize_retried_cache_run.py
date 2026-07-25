"""Finalize a cache batch after explicitly retried transient source failures.

The ordinary cache runner fails the whole terminal manifest if any acquisition
chunk fails.  A successful one-subject retry must not silently erase that
history or bless unrelated stale files.  This utility therefore accepts only a
complete failed batch from the current producer/runtime, verifies every
requested NPZ's embedded identity/status/schema/source digest, re-hashes every
file, and records the prior failure plus exact retry subjects in the replacement
terminal manifest.
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    atomic_json_dump,
    cache_code_sha256,
    file_sha256,
    npz_scalar_text,
    runtime_versions,
    write_run_manifest,
)


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def finalize(cache_dir, retried_subjects):
    cache_dir = os.path.abspath(cache_dir)
    manifest_path = os.path.join(cache_dir, "RUN_MANIFEST.json")
    prior_manifest_sha256 = file_sha256(manifest_path)
    with open(manifest_path) as handle:
        prior = json.load(handle)

    retried_subjects = list(retried_subjects)
    if len(retried_subjects) != len(set(retried_subjects)):
        raise RuntimeError("retried subjects contain duplicates")
    failed = list(prior.get("failed", []))
    failed_subjects = [
        value.get("subject") if isinstance(value, dict) else None
        for value in failed
    ]
    current_digest = cache_code_sha256(ROOT)
    if (
        prior.get("run_state") != "complete"
        or prior.get("pipeline") != "cache_lc_series"
        or prior.get("analysis_version") != ANALYSIS_VERSION
        or prior.get("cache_schema_version") != CACHE_SCHEMA_VERSION
        or prior.get("runtime_versions") != runtime_versions()
        or prior.get("config", {}).get("cache_code_sha256") != current_digest
        or set(failed_subjects) != set(retried_subjects)
        or any(subject is None for subject in failed_subjects)
    ):
        raise RuntimeError(
            "prior manifest is not the exact current failed cache batch "
            "described by --retried-subject")

    requested = list(prior.get("requested", []))
    skipped = list(prior.get("skipped", []))
    if any(
        not isinstance(value, dict)
        or not value.get("subject")
        or not value.get("reason")
        for value in skipped
    ):
        raise RuntimeError("prior manifest has malformed skip records")
    skipped_subjects = {value["subject"] for value in skipped}
    if (
        len(requested) != len(set(requested))
        or not set(retried_subjects) <= set(requested)
        or skipped_subjects & set(retried_subjects)
    ):
        raise RuntimeError("prior manifest subject partition is invalid")

    hashes = {}
    completed = []
    for subject in requested:
        path = os.path.join(cache_dir, f"{subject}.npz")
        if not os.path.isfile(path):
            raise RuntimeError(f"requested cache is missing: {path}")
        with np.load(path, allow_pickle=False) as cache:
            expected_status = "skip" if subject in skipped_subjects else "ok"
            if (
                npz_scalar_text(cache, "subject") != subject
                or npz_scalar_text(cache, "status") != expected_status
                or npz_scalar_text(
                    cache, "cache_schema_version") != CACHE_SCHEMA_VERSION
                or npz_scalar_text(
                    cache, "cache_code_sha256") != current_digest
            ):
                raise RuntimeError(
                    f"{path} embedded identity/status/schema/source digest is invalid")
        hashes[subject] = file_sha256(path)
        if subject not in skipped_subjects:
            completed.append(subject)

    manifest = write_run_manifest(
        cache_dir,
        pipeline="cache_lc_series",
        requested=requested,
        completed=completed,
        skipped=skipped,
        failed=[],
        config=prior["config"],
        run_state="complete",
        result_files_sha256=hashes,
    )
    manifest["recovery_provenance"] = {
        "prior_batch_run_id": prior["run_id"],
        "prior_batch_manifest_sha256": prior_manifest_sha256,
        "prior_failed_records": failed,
        "retried_subjects": retried_subjects,
        "retry_method": (
            "same current cache producer, seven-hour --force-equivalent "
            "single-subject streams run sequentially after transient chunk failures"),
        "retried_result_files_sha256": {
            subject: hashes[subject] for subject in retried_subjects
        },
    }
    atomic_json_dump(manifest, manifest_path)
    return manifest


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
