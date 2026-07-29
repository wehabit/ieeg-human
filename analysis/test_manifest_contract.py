"""Focused regressions for terminal run and completed-cache contracts."""
from __future__ import annotations

import json
import os
import tempfile

import numpy as np

from pipeline_version import (
    validate_completed_cache_failures,
    validate_full_interval_acquisition,
    validate_local_chunk_acquisition,
    validate_terminal_run_manifest,
    write_run_manifest,
)


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


def rejected(call, message):
    try:
        call()
    except (RuntimeError, ValueError) as exc:
        return message in str(exc)
    return False


def manifest(directory, **overrides):
    values = dict(
        pipeline="synthetic",
        requested=["A", "B"],
        completed=["A"],
        skipped=[{"subject": "B", "reason": "prespecified skip"}],
        failed=[],
        config={"purpose": "contract regression"},
        run_state="complete",
        result_files_sha256={"A": "a" * 64, "B": "b" * 64},
    )
    values.update(overrides)
    return write_run_manifest(directory, **values)


with tempfile.TemporaryDirectory() as directory:
    completed = manifest(
        directory,
        extra_fields={"recovery_provenance": {"method": "synthetic"}},
    )
    with open(os.path.join(directory, "RUN_MANIFEST.json")) as handle:
        on_disk = json.load(handle)
    check(
        "a complete disjoint partition is written once with extra provenance",
        completed == on_disk
        and on_disk["recovery_provenance"]["method"] == "synthetic",
    )
    duplicate_skipped = dict(completed)
    duplicate_skipped["skipped"] = completed["skipped"] * 2
    check(
        "read-side validation rejects duplicate skipped subjects in a tampered manifest",
        rejected(
            lambda: validate_terminal_run_manifest(
                duplicate_skipped, source="tampered"),
            "skipped contains duplicate subject identifiers",
        ),
    )

    check(
        "extra fields cannot replace reserved manifest identity",
        rejected(
            lambda: manifest(
                directory, extra_fields={"run_state": "forged"}),
            "reserved keys",
        ),
    )
    check(
        "terminal partitions must be disjoint",
        rejected(
            lambda: manifest(
                directory,
                completed=["A"],
                skipped=[],
                failed=[{"subject": "A", "error": "synthetic"}],
                run_state="failed",
                result_files_sha256={"A": "a" * 64},
            ),
            "disjoint",
        ),
    )
    check(
        "terminal partitions must account for every requested subject",
        rejected(
            lambda: manifest(
                directory,
                completed=["A"],
                skipped=[],
                result_files_sha256={"A": "a" * 64},
            ),
            "every requested subject",
        ),
    )
    check(
        "a complete manifest cannot hide failures",
        rejected(
            lambda: manifest(
                directory,
                completed=["A"],
                skipped=[],
                failed=[{"subject": "B", "error": "synthetic"}],
                result_files_sha256={"A": "a" * 64},
            ),
            "complete manifest cannot contain failed",
        ),
    )
    failed_manifest = manifest(
        directory,
        completed=["A"],
        skipped=[],
        failed=[{"subject": "B", "error": "synthetic"}],
        run_state="failed",
        result_files_sha256={"A": "a" * 64},
    )
    check(
        "a failed manifest records the exact completed-versus-failed partition",
        failed_manifest["run_state"] == "failed",
    )
    duplicate_failed = dict(failed_manifest)
    duplicate_failed["failed"] = failed_manifest["failed"] * 2
    check(
        "read-side validation rejects duplicate failed subjects in a tampered manifest",
        rejected(
            lambda: validate_terminal_run_manifest(
                duplicate_failed, source="tampered"),
            "failed contains duplicate subject identifiers",
        ),
    )
    check(
        "terminal hashes cover completed and skipped files exactly",
        rejected(
            lambda: manifest(
                directory,
                result_files_sha256={"A": "a" * 64}),
            "hashes must cover",
        ),
    )
    check(
        "an in-progress manifest cannot claim terminal results",
        rejected(
            lambda: manifest(
                directory,
                completed=["A"],
                skipped=[],
                failed=[],
                run_state="in_progress",
                result_files_sha256={"A": "a" * 64},
            ),
            "in-progress manifest cannot claim",
        ),
    )


class ArrayCache:
    def __init__(self, **values):
        self.values = {
            key: np.asarray(value) for key, value in values.items()
        }
        self.files = tuple(self.values)

    def __getitem__(self, key):
        return self.values[key]


clean = ArrayCache(
    failed_chunks_json=json.dumps([]),
    ecg_failures_json=json.dumps([]),
)
check(
    "completed cardiac caches require two explicit empty failure arrays",
    validate_completed_cache_failures(
        clean, source="synthetic", require_ecg=True)
    == {"failed_chunks": [], "ecg_failures": []},
)
check(
    "missing ECG failure metadata fails closed",
    rejected(
        lambda: validate_completed_cache_failures(
            ArrayCache(failed_chunks_json=json.dumps([])),
            source="synthetic",
            require_ecg=True,
        ),
        "lacks required cache field",
    ),
)
check(
    "nonempty ECG failure metadata fails closed",
    rejected(
        lambda: validate_completed_cache_failures(
            ArrayCache(
                failed_chunks_json=json.dumps([]),
                ecg_failures_json=json.dumps([{"chunk": 1}]),
            ),
            source="synthetic",
            require_ecg=True,
        ),
        "failed ECG-detector chunks",
    ),
)
check(
    "malformed failure metadata fails closed",
    rejected(
        lambda: validate_completed_cache_failures(
            ArrayCache(
                failed_chunks_json="not-json",
                ecg_failures_json=json.dumps([]),
            ),
            source="synthetic",
            require_ecg=True,
        ),
        "malformed JSON",
    ),
)


def portal_acquisition_cache(records):
    return ArrayCache(
        hours=1200.0 / 3600.0,
        night_s=10.0,
        sf=100.0,
        acquisition_sample_counts_json=json.dumps(records),
    )


def portal_subrequest(start_s, duration_s):
    sample_count = int(round(duration_s * 100.0))
    return {
        "purpose": "analysis_subrequest",
        "request_start_s": float(start_s),
        "request_duration_s": float(duration_s),
        "requested_sample_count": sample_count,
        "returned_sample_count": sample_count,
        "requested_channel_count": 7,
        "returned_channel_count": 7,
        "status": "ok",
    }


def portal_core(start_s, duration_s=600.0):
    sample_count = int(round(duration_s * 100.0))
    return {
        "purpose": "analysis_core",
        "analysis_start_s": float(start_s),
        "analysis_duration_s": float(duration_s),
        "requested_sample_count": sample_count,
        "returned_sample_count": sample_count,
        "status": "ok",
    }


valid_portal_records = [
    # Core [0, 600) requires the clipped padded pull [10, 640), split at
    # the portal's 600-second request boundary.
    portal_subrequest(10.0, 600.0),
    portal_subrequest(610.0, 30.0),
    portal_core(0.0),
    # Core [600, 1200) requires [580, 1210), likewise split 600 + 30.
    portal_subrequest(580.0, 600.0),
    portal_subrequest(1180.0, 30.0),
    portal_core(600.0),
]
valid_portal = validate_full_interval_acquisition(
    portal_acquisition_cache(valid_portal_records),
    source="synthetic portal cache",
    core_purpose="analysis_core",
    subrequest_purpose="analysis_subrequest",
    core_chunk_s=600.0,
    filter_edge_s=30.0,
)
check(
    "portal acquisition accepts exact padded pulls split across bounded subrequests",
    valid_portal["core_count"] == 2
    and valid_portal["subrequest_count"] == 4,
)

unpadded_portal_records = [
    portal_subrequest(10.0, 600.0),
    portal_core(0.0),
    portal_subrequest(610.0, 600.0),
    portal_core(600.0),
]
check(
    "portal acquisition rejects complete cores whose pulls omit filter padding",
    rejected(
        lambda: validate_full_interval_acquisition(
            portal_acquisition_cache(unpadded_portal_records),
            source="synthetic portal cache",
            core_purpose="analysis_core",
            subrequest_purpose="analysis_subrequest",
            core_chunk_s=600.0,
            filter_edge_s=30.0,
        ),
        "intended padded pull",
    ),
)


def local_acquisition_cache(records=None, *, include_ledger=True):
    values = {
        "hours": 2.0 / 3600.0,
        "sf": 100.0,
    }
    if include_ledger:
        values["acquisition_chunk_sample_counts_json"] = json.dumps(
            records)
    return ArrayCache(**values)


valid_local_records = [
    {
        "start_s": 0.0,
        "duration_s": 1.0,
        "requested_pull_samples": 200,
        "returned_pull_samples": 200,
    },
    {
        "start_s": 1.0,
        "duration_s": 1.0,
        "requested_pull_samples": 200,
        "returned_pull_samples": 200,
    },
]
check(
    "local acquisition ledger proves exact padded reads and complete core coverage",
    validate_local_chunk_acquisition(
        local_acquisition_cache(valid_local_records),
        source="synthetic local cache",
        filter_edge_s=30.0,
    )["chunk_count"] == 2,
)
check(
    "local acquisition validation rejects a missing ledger",
    rejected(
        lambda: validate_local_chunk_acquisition(
            local_acquisition_cache(include_ledger=False),
            source="synthetic local cache",
            filter_edge_s=30.0,
        ),
        "lacks required cache field",
    ),
)
short_pull_records = [dict(value) for value in valid_local_records]
short_pull_records[1]["requested_pull_samples"] = 100
short_pull_records[1]["returned_pull_samples"] = 100
check(
    "local acquisition validation rejects a short padded pull",
    rejected(
        lambda: validate_local_chunk_acquisition(
            local_acquisition_cache(short_pull_records),
            source="synthetic local cache",
            filter_edge_s=30.0,
        ),
        "exact padded pull sample count",
    ),
)
short_coverage_records = [dict(value) for value in valid_local_records]
short_coverage_records[1]["duration_s"] = 0.5
check(
    "local acquisition validation rejects short core coverage",
    rejected(
        lambda: validate_local_chunk_acquisition(
            local_acquisition_cache(short_coverage_records),
            source="synthetic local cache",
            filter_edge_s=30.0,
        ),
        "do not cover the complete cached interval",
    ),
)
gap_records = [dict(value) for value in valid_local_records]
gap_records[1].update(start_s=1.1, duration_s=0.9)
check(
    "local acquisition validation rejects a gap between cores",
    rejected(
        lambda: validate_local_chunk_acquisition(
            local_acquisition_cache(gap_records),
            source="synthetic local cache",
            filter_edge_s=30.0,
        ),
        "ordered gap-free",
    ),
)
overlap_records = [dict(value) for value in valid_local_records]
overlap_records[1].update(start_s=0.9, duration_s=1.1)
check(
    "local acquisition validation rejects overlapping cores",
    rejected(
        lambda: validate_local_chunk_acquisition(
            local_acquisition_cache(overlap_records),
            source="synthetic local cache",
            filter_edge_s=30.0,
        ),
        "ordered gap-free",
    ),
)

print("ALL MANIFEST CONTRACT CHECKS PASSED")
