"""Synthetic regressions for clean-checkout paired-evidence contracts."""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

from artifact_contracts import PAIRED_RESULT_SCHEMA
from paired_artifact_validation import (
    _load_json,
    _validate_result_metadata,
)
from pipeline_version import ANALYSIS_VERSION, CACHE_SCHEMA_VERSION


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


metadata = {
    "schema_version": PAIRED_RESULT_SCHEMA,
    "pipeline": "paired_scalp_ieeg_comparison",
    "analysis_version": ANALYSIS_VERSION,
    "cache_schema_version": CACHE_SCHEMA_VERSION,
    "source_tree_sha256": "a" * 64,
    "requested_subjects": ["S1", "S2"],
    "completed_subjects": ["S1", "S2"],
    "profile_id": "overlap11_endpoint_local",
}
saved_summary = {"example": {"n_pairs": 2}}
manifest = {
    **metadata,
    "run_state": "complete",
    "result_files_sha256": {},
}
subject_results = {
    **metadata,
    "subjects": [{"subject": "S1"}, {"subject": "S2"}],
    "group_summary": saved_summary,
}
group_summary = {
    **metadata,
    "group_summary": saved_summary,
}
check(
    "group summary intentionally needs no duplicate participant records",
    _validate_result_metadata(
        manifest, subject_results, group_summary) == ["S1", "S2"]
    and "subjects" not in group_summary,
)

bad_manifest = {
    **manifest,
    "completed_subjects": ["S1"],
}
try:
    _validate_result_metadata(
        bad_manifest, subject_results, group_summary)
except RuntimeError:
    incomplete_partition_rejected = True
else:
    incomplete_partition_rejected = False
check(
    "paired requested/completed/record partitions fail closed",
    incomplete_partition_rejected,
)

with tempfile.TemporaryDirectory() as directory:
    path = Path(os.path.join(directory, "invalid.json"))
    path.write_text('{"unused": NaN}\\n', encoding="utf-8")
    try:
        _load_json(path, "synthetic invalid JSON")
    except ValueError:
        nonstandard_json_rejected = True
    else:
        nonstandard_json_rejected = False
check(
    "paired artifact loader rejects NaN/Infinity JSON extensions",
    nonstandard_json_rejected,
)

print("ALL PAIRED-ARTIFACT CONTRACT CHECKS PASSED")
