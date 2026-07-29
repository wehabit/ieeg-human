"""Synthetic regressions for clean-checkout paired-evidence contracts."""
from __future__ import annotations

import copy
import hashlib
import json
import os
import tempfile
from pathlib import Path

from artifact_contracts import PAIRED_RESULT_SCHEMA
from paired_artifact_validation import (
    _validate_frozen_3a_evidence,
    _load_json,
    _validate_classifications,
    _validate_result_metadata,
    derive_inventory_classification,
    json_safe_inventory_activity,
)
from paired_reporting import _group_summary
from pipeline_version import ANALYSIS_VERSION, CACHE_SCHEMA_VERSION


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


def _identity(revision_id, data_check, *, duration_us=1_000_000.0):
    return {
        "revision_id": revision_id,
        "data_check": data_check,
        "start_time_us": 0,
        "end_time_us": 1_000_000,
        "duration_us": duration_us,
        "number_of_samples": 100,
        "sample_rate_hz": 100.0,
    }


def _record(*, nonflat=True, finite_count=100, dynamic_range=1.0):
    labels = ["A1", "C3"]
    reference = _identity("reference-revision", "reference-data-check")
    c3_identity = _identity("c3-revision", "c3-data-check")
    roles = {
        role: {
            "target_normalized_label": target,
            "matches": [],
        }
        for role, target in {
            "c3": "C3", "c4": "C4", "f3": "F3", "f4": "F4",
            "fz": "FZ", "a1": "A1", "a2": "A2", "m1": "M1",
            "m2": "M2",
        }.items()
    }
    roles["a1"]["matches"] = [{
        "label": "A1",
        "identity": reference,
        "same_geometry_as_pinned_reference": True,
    }]
    roles["c3"]["matches"] = [{
        "label": "C3",
        "identity": c3_identity,
        "same_geometry_as_pinned_reference": True,
    }]
    return {
        "subject": "S1",
        "query_status": "ok",
        "dataset_name": "S1",
        "snapshot_id": "snapshot-1",
        "ordered_channel_labels": labels,
        "ordered_channel_labels_sha256": hashlib.sha256(
            json.dumps(
                labels, ensure_ascii=False, separators=(",", ":")
            ).encode("utf-8")
        ).hexdigest(),
        "pinned_reference_channel": "A1",
        "pinned_reference_identity": reference,
        "roles": roles,
        "frozen_3a": {
            "record_available": True,
            "spectrum_available": True,
            "coherence_available": False,
            "cross_correlation_available": True,
        },
        "sidecar": {
            "path_relative": "data/derived/paired_scalp/S1.npz",
            "sha256": "b" * 64,
            "roles": {"c3": "C3"},
            "source_dataset": "S1",
            "source_snapshot_id": "snapshot-1",
            "source_identity_by_channel": {"C3": copy.deepcopy(c3_identity)},
            "channels": {
                "C3": {
                    "numerically_nonflat_full_interval": nonflat,
                    "finite_sample_count": finite_count,
                    "raw_dynamic_range": dynamic_range,
                },
            },
        },
    }


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


def _analysis_record(subject):
    unavailable_3a = {
        "peak": {"accepted": False, "peak_hz": None},
        "coherence": None,
        "cross_correlation": None,
        "endpoint_availability": {"spectrum": False},
        "negative_control": {},
    }
    return {
        "subject": subject,
        "ieeg": {
            "result_3a": copy.deepcopy(unavailable_3a),
            "result_3b": {"stages": {}},
        },
        "scalp": {
            "result_3a": copy.deepcopy(unavailable_3a),
            "result_3b": {},
        },
    }


analysis_records = [_analysis_record("S1"), _analysis_record("S2")]
saved_summary = _group_summary(analysis_records)
manifest = {
    **metadata,
    "run_state": "complete",
    "result_files_sha256": {},
}
subject_results = {
    **metadata,
    "subjects": analysis_records,
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

forged_subject_results = copy.deepcopy(subject_results)
forged_group_summary = copy.deepcopy(group_summary)
forged_subject_results[
    "group_summary"]["3A_C3"]["endpoint_counts"]["requested"] = 999
forged_group_summary[
    "group_summary"]["3A_C3"]["endpoint_counts"]["requested"] = 999
try:
    _validate_result_metadata(
        manifest, forged_subject_results, forged_group_summary)
except RuntimeError:
    jointly_forged_summary_rejected = True
else:
    jointly_forged_summary_rejected = False
check(
    "group summary is recomputed rather than trusting two matching copies",
    jointly_forged_summary_rejected,
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

for unresolved_classification in (
    "requires_full_interval_activity_sidecar",
    "c3_not_streamed_in_sidecar",
):
    pending_inventory = {
        "requested_subjects": ["S1"],
        "n_requested": 1,
        "subjects": [{
            "subject": "S1",
            "selection_classification": unresolved_classification,
        }],
        "query_error_subjects": [],
        "classification_counts": {unresolved_classification: 1},
        "classification_subjects": {unresolved_classification: ["S1"]},
    }
    try:
        _validate_classifications(pending_inventory)
    except RuntimeError:
        pending_sidecar_rejected = True
    else:
        pending_sidecar_rejected = False
    check(
        f"unresolved activity class {unresolved_classification} fails closed",
        pending_sidecar_rejected,
    )

forged_record = _record(nonflat=False, dynamic_range=0.0)
forged_record["selection_classification"] = "paired_3a_eligible"
mismatched_inventory = {
    "requested_subjects": ["S1"],
    "n_requested": 1,
    "subjects": [forged_record],
    "query_error_subjects": [],
    "classification_counts": {"paired_3a_eligible": 1},
    "classification_subjects": {"paired_3a_eligible": ["S1"]},
}
try:
    _validate_classifications(mismatched_inventory)
except RuntimeError:
    record_index_mismatch_rejected = True
else:
    record_index_mismatch_rejected = False
check(
    "classification is recomputed when the record and index are jointly forged",
    record_index_mismatch_rejected,
)

valid_record = _record()
check(
    "valid record evidence derives paired eligibility",
    derive_inventory_classification(valid_record) == "paired_3a_eligible",
)

for field in ("revision_id", "data_check"):
    tampered = copy.deepcopy(valid_record)
    tampered["sidecar"]["source_identity_by_channel"]["C3"][field] += "-changed"
    try:
        derive_inventory_classification(tampered)
    except RuntimeError:
        identity_tamper_rejected = True
    else:
        identity_tamper_rejected = False
    check(
        f"fresh portal evidence binds sidecar {field}",
        identity_tamper_rejected,
    )

tolerated_geometry = copy.deepcopy(valid_record)
tolerated_identity = tolerated_geometry[
    "sidecar"]["source_identity_by_channel"]["C3"]
tolerated_identity["duration_us"] += 0.5
tolerated_identity["sample_rate_hz"] += 0.5e-9
check(
    "sidecar identity uses the producer's geometry tolerances",
    derive_inventory_classification(tolerated_geometry)
    == "paired_3a_eligible",
)

changed_geometry = copy.deepcopy(valid_record)
changed_geometry[
    "sidecar"]["source_identity_by_channel"]["C3"]["duration_us"] += 1.01
try:
    derive_inventory_classification(changed_geometry)
except RuntimeError:
    geometry_tamper_rejected = True
else:
    geometry_tamper_rejected = False
check(
    "sidecar geometry outside the producer tolerance fails closed",
    geometry_tamper_rejected,
)

forged_roles = copy.deepcopy(valid_record)
forged_roles["roles"]["c3"]["matches"] = []
try:
    derive_inventory_classification(forged_roles)
except RuntimeError:
    role_tamper_rejected = True
else:
    role_tamper_rejected = False
check(
    "role matches are re-derived from the ordered portal labels",
    role_tamper_rejected,
)

empty_activity = json_safe_inventory_activity({
    "C3": {
        "numerically_nonflat_full_interval": False,
        "finite_sample_count": 0,
        "raw_dynamic_range": float("-inf"),
    },
})
empty_record = _record(
    nonflat=False, finite_count=0, dynamic_range=None)
check(
    "empty activity stays valid, flat, and standards-compliant JSON",
    empty_activity["C3"]["raw_dynamic_range"] is None
    and derive_inventory_classification(empty_record)
    == "c3_numerically_flat"
    and json.loads(json.dumps(empty_activity, allow_nan=False))
    == empty_activity,
)

locked = {
    "profiles": [{
        "subjects": [{
            "subject": "S1",
            "result_3a": {
                "endpoint_availability": {
                    "spectrum": True,
                    "fixed_0p02_coherence": False,
                    "cross_correlation": True,
                },
            },
        }],
    }],
}
_validate_frozen_3a_evidence({"subjects": [valid_record]}, locked)
forged_frozen = copy.deepcopy(valid_record)
forged_frozen["frozen_3a"]["spectrum_available"] = False
try:
    _validate_frozen_3a_evidence({"subjects": [forged_frozen]}, locked)
except RuntimeError:
    frozen_tamper_rejected = True
else:
    frozen_tamper_rejected = False
check(
    "frozen 3A availability is recomputed from locked QC evidence",
    frozen_tamper_rejected,
)

error_record = {
    "subject": "S1",
    "query_status": "error",
    "selection_classification": "inventory_error",
}
unindexed_error_inventory = {
    "requested_subjects": ["S1"],
    "n_requested": 1,
    "subjects": [error_record],
    "query_error_subjects": [],
    "classification_counts": {"inventory_error": 1},
    "classification_subjects": {"inventory_error": ["S1"]},
}
try:
    _validate_classifications(unindexed_error_inventory)
except RuntimeError:
    unindexed_error_rejected = True
else:
    unindexed_error_rejected = False
check(
    "terminal inventory cannot hide an error record from its error index",
    unindexed_error_rejected,
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
