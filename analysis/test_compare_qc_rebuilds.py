"""Regression checks for explicit before/after QC rebuild comparison."""

import copy
import json
import os
import tempfile

from compare_qc_rebuilds import (
    compare_profiles,
    validate_rebuild_comparison_summary,
)
from rebuild_comparison_evidence import (
    SNAPSHOT_RELATIVE,
    SUMMARY_RELATIVE,
    _legacy_v8_inputs,
    capture_v8_snapshot,
    validate_v8_snapshot,
)
from pipeline_version import file_sha256


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


def artifact(subjects, *, coherence_count, metric):
    return {
        "profiles": [{
            "qc_profile_id": "overlap11_endpoint_local",
            "summary": {
                "n_3a_coherence": coherence_count,
                "n_subjects": len(subjects),
            },
            "subjects": [
                {
                    "subject": subject,
                    "result_3a": {
                        "endpoint_availability": {
                            "spectrum": True,
                            "fixed_0p02_coherence": (
                                coherence_count > 0),
                            "cross_correlation": False,
                        },
                        "peak": {"accepted": False},
                        "coherence": (
                            {"at_0p02_hz": metric, "K": 5}
                            if coherence_count > 0
                            else None
                        ),
                    },
                    "result_3b": {"stages": {}},
                    "result_3d": {
                        "descriptive_effect_available": False,
                    },
                }
                for subject in subjects
            ],
        }],
    }


before = artifact(["S1"], coherence_count=0, metric=None)
after = artifact(["S1", "S2"], coherence_count=1, metric=0.4)
report = compare_profiles(before, after)

check(
    "comparison separates subject additions",
    report["subjects_added"] == ["S2"]
    and report["subjects_removed"] == [],
)
check(
    "comparison reports endpoint availability transitions",
    report["availability_changes"] == [{
        "subject": "S1",
        "field": "availability.3a_fixed_0p02_coherence",
        "before": False,
        "after": True,
    }],
)
check(
    "comparison reports newly available numerical values",
    any(
        value["field"] == "metric.3a_coherence_at_0p02_hz"
        and value["before"] is None
        and value["after"] == 0.4
        for value in report["numerical_or_support_changes"]
    ),
)
check(
    "comparison reports cohort-count changes",
    report["summary_changes"]["n_3a_coherence"][
        "after_minus_before"] == 1
    and report["summary_changes"]["n_subjects"][
        "after_minus_before"] == 1,
)

snapshot_path = os.path.join(ROOT, *SNAPSHOT_RELATIVE.split("/"))
with open(snapshot_path, encoding="utf-8") as handle:
    snapshot = json.load(handle)
validate_v8_snapshot(snapshot)
check(
    "tracked historical snapshot validates directly",
    snapshot["artifact_schema_version"]
    == "2026-07-v8-locked-profile-snapshot-v1",
)


def snapshot_tamper_is_rejected(mutator):
    attacked = copy.deepcopy(snapshot)
    mutator(attacked)
    try:
        validate_v8_snapshot(attacked)
    except RuntimeError:
        return True
    return False


check(
    "historical snapshot subject IDs must remain unique and sorted",
    snapshot_tamper_is_rejected(
        lambda value: value["event_count_sources"]["hup"][
            "subjects"][0].update({
                "subject": value["event_count_sources"]["hup"][
                    "subjects"][1]["subject"],
            })),
)
check(
    "historical source identities require valid SHA-256 digests",
    snapshot_tamper_is_rejected(
        lambda value: value["event_count_sources"]["hup"].update({
            "source_sha256": "not-a-digest",
        })),
)
check(
    "historical numeric field schema cannot omit a selected field",
    snapshot_tamper_is_rejected(
        lambda value: value["event_count_sources"]["hup"][
            "subjects"][0]["selected_numeric_or_support_fields"].pop(
                "metric.3a_peak_hz")),
)
check(
    "historical selected values must be finite numeric values or null",
    snapshot_tamper_is_rejected(
        lambda value: value["event_count_sources"]["hup"][
            "subjects"][0]["selected_numeric_or_support_fields"].update({
                "metric.3a_peak_hz": "0.02",
            })),
)

legacy_records = []
anchor_by_pair = {
    (value["grid_id"], value["cohort"]): value
    for value in snapshot["other_grid_identity_anchors"]
}
for grid_id in (
    "coverage_oat_v1",
    "staging_window_support_v1",
    "auxiliary_window_support_v1",
    "event_count_oat_v1",
):
    for cohort in ("hup", "respect"):
        source = (
            snapshot["event_count_sources"][cohort]
            if grid_id == "event_count_oat_v1"
            else anchor_by_pair[(grid_id, cohort)]
        )
        legacy_records.append({
            "after_sha256": "1" * 64,
            "before_sha256": source["source_sha256"],
            "cohort": cohort,
            "grid_id": grid_id,
            "profile_id": (
                "audit80"
                if grid_id == "coverage_oat_v1"
                else "overlap11_endpoint_local"
            ),
        })
legacy_index = _legacy_v8_inputs({
    "comparison_inputs": legacy_records,
})
check(
    "legacy capture ledger covers every historical grid/cohort identity",
    len(legacy_index) == 8
    and legacy_index[("event_count_oat_v1", "hup")][
        "before_sha256"
    ] == snapshot["event_count_sources"]["hup"]["source_sha256"],
)

with tempfile.TemporaryDirectory() as directory:
    hup_grid_path = os.path.join(directory, "hup.json")
    respect_grid_path = os.path.join(directory, "respect.json")
    legacy_path = os.path.join(directory, "legacy.json")
    with open(hup_grid_path, "w", encoding="utf-8") as handle:
        json.dump(
            artifact(["HUP_TEST"], coherence_count=0, metric=None),
            handle,
            allow_nan=False,
        )
    with open(respect_grid_path, "w", encoding="utf-8") as handle:
        json.dump(
            artifact(["RESP_TEST"], coherence_count=0, metric=None),
            handle,
            allow_nan=False,
        )
    fixture_records = copy.deepcopy(legacy_records)
    for record in fixture_records:
        if record["grid_id"] == "event_count_oat_v1":
            record["before_sha256"] = file_sha256(
                hup_grid_path
                if record["cohort"] == "hup"
                else respect_grid_path
            )
    with open(legacy_path, "w", encoding="utf-8") as handle:
        json.dump(
            {"comparison_inputs": fixture_records},
            handle,
            allow_nan=False,
        )
    captured = capture_v8_snapshot(
        hup_grid_path,
        respect_grid_path,
        legacy_summary_path=legacy_path,
    )
    check(
        "snapshot capture binds compact records to ledger-identified grids",
        captured["event_count_sources"]["hup"]["subjects"][0][
            "subject"
        ] == "HUP_TEST",
    )
    next(
        record for record in fixture_records
        if (
            record["grid_id"] == "event_count_oat_v1"
            and record["cohort"] == "hup"
        )
    )["before_sha256"] = "0" * 64
    with open(legacy_path, "w", encoding="utf-8") as handle:
        json.dump(
            {"comparison_inputs": fixture_records},
            handle,
            allow_nan=False,
        )
    try:
        capture_v8_snapshot(
            hup_grid_path,
            respect_grid_path,
            legacy_summary_path=legacy_path,
        )
    except RuntimeError:
        wrong_historical_grid_rejected = True
    else:
        wrong_historical_grid_rejected = False
    check(
        "snapshot capture rejects a grid outside the historical ledger",
        wrong_historical_grid_rejected,
    )

tracked = validate_rebuild_comparison_summary(ROOT)
check(
    "tracked v8-to-v9 summary validates against current public evidence",
    tracked["artifact_schema_version"]
    == "2026-07-v8-v9-comparison-summary-v3",
)

summary_path = os.path.join(
    ROOT, *SUMMARY_RELATIVE.split("/"))
with open(summary_path, encoding="utf-8") as handle:
    baseline_summary = json.load(handle)


def tamper_is_rejected(mutator):
    attacked = copy.deepcopy(baseline_summary)
    mutator(attacked)
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "V8_TO_V9_SUMMARY.json")
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(attacked, handle, allow_nan=False)
        try:
            validate_rebuild_comparison_summary(ROOT, path)
        except RuntimeError:
            return True
    return False


check(
    "manifest-rehashed public-evidence substitution is rejected",
    tamper_is_rejected(
        lambda value: value["inputs"].update({
            "v9_public_manifest_sha256": "0" * 64,
        })),
)
check(
    "historical event-count identity cannot be substituted",
    tamper_is_rejected(
        lambda value: value["inputs"]["v8_event_count_sources"][
            "hup"].update({
                "source_sha256": "0" * 64,
        })),
)
check(
    "fabricated locked v9 count is rejected",
    tamper_is_rejected(
        lambda value: value["cohorts"]["hup"][
            "counts_after"].update({"3a_spectrum": 21})),
)
check(
    "availability list cannot omit a stated endpoint-count transition",
    tamper_is_rejected(
        lambda value: value["cohorts"]["hup"].update({
            "availability_changes": [],
        })),
)
check(
    "numeric comparison totals are rebuilt rather than trusted",
    tamper_is_rejected(
        lambda value: value["cohorts"]["hup"][
            "selected_numeric_or_support_comparison"].update({
                "n_changed_selected_numeric_or_support_fields": 999999,
            })),
)
check(
    "changed-field ledger cannot be altered independently of evidence",
    tamper_is_rejected(
        lambda value: value["cohorts"]["hup"][
            "selected_numeric_or_support_comparison"][
                "changed_field_ledger"][0].update({
                    "direction": "fabricated",
                })),
)
check(
    "availability transitions are bound to real snapshot subjects",
    tamper_is_rejected(
        lambda value: value["cohorts"]["hup"][
            "availability_changes"][0].update({
                "subject": "NOT_A_REAL_SUBJECT",
            })),
)

print("ALL QC-REBUILD COMPARISON CHECKS PASSED")
