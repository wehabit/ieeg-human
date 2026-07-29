"""Synthetic checks for the bounded v8-to-v9 changed-field ledger."""

from rebuild_comparison_evidence import (
    CHANGED_FIELD_LEDGER_LIMIT,
    _selected_field_comparison,
)


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


def record(subject, fields):
    return {
        "subject": subject,
        "selected_numeric_or_support_fields": fields,
    }


before = [
    record("S02", {
        "metric.3a_peak_hz": 0.03,
        "support.3a_nrem_epochs": 2,
    }),
    record("S01", {
        "metric.3a_peak_accepted": False,
        "metric.3a_peak_hz": None,
        "support.3a_nrem_epochs": 10,
    }),
]
after = [
    record("S01", {
        "metric.3a_peak_accepted": True,
        "metric.3a_peak_hz": 0.02,
        "support.3a_nrem_epochs": 8,
    }),
    record("S02", {
        "metric.3a_peak_hz": None,
        "support.3a_nrem_epochs": 3,
    }),
]
comparison = _selected_field_comparison(before, after)
check(
    "ledger is ordered deterministically by subject then field",
    [
        (entry["subject"], entry["field"])
        for entry in comparison["changed_field_ledger"]
    ] == [
        ("S01", "metric.3a_peak_accepted"),
        ("S01", "metric.3a_peak_hz"),
        ("S01", "support.3a_nrem_epochs"),
        ("S02", "metric.3a_peak_hz"),
        ("S02", "support.3a_nrem_epochs"),
    ],
)
check(
    "ledger reports exact before/after values and neutral directions",
    [
        entry["direction"]
        for entry in comparison["changed_field_ledger"]
    ] == [
        "false_to_true",
        "null_to_value",
        "decreased",
        "value_to_null",
        "increased",
    ]
    and comparison["changed_field_ledger"][1]["before"] is None
    and comparison["changed_field_ledger"][1]["after"] == 0.02,
)
check(
    "comparison totals agree with its emitted ledger",
    comparison["n_changed_selected_numeric_or_support_fields"] == 5
    and comparison["n_unchanged_selected_numeric_or_support_fields"] == 0
    and comparison["n_changed_field_ledger_entries_emitted"] == 5
    and comparison["n_changed_field_ledger_entries_omitted"] == 0,
)
check(
    "canonical comparison is independent of input record order",
    comparison == _selected_field_comparison(
        list(reversed(before)), list(reversed(after))),
)

many_before = [
    record(f"S{index:03d}", {"support.3a_nrem_epochs": index})
    for index in range(CHANGED_FIELD_LEDGER_LIMIT + 7)
]
many_after = [
    record(f"S{index:03d}", {"support.3a_nrem_epochs": index + 1})
    for index in range(CHANGED_FIELD_LEDGER_LIMIT + 7)
]
bounded = _selected_field_comparison(many_before, many_after)
check(
    "ledger has an explicit deterministic size bound",
    bounded["changed_field_ledger_limit"]
    == CHANGED_FIELD_LEDGER_LIMIT
    and len(bounded["changed_field_ledger"])
    == CHANGED_FIELD_LEDGER_LIMIT
    and bounded["n_changed_field_ledger_entries_emitted"]
    == CHANGED_FIELD_LEDGER_LIMIT
    and bounded["n_changed_field_ledger_entries_omitted"] == 7,
)
check(
    "bounded ledger exposes the deterministic retained prefix",
    bounded["changed_field_ledger"][0]["subject"] == "S000"
    and bounded["changed_field_ledger"][-1]["subject"]
    == f"S{CHANGED_FIELD_LEDGER_LIMIT - 1:03d}",
)
check(
    "full-ledger digest covers entries beyond the retained prefix",
    bounded["full_changed_field_ledger_sha256"]
    != _selected_field_comparison(
        many_before[:-1], many_after[:-1]
    )["full_changed_field_ledger_sha256"],
)

try:
    _selected_field_comparison(before, after[:-1])
except RuntimeError:
    cohort_mismatch_rejected = True
else:
    cohort_mismatch_rejected = False
check(
    "comparison rejects unequal subject cohorts",
    cohort_mismatch_rejected,
)

partition_after = [
    record("S01", {
        "metric.3a_peak_accepted": True,
        "metric.3a_peak_hz": 0.02,
    }),
    after[1],
]
try:
    _selected_field_comparison(before, partition_after)
except RuntimeError:
    field_mismatch_rejected = True
else:
    field_mismatch_rejected = False
check(
    "comparison rejects unequal selected-field partitions",
    field_mismatch_rejected,
)

print("ALL REBUILD-COMPARISON LEDGER CHECKS PASSED")
