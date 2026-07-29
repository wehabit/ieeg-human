"""Regression checks for explicit before/after QC rebuild comparison."""

from compare_qc_rebuilds import compare_profiles


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

print("ALL QC-REBUILD COMPARISON CHECKS PASSED")
