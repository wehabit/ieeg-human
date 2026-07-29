"""Compare one locked QC profile before and after a pipeline rebuild.

The report separates endpoint availability changes from numerical changes.
It accepts either a full ``run_qc_grid.py`` artifact or the compact public
locked-profile snapshot, so the same command can be reused for future data.
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

from pipeline_version import atomic_json_dump, file_sha256, utc_now


DEFAULT_PROFILE = "overlap11_endpoint_local"


def _reject_json_constant(value):
    raise ValueError(f"non-standard JSON constant {value!r}")


def _load(path):
    with open(path, encoding="utf-8") as handle:
        payload = json.load(
            handle, parse_constant=_reject_json_constant)
    if not isinstance(payload, dict):
        raise RuntimeError(f"{path} must contain one JSON object")
    return payload


def validate_rebuild_comparison_summary(
        root=None, summary_path=None):
    """Validate v3 by rebuilding it from the tracked snapshot and public v9."""
    from rebuild_comparison_evidence import (
        SUMMARY_RELATIVE,
        validate_summary,
    )

    source_root = os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))
    root = os.path.abspath(root if root is not None else source_root)
    summary_path = os.path.abspath(
        summary_path
        if summary_path is not None
        else os.path.join(root, *SUMMARY_RELATIVE.split("/")))
    return validate_summary(root, summary_path)


def _profile(payload, profile_id, *, source):
    matches = [
        value for value in payload.get("profiles", [])
        if value.get("qc_profile_id") == profile_id
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"{source} must contain exactly one {profile_id!r} profile")
    profile = matches[0]
    subjects = profile.get("subjects")
    if not isinstance(subjects, list):
        raise RuntimeError(
            f"{source} profile {profile_id!r} lacks full subject records")
    return profile


def _nested(value, *path):
    for key in path:
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    return value


def _subject_snapshot(record):
    result_3a = record.get("result_3a") or {}
    result_3b = record.get("result_3b") or {}
    result_3d = record.get("result_3d") or {}
    snapshot = {
        "availability.3a_spectrum": bool(
            _nested(result_3a, "endpoint_availability", "spectrum")),
        "availability.3a_fixed_0p02_coherence": bool(
            _nested(
                result_3a, "endpoint_availability",
                "fixed_0p02_coherence")),
        "availability.3a_cross_correlation": bool(
            _nested(
                result_3a, "endpoint_availability",
                "cross_correlation")),
        "availability.3d_descriptive_effect": bool(
            result_3d.get("descriptive_effect_available", False)),
        "support.3a_selected_contacts": result_3a.get(
            "n_selected_contacts"),
        "support.3a_nrem_epochs": result_3a.get("n_nrem_epochs"),
        "support.3a_bout_seconds": result_3a.get("bout_seconds"),
        "metric.3a_peak_hz": _nested(
            result_3a, "peak", "peak_hz"),
        "metric.3a_peak_accepted": _nested(
            result_3a, "peak", "accepted"),
        "metric.3a_coherence_at_0p02_hz": _nested(
            result_3a, "coherence", "at_0p02_hz"),
        "metric.3a_coherence_K": _nested(
            result_3a, "coherence", "K"),
        "metric.3a_cross_correlation_peak_r": _nested(
            result_3a, "cross_correlation",
            "lecci_direction_peak_r"),
        "metric.3a_cross_correlation_peak_lag_s": _nested(
            result_3a, "cross_correlation",
            "lecci_direction_peak_lag_s"),
        "metric.3d_participant_R": _nested(
            result_3d, "descriptive_effect", "participant_R"),
        "metric.3d_participant_phase_deg": _nested(
            result_3d, "descriptive_effect",
            "participant_preferred_phase_deg"),
        "support.3d_paired_events": _nested(
            result_3d, "descriptive_effect",
            "n_paired_events_across_qualified_contacts"),
    }
    for stage in ("N2", "N3", "NREM"):
        stage_result = _nested(result_3b, "stages", stage) or {}
        estimate = stage_result.get("estimate") or {}
        prefix = f"3b_{stage}"
        snapshot[f"availability.{prefix}"] = bool(
            stage_result.get("available_under_profile", False))
        snapshot[f"support.{prefix}_stable_seconds"] = (
            stage_result.get("stable_seconds"))
        snapshot[f"metric.{prefix}_local_change_pct"] = (
            estimate.get("event_locked_local_change_pct"))
        snapshot[f"metric.{prefix}_pct_above_stage_mean"] = (
            estimate.get("pct_above_stage_mean"))
        snapshot[f"metric.{prefix}_peak_lag_s"] = (
            estimate.get("peak_lag_s"))
        snapshot[f"support.{prefix}_n_so"] = (
            estimate.get("n_so_total"))
    return snapshot


def _equal(left, right):
    if left is None or right is None:
        return left is right
    if (
        isinstance(left, (bool, str))
        or isinstance(right, (bool, str))
    ):
        return left == right
    try:
        return bool(np.isclose(
            float(left), float(right), rtol=1e-9, atol=1e-12))
    except (TypeError, ValueError, OverflowError):
        return left == right


def compare_profiles(before, after, *, profile_id=DEFAULT_PROFILE):
    before_profile = _profile(
        before, profile_id, source="before artifact")
    after_profile = _profile(
        after, profile_id, source="after artifact")
    before_subjects = {
        value.get("subject"): value
        for value in before_profile["subjects"]
    }
    after_subjects = {
        value.get("subject"): value
        for value in after_profile["subjects"]
    }
    if (
        None in before_subjects
        or None in after_subjects
        or len(before_subjects) != len(before_profile["subjects"])
        or len(after_subjects) != len(after_profile["subjects"])
    ):
        raise RuntimeError("subject records must have unique nonempty IDs")

    availability_changes = []
    numerical_changes = []
    unchanged_numerical_fields = 0
    common_subjects = sorted(set(before_subjects) & set(after_subjects))
    for subject in common_subjects:
        left = _subject_snapshot(before_subjects[subject])
        right = _subject_snapshot(after_subjects[subject])
        for path in sorted(set(left) | set(right)):
            old = left.get(path)
            new = right.get(path)
            if _equal(old, new):
                if path.startswith(("metric.", "support.")):
                    unchanged_numerical_fields += 1
                continue
            change = {
                "subject": subject,
                "field": path,
                "before": old,
                "after": new,
            }
            if path.startswith("availability."):
                availability_changes.append(change)
            else:
                if (
                    old is not None
                    and new is not None
                    and not isinstance(old, bool)
                    and not isinstance(new, bool)
                ):
                    try:
                        change["after_minus_before"] = (
                            float(new) - float(old))
                    except (TypeError, ValueError, OverflowError):
                        pass
                numerical_changes.append(change)

    before_summary = before_profile.get("summary") or {}
    after_summary = after_profile.get("summary") or {}
    summary_changes = {
        key: {
            "before": before_summary.get(key),
            "after": after_summary.get(key),
            "after_minus_before": (
                after_summary[key] - before_summary[key]
                if (
                    isinstance(before_summary.get(key), (int, float))
                    and not isinstance(before_summary.get(key), bool)
                    and isinstance(after_summary.get(key), (int, float))
                    and not isinstance(after_summary.get(key), bool)
                )
                else None
            ),
        }
        for key in sorted(set(before_summary) | set(after_summary))
        if not _equal(before_summary.get(key), after_summary.get(key))
    }
    return {
        "profile_id": profile_id,
        "subjects_before": sorted(before_subjects),
        "subjects_after": sorted(after_subjects),
        "subjects_added": sorted(set(after_subjects) - set(before_subjects)),
        "subjects_removed": sorted(
            set(before_subjects) - set(after_subjects)),
        "summary_before": before_summary,
        "summary_after": after_summary,
        "summary_changes": summary_changes,
        "availability_changes": availability_changes,
        "numerical_or_support_changes": numerical_changes,
        "n_unchanged_numerical_or_support_fields": (
            unchanged_numerical_fields),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--before")
    parser.add_argument("--after")
    parser.add_argument("--profile", default=DEFAULT_PROFILE)
    parser.add_argument("--output")
    writer = parser.add_mutually_exclusive_group()
    writer.add_argument(
        "--capture-v8-snapshot",
        nargs=2,
        metavar=("HUP_GRID", "RESPECT_GRID"),
    )
    parser.add_argument("--legacy-summary")
    writer.add_argument("--write-machine-summary", action="store_true")
    parser.add_argument("--root")
    args = parser.parse_args()
    source_root = os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))
    root = os.path.abspath(args.root or source_root)
    if args.capture_v8_snapshot:
        if not args.legacy_summary:
            parser.error(
                "--legacy-summary is required with --capture-v8-snapshot")
        if args.before or args.after:
            parser.error(
                "--before/--after cannot be combined with a writer mode")
        from rebuild_comparison_evidence import (
            SNAPSHOT_RELATIVE,
            capture_v8_snapshot,
            write_json,
        )
        output = args.output or os.path.join(
            root, *SNAPSHOT_RELATIVE.split("/"))
        write_json(capture_v8_snapshot(
            *args.capture_v8_snapshot,
            legacy_summary_path=args.legacy_summary,
        ), output)
        return
    if args.write_machine_summary:
        if args.before or args.after or args.legacy_summary:
            parser.error(
                "comparison/capture inputs cannot be combined with "
                "--write-machine-summary")
        from rebuild_comparison_evidence import (
            SNAPSHOT_RELATIVE,
            SUMMARY_RELATIVE,
            build_summary,
            write_json,
        )
        snapshot = os.path.join(root, *SNAPSHOT_RELATIVE.split("/"))
        output = args.output or os.path.join(
            root, *SUMMARY_RELATIVE.split("/"))
        write_json(build_summary(root, snapshot), output)
        return
    if args.legacy_summary:
        parser.error(
            "--legacy-summary is only valid with --capture-v8-snapshot")
    if not args.before or not args.after:
        parser.error(
            "--before and --after are required unless a writer mode is used")
    before_path = os.path.abspath(args.before)
    after_path = os.path.abspath(args.after)
    report = {
        "generated_at_utc": utc_now(),
        "before_path": before_path,
        "before_sha256": file_sha256(before_path),
        "after_path": after_path,
        "after_sha256": file_sha256(after_path),
        **compare_profiles(
            _load(before_path),
            _load(after_path),
            profile_id=args.profile,
        ),
    }
    if args.output:
        atomic_json_dump(report, os.path.abspath(args.output))
    else:
        print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
