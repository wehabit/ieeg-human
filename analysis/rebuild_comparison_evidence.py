"""Reproducible compact evidence for the historical-v8/current-v9 comparison."""
from __future__ import annotations

import hashlib
import json
import math
import os
import re

from pipeline_version import atomic_json_dump, file_sha256


PROFILE_ID = "overlap11_endpoint_local"
SNAPSHOT_SCHEMA = "2026-07-v8-locked-profile-snapshot-v1"
SUMMARY_SCHEMA = "2026-07-v8-v9-comparison-summary-v3"
SNAPSHOT_RELATIVE = (
    "outputs/rebuild_comparison/V8_LOCKED_PROFILE_SNAPSHOT.json")
SUMMARY_RELATIVE = "outputs/rebuild_comparison/V8_TO_V9_SUMMARY.json"
COHORTS = ("hup", "respect")
COUNT_ENDPOINTS = {
    "3a_spectrum": ("endpoint_available", "3a_spectrum"),
    "3a_coherence": (
        "endpoint_available", "3a_fixed_0p02_coherence"),
    "3a_cross_correlation": (
        "endpoint_available", "3a_cross_correlation"),
    "3b_N2": ("endpoint_available", "3b_N2"),
    "3b_N3": ("endpoint_available", "3b_N3"),
    "3b_pooled_NREM_exploratory": (
        "endpoint_available", "3b_NREM"),
    "3d_descriptive_effect": (
        "endpoint_available", "3d_descriptive_effect"),
    "3d_support_evaluable": ("3d_support_evaluable",),
}
OTHER_GRIDS = (
    "auxiliary_window_support_v1",
    "coverage_oat_v1",
    "staging_window_support_v1",
)
ALL_GRIDS = (*OTHER_GRIDS, "event_count_oat_v1")
SELECTED_FIELD_KEYS = frozenset({
    "metric.3a_coherence_K",
    "metric.3a_coherence_at_0p02_hz",
    "metric.3a_cross_correlation_peak_lag_s",
    "metric.3a_cross_correlation_peak_r",
    "metric.3a_peak_accepted",
    "metric.3a_peak_hz",
    "metric.3b_N2_local_change_pct",
    "metric.3b_N2_pct_above_stage_mean",
    "metric.3b_N2_peak_lag_s",
    "metric.3b_N3_local_change_pct",
    "metric.3b_N3_pct_above_stage_mean",
    "metric.3b_N3_peak_lag_s",
    "metric.3b_NREM_local_change_pct",
    "metric.3b_NREM_pct_above_stage_mean",
    "metric.3b_NREM_peak_lag_s",
    "metric.3d_participant_R",
    "metric.3d_participant_phase_deg",
    "support.3a_bout_seconds",
    "support.3a_nrem_epochs",
    "support.3a_selected_contacts",
    "support.3b_N2_n_so",
    "support.3b_N2_stable_seconds",
    "support.3b_N3_n_so",
    "support.3b_N3_stable_seconds",
    "support.3b_NREM_n_so",
    "support.3b_NREM_stable_seconds",
    "support.3d_paired_events",
})
# Serialization-only bound: totals still count every change and the full
# canonical-ledger digest binds any entries omitted from the readable prefix.
CHANGED_FIELD_LEDGER_LIMIT = 100
IDENTITY_ONLY = (
    "identity anchor only; historical full grid is not tracked and its "
    "contents are not independently proven"
)
SNAPSHOT_EVIDENCE = (
    "tracked compact snapshot captured from the byte-identified historical "
    "v8 full grid"
)
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _reject_json_constant(value):
    raise ValueError(f"non-standard JSON constant {value!r}")


def _load(path):
    with open(path, encoding="utf-8") as handle:
        payload = json.load(
            handle,
            parse_constant=_reject_json_constant,
        )
    if not isinstance(payload, dict):
        raise RuntimeError(f"{path} must contain one JSON object")
    return payload


def _exact_keys(value, keys, label):
    if not isinstance(value, dict) or set(value) != set(keys):
        raise RuntimeError(f"{label} fields differ")


def _sha256(value, label):
    if not isinstance(value, str) or not _SHA256.fullmatch(value):
        raise RuntimeError(f"{label} is not a lowercase SHA-256 digest")
    return value


def _nested(value, path):
    for key in path:
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    return value


def _profile(payload, source):
    from compare_qc_rebuilds import _profile
    return _profile(payload, PROFILE_ID, source=source)


def _compact_subject(record):
    from compact_qc_grid_artifacts import _subject_availability
    from compare_qc_rebuilds import _subject_snapshot

    availability = _subject_availability(record)
    selected = {
        key: value
        for key, value in _subject_snapshot(record).items()
        if key.startswith(("metric.", "support."))
    }
    return {
        "subject": availability["subject"],
        "endpoint_available": availability["endpoint_available"],
        "3d_support_evaluable": availability["3d_support_evaluable"],
        "selected_numeric_or_support_fields": selected,
    }


def _compact_profile(payload, source):
    records = sorted(
        (_compact_subject(record)
         for record in _profile(payload, source)["subjects"]),
        key=lambda value: value["subject"],
    )
    validate_subjects(records, source, require_selected=True)
    return records


def _legacy_v8_inputs(payload):
    """Validate and index the historical summary's frozen-v8 identities."""
    records = payload.get("comparison_inputs")
    expected_pairs = {
        (grid, cohort) for grid in ALL_GRIDS for cohort in COHORTS
    }
    if not isinstance(records, list) or len(records) != len(expected_pairs):
        raise RuntimeError(
            "legacy comparison summary has an incomplete input ledger")
    indexed = {}
    for index, record in enumerate(records):
        label = f"legacy comparison input {index}"
        _exact_keys(
            record,
            {
                "after_sha256",
                "before_sha256",
                "cohort",
                "grid_id",
                "profile_id",
            },
            label,
        )
        pair = (record["grid_id"], record["cohort"])
        expected_profile = (
            "audit80"
            if record["grid_id"] == "coverage_oat_v1"
            else PROFILE_ID
        )
        if (
            pair not in expected_pairs
            or pair in indexed
            or record["profile_id"] != expected_profile
        ):
            raise RuntimeError(f"{label} identity differs")
        _sha256(record["before_sha256"], f"{label} before")
        _sha256(record["after_sha256"], f"{label} after")
        indexed[pair] = record
    if set(indexed) != expected_pairs:
        raise RuntimeError(
            "legacy comparison summary input partition differs")
    return indexed


def validate_subjects(records, label, *, require_selected):
    endpoint_keys = {
        "3a_spectrum",
        "3a_fixed_0p02_coherence",
        "3a_cross_correlation",
        "3b_N2",
        "3b_N3",
        "3b_NREM",
        "3d_descriptive_effect",
    }
    record_keys = {
        "subject",
        "endpoint_available",
        "3d_support_evaluable",
        "selected_numeric_or_support_fields",
    }
    if not isinstance(records, list) or not records:
        raise RuntimeError(f"{label} has no subject evidence")
    subjects = []
    for index, record in enumerate(records):
        _exact_keys(record, record_keys, f"{label}[{index}]")
        subject = record["subject"]
        endpoints = record["endpoint_available"]
        selected = record["selected_numeric_or_support_fields"]
        if not isinstance(subject, str) or not subject:
            raise RuntimeError(f"{label}[{index}] subject is invalid")
        _exact_keys(endpoints, endpoint_keys, f"{label}[{index}] endpoints")
        if any(not isinstance(value, bool) for value in endpoints.values()):
            raise RuntimeError(f"{label}[{index}] endpoints are invalid")
        if not isinstance(record["3d_support_evaluable"], bool):
            raise RuntimeError(f"{label}[{index}] 3D support is invalid")
        expected_selected = SELECTED_FIELD_KEYS if require_selected else set()
        _exact_keys(
            selected, expected_selected,
            f"{label}[{index}] selected fields")
        for key, value in selected.items():
            if key == "metric.3a_peak_accepted":
                if not isinstance(value, bool):
                    raise RuntimeError(
                        f"{label}[{index}] {key} is not Boolean")
                continue
            if value is None:
                continue
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
            ):
                raise RuntimeError(
                    f"{label}[{index}] {key} is not numeric or null")
            try:
                finite = math.isfinite(value)
            except (OverflowError, TypeError, ValueError):
                finite = False
            if not finite:
                raise RuntimeError(
                    f"{label}[{index}] {key} is not finite")
        subjects.append(subject)
    if subjects != sorted(subjects) or len(subjects) != len(set(subjects)):
        raise RuntimeError(f"{label} subjects are not unique and sorted")
    return records


def capture_v8_snapshot(
        hup_path, respect_path, *, legacy_summary_path):
    """Capture the two historical event-count grids without workstation paths."""
    legacy_inputs = _legacy_v8_inputs(
        _load(os.path.abspath(legacy_summary_path)))
    sources = {}
    for cohort, path in (("hup", hup_path), ("respect", respect_path)):
        path = os.path.abspath(path)
        source_sha256 = file_sha256(path)
        if source_sha256 != legacy_inputs[
                ("event_count_oat_v1", cohort)]["before_sha256"]:
            raise RuntimeError(
                f"{cohort} event-count grid does not match the historical "
                "comparison ledger")
        sources[cohort] = {
            "logical_source": (
                f"outputs/qc_grid/event_count_oat_v1/"
                f"{cohort}_qc_grid.json"),
            "source_sha256": source_sha256,
            "evidence_status": SNAPSHOT_EVIDENCE,
            "subjects": _compact_profile(
                _load(path), f"historical v8 {cohort} grid"),
        }
    anchors = [{
        "grid_id": grid,
        "cohort": cohort,
        "source_sha256": legacy_inputs[(grid, cohort)]["before_sha256"],
        "evidence_status": IDENTITY_ONLY,
    } for grid in OTHER_GRIDS for cohort in COHORTS]
    anchors.sort(key=lambda value: (value["grid_id"], value["cohort"]))
    snapshot = {
        "artifact_schema_version": SNAPSHOT_SCHEMA,
        "artifact_role": (
            "compact historical v8 locked-profile event-count evidence"),
        "profile_id": PROFILE_ID,
        "event_count_sources": sources,
        "other_grid_identity_anchors": anchors,
    }
    return validate_v8_snapshot(snapshot)


def validate_v8_snapshot(snapshot):
    _exact_keys(
        snapshot,
        {
            "artifact_schema_version",
            "artifact_role",
            "profile_id",
            "event_count_sources",
            "other_grid_identity_anchors",
        },
        "v8 snapshot",
    )
    if (
        snapshot["artifact_schema_version"] != SNAPSHOT_SCHEMA
        or snapshot["artifact_role"] != (
            "compact historical v8 locked-profile event-count evidence")
        or snapshot["profile_id"] != PROFILE_ID
    ):
        raise RuntimeError("v8 snapshot header differs")
    sources = snapshot["event_count_sources"]
    _exact_keys(sources, COHORTS, "v8 event-count sources")
    for cohort in COHORTS:
        source = sources[cohort]
        _exact_keys(
            source,
            {"logical_source", "source_sha256", "evidence_status", "subjects"},
            f"v8 {cohort} source",
        )
        if (
            source["logical_source"] != (
                f"outputs/qc_grid/event_count_oat_v1/"
                f"{cohort}_qc_grid.json")
            or source["evidence_status"] != SNAPSHOT_EVIDENCE
        ):
            raise RuntimeError(f"v8 {cohort} source identity differs")
        _sha256(source["source_sha256"], f"v8 {cohort} source")
        validate_subjects(
            source["subjects"], f"v8 {cohort}", require_selected=True)
    anchors = snapshot["other_grid_identity_anchors"]
    expected_pairs = {
        (grid, cohort) for grid in OTHER_GRIDS for cohort in COHORTS
    }
    if not isinstance(anchors, list):
        raise RuntimeError("v8 identity anchors are not a list")
    pairs = []
    for index, anchor in enumerate(anchors):
        _exact_keys(
            anchor,
            {"grid_id", "cohort", "source_sha256", "evidence_status"},
            f"v8 identity anchor {index}",
        )
        pair = (anchor["grid_id"], anchor["cohort"])
        if pair not in expected_pairs or anchor["evidence_status"] != IDENTITY_ONLY:
            raise RuntimeError(f"v8 identity anchor {index} is invalid")
        _sha256(anchor["source_sha256"], f"v8 identity anchor {index}")
        pairs.append(pair)
    if pairs != sorted(pairs) or set(pairs) != expected_pairs:
        raise RuntimeError("v8 identity-anchor partition differs")
    return snapshot


def _counts(records):
    result = {"subjects": len(records)}
    for endpoint, path in COUNT_ENDPOINTS.items():
        result[endpoint] = sum(
            bool(_nested(record, path)) for record in records)
    return result


def _availability_changes(before, after):
    left = {record["subject"]: record for record in before}
    right = {record["subject"]: record for record in after}
    if set(left) != set(right):
        raise RuntimeError("v8/v9 subject cohorts differ")
    changes = []
    for subject in sorted(left):
        for endpoint, path in COUNT_ENDPOINTS.items():
            old = bool(_nested(left[subject], path))
            new = bool(_nested(right[subject], path))
            if old != new:
                changes.append({
                    "subject": subject,
                    "endpoint": endpoint,
                    "before": old,
                    "after": new,
                })
    return changes


def _change_direction(before, after):
    """Describe one selected-field transition without scientific inference."""
    if before is None:
        return "null_to_value"
    if after is None:
        return "value_to_null"
    if isinstance(before, bool):
        return "false_to_true" if after else "true_to_false"
    return "increased" if after > before else "decreased"


def _canonical_sha256(value):
    encoded = json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _selected_field_comparison(before, after):
    """Compare selected fields and retain an auditable bounded change ledger."""
    from compare_qc_rebuilds import _equal

    left = {record["subject"]: record for record in before}
    right = {record["subject"]: record for record in after}
    if set(left) != set(right):
        raise RuntimeError("v8/v9 subject cohorts differ")
    changed_field_ledger = []
    unchanged = 0
    for subject in sorted(left):
        old = left[subject]["selected_numeric_or_support_fields"]
        new = right[subject]["selected_numeric_or_support_fields"]
        if set(old) != set(new):
            raise RuntimeError(
                f"v8/v9 selected-field partitions differ for {subject}")
        for field in sorted(set(old) | set(new)):
            before_value = old.get(field)
            after_value = new.get(field)
            if _equal(before_value, after_value):
                unchanged += 1
            else:
                changed_field_ledger.append({
                    "subject": subject,
                    "field": field,
                    "before": before_value,
                    "after": after_value,
                    "direction": _change_direction(
                        before_value, after_value),
                })
    total_changed = len(changed_field_ledger)
    emitted_ledger = changed_field_ledger[:CHANGED_FIELD_LEDGER_LIMIT]
    return {
        "evidence_status": (
            "recomputed from the tracked v8 snapshot and current public v9 "
            "locked-profile subject records"
        ),
        "n_changed_selected_numeric_or_support_fields": total_changed,
        "n_unchanged_selected_numeric_or_support_fields": unchanged,
        "changed_field_ledger_order": "subject_then_field",
        "changed_field_ledger_limit": CHANGED_FIELD_LEDGER_LIMIT,
        "n_changed_field_ledger_entries_emitted": len(emitted_ledger),
        "n_changed_field_ledger_entries_omitted": (
            total_changed - len(emitted_ledger)),
        "full_changed_field_ledger_sha256": _canonical_sha256(
            changed_field_ledger),
        "changed_field_ledger": emitted_ledger,
    }


def _current_v9(public_root, cohort):
    if cohort == "hup":
        payload = _load(os.path.join(
            public_root, "locked",
            f"{PROFILE_ID}__hup.json"))
        return (
            _compact_profile(payload, "public v9 HUP locked profile"),
            payload["full_grid_sha256"],
        )
    payload = _load(os.path.join(
        public_root, "event_count_oat_v1",
        "respect_qc_grid_summary.json"))
    profiles = [
        value for value in payload.get("profiles", [])
        if value.get("qc_profile_id") == PROFILE_ID
    ]
    if len(profiles) != 1:
        raise RuntimeError("public v9 RESPect locked profile is ambiguous")
    records = [{
        **value,
        "selected_numeric_or_support_fields": {},
    } for value in profiles[0]["subject_availability"]]
    records.sort(key=lambda value: value["subject"])
    validate_subjects(records, "public v9 RESPect", require_selected=False)
    return records, payload["full_grid_sha256"]


def build_summary(root, snapshot_path):
    """Deterministically rebuild the v3 summary from validated public evidence."""
    from compact_qc_grid_artifacts import validate_public_artifacts

    root = os.path.abspath(root)
    public_root = os.path.join(root, "outputs", "qc_grid_public")
    validate_public_artifacts(public_root)
    snapshot_path = os.path.abspath(snapshot_path)
    snapshot = validate_v8_snapshot(_load(snapshot_path))
    manifest_path = os.path.join(public_root, "RUN_MANIFEST.json")
    cohorts = {}
    v9_sources = {}
    for cohort in COHORTS:
        before = snapshot["event_count_sources"][cohort]["subjects"]
        after, source_sha256 = _current_v9(public_root, cohort)
        before_counts = _counts(before)
        after_counts = _counts(after)
        cohorts[cohort] = {
            "subjects_before": [value["subject"] for value in before],
            "subjects_after": [value["subject"] for value in after],
            "counts_before": before_counts,
            "counts_after": after_counts,
            "endpoint_count_changes": {
                endpoint: after_counts[endpoint] - before_counts[endpoint]
                for endpoint in sorted(COUNT_ENDPOINTS)
                if after_counts[endpoint] != before_counts[endpoint]
            },
            "availability_changes": _availability_changes(before, after),
            "selected_numeric_or_support_comparison": (
                _selected_field_comparison(before, after)
                if cohort == "hup"
                else {
                    "evidence_status": (
                        "not in the machine contract because the public v9 "
                        "RESPect artifact intentionally omits numeric vectors"
                    ),
                }
            ),
        }
        v9_sources[cohort] = {
            "logical_source": (
                f"outputs/qc_grid/event_count_oat_v1/"
                f"{cohort}_qc_grid.json"),
            "source_sha256": source_sha256,
        }
    return {
        "artifact_schema_version": SUMMARY_SCHEMA,
        "artifact_role": (
            "reproducible machine comparison of compact historical v8 "
            "evidence with current validated public v9 artifacts"
        ),
        "profile_id": PROFILE_ID,
        "inputs": {
            "v8_snapshot_path_relative": os.path.relpath(
                snapshot_path, root).replace(os.sep, "/"),
            "v8_snapshot_sha256": file_sha256(snapshot_path),
            "v9_public_manifest_path_relative": (
                "outputs/qc_grid_public/RUN_MANIFEST.json"),
            "v9_public_manifest_sha256": file_sha256(manifest_path),
            "v8_event_count_sources": {
                cohort: {
                    key: snapshot["event_count_sources"][cohort][key]
                    for key in (
                        "logical_source", "source_sha256", "evidence_status")
                }
                for cohort in COHORTS
            },
            "v9_event_count_sources": v9_sources,
            "other_v8_grid_identity_anchors": snapshot[
                "other_grid_identity_anchors"],
        },
        "cohorts": cohorts,
    }


def validate_summary(root, summary_path):
    root = os.path.abspath(root)
    snapshot_path = os.path.join(root, *SNAPSHOT_RELATIVE.split("/"))
    actual = _load(os.path.abspath(summary_path))
    expected = build_summary(root, snapshot_path)
    if actual != expected:
        raise RuntimeError(
            "v8-to-v9 machine summary differs from a fresh evidence rebuild")
    return actual


def write_json(payload, path):
    atomic_json_dump(payload, os.path.abspath(path))
