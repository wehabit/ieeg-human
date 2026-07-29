"""Build and validate compact, reviewable QC-grid evidence.

Full grids contain repeated spectra, surrogate diagnostics, and event details
for every profile. They remain ignored local build products. This module
publishes:

* one compact summary per cohort/grid, retaining every profile definition,
  cohort count, and per-subject endpoint availability;
* one complete locked HUP profile used by the paired scalp–iEEG comparison;
* one terminal manifest hashing the exact nine-file public evidence set.

The terminal manifest is written last. An interrupted build therefore cannot
leave a mixed artifact set that validates as complete.
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import re
import uuid

from artifact_contracts import (
    QC_LOCKED_SNAPSHOT_SCHEMA,
    QC_PUBLIC_MANIFEST_SCHEMA,
    QC_PUBLIC_PIPELINE,
    QC_PUBLIC_SUMMARY_SCHEMA,
)
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    atomic_json_dump,
    file_sha256,
    git_is_dirty,
    git_revision,
    runtime_versions,
    source_tree_sha256,
    utc_now,
)
from qc_profiles import (
    PROFILE_SET_SCHEMA,
    expand_qc_grid,
    profile_file_sha256,
    qc_profile_sha256,
)


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRID_IDS = (
    "coverage_oat_v1",
    "staging_window_support_v1",
    "auxiliary_window_support_v1",
    "event_count_oat_v1",
)
COHORTS = ("hup", "respect")
LOCKED_PROFILE_ID = "overlap11_endpoint_local"
SUMMARY_SCHEMA = QC_PUBLIC_SUMMARY_SCHEMA
LOCKED_SCHEMA = QC_LOCKED_SNAPSHOT_SCHEMA
MANIFEST_SCHEMA = QC_PUBLIC_MANIFEST_SCHEMA
DEFAULT_INPUT = os.path.join(ROOT, "outputs", "qc_grid")
DEFAULT_OUTPUT = os.path.join(ROOT, "outputs", "qc_grid_public")
LOCKED_RELATIVE_PATH = os.path.join(
    "locked", f"{LOCKED_PROFILE_ID}__hup.json")
DEFAULT_LOCKED = os.path.join(DEFAULT_OUTPUT, LOCKED_RELATIVE_PATH)
MANIFEST_NAME = "RUN_MANIFEST.json"
_SHA256 = re.compile(r"^[0-9a-f]{64}$")

FULL_GRID_FIELDS = frozenset({
    "analysed_subjects",
    "analysis_version",
    "cache_directory_relative",
    "cache_files_sha256",
    "cache_manifest_run_id",
    "cache_manifest_sha256",
    "cache_pipeline",
    "cache_runtime_versions",
    "cache_schema_version",
    "explicitly_skipped_subjects",
    "generated_at_utc",
    "grid_id",
    "grid_role",
    "limitations",
    "profiles",
    "requested_subjects",
    "runtime_versions",
    "source_tree_sha256",
    "staging_calibration_application_to_this_cache",
    "staging_calibration_provenance",
})
PUBLIC_ARTIFACT_FIELDS = frozenset({
    *(FULL_GRID_FIELDS - {"profiles"}),
    "artifact_role",
    "artifact_schema_version",
    "full_grid_path_relative",
    "full_grid_sha256",
    "profiles",
})
FULL_PROFILE_FIELDS = frozenset({
    "qc_profile",
    "qc_profile_file_sha256",
    "qc_profile_id",
    "qc_profile_set_schema",
    "qc_profile_sha256",
    "subjects",
    "summary",
})
COMPACT_PROFILE_FIELDS = frozenset({
    *(FULL_PROFILE_FIELDS - {"subjects"}),
    "subject_availability",
})
SUBJECT_AVAILABILITY_FIELDS = frozenset({
    "3d_support_evaluable",
    "endpoint_available",
    "subject",
})
LOCKED_SUBJECT_FIELDS = frozenset({
    "result_3a",
    "result_3b",
    "result_3d",
    "subject",
    "support_diagnostics",
})
SUMMARY_COUNT_TO_ENDPOINT = {
    "n_3a_spectrum": "3a_spectrum",
    "n_3a_coherence": "3a_coherence",
    "n_3a_cross_correlation": "3a_cross_correlation",
    "n_3b_N2": "3b_N2",
    "n_3b_N3": "3b_N3",
    "n_3b_pooled_NREM_exploratory": "3b_pooled_NREM_exploratory",
    "n_3d_descriptive_effect": "3d_descriptive_effect",
    "n_3d_support_evaluable": "3d_support_evaluable",
}
SUMMARY_FIELDS = frozenset({
    "n_subjects",
    "minimum_participants",
    "meets_minimum_participants",
    *SUMMARY_COUNT_TO_ENDPOINT,
})
SUMMARY_ARTIFACT_ROLE = (
    "public QC sensitivity evidence: all profile definitions, "
    "cohort counts, and per-subject endpoint availability; numeric "
    "endpoint vectors remain in ignored regenerable full grids"
)
LOCKED_ARTIFACT_ROLE = (
    "full frozen HUP subject records for the locked profile used by "
    "the paired scalp–iEEG exact-match analysis"
)


def public_artifact_relative_paths():
    """Return the exact terminal public-evidence file set."""
    paths = [
        os.path.join(
            grid_id, f"{cohort}_qc_grid_summary.json"
        ).replace(os.sep, "/")
        for grid_id in GRID_IDS
        for cohort in COHORTS
    ]
    paths.append(LOCKED_RELATIVE_PATH.replace(os.sep, "/"))
    return tuple(paths)


PUBLIC_ARTIFACT_RELATIVE_PATHS = public_artifact_relative_paths()


def _reject_json_constant(value):
    raise ValueError(f"non-standard JSON constant {value!r}")


def _load_json(path):
    with open(path, encoding="utf-8") as handle:
        value = json.load(handle, parse_constant=_reject_json_constant)
    if not isinstance(value, dict):
        raise RuntimeError(f"{path} does not contain one JSON object")
    return value


def _require_sha256(value, label):
    if not isinstance(value, str) or not _SHA256.fullmatch(value):
        raise RuntimeError(f"{label} is not a lowercase SHA-256 digest")
    return value


def _require_exact_keys(value, expected, label):
    if not isinstance(value, dict):
        raise RuntimeError(f"{label} is not an object")
    actual = set(value)
    if actual != set(expected):
        raise RuntimeError(
            f"{label} fields differ: "
            f"missing={sorted(set(expected) - actual)}, "
            f"extra={sorted(actual - set(expected))}")


def _validate_locked_subject_records(records, path):
    if not isinstance(records, list):
        raise RuntimeError(f"{path} lacks full subject records")
    for index, record in enumerate(records):
        label = f"{path} subject record {index}"
        _require_exact_keys(record, LOCKED_SUBJECT_FIELDS, label)
        if not isinstance(record["subject"], str) or not record["subject"]:
            raise RuntimeError(f"{label} has an invalid subject")
        for key in (
            "result_3a",
            "result_3b",
            "result_3d",
            "support_diagnostics",
        ):
            if not isinstance(record[key], dict):
                raise RuntimeError(f"{label} {key} is not an object")


def _validate_summary_contract(profile, path):
    summary = profile.get("summary")
    _require_exact_keys(summary, SUMMARY_FIELDS, f"{path} summary")
    for key in ("n_subjects", *SUMMARY_COUNT_TO_ENDPOINT):
        value = summary[key]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise RuntimeError(
                f"{path} summary {key} is not a nonnegative integer")
    expected_minimum = (
        profile.get("qc_profile", {})
        .get("cohort", {})
        .get("minimum_participants")
    )
    if (
        isinstance(expected_minimum, bool)
        or not isinstance(expected_minimum, int)
        or expected_minimum < 1
    ):
        raise RuntimeError(
            f"{path} QC profile has an invalid cohort minimum")
    actual_minimum = summary["minimum_participants"]
    if (
        isinstance(actual_minimum, bool)
        or not isinstance(actual_minimum, int)
        or actual_minimum != expected_minimum
    ):
        raise RuntimeError(
            f"{path} summary minimum_participants differs from its QC profile")
    expected_meets = {
        endpoint: bool(summary[count] >= expected_minimum)
        for count, endpoint in SUMMARY_COUNT_TO_ENDPOINT.items()
    }
    actual_meets = summary["meets_minimum_participants"]
    if (
        not isinstance(actual_meets, dict)
        or set(actual_meets) != set(expected_meets)
        or any(not isinstance(value, bool) for value in actual_meets.values())
        or actual_meets != expected_meets
    ):
        raise RuntimeError(
            f"{path} summary meets_minimum_participants differs from its counts")


def _logical_full_grid_path(grid_id, cohort):
    return (
        f"outputs/qc_grid/{grid_id}/{cohort}_qc_grid.json")


def _subject_availability(record):
    """Keep endpoint availability only, never numeric vectors/contact records."""
    result_3a = record.get("result_3a", {})
    result_3b = record.get("result_3b", {})
    result_3d = record.get("result_3d", {})
    stages = result_3b.get("stages") or {}
    endpoints_3a = result_3a.get("endpoint_availability") or {}
    return {
        "subject": record.get("subject"),
        "endpoint_available": {
            "3a_spectrum": bool(endpoints_3a.get("spectrum", False)),
            "3a_fixed_0p02_coherence": bool(
                endpoints_3a.get("fixed_0p02_coherence", False)),
            "3a_cross_correlation": bool(
                endpoints_3a.get("cross_correlation", False)),
            "3b_N2": bool(
                (stages.get("N2") or {}).get(
                    "available_under_profile", False)),
            "3b_N3": bool(
                (stages.get("N3") or {}).get(
                    "available_under_profile", False)),
            "3b_NREM": bool(
                (stages.get("NREM") or {}).get(
                    "available_under_profile", False)),
            "3d_descriptive_effect": bool(
                result_3d.get("descriptive_effect_available", False)),
        },
        "3d_support_evaluable": (
            result_3d.get("support_passes_profile") is not None),
    }


def compact_summary(full_grid, *, source_path, source_sha256):
    _require_exact_keys(full_grid, FULL_GRID_FIELDS, "full QC grid")
    compact_profiles = []
    for profile in full_grid.get("profiles", []):
        _require_exact_keys(
            profile, FULL_PROFILE_FIELDS, "full QC-grid profile")
        compact = {
            key: copy.deepcopy(value)
            for key, value in profile.items()
            if key in COMPACT_PROFILE_FIELDS
        }
        compact["subject_availability"] = [
            _subject_availability(record)
            for record in profile.get("subjects", [])
        ]
        compact_profiles.append(compact)

    compact = {
        key: copy.deepcopy(full_grid[key])
        for key in FULL_GRID_FIELDS
        if key != "profiles"
    }
    compact.update({
        "artifact_schema_version": SUMMARY_SCHEMA,
        "artifact_role": SUMMARY_ARTIFACT_ROLE,
        "full_grid_path_relative": source_path,
        "full_grid_sha256": source_sha256,
        "profiles": compact_profiles,
    })
    return compact


def locked_profile_snapshot(full_grid, *, source_path, source_sha256):
    _require_exact_keys(full_grid, FULL_GRID_FIELDS, "full QC grid")
    matches = [
        profile for profile in full_grid.get("profiles", [])
        if profile.get("qc_profile_id") == LOCKED_PROFILE_ID
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one {LOCKED_PROFILE_ID!r} profile, found {len(matches)}")
    _require_exact_keys(
        matches[0], FULL_PROFILE_FIELDS, "locked full QC-grid profile")
    snapshot = {
        key: copy.deepcopy(full_grid[key])
        for key in FULL_GRID_FIELDS
        if key != "profiles"
    }
    snapshot.update({
        "artifact_schema_version": LOCKED_SCHEMA,
        "artifact_role": LOCKED_ARTIFACT_ROLE,
        "full_grid_path_relative": source_path,
        "full_grid_sha256": source_sha256,
        "profiles": [copy.deepcopy(matches[0])],
    })
    return snapshot


def _expected_profile_records(grid_id):
    return expand_qc_grid(grid_id)


def _validate_profile_contracts(profiles, grid_id, *, compact):
    expected = _expected_profile_records(grid_id)
    expected_ids = [profile["profile_id"] for profile in expected]
    actual_ids = [profile.get("qc_profile_id") for profile in profiles]
    if actual_ids != expected_ids:
        raise RuntimeError(
            f"{grid_id} profile set/order differs: "
            f"expected={expected_ids}, actual={actual_ids}")
    current_profile_file = profile_file_sha256()
    for actual, expected_profile in zip(profiles, expected):
        profile_id = expected_profile["profile_id"]
        expected_fields = (
            COMPACT_PROFILE_FIELDS if compact else FULL_PROFILE_FIELDS)
        _require_exact_keys(
            actual,
            expected_fields,
            f"{grid_id}/{profile_id} profile",
        )
        if actual.get("qc_profile") != expected_profile:
            raise RuntimeError(
                f"{grid_id}/{profile_id} QC profile definition differs")
        if actual.get("qc_profile_id") != profile_id:
            raise RuntimeError(
                f"{grid_id}/{profile_id} profile ID differs")
        if actual.get("qc_profile_sha256") != qc_profile_sha256(
                expected_profile):
            raise RuntimeError(
                f"{grid_id}/{profile_id} profile digest differs")
        if actual.get("qc_profile_file_sha256") != current_profile_file:
            raise RuntimeError(
                f"{grid_id}/{profile_id} profile-file digest is stale")
        if actual.get("qc_profile_set_schema") != PROFILE_SET_SCHEMA:
            raise RuntimeError(
                f"{grid_id}/{profile_id} profile-set schema differs")
        records_key = "subject_availability" if compact else "subjects"
        if not isinstance(actual.get(records_key), list):
            raise RuntimeError(
                f"{grid_id}/{profile_id} lacks {records_key}")
        _validate_summary_contract(actual, f"{grid_id}/{profile_id}")
        if not compact:
            _validate_locked_subject_records(
                actual["subjects"], f"{grid_id}/{profile_id}")


def _validate_full_grid(full_grid, grid_id, cohort):
    _require_exact_keys(
        full_grid, FULL_GRID_FIELDS, f"{grid_id}/{cohort} full grid")
    if full_grid.get("grid_id") != grid_id:
        raise RuntimeError(f"{grid_id}/{cohort} full grid has a different ID")
    for key, expected in (
        ("analysis_version", ANALYSIS_VERSION),
        ("cache_schema_version", CACHE_SCHEMA_VERSION),
        ("source_tree_sha256", source_tree_sha256(ROOT)),
        ("runtime_versions", runtime_versions()),
        ("cache_runtime_versions", runtime_versions()),
    ):
        if full_grid.get(key) != expected:
            raise RuntimeError(
                f"{grid_id}/{cohort} full grid differs at {key}")
    profiles = full_grid.get("profiles")
    if not isinstance(profiles, list) or not profiles:
        raise RuntimeError(f"{grid_id}/{cohort} full grid has no profiles")
    _validate_profile_contracts(profiles, grid_id, compact=False)
    _validate_grid_lineage(
        full_grid, f"{grid_id}/{cohort} full grid")
    subjects = list(full_grid.get("analysed_subjects", []))
    if not subjects or len(subjects) != len(set(subjects)):
        raise RuntimeError(
            f"{grid_id}/{cohort} full grid has an invalid analysed cohort")
    for profile in profiles:
        if [record.get("subject") for record in profile["subjects"]] != subjects:
            raise RuntimeError(
                f"{grid_id}/{cohort}/{profile['qc_profile_id']} "
                "subjects differ from the analysed cohort")


def _validate_grid_lineage(payload, path):
    requested = payload.get("requested_subjects")
    analysed = payload.get("analysed_subjects")
    skipped_records = payload.get("explicitly_skipped_subjects")
    if (
        not isinstance(requested, list)
        or not isinstance(analysed, list)
        or not isinstance(skipped_records, list)
        or any(
            not isinstance(subject, str) or not subject
            for subject in requested + analysed
        )
        or len(requested) != len(set(requested))
        or len(analysed) != len(set(analysed))
        or any(
            not isinstance(record, dict)
            or not isinstance(record.get("subject"), str)
            or not record.get("subject")
            or not record.get("reason")
            for record in skipped_records
        )
    ):
        raise RuntimeError(f"{path} has malformed cohort lineage")
    skipped = [record["subject"] for record in skipped_records]
    if (
        len(skipped) != len(set(skipped))
        or set(analysed) & set(skipped)
        or set(analysed) | set(skipped) != set(requested)
    ):
        raise RuntimeError(
            f"{path} analysed/skipped subjects do not partition requested")
    cache_hashes = payload.get("cache_files_sha256")
    if (
        not isinstance(cache_hashes, dict)
        or set(cache_hashes) != set(requested)
        or any(
            not _SHA256.fullmatch(value)
            for value in cache_hashes.values()
            if isinstance(value, str)
        )
        or any(not isinstance(value, str) for value in cache_hashes.values())
    ):
        raise RuntimeError(f"{path} has invalid cache-file hashes")
    _require_sha256(
        payload.get("cache_manifest_sha256"),
        f"{path} cache_manifest_sha256",
    )


def _validate_compaction(full_grid, compact):
    if len(compact["profiles"]) != len(full_grid.get("profiles", [])):
        raise RuntimeError("compact output lost a QC profile")
    for full_profile, compact_profile in zip(
            full_grid.get("profiles", []), compact["profiles"]):
        for key in (
            "qc_profile",
            "qc_profile_file_sha256",
            "qc_profile_id",
            "qc_profile_set_schema",
            "qc_profile_sha256",
            "summary",
        ):
            if compact_profile.get(key) != full_profile.get(key):
                raise RuntimeError(
                    f"compact profile differs at {key}: "
                    f"{full_profile.get('qc_profile_id')}")
        expected_availability = [
            _subject_availability(record)
            for record in full_profile.get("subjects", [])
        ]
        if compact_profile.get(
                "subject_availability") != expected_availability:
            raise RuntimeError(
                "compact per-subject availability differs from the full grid")


def _manifest_metadata(*, run_id, run_state, generated_at_utc):
    return {
        "schema_version": MANIFEST_SCHEMA,
        "pipeline": QC_PUBLIC_PIPELINE,
        "run_id": run_id,
        "run_state": run_state,
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "generated_at_utc": generated_at_utc,
        "code_revision": git_revision(ROOT),
        "code_dirty_at_start": git_is_dirty(ROOT),
        "source_tree_sha256": source_tree_sha256(ROOT),
        "profile_file_sha256": profile_file_sha256(),
        "runtime_versions": runtime_versions(),
    }


def _validate_summary_counts(path, profile):
    availability = profile.get("subject_availability")
    summary = profile.get("summary")
    if not isinstance(availability, list) or not isinstance(summary, dict):
        raise RuntimeError(f"{path} lacks availability/summary evidence")
    subjects = [value.get("subject") for value in availability]
    if (
        any(not isinstance(subject, str) or not subject for subject in subjects)
        or len(subjects) != len(set(subjects))
        or len(subjects) != int(summary.get("n_subjects", -1))
    ):
        raise RuntimeError(f"{path} has inconsistent subject availability")
    endpoint_keys = {
        "3a_spectrum",
        "3a_fixed_0p02_coherence",
        "3a_cross_correlation",
        "3b_N2",
        "3b_N3",
        "3b_NREM",
        "3d_descriptive_effect",
    }
    for value in availability:
        _require_exact_keys(
            value, SUBJECT_AVAILABILITY_FIELDS,
            f"{path} subject availability")
        endpoint_available = value.get("endpoint_available")
        if (
            not isinstance(endpoint_available, dict)
            or set(endpoint_available) != endpoint_keys
            or any(
                not isinstance(flag, bool)
                for flag in endpoint_available.values()
            )
            or not isinstance(value.get("3d_support_evaluable"), bool)
        ):
            raise RuntimeError(
                f"{path} has malformed endpoint-availability evidence")
    _validate_summary_contract(profile, path)
    counts = {
        "n_3a_spectrum": "3a_spectrum",
        "n_3a_coherence": "3a_fixed_0p02_coherence",
        "n_3a_cross_correlation": "3a_cross_correlation",
        "n_3b_N2": "3b_N2",
        "n_3b_N3": "3b_N3",
        "n_3b_pooled_NREM_exploratory": "3b_NREM",
        "n_3d_descriptive_effect": "3d_descriptive_effect",
    }
    for summary_key, availability_key in counts.items():
        actual = sum(
            value["endpoint_available"][availability_key]
            for value in availability
        )
        if actual != int(summary.get(summary_key, -1)):
            raise RuntimeError(f"{path} count differs at {summary_key}")
    evaluable = sum(
        value["3d_support_evaluable"] for value in availability)
    if evaluable != int(summary.get("n_3d_support_evaluable", -1)):
        raise RuntimeError(
            f"{path} count differs at n_3d_support_evaluable")


def _artifact_without_profile(payload):
    return {
        key: value for key, value in payload.items()
        if key not in {"artifact_schema_version", "artifact_role", "profiles"}
    }


def validate_public_artifacts(
    output_root=DEFAULT_OUTPUT,
    locked_path=None,
):
    """Fail closed unless the exact current nine-file evidence set agrees."""
    output_root = os.path.abspath(output_root)
    canonical_locked = os.path.join(output_root, LOCKED_RELATIVE_PATH)
    locked_path = (
        canonical_locked if locked_path is None
        else os.path.abspath(locked_path)
    )
    if locked_path != canonical_locked:
        raise RuntimeError(
            "locked snapshot must live under the selected public output root")

    manifest_path = os.path.join(output_root, MANIFEST_NAME)
    manifest = _load_json(manifest_path)
    expected_manifest_values = {
        "schema_version": MANIFEST_SCHEMA,
        "pipeline": QC_PUBLIC_PIPELINE,
        "run_state": "complete",
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "source_tree_sha256": source_tree_sha256(ROOT),
        "profile_file_sha256": profile_file_sha256(),
        "runtime_versions": runtime_versions(),
    }
    for key, expected in expected_manifest_values.items():
        if manifest.get(key) != expected:
            raise RuntimeError(
                f"{manifest_path} differs at {key}")
    result_hashes = manifest.get("result_files_sha256")
    if (
        not isinstance(result_hashes, dict)
        or set(result_hashes) != set(PUBLIC_ARTIFACT_RELATIVE_PATHS)
    ):
        raise RuntimeError(
            f"{manifest_path} does not name the exact public artifact set")
    for relative in PUBLIC_ARTIFACT_RELATIVE_PATHS:
        expected = _require_sha256(
            result_hashes[relative],
            f"result_files_sha256[{relative!r}]",
        )
        path = os.path.join(output_root, *relative.split("/"))
        if not os.path.isfile(path) or file_sha256(path) != expected:
            raise RuntimeError(
                f"public artifact bytes differ or are missing: {relative}")

    validated = []
    summaries = {}
    cohort_lineage = {}
    current_source = source_tree_sha256(ROOT)
    current_profile_file = profile_file_sha256()
    for grid_id in GRID_IDS:
        for cohort in COHORTS:
            relative = (
                f"{grid_id}/{cohort}_qc_grid_summary.json")
            path = os.path.join(output_root, *relative.split("/"))
            payload = _load_json(path)
            _require_exact_keys(
                payload, PUBLIC_ARTIFACT_FIELDS, f"{path} public summary")
            for key, expected in (
                ("artifact_schema_version", SUMMARY_SCHEMA),
                ("artifact_role", SUMMARY_ARTIFACT_ROLE),
                ("grid_id", grid_id),
                ("analysis_version", ANALYSIS_VERSION),
                ("cache_schema_version", CACHE_SCHEMA_VERSION),
                ("source_tree_sha256", current_source),
                ("runtime_versions", runtime_versions()),
                ("cache_runtime_versions", runtime_versions()),
                (
                    "full_grid_path_relative",
                    _logical_full_grid_path(grid_id, cohort),
                ),
            ):
                if payload.get(key) != expected:
                    raise RuntimeError(f"{path} differs at {key}")
            _validate_grid_lineage(payload, path)
            source_sha = _require_sha256(
                payload.get("full_grid_sha256"),
                f"{path} full_grid_sha256",
            )
            profiles = payload.get("profiles")
            if not isinstance(profiles, list) or not profiles:
                raise RuntimeError(f"{path} has no profiles")
            _validate_profile_contracts(profiles, grid_id, compact=True)
            subjects = list(payload.get("analysed_subjects", []))
            if (
                not subjects
                or len(subjects) != len(set(subjects))
            ):
                raise RuntimeError(f"{path} has an invalid analysed cohort")
            for profile in profiles:
                _validate_summary_counts(path, profile)
                if [
                    value.get("subject")
                    for value in profile["subject_availability"]
                ] != subjects:
                    raise RuntimeError(
                        f"{path} profile subjects differ from the cohort")
                if profile.get(
                        "qc_profile_file_sha256") != current_profile_file:
                    raise RuntimeError(
                        f"{path} has a stale profile-file digest")

            lineage = {
                key: copy.deepcopy(payload.get(key))
                for key in (
                    "analysis_version",
                    "cache_schema_version",
                    "source_tree_sha256",
                    "cache_pipeline",
                    "cache_manifest_sha256",
                    "cache_manifest_run_id",
                    "cache_files_sha256",
                    "requested_subjects",
                    "analysed_subjects",
                    "explicitly_skipped_subjects",
                    "staging_calibration_provenance",
                    "staging_calibration_application_to_this_cache",
                )
            }
            previous = cohort_lineage.setdefault(cohort, lineage)
            if previous != lineage:
                raise RuntimeError(
                    f"{path} does not share one {cohort} cache lineage")
            summaries[(grid_id, cohort)] = payload
            validated.append(path)

    source_hashes = manifest.get("source_full_grid_sha256")
    expected_source_keys = {
        _logical_full_grid_path(grid_id, cohort)
        for grid_id in GRID_IDS for cohort in COHORTS
    }
    if (
        not isinstance(source_hashes, dict)
        or set(source_hashes) != expected_source_keys
    ):
        raise RuntimeError(
            f"{manifest_path} lacks the exact full-grid source hash set")
    for (grid_id, cohort), payload in summaries.items():
        logical_path = _logical_full_grid_path(grid_id, cohort)
        expected = _require_sha256(
            source_hashes[logical_path],
            f"source_full_grid_sha256[{logical_path!r}]",
        )
        if payload["full_grid_sha256"] != expected:
            raise RuntimeError(
                f"{logical_path} source hash differs across public evidence")

    locked = _load_json(locked_path)
    _require_exact_keys(
        locked, PUBLIC_ARTIFACT_FIELDS, "locked-profile snapshot")
    event_summary = summaries[("event_count_oat_v1", "hup")]
    if (
        locked.get("artifact_schema_version") != LOCKED_SCHEMA
        or locked.get("artifact_role") != LOCKED_ARTIFACT_ROLE
        or _artifact_without_profile(locked)
        != _artifact_without_profile(event_summary)
    ):
        raise RuntimeError(
            "locked-profile snapshot and event-count summary lineage differ")
    locked_profiles = locked.get("profiles")
    if not isinstance(locked_profiles, list) or len(locked_profiles) != 1:
        raise RuntimeError("locked-profile snapshot must contain one profile")
    locked_profile = locked_profiles[0]
    _require_exact_keys(
        locked_profile, FULL_PROFILE_FIELDS, "locked-profile definition")
    _validate_locked_subject_records(
        locked_profile.get("subjects"), "locked-profile snapshot")
    _validate_summary_contract(locked_profile, "locked-profile snapshot")
    expected_locked = _expected_profile_records("event_count_oat_v1")[0]
    if (
        locked_profile.get("qc_profile_id") != LOCKED_PROFILE_ID
        or locked_profile.get("qc_profile") != expected_locked
        or locked_profile.get("qc_profile_sha256")
        != qc_profile_sha256(expected_locked)
        or locked_profile.get("qc_profile_file_sha256")
        != current_profile_file
        or locked_profile.get("qc_profile_set_schema") != PROFILE_SET_SCHEMA
        or not isinstance(locked_profile.get("subjects"), list)
    ):
        raise RuntimeError("locked-profile definition differs from current")
    compact_locked = next(
        profile for profile in event_summary["profiles"]
        if profile["qc_profile_id"] == LOCKED_PROFILE_ID
    )
    if (
        locked_profile.get("summary") != compact_locked.get("summary")
        or [
            _subject_availability(record)
            for record in locked_profile.get("subjects", [])
        ] != compact_locked.get("subject_availability")
    ):
        raise RuntimeError(
            "locked full records do not derive the compact locked evidence")
    if [
        value.get("subject") for value in locked_profile["subjects"]
    ] != list(locked.get("analysed_subjects", [])):
        raise RuntimeError(
            "locked-profile subjects differ from the analysed cohort")
    validated.append(locked_path)
    return validated


def build_public_artifacts(
    input_root=DEFAULT_INPUT,
    output_root=DEFAULT_OUTPUT,
    locked_path=None,
):
    """Build all public artifacts and publish one terminal manifest last."""
    input_root = os.path.abspath(input_root)
    output_root = os.path.abspath(output_root)
    canonical_locked = os.path.join(output_root, LOCKED_RELATIVE_PATH)
    locked_path = (
        canonical_locked if locked_path is None
        else os.path.abspath(locked_path)
    )
    if locked_path != canonical_locked:
        raise ValueError(
            "locked_path must be the canonical locked file inside output_root")

    os.makedirs(output_root, exist_ok=True)
    run_id = str(uuid.uuid4())
    generated_at = utc_now()
    manifest_path = os.path.join(output_root, MANIFEST_NAME)
    start_manifest = _manifest_metadata(
        run_id=run_id,
        run_state="in_progress",
        generated_at_utc=generated_at,
    )
    start_manifest.update({
        "result_files_sha256": {},
        "source_full_grid_sha256": {},
    })
    atomic_json_dump(start_manifest, manifest_path)

    written = []
    source_hashes = {}
    for grid_id in GRID_IDS:
        for cohort in COHORTS:
            source = os.path.join(
                input_root, grid_id, f"{cohort}_qc_grid.json")
            full_grid = _load_json(source)
            _validate_full_grid(full_grid, grid_id, cohort)
            source_sha256 = file_sha256(source)
            source_relative = _logical_full_grid_path(grid_id, cohort)
            source_hashes[source_relative] = source_sha256
            compact = compact_summary(
                full_grid,
                source_path=source_relative,
                source_sha256=source_sha256,
            )
            _validate_compaction(full_grid, compact)
            destination = os.path.join(
                output_root, grid_id, f"{cohort}_qc_grid_summary.json")
            atomic_json_dump(compact, destination)
            written.append(destination)

            if grid_id == "event_count_oat_v1" and cohort == "hup":
                locked = locked_profile_snapshot(
                    full_grid,
                    source_path=source_relative,
                    source_sha256=source_sha256,
                )
                atomic_json_dump(locked, locked_path)
                written.append(locked_path)

    if set(
        os.path.relpath(path, output_root).replace(os.sep, "/")
        for path in written
    ) != set(PUBLIC_ARTIFACT_RELATIVE_PATHS):
        raise RuntimeError("builder did not produce the exact public artifact set")
    complete_manifest = {
        **start_manifest,
        "run_state": "complete",
        "result_files_sha256": {
            relative: file_sha256(os.path.join(
                output_root, *relative.split("/")))
            for relative in PUBLIC_ARTIFACT_RELATIVE_PATHS
        },
        "source_full_grid_sha256": source_hashes,
    }
    atomic_json_dump(complete_manifest, manifest_path)
    try:
        validate_public_artifacts(output_root, locked_path)
    except Exception as error:
        atomic_json_dump(
            {
                **complete_manifest,
                "run_state": "failed",
                "validation_error": f"{type(error).__name__}: {error}",
            },
            manifest_path,
        )
        raise
    return written


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", default=DEFAULT_INPUT)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--locked-output",
        default=None,
        help=(
            "optional compatibility argument; when supplied it must equal "
            "<output-root>/locked/overlap11_endpoint_local__hup.json"
        ),
    )
    args = parser.parse_args()
    output_root = os.path.abspath(args.output_root)
    locked_output = (
        None if args.locked_output is None
        else os.path.abspath(args.locked_output)
    )
    paths = build_public_artifacts(
        input_root=os.path.abspath(args.input_root),
        output_root=output_root,
        locked_path=locked_output,
    )
    for path in paths:
        print(
            f"wrote {os.path.relpath(path, ROOT)} "
            f"sha256={file_sha256(path)}",
            flush=True,
        )
    manifest_path = os.path.join(output_root, MANIFEST_NAME)
    print(
        f"wrote {os.path.relpath(manifest_path, ROOT)} "
        f"sha256={file_sha256(manifest_path)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
