"""Strict clean-checkout validation for checked-in paired-EEG evidence.

Raw EEG and derived NPZ caches are intentionally absent from Git. A clean
checkout can still prove that every public inventory/result file is present,
byte-pinned, mutually consistent, and tied to the exact current source,
profile, schemas, and locked QC evidence. Publication mode separately requires
the ignored private cache manifests to be present.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

from artifact_contracts import (
    PAIRED_RESULT_SCHEMA,
    QC_LOCKED_SNAPSHOT_SCHEMA,
    SCALP_INVENTORY_SCHEMA,
)
from paired_reporting import (
    _group_summary,
    _normalized_csv_rows,
    _role_pair_csv_rows,
    _validate_csv_artifact,
)
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    file_sha256,
    source_tree_sha256,
    validate_recorded_release_provenance,
)
from qc_profiles import (
    load_qc_profile,
    profile_file_sha256,
    qc_profile_sha256,
)


INVENTORY_DIRECTORY = "outputs/paired_scalp_inventory"
INVENTORY_FILE = "hup_scalp_channel_inventory.json"
RESULT_DIRECTORY = "outputs/paired_scalp_ieeg"
LOCKED_PROFILE_ID = "overlap11_endpoint_local"
RESULT_FILES = frozenset({
    "group_summary.json",
    "paired_3A_C3.png",
    "paired_3A_C3.svg",
    "paired_3B_F3_Fz.png",
    "paired_3B_F3_Fz.svg",
    "paired_metrics.csv",
    "role_pair_metrics.csv",
    "subject_results.json",
})
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_INVENTORY_ROLE_TARGETS = {
    "c3": "C3",
    "c4": "C4",
    "f3": "F3",
    "f4": "F4",
    "fz": "FZ",
    "a1": "A1",
    "a2": "A2",
    "m1": "M1",
    "m2": "M2",
}
_CHANNEL_IDENTITY_FIELDS = {
    "revision_id",
    "data_check",
    "start_time_us",
    "end_time_us",
    "duration_us",
    "number_of_samples",
    "sample_rate_hz",
}


def _reject_json_constant(value):
    raise ValueError(f"non-standard JSON constant {value!r}")


def _require_file(path: Path, label: str) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"required {label} is missing: {path}")
    return path


def _load_json(path: Path, label: str) -> dict:
    _require_file(path, label)
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle, parse_constant=_reject_json_constant)
    if not isinstance(value, dict):
        raise RuntimeError(f"{label} must contain one JSON object")
    return value


def _require_sha256(value, label: str) -> str:
    if not isinstance(value, str) or not _SHA256.fullmatch(value):
        raise RuntimeError(f"{label} is not a lowercase SHA-256 digest")
    return value


def _repository_path(root: Path, relative: str, label: str) -> Path:
    """Resolve one repository-relative path without permitting traversal."""
    if (
        not isinstance(relative, str)
        or not relative
        or Path(relative).is_absolute()
    ):
        raise RuntimeError(f"{label} is not a repository-relative path")
    candidate = (root / relative).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as error:
        raise RuntimeError(f"{label} escapes the repository") from error
    return candidate


def _require_equal(left, right, label: str) -> None:
    if left != right:
        raise RuntimeError(f"paired evidence differs at {label}")


def _subject_list(value, label, *, allow_empty=False):
    if (
        not isinstance(value, list)
        or (not allow_empty and not value)
        or any(not isinstance(item, str) or not item for item in value)
        or len(value) != len(set(value))
    ):
        raise RuntimeError(f"{label} is not a unique subject list")
    return value


def _require_current_header(
    payload,
    *,
    label,
    pipeline,
    schema,
    terminal,
):
    expected = {
        "pipeline": pipeline,
        "schema_version": schema,
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
    }
    if terminal:
        expected["run_state"] = "complete"
    for key, value in expected.items():
        if payload.get(key) != value:
            raise RuntimeError(f"{label} differs at {key}")


def _normalized_inventory_label(label):
    if not isinstance(label, str) or not label:
        raise RuntimeError("inventory channel labels must be nonempty strings")
    value = label.strip().upper()
    head = value.rstrip("0123456789")
    tail = value[len(head):]
    if tail:
        tail = str(int(tail))
    return head + tail


def _validated_channel_identity(identity, label):
    if (
        not isinstance(identity, dict)
        or set(identity) != _CHANNEL_IDENTITY_FIELDS
        or not isinstance(identity.get("revision_id"), str)
        or not identity["revision_id"]
        or not isinstance(identity.get("data_check"), str)
        or not identity["data_check"]
    ):
        raise RuntimeError(f"{label} lacks a complete portal channel identity")
    for field in ("start_time_us", "end_time_us", "number_of_samples"):
        value = identity[field]
        if isinstance(value, bool) or not isinstance(value, int):
            raise RuntimeError(f"{label} has invalid identity field {field}")
    for field in ("duration_us", "sample_rate_hz"):
        value = identity[field]
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
        ):
            raise RuntimeError(f"{label} has invalid identity field {field}")
    if (
        identity["end_time_us"] <= identity["start_time_us"]
        or identity["number_of_samples"] <= 0
        or identity["duration_us"] <= 0
        or identity["sample_rate_hz"] <= 0
    ):
        raise RuntimeError(f"{label} has nonpositive portal series geometry")
    return identity


def _same_inventory_geometry(left, right):
    return (
        all(
            left[field] == right[field]
            for field in ("start_time_us", "end_time_us", "number_of_samples")
        )
        and abs(float(left["duration_us"]) - float(right["duration_us"])) <= 1
        and abs(
            float(left["sample_rate_hz"])
            - float(right["sample_rate_hz"])
        ) <= 1e-9
    )


def _same_inventory_channel_identity(left, right):
    return (
        left["revision_id"] == right["revision_id"]
        and left["data_check"] == right["data_check"]
        and _same_inventory_geometry(left, right)
    )


def _validated_json_activity(activity, label):
    if not isinstance(activity, dict):
        raise RuntimeError(f"{label} lacks activity evidence")
    nonflat = activity.get("numerically_nonflat_full_interval")
    finite_count = activity.get("finite_sample_count")
    dynamic_range = activity.get("raw_dynamic_range")
    if (
        not isinstance(nonflat, bool)
        or isinstance(finite_count, bool)
        or not isinstance(finite_count, int)
        or finite_count < 0
    ):
        raise RuntimeError(f"{label} has malformed activity evidence")
    if finite_count == 0:
        if nonflat or dynamic_range is not None:
            raise RuntimeError(
                f"{label} empty activity must be flat with a null range")
    elif (
        isinstance(dynamic_range, bool)
        or not isinstance(dynamic_range, (int, float))
        or not math.isfinite(float(dynamic_range))
        or dynamic_range < 0
        or (nonflat and (finite_count < 2 or dynamic_range <= 0))
    ):
        raise RuntimeError(f"{label} has invalid finite activity evidence")
    return activity


def json_safe_inventory_activity(activity_by_channel):
    """Convert validated NPZ activity summaries to standards-compliant JSON."""
    if not isinstance(activity_by_channel, dict):
        raise RuntimeError("sidecar activity evidence must be a channel map")
    safe = {}
    for channel, activity in activity_by_channel.items():
        if not isinstance(channel, str) or not isinstance(activity, dict):
            raise RuntimeError("sidecar activity evidence is malformed")
        finite_count = activity.get("finite_sample_count")
        nonflat = activity.get("numerically_nonflat_full_interval")
        dynamic_range = activity.get("raw_dynamic_range")
        if (
            isinstance(finite_count, bool)
            or not isinstance(finite_count, int)
            or finite_count < 0
            or not isinstance(nonflat, bool)
        ):
            raise RuntimeError(
                f"sidecar activity evidence is malformed for {channel}")
        if finite_count == 0:
            if nonflat or not (
                isinstance(dynamic_range, (int, float))
                and not isinstance(dynamic_range, bool)
                and math.isinf(float(dynamic_range))
                and float(dynamic_range) < 0
            ):
                raise RuntimeError(
                    f"empty sidecar activity is inconsistent for {channel}")
            dynamic_range = None
        elif (
            isinstance(dynamic_range, bool)
            or not isinstance(dynamic_range, (int, float))
            or not math.isfinite(float(dynamic_range))
            or dynamic_range < 0
        ):
            raise RuntimeError(
                f"finite sidecar activity is invalid for {channel}")
        safe[channel] = {
            "numerically_nonflat_full_interval": nonflat,
            "finite_sample_count": finite_count,
            "raw_dynamic_range": dynamic_range,
        }
        _validated_json_activity(safe[channel], f"{channel} sidecar")
    return safe


def derive_inventory_classification(record):
    """Validate one inventory record and derive its selection class."""
    if not isinstance(record, dict):
        raise RuntimeError("inventory subject record must be an object")
    if record.get("query_status") != "ok":
        return "inventory_error"
    subject = record.get("subject")
    labels = record.get("ordered_channel_labels")
    if (
        not isinstance(subject, str)
        or not subject
        or record.get("dataset_name") != subject
        or not isinstance(record.get("snapshot_id"), str)
        or not record["snapshot_id"]
        or not isinstance(labels, list)
        or any(not isinstance(label, str) or not label for label in labels)
    ):
        raise RuntimeError("inventory record has malformed portal metadata")
    labels_digest = hashlib.sha256(
        json.dumps(
            labels, ensure_ascii=False, separators=(",", ":")
        ).encode("utf-8")
    ).hexdigest()
    if record.get("ordered_channel_labels_sha256") != labels_digest:
        raise RuntimeError(
            f"{subject} ordered channel-label evidence differs from its hash")
    reference_label = record.get("pinned_reference_channel")
    if not isinstance(reference_label, str) or reference_label not in labels:
        raise RuntimeError(f"{subject} lacks its pinned reference channel")
    reference = _validated_channel_identity(
        record.get("pinned_reference_identity"),
        f"{subject} pinned reference",
    )

    roles = record.get("roles")
    if not isinstance(roles, dict) or set(roles) != set(
            _INVENTORY_ROLE_TARGETS):
        raise RuntimeError(f"{subject} has malformed role evidence")
    for role, target in _INVENTORY_ROLE_TARGETS.items():
        role_record = roles[role]
        expected_labels = [
            label for label in labels
            if _normalized_inventory_label(label) == target
        ]
        if (
            not isinstance(role_record, dict)
            or role_record.get("target_normalized_label") != target
            or not isinstance(role_record.get("matches"), list)
            or [
                match.get("label")
                for match in role_record["matches"]
                if isinstance(match, dict)
            ] != expected_labels
            or len(role_record["matches"]) != len(expected_labels)
        ):
            raise RuntimeError(f"{subject} role evidence differs for {role}")
        for match in role_record["matches"]:
            identity = _validated_channel_identity(
                match.get("identity"),
                f"{subject} {match.get('label')} inventory",
            )
            same_geometry = _same_inventory_geometry(identity, reference)
            if match.get("same_geometry_as_pinned_reference") is not same_geometry:
                raise RuntimeError(
                    f"{subject} has stale geometry evidence for "
                    f"{match.get('label')}")

    frozen = record.get("frozen_3a")
    frozen_fields = {
        "record_available",
        "spectrum_available",
        "coherence_available",
        "cross_correlation_available",
    }
    if (
        not isinstance(frozen, dict)
        or set(frozen) != frozen_fields
        or any(not isinstance(frozen[field], bool) for field in frozen_fields)
    ):
        raise RuntimeError(f"{subject} has malformed frozen 3A evidence")

    sidecar = record.get("sidecar")
    if sidecar is not None:
        sidecar_roles = sidecar.get("roles") if isinstance(sidecar, dict) else None
        identities = (
            sidecar.get("source_identity_by_channel")
            if isinstance(sidecar, dict) else None
        )
        activity = sidecar.get("channels") if isinstance(sidecar, dict) else None
        if (
            not isinstance(sidecar_roles, dict)
            or not sidecar_roles
            or not set(sidecar_roles).issubset(_INVENTORY_ROLE_TARGETS)
            or not isinstance(identities, dict)
            or not isinstance(activity, dict)
            or set(identities) != set(sidecar_roles.values())
            or set(activity) != set(sidecar_roles.values())
            or sidecar.get("source_dataset") != subject
            or sidecar.get("source_snapshot_id") != record["snapshot_id"]
            or not isinstance(sidecar.get("path_relative"), str)
            or not _SHA256.fullmatch(str(sidecar.get("sha256", "")))
        ):
            raise RuntimeError(f"{subject} has malformed sidecar evidence")
        for role, channel in sidecar_roles.items():
            matches = roles[role]["matches"]
            if len(matches) != 1 or matches[0]["label"] != channel:
                raise RuntimeError(
                    f"{subject} sidecar role {role} is not freshly resolved")
            current_identity = matches[0]["identity"]
            sidecar_identity = _validated_channel_identity(
                identities[channel],
                f"{subject} {channel} sidecar",
            )
            if not _same_inventory_channel_identity(
                    sidecar_identity, current_identity):
                raise RuntimeError(
                    f"{subject} sidecar source changed for {channel}")
            _validated_json_activity(
                activity[channel], f"{subject} {channel}")

    c3 = roles["c3"]
    if len(c3["matches"]) == 0:
        return "no_c3_or_c03_label"
    if len(c3["matches"]) != 1:
        return "ambiguous_c3_or_c03_labels"
    if not c3["matches"][0]["same_geometry_as_pinned_reference"]:
        return "c3_geometry_mismatch"
    if not frozen["spectrum_available"]:
        return "c3_present_but_frozen_3a_spectrum_unavailable"
    if sidecar is None:
        return "requires_full_interval_activity_sidecar"
    label = c3["matches"][0]["label"]
    c3_activity = sidecar["channels"].get(label)
    if c3_activity is None:
        return "c3_not_streamed_in_sidecar"
    if not c3_activity["numerically_nonflat_full_interval"]:
        return "c3_numerically_flat"
    return "paired_3a_eligible"


def _validate_classifications(inventory):
    requested = _subject_list(
        inventory.get("requested_subjects"),
        "inventory requested_subjects",
    )
    records = inventory.get("subjects")
    if (
        not isinstance(records, list)
        or [record.get("subject") for record in records] != requested
        or inventory.get("n_requested") != len(requested)
    ):
        raise RuntimeError(
            "inventory records do not exactly preserve the requested cohort")
    query_errors = [
        record["subject"]
        for record in records
        if record.get("query_status") != "ok"
    ]
    if inventory.get("query_error_subjects") != query_errors:
        raise RuntimeError(
            "inventory query-error index differs from its subject records")
    if query_errors:
        raise RuntimeError("inventory contains unresolved query errors")
    classes = inventory.get("classification_subjects")
    counts = inventory.get("classification_counts")
    if not isinstance(classes, dict) or not isinstance(counts, dict):
        raise RuntimeError("inventory lacks classification evidence")
    flattened = []
    for key, subjects in classes.items():
        _subject_list(
            subjects, f"inventory classification {key}", allow_empty=True)
        if counts.get(key) != len(subjects):
            raise RuntimeError(
                f"inventory classification count differs for {key}")
        flattened.extend(subjects)
    if (
        len(flattened) != len(set(flattened))
        or set(flattened) != set(requested)
        or set(counts) != set(classes)
    ):
        raise RuntimeError(
            "inventory classifications do not partition the requested cohort")
    declared_class = {
        subject: classification
        for classification, subjects in classes.items()
        for subject in subjects
    }
    for record in records:
        subject = record["subject"]
        derived = derive_inventory_classification(record)
        if (
            record.get("selection_classification") != derived
            or declared_class[subject] != derived
        ):
            raise RuntimeError(
                f"inventory classification is not supported for {subject}")
    unresolved_activity = [
        subject
        for classification in (
            "requires_full_interval_activity_sidecar",
            "c3_not_streamed_in_sidecar",
        )
        for subject in classes.get(classification, [])
    ]
    if unresolved_activity:
        raise RuntimeError(
            "inventory has spectrum-eligible C3 candidates without "
            "usable full-interval activity sidecars")
    return requested


def _validate_result_metadata(
    result_manifest,
    subject_results,
    group_summary,
):
    for payload, label, terminal in (
        (result_manifest, "paired-result manifest", True),
        (subject_results, "paired subject results", False),
        (group_summary, "paired group summary", False),
    ):
        _require_current_header(
            payload,
            label=label,
            pipeline="paired_scalp_ieeg_comparison",
            schema=PAIRED_RESULT_SCHEMA,
            terminal=terminal,
        )

    # The producer writes one metadata object into all three artifacts. Compare
    # that complete object, not a hand-picked subset that could drift.
    metadata = {
        key: value for key, value in group_summary.items()
        if key != "group_summary"
    }
    subject_metadata = {
        key: value for key, value in subject_results.items()
        if key not in {"group_summary", "subjects"}
    }
    manifest_metadata = {
        key: value for key, value in result_manifest.items()
        if key not in {"run_state", "result_files_sha256"}
    }
    _require_equal(subject_metadata, metadata, "result metadata")
    _require_equal(manifest_metadata, metadata, "manifest metadata")
    _require_equal(
        subject_results.get("group_summary"),
        group_summary.get("group_summary"),
        "saved group summary",
    )

    requested = _subject_list(
        result_manifest.get("requested_subjects"),
        "paired requested_subjects",
    )
    completed = _subject_list(
        result_manifest.get("completed_subjects"),
        "paired completed_subjects",
    )
    records = subject_results.get("subjects")
    if (
        requested != completed
        or not isinstance(records, list)
        or [record.get("subject") for record in records] != completed
    ):
        raise RuntimeError(
            "paired subject records do not exactly match requested/completed")
    recomputed_group_summary = _group_summary(records)
    _require_equal(
        subject_results.get("group_summary"),
        recomputed_group_summary,
        "group summary recomputed from subject records",
    )
    _require_equal(
        group_summary.get("group_summary"),
        recomputed_group_summary,
        "saved group summary recomputed from subject records",
    )
    return requested


def _validate_frozen_3a_evidence(inventory, locked):
    """Recompute every recorded 3A availability flag from locked QC evidence."""
    profiles = locked.get("profiles") if isinstance(locked, dict) else None
    if not isinstance(profiles, list) or len(profiles) != 1:
        raise RuntimeError("locked QC evidence lacks one profile")
    locked_records = profiles[0].get("subjects")
    if not isinstance(locked_records, list):
        raise RuntimeError("locked QC evidence lacks subject records")
    locked_by_subject = {}
    for value in locked_records:
        subject = value.get("subject") if isinstance(value, dict) else None
        if not isinstance(subject, str) or subject in locked_by_subject:
            raise RuntimeError("locked QC subject evidence is malformed")
        locked_by_subject[subject] = value
    for record in inventory["subjects"]:
        subject = record["subject"]
        locked_record = locked_by_subject.get(subject)
        result_3a = (
            locked_record.get("result_3a", {})
            if isinstance(locked_record, dict) else {}
        )
        availability = result_3a.get("endpoint_availability", {})
        availability_fields = (
            "spectrum",
            "fixed_0p02_coherence",
            "cross_correlation",
        )
        if locked_record is not None and (
            not isinstance(result_3a, dict)
            or not isinstance(availability, dict)
            or any(
                not isinstance(availability.get(field), bool)
                for field in availability_fields
            )
        ):
            raise RuntimeError(
                f"{subject} locked 3A availability is malformed")
        expected = {
            "record_available": subject in locked_by_subject,
            "spectrum_available": availability.get("spectrum", False),
            "coherence_available": availability.get(
                "fixed_0p02_coherence", False),
            "cross_correlation_available": availability.get(
                "cross_correlation", False),
        }
        if record.get("frozen_3a") != expected:
            raise RuntimeError(
                f"{subject} frozen 3A evidence differs from the locked grid")


def _validate_tracked_lineage(
    root,
    inventory,
    inventory_manifest,
    result_manifest,
):
    source_pin_relative = inventory.get("source_pin_path_relative")
    qc_grid_relative = inventory.get("qc_grid_path_relative")
    source_pin = _repository_path(
        root, source_pin_relative, "source-pin path")
    qc_grid = _repository_path(root, qc_grid_relative, "QC-grid path")
    _require_file(source_pin, "source pin")
    _require_file(qc_grid, "locked QC grid")

    source_pin_sha = file_sha256(source_pin)
    qc_grid_sha = file_sha256(qc_grid)
    _require_equal(
        source_pin_sha,
        inventory.get("source_pin_sha256"),
        "inventory source_pin_sha256",
    )
    _require_equal(
        source_pin_sha,
        inventory_manifest.get("source_pin_sha256"),
        "inventory-manifest source_pin_sha256",
    )
    for payload, label in (
        (inventory, "inventory"),
        (inventory_manifest, "inventory manifest"),
    ):
        _require_equal(
            payload.get("qc_grid_sha256"),
            qc_grid_sha,
            f"{label} qc_grid_sha256",
        )
    _require_equal(
        result_manifest.get("frozen_qc_grid_path_relative"),
        qc_grid_relative,
        "result frozen_qc_grid_path_relative",
    )
    _require_equal(
        result_manifest.get("frozen_qc_grid_sha256"),
        qc_grid_sha,
        "result frozen_qc_grid_sha256",
    )

    profile = load_qc_profile(LOCKED_PROFILE_ID)
    profile_sha = qc_profile_sha256(profile)
    profile_file_sha = profile_file_sha256()
    for payload, label in (
        (inventory, "inventory"),
        (inventory_manifest, "inventory manifest"),
    ):
        _require_equal(
            payload.get("qc_profile_id"),
            LOCKED_PROFILE_ID,
            f"{label} qc_profile_id",
        )
        _require_equal(
            payload.get("qc_profile_sha256"),
            profile_sha,
            f"{label} qc_profile_sha256",
        )
        _require_equal(
            payload.get("current_locked_qc_profile_sha256"),
            profile_sha,
            f"{label} current_locked_qc_profile_sha256",
        )
    _require_equal(
        result_manifest.get("profile_id"),
        LOCKED_PROFILE_ID,
        "result profile_id",
    )
    _require_equal(
        result_manifest.get("profile_sha256"),
        profile_sha,
        "result profile_sha256",
    )
    _require_equal(
        result_manifest.get("profile_file_sha256"),
        profile_file_sha,
        "result profile_file_sha256",
    )

    locked = _load_json(qc_grid, "locked QC grid")
    expected_locked_header = {
        "artifact_schema_version": QC_LOCKED_SNAPSHOT_SCHEMA,
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "source_tree_sha256": source_tree_sha256(root),
        "grid_id": "event_count_oat_v1",
        "full_grid_path_relative": (
            "outputs/qc_grid/event_count_oat_v1/hup_qc_grid.json"),
    }
    for key, expected in expected_locked_header.items():
        if locked.get(key) != expected:
            raise RuntimeError(f"locked QC grid differs at {key}")
    profiles = locked.get("profiles")
    if not isinstance(profiles, list) or len(profiles) != 1:
        raise RuntimeError("locked QC grid must contain exactly one profile")
    locked_profile = profiles[0]
    if (
        locked_profile.get("qc_profile_id") != LOCKED_PROFILE_ID
        or locked_profile.get("qc_profile") != profile
        or locked_profile.get("qc_profile_sha256") != profile_sha
        or locked_profile.get("qc_profile_file_sha256") != profile_file_sha
    ):
        raise RuntimeError("locked QC grid profile differs from current")
    _require_equal(
        result_manifest.get("frozen_qc_grid_profile_sha256"),
        profile_sha,
        "result frozen_qc_grid_profile_sha256",
    )
    return locked


def validate_checked_in_paired_evidence(root) -> dict:
    """Validate all public paired evidence without reading private NPZ caches."""
    root = Path(root).resolve()
    inventory_dir = root / INVENTORY_DIRECTORY
    result_dir = root / RESULT_DIRECTORY
    inventory_path = inventory_dir / INVENTORY_FILE
    inventory_manifest_path = inventory_dir / "RUN_MANIFEST.json"
    result_manifest_path = result_dir / "RUN_MANIFEST.json"

    inventory = _load_json(inventory_path, "scalp inventory")
    inventory_manifest = _load_json(
        inventory_manifest_path, "scalp-inventory manifest")
    result_manifest = _load_json(
        result_manifest_path, "paired-result manifest")
    subject_results = _load_json(
        result_dir / "subject_results.json", "paired subject results")
    group_summary = _load_json(
        result_dir / "group_summary.json", "paired group summary")

    for payload, label in (
        (inventory, "scalp inventory"),
        (inventory_manifest, "scalp-inventory manifest"),
    ):
        _require_current_header(
            payload,
            label=label,
            pipeline="audit_hup_scalp_inventory",
            schema=SCALP_INVENTORY_SCHEMA,
            terminal=True,
        )
    inventory_provenance = validate_recorded_release_provenance(
        inventory, root, label="scalp inventory")
    inventory_manifest_provenance = validate_recorded_release_provenance(
        inventory_manifest, root, label="scalp-inventory manifest")
    _require_equal(
        inventory_provenance,
        inventory_manifest_provenance,
        "inventory release provenance",
    )
    inventory_subjects = _validate_classifications(inventory)
    if inventory_manifest.get("inventory_file") != INVENTORY_FILE:
        raise RuntimeError("scalp-inventory manifest names a different file")
    inventory_sha = file_sha256(inventory_path)
    _require_equal(
        inventory_sha,
        inventory_manifest.get("inventory_file_sha256"),
        "inventory_file_sha256",
    )
    for key in (
        "requested_subjects",
        "n_requested",
        "classification_counts",
        "classification_subjects",
        "query_error_subjects",
        "source_tree_sha256",
        "frozen_cache_manifest_sha256",
        "source_pin_sha256",
        "qc_grid_sha256",
        "qc_profile_id",
        "qc_profile_sha256",
        "current_locked_qc_profile_sha256",
        "sidecar_manifest_sha256",
    ):
        _require_equal(
            inventory.get(key),
            inventory_manifest.get(key),
            f"inventory manifest {key}",
        )

    paired_subjects = _validate_result_metadata(
        result_manifest, subject_results, group_summary)
    paired_provenance = validate_recorded_release_provenance(
        result_manifest, root, label="paired-result manifest")
    for payload, label in (
        (subject_results, "paired subject results"),
        (group_summary, "paired group summary"),
    ):
        _require_equal(
            validate_recorded_release_provenance(
                payload, root, label=label),
            paired_provenance,
            f"{label} release provenance",
        )
    declared_results = result_manifest.get("result_files_sha256")
    if (
        not isinstance(declared_results, dict)
        or set(declared_results) != RESULT_FILES
    ):
        missing = sorted(RESULT_FILES - set(declared_results or {}))
        extra = sorted(set(declared_results or {}) - RESULT_FILES)
        raise RuntimeError(
            f"paired-result evidence differs: missing={missing}, extra={extra}")
    for name in sorted(RESULT_FILES):
        expected = _require_sha256(
            declared_results[name], f"result_files_sha256[{name!r}]")
        actual = file_sha256(_require_file(
            result_dir / name, f"paired result {name}"))
        _require_equal(actual, expected, f"result_files_sha256[{name!r}]")

    role_pair_rows = _role_pair_csv_rows(subject_results["subjects"])
    normalized_rows = _normalized_csv_rows(
        subject_results["subjects"], role_pair_rows)
    _validate_csv_artifact(
        result_dir / "paired_metrics.csv",
        normalized_rows,
    )
    _validate_csv_artifact(
        result_dir / "role_pair_metrics.csv",
        role_pair_rows,
    )

    current_source = source_tree_sha256(root)
    for payload, label in (
        (inventory, "inventory"),
        (inventory_manifest, "inventory manifest"),
        (result_manifest, "paired-result manifest"),
        (subject_results, "subject results"),
        (group_summary, "group summary"),
    ):
        _require_equal(
            payload.get("source_tree_sha256"),
            current_source,
            f"{label} source_tree_sha256",
        )

    inventory_manifest_sha = file_sha256(inventory_manifest_path)
    unique_f3_f4 = [
        record["subject"]
        for record in inventory["subjects"]
        if (
            len(
                record.get("roles", {}).get(
                    "f3", {}).get("matches", [])
            ) == 1
            and len(
                record.get("roles", {}).get(
                    "f4", {}).get("matches", [])
            ) == 1
        )
    ]
    for payload, label in (
        (result_manifest, "paired-result manifest"),
        (subject_results, "paired subject results"),
        (group_summary, "paired group summary"),
    ):
        pinned = payload.get("scalp_inventory")
        if not isinstance(pinned, dict):
            raise RuntimeError(f"{label} does not pin the scalp inventory")
        expected_pin = {
            "path_relative": (
                f"{INVENTORY_DIRECTORY}/{INVENTORY_FILE}"),
            "sha256": inventory_sha,
            "manifest_path_relative": (
                f"{INVENTORY_DIRECTORY}/RUN_MANIFEST.json"),
            "manifest_sha256": inventory_manifest_sha,
            "n_frozen_hup_participants": len(inventory_subjects),
            "classification_counts": inventory["classification_counts"],
            "classification_subjects": inventory[
                "classification_subjects"],
            "paired_3a_eligible_subjects": inventory[
                "classification_subjects"].get("paired_3a_eligible", []),
            "participants_with_unique_f3_and_f4_labels": unique_f3_f4,
            "reference_warning": inventory["reference_warning"],
            "selection_was_blind_to_scalp_endpoint_values": inventory[
                "selection_was_blind_to_scalp_endpoint_values"],
        }
        _require_equal(pinned, expected_pin, f"{label} inventory pin")
        _require_equal(
            paired_subjects,
            pinned["paired_3a_eligible_subjects"],
            f"{label} paired subject intersection",
        )

    locked = _validate_tracked_lineage(
        root, inventory, inventory_manifest, result_manifest)
    _validate_frozen_3a_evidence(inventory, locked)

    cache_sha = _require_sha256(
        inventory.get("frozen_cache_manifest_sha256"),
        "frozen_cache_manifest_sha256",
    )
    sidecar_sha = _require_sha256(
        inventory.get("sidecar_manifest_sha256"),
        "sidecar_manifest_sha256",
    )
    _require_equal(
        locked.get("cache_manifest_sha256"),
        cache_sha,
        "locked-grid cache-manifest lineage",
    )
    for payload, label in (
        (result_manifest, "paired-result manifest"),
        (subject_results, "paired subject results"),
        (group_summary, "paired group summary"),
    ):
        _require_equal(
            payload.get("frozen_qc_grid_cache_manifest_sha256"),
            cache_sha,
            f"{label} cache-manifest lineage",
        )

    skipped = []
    for key, expected_sha in (
        ("frozen_cache_manifest_path_relative", cache_sha),
        ("sidecar_manifest_path_relative", sidecar_sha),
    ):
        relative = inventory.get(key)
        path = _repository_path(root, relative, key)
        if path.is_file():
            _require_equal(file_sha256(path), expected_sha, key)
        else:
            if not relative.startswith("data/derived/"):
                raise RuntimeError(
                    "only ignored data/derived lineage may be absent: "
                    f"{relative}")
            skipped.append(relative)

    return {
        "inventory": inventory,
        "inventory_manifest": inventory_manifest,
        "result_manifest": result_manifest,
        "subject_results": subject_results,
        "group_summary": group_summary,
        "skipped_derived_lineage": skipped,
    }
