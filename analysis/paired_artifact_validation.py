"""Strict clean-checkout validation for checked-in paired-EEG evidence.

Raw EEG and derived NPZ caches are intentionally absent from Git. A clean
checkout can still prove that every public inventory/result file is present,
byte-pinned, mutually consistent, and tied to the exact current source,
profile, schemas, and locked QC evidence. Publication mode separately requires
the ignored private cache manifests to be present.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from artifact_contracts import (
    PAIRED_RESULT_SCHEMA,
    QC_LOCKED_SNAPSHOT_SCHEMA,
    SCALP_INVENTORY_SCHEMA,
)
from paired_reporting import (
    _normalized_csv_rows,
    _role_pair_csv_rows,
    _validate_csv_artifact,
)
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    file_sha256,
    source_tree_sha256,
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
    if inventory.get("query_error_subjects"):
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
    return requested


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
