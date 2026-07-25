"""Validated, hashable QC profiles for LC-proxy sensitivity analyses.

Scientific method constants and data-availability/QC choices are deliberately separate.  The
coverage profiles in ``qc_profiles_v1.json`` are post-audit sensitivity specifications, not
thresholds supplied by the cited papers and not prospective preregistration.
"""
from __future__ import annotations

import copy
import hashlib
import json
import os


HERE = os.path.dirname(os.path.abspath(__file__))
PROFILE_PATH = os.path.join(HERE, "qc_profiles_v1.json")
PROFILE_SET_SCHEMA = "qc-profile-set-v1"
ROOT = os.path.dirname(HERE)

_FRACTION_PATHS = {
    "hr.minimum_coverage",
    "power.minimum_aggregate_coverage",
    "power.minimum_clean_fraction_per_second",
    "power.minimum_contact_coverage",
    "power.minimum_contact_fraction_per_bin",
    "staging.minimum_contact_feature_coverage",
    "staging.minimum_contact_fraction_per_epoch",
    "staging.minimum_auxiliary_channel_coverage",
    "staging.minimum_auxiliary_channel_fraction_per_epoch",
    "endpoint_3d.minimum_contact_acquisition_fraction",
    "endpoint_3d.minimum_contact_event_valid_fraction",
    "endpoint_3d.minimum_valid_nrem_fraction_per_contact",
}
_POSITIVE_INTEGER_PATHS = {
    "cohort.minimum_participants",
    "endpoint_3a.minimum_cross_correlation_windows",
    "endpoint_3a.minimum_nrem_epochs",
    "endpoint_3b.minimum_contacts",
    "endpoint_3b.minimum_finite_rr_samples",
    "endpoint_3b.minimum_finite_stage_samples",
    "endpoint_3b.minimum_so_per_contact",
    "endpoint_3d.minimum_contacts",
    "endpoint_3d.minimum_events_per_contact",
    "endpoint_3d.minimum_paired_events",
    "endpoint_3d.minimum_valid_nrem_seconds_per_contact",
    "power.minimum_contacts",
    "staging.minimum_contacts",
    "staging.minimum_valid_auxiliary_windows",
    "staging.minimum_valid_welch_windows",
}
_PERCENTILE_PATHS = {
    "endpoint_3b.so_amplitude_percentile",
}
_REQUIRED_SECTIONS = {
    "cohort",
    "endpoint_3a",
    "endpoint_3b",
    "endpoint_3d",
    "hr",
    "power",
    "staging",
}


def _canonical_bytes(value):
    return json.dumps(
        value, allow_nan=False, separators=(",", ":"), sort_keys=True
    ).encode("utf-8")


def qc_profile_sha256(profile):
    """Canonical SHA-256 for one fully expanded profile."""
    return hashlib.sha256(_canonical_bytes(profile)).hexdigest()


def profile_file_sha256(path=PROFILE_PATH):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validated_staging_calibration(path=PROFILE_PATH):
    """Verify the exact outcome-blind artifact used to motivate profile 11.

    The artifact is pinned separately from the profile hash because generated
    outputs are intentionally outside ``source_tree_sha256``.  This prevents a
    replaced calibration JSON from silently inheriting the word "calibrated."
    """
    with open(path) as handle:
        profile_set = json.load(handle)
    specification = (
        profile_set.get("method_config", {}).get("staging_calibration"))
    if not isinstance(specification, dict):
        raise RuntimeError("QC profile set lacks staging calibration provenance")
    required = {
        "analysis_version",
        "artifact_relative_path",
        "artifact_sha256",
        "cache_manifest_sha256",
        "cache_pipeline",
        "cache_run_id",
        "cache_schema_version",
        "calibration_scope",
        "hup_transport_status",
        "not_a_paper_requirement",
        "recommended_minimum_valid_windows",
    }
    missing = sorted(required - set(specification))
    if missing:
        raise RuntimeError(
            f"staging calibration provenance lacks fields: {missing}")
    artifact_path = os.path.join(
        ROOT, *specification["artifact_relative_path"].split("/"))
    if not os.path.isfile(artifact_path):
        raise RuntimeError(
            f"pinned staging calibration artifact is missing: {artifact_path}")
    artifact_hash = profile_file_sha256(artifact_path)
    if artifact_hash != specification["artifact_sha256"]:
        raise RuntimeError(
            "staging calibration artifact bytes differ from the QC profile pin")
    with open(artifact_path) as handle:
        artifact = json.load(handle)
    calibration_source_files = artifact.get(
        "calibration_source_files_sha256")
    required_calibration_sources = {
        "analysis/calibrate_staging_windows.py",
        "analysis/pipeline_version.py",
    }
    if (
        not isinstance(calibration_source_files, dict)
        or set(calibration_source_files) != required_calibration_sources
    ):
        raise RuntimeError(
            "staging calibration lacks exact purpose-specific source hashes")
    for relative, expected_hash in calibration_source_files.items():
        current_path = os.path.join(ROOT, *relative.split("/"))
        if (
            not os.path.isfile(current_path)
            or profile_file_sha256(current_path) != expected_hash
        ):
            raise RuntimeError(
                f"staging calibration source changed after calibration: {relative}")
    for key in (
        "analysis_version",
        "cache_schema_version",
        "cache_pipeline",
        "cache_run_id",
        "cache_manifest_sha256",
        "recommended_minimum_valid_windows",
    ):
        if artifact.get(key) != specification[key]:
            raise RuntimeError(
                f"staging calibration artifact {key} differs from its profile pin")
    if specification["not_a_paper_requirement"] is not True:
        raise RuntimeError(
            "staging calibration must be explicitly labeled non-paper-derived")
    recommended = int(specification["recommended_minimum_valid_windows"])
    exact_row = next(
        (
            row for row in artifact.get("support_results", [])
            if row.get("minimum_valid_windows") == recommended
        ),
        None,
    )
    if (
        exact_row is None
        or not exact_row.get("exact_support", {}).get(
            "meets_calibration_targets")
    ):
        raise RuntimeError(
            "pinned staging recommendation does not meet its exact-support targets")
    return {
        **copy.deepcopy(specification),
        "artifact_path": artifact_path,
        "cache_files_sha256": copy.deepcopy(
            artifact.get("cache_files_sha256", {})),
        "runtime_versions": copy.deepcopy(artifact.get("runtime_versions", {})),
        "calibration_source_files_sha256": copy.deepcopy(
            calibration_source_files),
        "n_calibration_records": artifact.get("n_calibration_records"),
        "selection_blinding": artifact.get("selection_blinding"),
    }


def load_profile_set(path=PROFILE_PATH):
    with open(path) as handle:
        payload = json.load(handle)
    if payload.get("schema_version") != PROFILE_SET_SCHEMA:
        raise ValueError(
            f"QC profile set has schema {payload.get('schema_version')!r}; "
            f"expected {PROFILE_SET_SCHEMA!r}"
        )
    if not isinstance(payload.get("profiles"), dict) or not payload["profiles"]:
        raise ValueError("QC profile set must contain at least one profile")
    if not isinstance(payload.get("grids"), dict):
        raise ValueError("QC profile set grids must be an object")
    for profile_id, profile in payload["profiles"].items():
        validate_qc_profile(profile, expected_id=profile_id)
    return payload


def _path_get(payload, path):
    value = payload
    for key in path.split("."):
        if not isinstance(value, dict) or key not in value:
            raise ValueError(f"QC profile is missing {path!r}")
        value = value[key]
    return value


def _path_set(payload, path, value):
    target = payload
    keys = path.split(".")
    for key in keys[:-1]:
        if not isinstance(target, dict) or key not in target:
            raise ValueError(f"QC grid path does not exist: {path!r}")
        target = target[key]
    if not isinstance(target, dict) or keys[-1] not in target:
        raise ValueError(f"QC grid path does not exist: {path!r}")
    target[keys[-1]] = value


def validate_qc_profile(profile, expected_id=None):
    if not isinstance(profile, dict):
        raise ValueError("QC profile must be an object")
    missing = sorted(_REQUIRED_SECTIONS - set(profile))
    if missing:
        raise ValueError(f"QC profile is missing sections: {missing}")
    profile_id = profile.get("profile_id")
    if not isinstance(profile_id, str) or not profile_id:
        raise ValueError("QC profile_id must be a nonempty string")
    if expected_id is not None and profile_id != expected_id:
        raise ValueError(
            f"QC profile key {expected_id!r} disagrees with profile_id {profile_id!r}"
        )
    if not isinstance(profile.get("role"), str) or not profile["role"]:
        raise ValueError("QC profile role must be a nonempty string")
    if not isinstance(profile.get("selection_blinded_to_endpoint_values"), bool):
        raise ValueError(
            "QC profile selection_blinded_to_endpoint_values must be boolean"
        )
    for path in sorted(_FRACTION_PATHS):
        value = _path_get(profile, path)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"{path} must be numeric")
        if not 0.0 <= float(value) <= 1.0:
            raise ValueError(f"{path} must lie in [0, 1]")
    for path in sorted(_POSITIVE_INTEGER_PATHS):
        value = _path_get(profile, path)
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{path} must be a positive integer")
    for path in sorted(_PERCENTILE_PATHS):
        value = _path_get(profile, path)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"{path} must be numeric")
        if not 0.0 <= float(value) <= 100.0:
            raise ValueError(f"{path} must lie in [0, 100]")
    if profile["power"].get("aggregation") not in {
        "fixed_contact_mean",
        "overlap_connected_median_polish",
    }:
        raise ValueError("unknown power aggregation")
    if profile["staging"].get("aggregation") not in {
        "fixed_contact_median",
        "overlap_connected_median_polish",
    }:
        raise ValueError("unknown staging aggregation")
    if profile["staging"].get("auxiliary_aggregation") not in {
        "fixed_channel_geometric_median",
        "overlap_connected_median_polish",
    }:
        raise ValueError("unknown staging auxiliary aggregation")
    return profile


def load_qc_profile(profile_id, path=PROFILE_PATH):
    payload = load_profile_set(path)
    try:
        profile = copy.deepcopy(payload["profiles"][profile_id])
    except KeyError as exc:
        raise KeyError(f"unknown QC profile {profile_id!r}") from exc
    return validate_qc_profile(profile, expected_id=profile_id)


def expand_qc_grid(grid_id, path=PROFILE_PATH):
    """Expand a deterministic one-factor-at-a-time grid from its validated base."""
    payload = load_profile_set(path)
    try:
        grid = payload["grids"][grid_id]
    except KeyError as exc:
        raise KeyError(f"unknown QC grid {grid_id!r}") from exc
    base = load_qc_profile(grid["base_profile"], path)
    axes = list(grid.get("axes", []))
    global_levels = list(grid.get("levels", []))
    if not axes:
        raise ValueError(f"QC grid {grid_id!r} requires axes")
    expanded = [copy.deepcopy(base)]
    for axis in axes:
        axis_id = axis.get("id")
        axis_path = axis.get("path")
        if not isinstance(axis_id, str) or not axis_id:
            raise ValueError("QC grid axis id must be a nonempty string")
        levels = list(axis.get("levels", global_levels))
        if not levels:
            raise ValueError(f"QC grid axis {axis_id!r} requires levels")
        baseline = _path_get(base, axis_path)
        for level in levels:
            if level == baseline:
                continue
            profile = copy.deepcopy(base)
            _path_set(profile, axis_path, level)
            suffix = (
                str(int(round(level * 100)))
                if isinstance(level, float)
                else str(level)
            )
            profile["profile_id"] = f"{axis_id}_{suffix}"
            profile["role"] = "post_audit_sensitivity"
            profile["selection_blinded_to_endpoint_values"] = True
            profile["grid_id"] = grid_id
            profile["changed_axis"] = axis_id
            profile["changed_path"] = axis_path
            profile["changed_value"] = level
            validate_qc_profile(profile)
            expanded.append(profile)
    ids = [value["profile_id"] for value in expanded]
    if len(ids) != len(set(ids)):
        raise ValueError(f"QC grid {grid_id!r} produced duplicate profile ids")
    return expanded


def manifest_qc_config(profile, path=PROFILE_PATH):
    profile = copy.deepcopy(validate_qc_profile(profile))
    return {
        "qc_profile": profile,
        "qc_profile_id": profile["profile_id"],
        "qc_profile_sha256": qc_profile_sha256(profile),
        "qc_profile_file_sha256": profile_file_sha256(path),
        "qc_profile_set_schema": PROFILE_SET_SCHEMA,
    }
