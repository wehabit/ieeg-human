"""Validation contracts for paired-scalp terminal manifests and NPZ sidecars.

This module deliberately does not import :mod:`cache_paired_scalp`.  The cache
producer owns every constant that can affect emitted bytes and passes those
values in through :class:`SidecarValidationContract`.  Keeping that direction
of dependency makes validator maintenance independent of the producer digest
while avoiding a producer/validator import cycle.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta
import hashlib
import json
import os
import subprocess
import uuid
from typing import Callable, Iterable, Mapping

import numpy as np

from cache_lc_series import power_from_binned_support
from hup_portal import same_series_geometry
from pipeline_version import (
    file_sha256,
    npz_scalar_text,
    validate_completed_cache_failures,
    validate_full_interval_acquisition,
)


TERMINAL_SIDECAR_MANIFEST_FIELDS = frozenset({
    "run_id",
    "schema_version",
    "pipeline",
    "run_state",
    "analysis_version",
    "cache_schema_version",
    "generated_at_utc",
    "code_revision",
    "code_dirty",
    "runtime_versions",
    "cache_dependency_sha256",
    "requested",
    "completed",
    "reused",
    "failed",
    "channel_plan",
    "config",
    "result_files_sha256",
})


@dataclass(frozen=True)
class SidecarValidationContract:
    """Producer-owned values needed to validate one sidecar schema."""

    analysis_version: str
    ieeg_cache_schema: str
    scalp_cache_schema: str
    pipeline: str
    channel_plan: Mapping[str, Mapping[str, str]]
    filter_edge_s: float
    chunk_s: float
    historical_power_support: float
    sigma_band_hz: tuple[float, float]
    swa_band_hz: tuple[float, float]
    so_band_hz: tuple[float, float]


def _is_lower_hex(value, *, length):
    return (
        isinstance(value, str)
        and len(value) == length
        and all(character in "0123456789abcdef" for character in value)
    )


def _valid_manifest_run_id(value):
    if not isinstance(value, str) or not value:
        return False
    try:
        parsed = uuid.UUID(value)
    except (AttributeError, ValueError):
        return False
    return str(parsed) == value


def recorded_commit_exists(value, *, root):
    if not _is_lower_hex(value, length=40):
        return False
    try:
        completed = subprocess.run(
            ["git", "-C", root, "cat-file", "-e", f"{value}^{{commit}}"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    except OSError:
        return False
    return completed.returncode == 0


def dependency_sha256_at_revision(
        revision, *, root, dependency_files: Iterable[str],
        commit_exists: Callable[[str], bool] | None = None):
    """Hash producer dependency bytes stored by one recorded Git commit."""
    commit_exists = (
        (lambda value: recorded_commit_exists(value, root=root))
        if commit_exists is None
        else commit_exists
    )
    if not commit_exists(revision):
        return None
    digest = hashlib.sha256()
    for relative in sorted(dependency_files):
        try:
            completed = subprocess.run(
                ["git", "-C", root, "show", f"{revision}:{relative}"],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                check=False,
            )
        except OSError:
            return None
        if completed.returncode != 0:
            return None
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(completed.stdout)
        digest.update(b"\0")
    return digest.hexdigest()


def _valid_utc_timestamp(value):
    if not isinstance(value, str) or not value:
        return False
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return False
    return (
        parsed.tzinfo is not None
        and parsed.utcoffset() == timedelta(0)
    )


def _terminal_sidecar_partitions(manifest, *, source):
    requested = manifest["requested"]
    completed = manifest["completed"]
    reused = manifest["reused"]
    failed = manifest["failed"]
    if any(
        not isinstance(values, list)
        for values in (requested, completed, reused, failed)
    ):
        raise RuntimeError(f"{source} has malformed terminal partitions")
    if any(
        not isinstance(subject, str) or not subject
        for values in (requested, completed, reused)
        for subject in values
    ):
        raise RuntimeError(f"{source} has malformed terminal subject IDs")
    failed_subjects = [
        value.get("subject") if isinstance(value, dict) else None
        for value in failed
    ]
    if any(
        not isinstance(subject, str) or not subject
        for subject in failed_subjects
    ):
        raise RuntimeError(f"{source} has malformed failed subject IDs")
    if any(
        len(values) != len(set(values))
        for values in (requested, completed, reused, failed_subjects)
    ):
        raise RuntimeError(f"{source} terminal partitions contain duplicates")
    partitions = (set(completed), set(reused), set(failed_subjects))
    if (
        not requested
        or failed
        or any(value is None for value in failed_subjects)
        or any(
            partitions[left] & partitions[right]
            for left, right in ((0, 1), (0, 2), (1, 2))
        )
        or set.union(*partitions) != set(requested)
    ):
        raise RuntimeError(
            f"{source} does not contain one exact successful partition")
    return list(requested), list(completed), list(reused)


def validate_terminal_sidecar_manifest(
        manifest, *, source, output_dir, cache_dir,
        expected_requested=None, contract,
        dependency_sha256: Callable[[], str],
        manifest_config: Callable[[str], dict],
        current_runtime_versions: Callable[[], dict],
        commit_exists: Callable[[str], bool],
        dependency_sha256_for_revision: Callable[[str], str | None]):
    """Prove one exact, complete current paired-scalp terminal manifest."""
    if not isinstance(manifest, dict):
        raise RuntimeError(f"{source} must contain one JSON object")
    if set(manifest) != TERMINAL_SIDECAR_MANIFEST_FIELDS:
        missing = sorted(TERMINAL_SIDECAR_MANIFEST_FIELDS - set(manifest))
        extra = sorted(set(manifest) - TERMINAL_SIDECAR_MANIFEST_FIELDS)
        raise RuntimeError(
            f"{source} has a noncanonical terminal schema "
            f"(missing={missing}, extra={extra})")
    requested, completed, reused = _terminal_sidecar_partitions(
        manifest, source=source)
    dependency_digest = dependency_sha256()
    recorded_hashes = manifest["result_files_sha256"]
    expected_plan = {
        subject: contract.channel_plan[subject]
        for subject in requested
        if subject in contract.channel_plan
    }
    if (
        manifest["pipeline"] != contract.pipeline
        or manifest["run_state"] != "complete"
        or manifest["schema_version"] != contract.scalp_cache_schema
        or manifest["analysis_version"] != contract.analysis_version
        or manifest["cache_schema_version"] != contract.ieeg_cache_schema
        or not _valid_manifest_run_id(manifest["run_id"])
        or not _valid_utc_timestamp(manifest["generated_at_utc"])
        or not commit_exists(manifest["code_revision"])
        or dependency_sha256_for_revision(manifest["code_revision"])
        != dependency_digest
        or manifest["code_dirty"] is not False
        or manifest["runtime_versions"] != current_runtime_versions()
        or manifest["cache_dependency_sha256"] != dependency_digest
        or set(expected_plan) != set(requested)
        or manifest["channel_plan"] != expected_plan
        or manifest["config"] != manifest_config(cache_dir)
        or not isinstance(recorded_hashes, dict)
        or set(recorded_hashes) != set(requested)
        or any(
            not _is_lower_hex(digest, length=64)
            for digest in recorded_hashes.values()
        )
    ):
        raise RuntimeError(
            f"{source} is not an exact complete current paired-scalp run")
    if (
        expected_requested is not None
        and requested != list(expected_requested)
    ):
        raise RuntimeError(
            f"{source} requested partition differs from the required scope")
    output_dir = os.path.abspath(output_dir)
    for subject, expected_sha256 in recorded_hashes.items():
        path = os.path.join(output_dir, f"{subject}.npz")
        if not os.path.isfile(path) or file_sha256(path) != expected_sha256:
            raise RuntimeError(
                f"{source} result bytes differ for {subject}")
    return {
        "requested": list(requested),
        "completed": list(completed),
        "reused": list(reused),
        "result_files_sha256": dict(recorded_hashes),
        "cache_dependency_sha256": dependency_digest,
    }


def _required_array(cache, key, shape, *, path, dtype=float):
    if key not in cache.files:
        raise RuntimeError(f"{path} lacks required sidecar array {key!r}")
    try:
        value = np.asarray(cache[key], dtype=dtype)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(
            f"{path} has a non-numeric sidecar array {key!r}") from exc
    if value.shape != tuple(shape):
        raise RuntimeError(
            f"{path} sidecar array {key!r} has shape {value.shape}, "
            f"expected {tuple(shape)}")
    return value


def _validate_power_support(
        cache, *, path, n_channels, total_s, sf,
        historical_power_support):
    try:
        samples_per_second = int(np.asarray(
            cache["power_samples_per_second"]).item())
        minimum_support = float(np.asarray(
            cache[
                "power_historical_minimum_clean_fraction_per_second"
            ]).item())
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError(f"{path} has malformed power-support metadata") from exc
    if (
        samples_per_second != int(round(sf))
        or not np.isclose(
            minimum_support,
            historical_power_support,
            rtol=0,
            atol=1e-12,
        )
    ):
        raise RuntimeError(f"{path} has stale power-support metadata")

    shape = (n_channels, total_s)
    materialized = {}
    clean_counts = {}
    for prefix, values_key in (
        ("sigma_fixed", "sigma_fixed_by_channel"),
        ("swa", "swa_by_channel"),
    ):
        values = _required_array(
            cache, values_key, shape, path=path, dtype=float)
        numerator = _required_array(
            cache, f"{prefix}_power_numerator_by_channel",
            shape, path=path, dtype=float)
        denominator = _required_array(
            cache, f"{prefix}_clean_sample_count_by_channel",
            shape, path=path, dtype=float)
        if (
            not np.isfinite(numerator).all()
            or np.any(numerator < 0)
            or not np.isfinite(denominator).all()
            or np.any(denominator < 0)
            or np.any(denominator > samples_per_second)
            or not np.equal(denominator, np.floor(denominator)).all()
        ):
            raise RuntimeError(
                f"{path} has invalid {prefix} clean-sample counts")
        reconstructed = power_from_binned_support(
            numerator,
            denominator,
            samples_per_second,
            minimum_support,
        )
        if not np.array_equal(values, reconstructed, equal_nan=True):
            raise RuntimeError(
                f"{path} {values_key} does not match its reversible support")
        materialized[prefix] = values
        clean_counts[prefix] = denominator
    if not np.array_equal(
            clean_counts["sigma_fixed"], clean_counts["swa"]):
        raise RuntimeError(
            f"{path} sigma and SWA clean-sample support differs")
    materialized["clean_sample_count"] = clean_counts["sigma_fixed"]
    return materialized


def _validate_activity(cache, *, path, channels, total_samples):
    n_channels = len(channels)
    minimum = _required_array(
        cache, "scalp_signal_raw_minimum", (n_channels,),
        path=path, dtype=float)
    maximum = _required_array(
        cache, "scalp_signal_raw_maximum", (n_channels,),
        path=path, dtype=float)
    dynamic = _required_array(
        cache, "scalp_signal_raw_dynamic_range", (n_channels,),
        path=path, dtype=float)
    tolerance = _required_array(
        cache, "scalp_signal_numerical_flat_tolerance", (n_channels,),
        path=path, dtype=float)
    finite_count = _required_array(
        cache, "scalp_signal_finite_sample_count", (n_channels,),
        path=path, dtype=float)
    raw_nonflat = np.asarray(
        cache["scalp_signal_nonflat_mask"]
        if "scalp_signal_nonflat_mask" in cache.files
        else [])
    if raw_nonflat.shape != (n_channels,) or raw_nonflat.dtype.kind != "b":
        raise RuntimeError(f"{path} has invalid activity nonflat flags")
    nonflat = raw_nonflat.astype(bool, copy=False)
    if (
        not np.isfinite(tolerance).all()
        or np.any(tolerance < 0)
        or not np.isfinite(finite_count).all()
        or np.any(finite_count < 0)
        or np.any(finite_count > total_samples)
        or not np.equal(finite_count, np.floor(finite_count)).all()
    ):
        raise RuntimeError(f"{path} has invalid full-interval activity support")
    empty = finite_count == 0
    if (
        np.any(empty & ~(
            np.isposinf(minimum) & np.isneginf(maximum)
        ))
        or np.any(~empty & ~(
            np.isfinite(minimum)
            & np.isfinite(maximum)
            & (maximum >= minimum)
        ))
    ):
        raise RuntimeError(f"{path} has invalid full-interval activity extrema")
    expected_dynamic = maximum - minimum
    if not np.array_equal(dynamic, expected_dynamic, equal_nan=True):
        raise RuntimeError(
            f"{path} activity dynamic range does not match its extrema")
    maximum_absolute = np.where(
        empty, 0.0, np.maximum(np.abs(minimum), np.abs(maximum)))
    expected_tolerance = (
        64.0 * np.finfo(float).eps
        * np.maximum(1.0, maximum_absolute)
    )
    if not np.array_equal(tolerance, expected_tolerance):
        raise RuntimeError(
            f"{path} activity tolerances do not match their extrema")
    expected_nonflat = (
        (finite_count >= 2)
        & np.isfinite(dynamic)
        & (dynamic > tolerance)
    )
    if not np.array_equal(nonflat, expected_nonflat):
        raise RuntimeError(
            f"{path} activity nonflat flags do not match their evidence")
    return {
        channel: {
            "numerically_nonflat_full_interval": bool(nonflat[index]),
            "finite_sample_count": int(finite_count[index]),
            "raw_dynamic_range": float(dynamic[index]),
        }
        for index, channel in enumerate(channels)
    }


def _validate_so_candidates(cache, *, path, channels, total_s):
    for channel in channels:
        keys = [
            f"so_candidate_{field}_{channel}"
            for field in ("t", "down", "up", "p2p")
        ]
        missing = [key for key in keys if key not in cache.files]
        if missing:
            raise RuntimeError(
                f"{path} lacks slow-oscillation candidate arrays {missing}")
        values = [np.asarray(cache[key], float) for key in keys]
        if any(value.ndim != 1 for value in values) or len({
                value.shape for value in values}) != 1:
            raise RuntimeError(
                f"{path} slow-oscillation candidate arrays do not align "
                f"for {channel}")
        times, down, up, peak_to_peak = values
        if (
            not all(np.isfinite(value).all() for value in values)
            or np.any(times < 0)
            or np.any(times >= total_s)
            or np.any(np.diff(times) <= 0)
            or np.any(down <= 0)
            or np.any(up <= 0)
            or np.any(peak_to_peak <= 0)
            or not np.allclose(
                peak_to_peak, down + up, rtol=1e-12, atol=1e-12)
        ):
            raise RuntimeError(
                f"{path} has invalid slow-oscillation candidates for {channel}")


def validate_scalp_sidecar_payload(
        cache, *, path, subject, pinned, dependency_digest, contract,
        normalize_scalp_label: Callable[[str], str]):
    """Validate one exact result-producing scalp sidecar payload."""
    required_scalars = {
        "status": "ok",
        "cache_schema_version": contract.scalp_cache_schema,
        "cache_dependency_sha256": dependency_digest,
        "subject": subject,
        "source_dataset": subject,
        "source_kind": "iEEG.org API",
        "ieeg_cache_sha256": pinned["sha256"],
        "ieeg_cache_manifest_sha256": pinned["manifest_sha256"],
        "ieeg_cache_manifest_run_id": str(
            pinned["manifest_run_id"] or ""),
        "ieeg_cache_schema_version": contract.ieeg_cache_schema,
    }
    for key, expected in required_scalars.items():
        if npz_scalar_text(cache, key) != str(expected):
            raise RuntimeError(
                f"{path} has stale or invalid embedded field {key!r}")
    validate_completed_cache_failures(
        cache, source=path, require_ecg=False)
    acquisition = validate_full_interval_acquisition(
        cache,
        source=path,
        core_purpose="paired_scalp_core",
        subrequest_purpose="paired_scalp_subrequest",
        core_chunk_s=contract.chunk_s,
        filter_edge_s=contract.filter_edge_s,
    )
    for key, expected in (
        ("night_s", pinned["night_s"]),
        ("hours", pinned["hours"]),
        ("sf", pinned["sf"]),
        ("filter_edge_seconds", contract.filter_edge_s),
        ("chunk_seconds", contract.chunk_s),
    ):
        if key not in cache.files or not np.isclose(
                float(np.asarray(cache[key]).item()),
                float(expected), rtol=0, atol=1e-9):
            raise RuntimeError(f"{path} differs from the frozen {key}")
    try:
        roles = json.loads(npz_scalar_text(cache, "channel_roles_json"))
        source_identity = json.loads(
            npz_scalar_text(cache, "source_identity_json"))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"{path} has malformed channel/source identity metadata") from exc
    if not isinstance(roles, dict) or not isinstance(source_identity, dict):
        raise RuntimeError(
            f"{path} channel/source identity metadata must be JSON objects")
    if subject not in contract.channel_plan:
        raise RuntimeError(f"{subject} is outside the frozen scalp plan")
    plan = contract.channel_plan[subject]
    if acquisition["channel_count"] != len(plan):
        raise RuntimeError(
            f"{path} acquisition ledger channel count differs from its plan")
    if (
        set(roles) != set(plan)
        or any(
            normalize_scalp_label(roles[role])
            != normalize_scalp_label(plan[role])
            for role in plan
        )
        or source_identity.get("dataset_name") != subject
        or source_identity.get("snapshot_id")
        != pinned["source_identity"].get("snapshot_id")
        or source_identity.get("scalp_role_to_channel") != roles
    ):
        raise RuntimeError(f"{path} has stale scalp source/channel identity")
    channels = [
        str(value) for value in _required_array(
            cache, "scalp_chans", (len(plan),), path=path, dtype=str)
    ]
    if (
        channels != [roles[role] for role in plan]
        or len(set(channels)) != len(channels)
    ):
        raise RuntimeError(f"{path} scalp channel order/identity is invalid")
    source_channels = source_identity.get("scalp_channels")
    reference_label = source_identity.get("reference_ieeg_channel")
    reference_identity = source_identity.get("reference_ieeg_identity")
    pinned_channels = pinned["source_identity"].get("channels", {})
    identity_fields = {
        "revision_id",
        "data_check",
        "start_time_us",
        "end_time_us",
        "duration_us",
        "number_of_samples",
        "sample_rate_hz",
    }
    pinned_reference = pinned_channels.get(reference_label)
    if (
        not isinstance(source_channels, dict)
        or set(source_channels) != set(channels)
        or any(
            not isinstance(source_channels[channel], dict)
            or set(source_channels[channel]) != identity_fields
            for channel in channels
        )
        or not isinstance(pinned_reference, dict)
        or not isinstance(reference_identity, dict)
        or set(reference_identity) != identity_fields
        or reference_identity.get("revision_id")
        != pinned_reference.get("revision_id")
        or reference_identity.get("data_check")
        != pinned_reference.get("data_check")
        or not same_series_geometry(
            reference_identity, pinned_reference)
        or any(
            not same_series_geometry(
                source_channels[channel], reference_identity)
            for channel in channels
        )
    ):
        raise RuntimeError(
            f"{path} lacks exact scalp-channel source identities")
    for key, expected in (
        ("sigma_band_hz", contract.sigma_band_hz),
        ("swa_band_hz", contract.swa_band_hz),
        ("so_band_hz", contract.so_band_hz),
    ):
        value = _required_array(
            cache, key, (2,), path=path, dtype=float)
        if not np.array_equal(value, np.asarray(expected, float)):
            raise RuntimeError(f"{path} has stale band metadata {key!r}")
    for key in ("ecg_reused_not_redetected", "staging_reused_not_recomputed"):
        try:
            flag = np.asarray(cache[key])
        except (KeyError, ValueError) as exc:
            raise RuntimeError(f"{path} lacks required flag {key!r}") from exc
        if flag.shape != () or flag.dtype.kind != "b" or not bool(flag.item()):
            raise RuntimeError(f"{path} does not assert {key!r}")
    sf = float(pinned["sf"])
    total_s_float = float(pinned["hours"]) * 3600.0
    total_s = int(round(total_s_float))
    if not np.isclose(total_s_float, total_s, rtol=0, atol=1e-9):
        raise RuntimeError(f"{path} frozen duration is not whole-second binned")
    power = _validate_power_support(
        cache,
        path=path,
        n_channels=len(channels),
        total_s=total_s,
        sf=sf,
        historical_power_support=contract.historical_power_support,
    )
    activity = _validate_activity(
        cache,
        path=path,
        channels=channels,
        total_samples=total_s * int(round(sf)),
    )
    finite_count = np.asarray([
        activity[channel]["finite_sample_count"]
        for channel in channels
    ])
    if np.any(np.sum(
            power["clean_sample_count"], axis=1) > finite_count):
        raise RuntimeError(
            f"{path} power support exceeds raw finite-sample support")
    _validate_so_candidates(
        cache, path=path, channels=channels, total_s=total_s)
    return {
        "roles": roles,
        "channels": channels,
        "activity_by_channel": activity,
        "power": power,
    }


def validate_reusable_sidecar_run(
        output_dir, requested, *, cache_dir, frozen_inputs,
        freeze_pinned_inputs: Callable,
        assert_pinned_inputs_unchanged: Callable,
        validate_terminal_manifest: Callable,
        validate_sidecar_payload: Callable):
    """Prove a complete prior sidecar run before reusing existing bytes."""
    output_dir = os.path.abspath(output_dir)
    cache_dir = os.path.abspath(cache_dir)
    requested = list(requested)
    if len(requested) != len(set(requested)):
        raise RuntimeError("paired-scalp reuse partition contains duplicates")
    if frozen_inputs is None:
        frozen_inputs = freeze_pinned_inputs(requested, cache_dir)
    if (
        frozen_inputs.get("cache_dir") != cache_dir
        or set(frozen_inputs.get("subjects", {})) != set(requested)
    ):
        raise RuntimeError(
            "frozen iEEG inputs do not exactly match paired-scalp reuse scope")
    assert_pinned_inputs_unchanged(frozen_inputs, requested)
    manifest_path = os.path.join(output_dir, "RUN_MANIFEST.json")
    if not os.path.isfile(manifest_path):
        raise RuntimeError(
            "existing paired-scalp outputs have no terminal manifest; "
            "rerun with --force")
    try:
        with open(manifest_path) as handle:
            manifest = json.load(handle)
    except Exception as exc:
        raise RuntimeError(
            "prior paired-scalp manifest cannot be read; rerun with --force"
        ) from exc

    manifest_pin = validate_terminal_manifest(
        manifest,
        source=manifest_path,
        output_dir=output_dir,
        cache_dir=cache_dir,
        expected_requested=requested,
    )
    dependency_digest = manifest_pin["cache_dependency_sha256"]
    recorded_hashes = manifest_pin["result_files_sha256"]

    validated = {}
    for subject in requested:
        path = os.path.join(output_dir, f"{subject}.npz")
        if (
            not os.path.isfile(path)
            or file_sha256(path) != recorded_hashes[subject]
        ):
            raise RuntimeError(
                f"{subject} bytes differ from the prior terminal manifest; "
                "rerun with --force")
        pinned = frozen_inputs["subjects"][subject]
        with np.load(path, allow_pickle=False) as cache:
            validate_sidecar_payload(
                cache,
                path=path,
                subject=subject,
                pinned=pinned,
                dependency_digest=dependency_digest,
            )
        validated[subject] = {
            "path": path,
            "sha256": recorded_hashes[subject],
        }
    assert_pinned_inputs_unchanged(frozen_inputs, requested)
    return {
        "manifest_path": manifest_path,
        "manifest_sha256": file_sha256(manifest_path),
        "requested": requested,
        "subjects": validated,
    }
