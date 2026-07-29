"""Validate HUP138 bilateral scalp-only artifacts offline or for publication."""
from __future__ import annotations

import argparse
import json
import os
import re
import tempfile
from pathlib import Path

import numpy as np

from artifact_contracts import SCALP_F3_F4_EXPLORATORY_SCHEMA
from event_3b_estimators import FS_RR, HALF_WIN
from materialize_qc_cache import materialize
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    file_sha256,
    npz_scalar_text,
    runtime_versions,
    source_tree_sha256,
)
from qc_profiles import (
    load_qc_profile,
    profile_file_sha256,
    qc_profile_sha256,
)
from scalp_f3_f4_exploratory import (
    CSV_FIELDS,
    DEFAULT_OUTPUT,
    INPUT_LINEAGE_FIELDS,
    MANIFEST_FIELDS,
    N_SURROGATES,
    PIPELINE,
    PROFILE_ID,
    RESULT_CSV,
    RESULT_FILES,
    RESULT_JSON,
    RESULT_PNG,
    RESULT_SVG,
    RESULT_FIELDS,
    RESULT_SCHEMA,
    ROLES,
    STAGES,
    SUBJECT,
    _analyse_bilateral,
    _array_sha256,
    _canonical_json,
    _csv_rows,
    _make_figure,
    _sha256_bytes,
    _stage_interpretation,
    _write_csv,
)
from staging_helpers import EPOCH


ROOT = Path(__file__).resolve().parent.parent
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_RESULT_VALUE_FIELDS = frozenset({
    "subject",
    "analysis_scope",
    "exploratory",
    "roles",
    "channels_by_role",
    "sensor_nonflat_by_role",
    "tachogram_domain",
    "minimum_channels",
    "minimum_events_per_channel",
    "stages",
})
_STAGE_VALUE_FIELDS = frozenset({
    "stage",
    "stage_interpretation",
    "pooled_nrem_exploratory",
    "raw_epochs",
    "stable_epochs",
    "stable_seconds",
    "finite_stage_rr_samples",
    "eligible_so_count_by_role",
    "so_amplitude_percentile",
    "so_amplitude_rule_provenance",
    "available",
    "support_reasons",
    "estimate",
})
_ESTIMATE_FIELDS = frozenset({
    "n_channels",
    "event_contact_ids",
    "n_so_total",
    "pct_above_stage_mean",
    "peak_lag_s",
    "channel_peak_lag_s",
    "participant_curve_peak_lag_s",
    "event_locked_local_change_pct",
    "event_curve_peak_to_peak_pct",
    "local_pre_event_mean_hr",
    "z",
    "p_upper",
    "stage_shift_z_diagnostic",
    "stage_shift_p_upper_diagnostic",
    "null_mean_pct",
    "null_sd_pct",
    "n_surrogates",
    "n_surrogate_pool",
    "null_method",
    "inference_status",
    "tachogram_domain",
    "curve",
    "rr_curve",
    "lag_s",
    "stage_mean_hr",
    "channel_peak_lag_by_role_s",
})


def _reject_json_constant(value):
    raise ValueError(f"non-standard JSON constant {value!r}")


def _load_json(path, label):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"required {label} is missing: {path}")
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle, parse_constant=_reject_json_constant)
    if not isinstance(value, dict):
        raise RuntimeError(f"{label} must contain one JSON object")
    return value


def _require_sha256(value, label):
    if not isinstance(value, str) or not _SHA256.fullmatch(value):
        raise RuntimeError(f"{label} is not a lowercase SHA-256 digest")
    return value


def _repository_path(root, relative, label):
    root = Path(root).resolve()
    if (
        not isinstance(relative, str)
        or not relative
        or Path(relative).is_absolute()
    ):
        raise RuntimeError(f"{label} is not repository-relative")
    path = (root / relative).resolve()
    try:
        path.relative_to(root)
    except ValueError as error:
        raise RuntimeError(f"{label} escapes the repository") from error
    return path


def _require_header(payload, label, *, manifest=False):
    expected = {
        "schema_version": SCALP_F3_F4_EXPLORATORY_SCHEMA,
        "pipeline": PIPELINE,
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "subject": SUBJECT,
        "analysis_scope": "participant_level_scalp_only",
        "paired_comparison_performed": False,
        "group_inference_performed": False,
    }
    if manifest:
        expected["run_state"] = "complete"
    for key, value in expected.items():
        if payload.get(key) != value:
            raise RuntimeError(f"{label} differs at {key}")


def _is_nonnegative_int(value):
    return (
        isinstance(value, int)
        and not isinstance(value, bool)
        and value >= 0
    )


def _finite_float(value, label, *, positive=False):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not np.isfinite(float(value))
        or (positive and float(value) <= 0)
    ):
        raise RuntimeError(f"{label} must be a finite numeric value")
    return float(value)


def _require_close(actual, expected, label):
    actual = _finite_float(actual, label)
    expected = _finite_float(expected, f"expected {label}")
    if not np.isclose(actual, expected, rtol=0, atol=1e-12):
        raise RuntimeError(f"{label} is not derived from the saved curve")


def _support_reason_contract(
        stage, value, *, nonflat, endpoint, counts):
    reasons = value["support_reasons"]
    if (
        not isinstance(reasons, list)
        or len(reasons) != len(set(reasons))
        or any(not isinstance(reason, str) or not reason for reason in reasons)
    ):
        raise RuntimeError(f"{stage} support reasons are malformed")
    expected = []
    for role in ROLES:
        reason = f"{role.upper()} is numerically flat"
        if not nonflat[role]:
            expected.append(reason)
        elif reason in reasons:
            raise RuntimeError(f"{stage} has a spurious {role} flat reason")
    if value["stable_epochs"] < 6:
        expected.append("no uninterrupted 180-s stage run")
    elif "no uninterrupted 180-s stage run" in reasons:
        raise RuntimeError(f"{stage} has a spurious stable-stage reason")
    stage_reason = (
        f"{value['finite_stage_rr_samples']} finite stage samples < "
        f"{endpoint['minimum_finite_stage_samples']}"
    )
    if (
        value["finite_stage_rr_samples"]
        < endpoint["minimum_finite_stage_samples"]
    ):
        expected.append(stage_reason)
    elif any(
        reason.endswith(
            f"finite stage samples < "
            f"{endpoint['minimum_finite_stage_samples']}")
        for reason in reasons
    ):
        raise RuntimeError(f"{stage} has a spurious finite-stage reason")
    for role in ROLES:
        reason = (
            f"{role.upper()} {counts[role]} eligible SOs < "
            f"{endpoint['minimum_so_per_contact']}"
        )
        if counts[role] < endpoint["minimum_so_per_contact"]:
            expected.append(reason)
        elif any(
            item.startswith(f"{role.upper()} ")
            and item.endswith(
                f"eligible SOs < {endpoint['minimum_so_per_contact']}")
            for item in reasons
        ):
            raise RuntimeError(f"{stage} has a spurious {role} event reason")
    missing = [reason for reason in expected if reason not in reasons]
    if missing:
        raise RuntimeError(
            f"{stage} support reasons omit locked gate failures {missing}")

    allowed_unreconstructable = {
        "cardiac coverage below locked profile",
        (
            "complete in-stage RR windows leave insufficient "
            "two-channel SO support"
        ),
    }
    complete_window_reason = (
        "complete in-stage RR windows leave insufficient "
        "two-channel SO support"
    )
    if complete_window_reason in reasons and reasons != [
            complete_window_reason]:
        raise RuntimeError(
            f"{stage} complete-window failure coexists with earlier gates")
    for reason in reasons:
        if reason in expected or reason in allowed_unreconstructable:
            continue
        match = re.fullmatch(r"(\d+) finite RR samples < (\d+)", reason)
        if (
            match is None
            or int(match.group(2))
            != int(endpoint["minimum_finite_rr_samples"])
            or int(match.group(1))
            >= int(endpoint["minimum_finite_rr_samples"])
        ):
            raise RuntimeError(f"{stage} has an unrecognized support reason")
    return reasons


def _validate_estimate(
        stage, estimate, stage_value, *, channels, endpoint, total_rr):
    if not isinstance(estimate, dict) or set(estimate) != _ESTIMATE_FIELDS:
        raise RuntimeError(f"{stage} estimate schema is nonexact")
    if (
        estimate.get("n_channels") != endpoint["minimum_contacts"]
        or estimate.get("event_contact_ids")
        != [channels[role] for role in ROLES]
        or len(set(estimate["event_contact_ids"])) != len(ROLES)
        or estimate.get("tachogram_domain") != "rr"
        or estimate.get("z") is not None
        or estimate.get("p_upper") is not None
        or estimate.get("null_method")
        != "one shared circular shift in eligible stage-time across all channels"
        or not isinstance(estimate.get("inference_status"), str)
        or not estimate["inference_status"].startswith("disabled:")
    ):
        raise RuntimeError(
            f"{stage} estimate violates the joint descriptive identity")

    counts = stage_value["eligible_so_count_by_role"]
    n_so = estimate.get("n_so_total")
    if (
        not _is_nonnegative_int(n_so)
        or n_so < len(ROLES) * endpoint["minimum_so_per_contact"]
        or n_so > sum(counts.values())
    ):
        raise RuntimeError(f"{stage} event counts do not reconcile")
    n_pool = estimate.get("n_surrogate_pool")
    if (
        not _is_nonnegative_int(n_pool)
        or n_pool < endpoint["minimum_finite_stage_samples"]
        or n_pool > stage_value["finite_stage_rr_samples"]
        or n_pool > total_rr
        or estimate.get("n_surrogates") != N_SURROGATES
    ):
        raise RuntimeError(f"{stage} surrogate support violates locked minima")
    if n_so > len(ROLES) * n_pool:
        raise RuntimeError(
            f"{stage} retained events exceed two-channel stage support")

    expected_n = int(2 * HALF_WIN * FS_RR)
    expected_lag = (
        np.arange(expected_n, dtype=float)
        - int(HALF_WIN * FS_RR)
    ) / FS_RR
    lag = np.asarray(estimate.get("lag_s"), dtype=float)
    curve = np.asarray(estimate.get("curve"), dtype=float)
    rr_curve = np.asarray(estimate.get("rr_curve"), dtype=float)
    if (
        lag.shape != (expected_n,)
        or curve.shape != (expected_n,)
        or rr_curve.shape != (expected_n,)
        or not np.isfinite(lag).all()
        or not np.isfinite(curve).all()
        or not np.isfinite(rr_curve).all()
        or np.any(curve <= 0)
        or np.any(rr_curve <= 0)
        or not np.array_equal(lag, expected_lag)
        or not np.all(np.diff(lag) > 0)
        or not np.allclose(curve, 60.0 / rr_curve, rtol=0, atol=1e-12)
    ):
        raise RuntimeError(f"{stage} has an invalid fixed RR/HR curve grid")
    post = lag >= 0
    pre = lag < 0
    post_grid = lag[post]
    channel_lags = estimate.get("channel_peak_lag_s")
    by_role = estimate.get("channel_peak_lag_by_role_s")
    if (
        not isinstance(channel_lags, list)
        or len(channel_lags) != len(ROLES)
        or not isinstance(by_role, dict)
        or list(by_role) != list(ROLES)
    ):
        raise RuntimeError(f"{stage} channel peak identities are malformed")
    channel_lags = np.asarray([
        _finite_float(value, f"{stage} channel peak lag")
        for value in channel_lags
    ])
    if any(
        not np.any(np.isclose(post_grid, value, rtol=0, atol=1e-12))
        for value in channel_lags
    ):
        raise RuntimeError(f"{stage} channel peaks are not on the post-SO grid")
    for role, lag_value in zip(ROLES, channel_lags):
        _require_close(
            by_role.get(role), lag_value, f"{stage} {role} peak lag")
    _require_close(
        estimate.get("peak_lag_s"),
        float(np.mean(channel_lags)),
        f"{stage} bilateral peak lag",
    )

    post_relative_index = int(np.argmin(rr_curve[post]))
    participant_lag = float(post_grid[post_relative_index])
    _require_close(
        estimate.get("participant_curve_peak_lag_s"),
        participant_lag,
        f"{stage} participant-curve peak lag",
    )
    stage_mean_hr = _finite_float(
        estimate.get("stage_mean_hr"),
        f"{stage} stage mean HR",
        positive=True,
    )
    peak_hr = float(curve[post][post_relative_index])
    expected_pct = 100.0 * (
        peak_hr - stage_mean_hr) / stage_mean_hr
    _require_close(
        estimate.get("pct_above_stage_mean"),
        expected_pct,
        f"{stage} percent above stage mean",
    )
    local_baseline = float(np.mean(curve[pre]))
    _require_close(
        estimate.get("local_pre_event_mean_hr"),
        local_baseline,
        f"{stage} local pre-event mean HR",
    )
    local_change = 100.0 * (
        float(np.max(curve[post])) - local_baseline) / local_baseline
    _require_close(
        estimate.get("event_locked_local_change_pct"),
        local_change,
        f"{stage} local HR change",
    )
    peak_to_peak = 100.0 * (
        float(np.max(curve)) - float(np.min(curve))
    ) / float(np.mean(curve))
    _require_close(
        estimate.get("event_curve_peak_to_peak_pct"),
        peak_to_peak,
        f"{stage} curve peak-to-peak magnitude",
    )

    null_mean = _finite_float(
        estimate.get("null_mean_pct"), f"{stage} null mean")
    null_sd = _finite_float(
        estimate.get("null_sd_pct"), f"{stage} null SD")
    if null_sd < 0:
        raise RuntimeError(f"{stage} null SD is negative")
    expected_z = (
        expected_pct - null_mean) / (null_sd + 1e-12)
    _require_close(
        estimate.get("stage_shift_z_diagnostic"),
        expected_z,
        f"{stage} diagnostic z",
    )
    diagnostic_p = _finite_float(
        estimate.get("stage_shift_p_upper_diagnostic"),
        f"{stage} diagnostic p",
    )
    if (
        diagnostic_p <= 0
        or diagnostic_p > 1
        or not np.isclose(
            diagnostic_p * (N_SURROGATES + 1),
            round(diagnostic_p * (N_SURROGATES + 1)),
            rtol=0,
            atol=1e-12,
        )
    ):
        raise RuntimeError(f"{stage} diagnostic p is not a valid rank p-value")


def _validate_result_contract(payload):
    if set(payload) != RESULT_FIELDS:
        missing = sorted(RESULT_FIELDS - set(payload))
        extra = sorted(set(payload) - RESULT_FIELDS)
        raise RuntimeError(
            "bilateral result schema is nonexact "
            f"(missing={missing}, extra={extra})")
    _require_header(payload, "bilateral result")
    if (
        payload.get("event_locking_inference_available") is not False
        or payload.get("exploratory") is not True
        or "group_summary" in payload
        or "paired_results" in payload
    ):
        raise RuntimeError("bilateral result exceeds its descriptive scope")
    profile = load_qc_profile(PROFILE_ID)
    endpoint = profile["endpoint_3b"]
    result = payload.get("result")
    if (
        not isinstance(result, dict)
        or set(result) != _RESULT_VALUE_FIELDS
        or result.get("subject") != SUBJECT
        or result.get("analysis_scope") != "participant_level_scalp_only"
        or result.get("exploratory") is not True
        or result.get("roles") != list(ROLES)
        or result.get("minimum_channels") != endpoint["minimum_contacts"]
        or result.get("minimum_events_per_channel")
        != endpoint["minimum_so_per_contact"]
        or result.get("tachogram_domain") != "rr"
        or result.get("channels_by_role") != {"f3": "F3", "f4": "F4"}
        or not isinstance(result.get("sensor_nonflat_by_role"), dict)
        or list(result["sensor_nonflat_by_role"]) != list(ROLES)
        or any(
            type(result["sensor_nonflat_by_role"][role]) is not bool
            for role in ROLES
        )
        or set(result.get("stages", {})) != set(STAGES)
    ):
        raise RuntimeError("bilateral result has malformed method geometry")
    method = payload.get("method")
    if (
        not isinstance(method, dict)
        or method.get("estimator")
        != "event_3b_estimators.subject_so_triggered"
        or method.get("tachogram_domain") != "rr"
        or method.get("channels") != ["F3", "F4"]
        or method.get("minimum_channels") != endpoint["minimum_contacts"]
        or method.get("n_surrogates_for_diagnostic_only") != N_SURROGATES
        or method.get("event_locking_z_p") is not None
    ):
        raise RuntimeError("bilateral method identity is malformed")
    locked_context = payload.get("locked_ieeg_context")
    if (
        not isinstance(locked_context, dict)
        or set(locked_context) != {
            "available_3b_stages",
            "used_as_comparison_arm",
            "note",
        }
        or not isinstance(
            locked_context.get("available_3b_stages"), list)
        or len(locked_context["available_3b_stages"])
        != len(set(locked_context["available_3b_stages"]))
        or any(
            stage not in STAGES
            for stage in locked_context["available_3b_stages"]
        )
        or locked_context.get("used_as_comparison_arm") is not False
        or not isinstance(locked_context.get("note"), str)
        or not locked_context["note"]
    ):
        raise RuntimeError("locked iEEG context exceeds scalp-only scope")
    total_rr = payload.get("input_lineage", {}).get("rr_4_n_samples")
    total_epochs = payload.get("input_lineage", {}).get(
        "stage_lab_n_epochs")
    if (
        not _is_nonnegative_int(total_rr)
        or not _is_nonnegative_int(total_epochs)
    ):
        raise RuntimeError("input RR/stage counts are malformed")
    nonflat = result["sensor_nonflat_by_role"]
    for stage in STAGES:
        value = result["stages"][stage]
        if not isinstance(value, dict) or set(value) != _STAGE_VALUE_FIELDS:
            raise RuntimeError(f"{stage} stage schema is nonexact")
        estimate = value.get("estimate")
        available = value.get("available")
        if (
            value.get("stage") != stage
            or value.get("stage_interpretation")
            != _stage_interpretation(stage)
            or available is not (estimate is not None)
        ):
            raise RuntimeError(f"{stage} availability differs from estimate")
        if (
            value.get("pooled_nrem_exploratory")
            is not (stage == "NREM")
        ):
            raise RuntimeError(f"{stage} pooled-NREM flag is invalid")
        raw_epochs = value.get("raw_epochs")
        stable_epochs = value.get("stable_epochs")
        stable_seconds = value.get("stable_seconds")
        finite_stage = value.get("finite_stage_rr_samples")
        counts = value.get("eligible_so_count_by_role")
        if (
            not _is_nonnegative_int(raw_epochs)
            or not _is_nonnegative_int(stable_epochs)
            or stable_epochs > raw_epochs
            or raw_epochs > total_epochs
            or not _is_nonnegative_int(stable_seconds)
            or stable_seconds != stable_epochs * EPOCH
            or not _is_nonnegative_int(finite_stage)
            or finite_stage > stable_seconds * FS_RR
            or finite_stage > total_rr
            or not isinstance(counts, dict)
            or list(counts) != list(ROLES)
            or any(
                not _is_nonnegative_int(counts[role]) for role in ROLES)
            or value.get("so_amplitude_percentile")
            != float(endpoint["so_amplitude_percentile"])
            or value.get("so_amplitude_rule_provenance")
            != (
                "within-sensor/stage percentile adaptation; not Naji's "
                "fixed referenced-scalp voltage criteria"
            )
        ):
            raise RuntimeError(f"{stage} support/count geometry is malformed")
        reasons = _support_reason_contract(
            stage,
            value,
            nonflat=nonflat,
            endpoint=endpoint,
            counts=counts,
        )
        visible_gate_failed = (
            not all(nonflat.values())
            or stable_epochs < 6
            or finite_stage < endpoint["minimum_finite_stage_samples"]
            or any(
                counts[role] < endpoint["minimum_so_per_contact"]
                for role in ROLES
            )
        )
        if estimate is not None and (reasons or visible_gate_failed):
            raise RuntimeError(
                f"{stage} estimate bypasses locked support gates")
        if estimate is None:
            if not reasons:
                raise RuntimeError(f"{stage} unavailable without a reason")
            continue
        _validate_estimate(
            stage,
            estimate,
            value,
            channels=result["channels_by_role"],
            endpoint=endpoint,
            total_rr=total_rr,
        )
    if (
        result["stages"]["N2"]["raw_epochs"]
        + result["stages"]["N3"]["raw_epochs"]
        > total_epochs
        or result["stages"]["NREM"]["raw_epochs"]
        < (
            result["stages"]["N2"]["raw_epochs"]
            + result["stages"]["N3"]["raw_epochs"]
        )
    ):
        raise RuntimeError("stage epoch counts do not reconcile")
    return payload


def _validate_current_provenance(payload, root, *, require_current_runtime):
    provenance = payload.get("provenance")
    method = payload.get("method")
    profile = load_qc_profile(PROFILE_ID)
    if not isinstance(provenance, dict):
        raise RuntimeError("bilateral result lacks source provenance")
    if provenance.get("source_tree_sha256") != source_tree_sha256(str(root)):
        raise RuntimeError(
            "bilateral result provenance is stale at source_tree_sha256")
    recorded_runtime = provenance.get("runtime_versions")
    current_runtime = runtime_versions()
    if (
        not isinstance(recorded_runtime, dict)
        or set(recorded_runtime) != set(current_runtime)
        or any(
            not isinstance(key, str)
            or not key
            or not isinstance(value, str)
            or not value
            for key, value in recorded_runtime.items()
        )
    ):
        raise RuntimeError(
            "bilateral result has invalid recorded runtime_versions")
    if require_current_runtime and recorded_runtime != current_runtime:
        raise RuntimeError(
            "publication runtime_versions differ from the producing runtime")
    if (
        not isinstance(method, dict)
        or method.get("source_profile_id") != PROFILE_ID
        or method.get("source_profile_sha256")
        != qc_profile_sha256(profile)
        or method.get("source_profile_file_sha256")
        != profile_file_sha256()
    ):
        raise RuntimeError("bilateral result profile provenance is stale")


def _validate_public_pins(payload, root):
    lineage = payload.get("input_lineage")
    if (
        not isinstance(lineage, dict)
        or set(lineage) != INPUT_LINEAGE_FIELDS
    ):
        raise RuntimeError("bilateral result has nonexact input lineage")
    for prefix in ("inventory", "inventory_manifest", "qc_grid"):
        path = _repository_path(
            root,
            lineage.get(f"{prefix}_path_relative"),
            f"{prefix} path",
        )
        expected = _require_sha256(
            lineage.get(f"{prefix}_sha256"),
            f"{prefix}_sha256",
        )
        if not path.is_file() or file_sha256(path) != expected:
            raise RuntimeError(f"{prefix} bytes differ from the result pin")
    return lineage


def _validate_csv_bytes(payload, result_dir):
    expected_rows = _csv_rows(payload)
    if len(expected_rows) != len(STAGES):
        raise RuntimeError("bilateral CSV row count is not stage-complete")
    with tempfile.TemporaryDirectory() as directory:
        expected_path = os.path.join(directory, RESULT_CSV)
        _write_csv(expected_path, expected_rows)
        with open(expected_path, "rb") as handle:
            expected = handle.read()
    with open(Path(result_dir) / RESULT_CSV, "rb") as handle:
        actual = handle.read()
    if actual != expected:
        raise RuntimeError(
            f"{RESULT_CSV} does not match the canonical {CSV_FIELDS!r} export")


def _validate_figure_bytes(payload, result_dir):
    """Re-render figures only on the publication runtime and compare bytes."""
    with tempfile.TemporaryDirectory() as directory:
        expected_png = os.path.join(directory, RESULT_PNG)
        expected_svg = os.path.join(directory, RESULT_SVG)
        _make_figure(payload, expected_png, expected_svg)
        for name, expected_path in (
            (RESULT_PNG, expected_png),
            (RESULT_SVG, expected_svg),
        ):
            with open(expected_path, "rb") as handle:
                expected = handle.read()
            with open(Path(result_dir) / name, "rb") as handle:
                actual = handle.read()
            if actual != expected:
                raise RuntimeError(
                    f"{name} does not match the canonical result rendering")


def _validate_private_lineage(payload, root):
    from cache_paired_scalp import validate_pinned_ieeg_cache
    from paired_scalp_ieeg_comparison import _validate_scalp_cache

    lineage = payload["input_lineage"]
    for key in (
        "ieeg_cache_sha256",
        "ieeg_cache_manifest_sha256",
        "scalp_sidecar_sha256",
        "scalp_sidecar_manifest_sha256",
        "rr_4_sha256",
        "stage_lab_sha256",
        "channel_roles_json_sha256",
        "source_identity_json_sha256",
    ):
        _require_sha256(lineage.get(key), key)

    ieeg_path = _repository_path(
        root, lineage["ieeg_cache_path_relative"], "iEEG cache path")
    ieeg_dir = ieeg_path.parent
    sidecar_path = _repository_path(
        root, lineage["scalp_sidecar_path_relative"], "sidecar path")
    sidecar_dir = sidecar_path.parent
    for path, key in (
        (ieeg_path, "ieeg_cache_sha256"),
        (
            _repository_path(
                root,
                lineage["ieeg_cache_manifest_path_relative"],
                "iEEG manifest path",
            ),
            "ieeg_cache_manifest_sha256",
        ),
        (sidecar_path, "scalp_sidecar_sha256"),
        (
            _repository_path(
                root,
                lineage["scalp_sidecar_manifest_path_relative"],
                "sidecar manifest path",
            ),
            "scalp_sidecar_manifest_sha256",
        ),
    ):
        if not path.is_file() or file_sha256(path) != lineage[key]:
            raise RuntimeError(f"publication lineage differs at {key}")

    pinned = validate_pinned_ieeg_cache(SUBJECT, str(ieeg_dir))
    sidecar = _validate_scalp_cache(SUBJECT, str(sidecar_dir), pinned)
    if (
        sidecar["sha256"] != lineage["scalp_sidecar_sha256"]
        or sidecar["manifest_sha256"]
        != lineage["scalp_sidecar_manifest_sha256"]
        or sidecar["roles"] != lineage["validated_sidecar_roles"]
        or set(sidecar["roles"]) < set(ROLES)
    ):
        raise RuntimeError("validated sidecar differs from recorded lineage")

    profile = load_qc_profile(PROFILE_ID)
    with (
        np.load(pinned["path"], allow_pickle=False) as ieeg,
        np.load(sidecar["path"], allow_pickle=False) as scalp,
    ):
        base = materialize(ieeg, profile)
        source_identity_text = npz_scalar_text(
            scalp, "source_identity_json")
        if (
            _array_sha256(base["rr_4"]) != lineage["rr_4_sha256"]
            or _array_sha256(base["stage_lab"])
            != lineage["stage_lab_sha256"]
            or _sha256_bytes(
                npz_scalar_text(
                    scalp, "channel_roles_json").encode("utf-8")
            ) != lineage["channel_roles_json_sha256"]
            or _sha256_bytes(
                source_identity_text.encode("utf-8")
            ) != lineage["source_identity_json_sha256"]
            or json.loads(source_identity_text)
            != lineage.get("portal_source_identity")
        ):
            raise RuntimeError(
                "RR/stage/source arrays differ from recorded hashes")
        recomputed = _analyse_bilateral(scalp, base, profile)
    if _canonical_json(recomputed) != _canonical_json(payload["result"]):
        raise RuntimeError(
            "publication result differs from the validated private inputs")


def validate_exploratory_artifacts(
        root=ROOT, output_dir=DEFAULT_OUTPUT, *,
        mode="offline", allow_absent=True):
    """Validate public bytes; publication mode also recomputes private lineage."""
    if mode not in ("offline", "publication"):
        raise ValueError("mode must be 'offline' or 'publication'")
    root = Path(root).resolve()
    result_dir = Path(output_dir).resolve()
    manifest_path = result_dir / "RUN_MANIFEST.json"
    if not manifest_path.is_file():
        present = [
            name for name in RESULT_FILES
            if (result_dir / name).exists()
        ]
        if allow_absent and not present:
            return {
                "status": "absent",
                "mode": mode,
                "output_dir": str(result_dir),
                "note": (
                    "no bilateral artifact has been published; source tests "
                    "and contracts remain enforceable"),
            }
        raise FileNotFoundError(
            f"bilateral manifest is missing from {result_dir}")

    manifest = _load_json(manifest_path, "bilateral manifest")
    payload = _load_json(result_dir / RESULT_JSON, "bilateral result")
    if set(manifest) != MANIFEST_FIELDS:
        missing = sorted(MANIFEST_FIELDS - set(manifest))
        extra = sorted(set(manifest) - MANIFEST_FIELDS)
        raise RuntimeError(
            "bilateral manifest schema is nonexact "
            f"(missing={missing}, extra={extra})")
    _require_header(manifest, "bilateral manifest", manifest=True)
    _validate_result_contract(payload)
    _validate_current_provenance(
        payload,
        root,
        require_current_runtime=mode == "publication",
    )
    declared = manifest.get("result_files_sha256")
    if not isinstance(declared, dict) or set(declared) != set(RESULT_FILES):
        raise RuntimeError("bilateral manifest has a nonexact result-file set")
    for name in RESULT_FILES:
        expected = _require_sha256(declared[name], f"hash for {name}")
        path = result_dir / name
        if not path.is_file() or file_sha256(path) != expected:
            raise RuntimeError(f"bilateral artifact bytes differ for {name}")
    for key in (
        "code_revision",
        "code_dirty",
        "source_tree_sha256",
        "runtime_versions",
    ):
        if manifest.get(key) != payload.get("provenance", {}).get(key):
            raise RuntimeError(f"manifest/result provenance differs at {key}")
    if manifest.get("code_dirty") is not False:
        raise RuntimeError("published bilateral manifest claims a dirty tree")
    if manifest.get("input_lineage") != payload.get("input_lineage"):
        raise RuntimeError("manifest/result input lineage differs")
    _validate_csv_bytes(payload, result_dir)
    lineage = _validate_public_pins(payload, root)

    skipped = []
    sidecar = _repository_path(
        root,
        lineage["scalp_sidecar_path_relative"],
        "sidecar path",
    )
    if mode == "publication":
        _validate_private_lineage(payload, root)
        _validate_figure_bytes(payload, result_dir)
    else:
        skipped.append(str(sidecar.relative_to(root)))
    return {
        "status": "complete",
        "mode": mode,
        "payload": payload,
        "manifest": manifest,
        "skipped_private_lineage": skipped,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--artifact-validation",
        choices=("offline", "publication"),
        default="offline",
    )
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = validate_exploratory_artifacts(
        ROOT,
        args.output_dir,
        mode=args.artifact_validation,
        allow_absent=True,
    )
    print(
        f"HUP138 bilateral artifact validation: {result['status']} "
        f"({args.artifact_validation})",
        flush=True,
    )


if __name__ == "__main__":
    main()
