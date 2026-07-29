"""Summaries, tabular exports, and figures for paired scalp/iEEG results.

This module is intentionally limited to reporting over completed subject records.
It does not load caches, materialize signals, or compute subject-level endpoints.
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import os
from collections import defaultdict

import numpy as np
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


STAGES_3B = ("N2", "N3", "NREM")
TARGET_COMPATIBLE_PEAK_BAND_HZ = (0.015, 0.025)
EXACT_SIGNFLIP_MAX_PAIRS = 15
MONTE_CARLO_SIGNFLIP_DRAWS = 20_000

_PAIRED_OBSERVATION_CSV_FIELDS = (
    "subject",
    "question",
    "comparison_role",
    "sensor",
    "stage",
    "modality",
    "sensor_available",
    "available",
    "unavailability_reason",
    "n_events_or_bouts",
    "pre_window_candidate_count",
    "peak_hz",
    "coherence_at_0p02_hz",
    "coherence_analytic_threshold",
    "coherence_above_threshold",
    "coherence_K",
    "shared_eeg_support_fraction",
    "shared_cardiac_support_fraction",
    "sigma_peak_window_mean",
    "swa_same_window_mean",
    "sigma_to_swa_same_peak_window_ratio",
    "xcorr_lecci_direction_peak_r",
    "xcorr_lecci_direction_peak_lag_s",
    "local_hr_change_pct",
    "pct_above_stage_mean",
    "so_hr_peak_lag_s",
    "exploratory",
)
# These are separate artifact contracts even though their current ordered
# columns are identical. Keeping both names prevents either consumer from
# treating an observed file header as its schema.
PAIRED_METRICS_CSV_FIELDS = _PAIRED_OBSERVATION_CSV_FIELDS
ROLE_PAIR_METRICS_CSV_FIELDS = _PAIRED_OBSERVATION_CSV_FIELDS
_PAIRED_CSV_FIELDS_BY_NAME = {
    "paired_metrics.csv": PAIRED_METRICS_CSV_FIELDS,
    "role_pair_metrics.csv": ROLE_PAIR_METRICS_CSV_FIELDS,
}


def _metric(result, path):
    value = result
    for key in path:
        if not isinstance(value, dict) or value.get(key) is None:
            return None
        value = value[key]
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def _negative_control_ratio(result):
    sigma = _metric(
        result, ("negative_control", "sigma_window_mean"))
    swa = _metric(
        result, ("negative_control", "swa_same_window_mean"))
    if sigma is None or swa is None or swa <= 0:
        return None
    return float(sigma / swa)


def _pair_summary(pairs):
    pairs = [
        (float(left), float(right))
        for left, right in pairs
        if left is not None and right is not None
        and np.isfinite(left) and np.isfinite(right)
    ]
    if not pairs:
        return {
            "n_pairs": 0,
            "ieeg_median": None,
            "scalp_median": None,
            "ieeg_mean": None,
            "scalp_mean": None,
            "scalp_minus_ieeg_median": None,
            "scalp_minus_ieeg_mean": None,
            "exploratory_exact_signflip_p": None,
            "exploratory_signflip_p": None,
            "signflip_method": None,
            "signflip_draws": 0,
            "signflip_seed": None,
            "signflip_statistic": "absolute mean paired difference",
        }
    array = np.asarray(pairs, float)
    differences = array[:, 1] - array[:, 0]
    observed = abs(float(np.mean(differences)))
    if len(differences) <= EXACT_SIGNFLIP_MAX_PAIRS:
        signed = np.asarray([
            abs(float(np.mean(differences * np.asarray(signs))))
            for signs in itertools.product((-1.0, 1.0), repeat=len(differences))
        ])
        p_value = float(np.mean(signed >= observed - 1e-15))
        exact_p_value = p_value
        signflip_method = "exact enumeration"
        signflip_draws = int(2 ** len(differences))
        signflip_seed = None
    else:
        seed = int.from_bytes(
            hashlib.sha256(
                np.asarray(differences, dtype="<f8").tobytes()
            ).digest()[:4],
            "little",
        )
        rng = np.random.RandomState(seed)
        exceedances = 0
        generated = 0
        while generated < MONTE_CARLO_SIGNFLIP_DRAWS:
            batch = min(
                1024,
                MONTE_CARLO_SIGNFLIP_DRAWS - generated,
            )
            signs = rng.choice(
                (-1.0, 1.0),
                size=(batch, len(differences)),
            )
            statistics = np.abs(
                np.mean(signs * differences[None, :], axis=1))
            exceedances += int(
                np.sum(statistics >= observed - 1e-15))
            generated += batch
        p_value = float(
            (1 + exceedances)
            / (MONTE_CARLO_SIGNFLIP_DRAWS + 1)
        )
        exact_p_value = None
        signflip_method = "deterministic Monte Carlo sign flips"
        signflip_draws = MONTE_CARLO_SIGNFLIP_DRAWS
        signflip_seed = seed
    try:
        wilcoxon = stats.wilcoxon(
            differences,
            alternative="two-sided",
            method="auto",
        )
        wilcoxon_p = float(wilcoxon.pvalue)
    except ValueError:
        wilcoxon_p = None
    return {
        "n_pairs": int(len(array)),
        "ieeg_median": float(np.median(array[:, 0])),
        "scalp_median": float(np.median(array[:, 1])),
        "ieeg_mean": float(np.mean(array[:, 0])),
        "scalp_mean": float(np.mean(array[:, 1])),
        "scalp_minus_ieeg_median": float(np.median(differences)),
        "scalp_minus_ieeg_mean": float(np.mean(differences)),
        "scalp_minus_ieeg_values": differences.tolist(),
        "exploratory_exact_signflip_p": exact_p_value,
        "exploratory_signflip_p": p_value,
        "signflip_method": signflip_method,
        "signflip_draws": signflip_draws,
        "signflip_seed": signflip_seed,
        "exploratory_wilcoxon_p": wilcoxon_p,
        "signflip_statistic": "absolute mean paired difference",
        "inference_warning": (
            "post-audit exploratory, uncorrected, small-n comparison; a "
            "non-significant difference is not evidence of equivalence"),
    }


def _group_summary(records):
    result = {
        "3A_C3": {},
        "3B": {},
    }
    metrics_3a = {
        "accepted_peak_hz": ("peak", "peak_hz"),
        "coherence_at_0p02_hz": ("coherence", "at_0p02_hz"),
        "lecci_direction_peak_r": (
            "cross_correlation", "lecci_direction_peak_r"),
        "lecci_direction_peak_lag_s": (
            "cross_correlation", "lecci_direction_peak_lag_s"),
    }
    for name, path in metrics_3a.items():
        pairs = []
        participants = []
        for record in records:
            left = _metric(record["ieeg"]["result_3a"], path)
            right = _metric(record["scalp"]["result_3a"], path)
            if left is not None and right is not None:
                pairs.append((left, right))
                participants.append(record["subject"])
        result["3A_C3"][name] = {
            **_pair_summary(pairs),
            "participants": participants,
        }
    negative_control_pairs = []
    negative_control_participants = []
    for record in records:
        left = _negative_control_ratio(record["ieeg"]["result_3a"])
        right = _negative_control_ratio(record["scalp"]["result_3a"])
        if left is not None and right is not None:
            negative_control_pairs.append((left, right))
            negative_control_participants.append(record["subject"])
    result["3A_C3"]["sigma_to_swa_same_peak_window_ratio"] = {
        **_pair_summary(negative_control_pairs),
        "participants": negative_control_participants,
        "interpretation": (
            "descriptive frequency-specificity control at each arm's fitted "
            "sigma-peak window; ratio >1 means normalized sigma-window power "
            "exceeded normalized SWA power, not a significance test"),
    }
    ieeg_peak = [
        value["subject"] for value in records
        if bool(value["ieeg"]["result_3a"]["peak"].get("accepted"))
    ]
    scalp_peak = [
        value["subject"] for value in records
        if bool(value["scalp"]["result_3a"]["peak"].get("accepted"))
    ]
    ieeg_coherence = [
        value["subject"] for value in records
        if bool(
            (value["ieeg"]["result_3a"].get("coherence") or {})
            .get("above_analytic_threshold"))
    ]
    scalp_coherence = [
        value["subject"] for value in records
        if bool(
            (value["scalp"]["result_3a"].get("coherence") or {})
            .get("above_analytic_threshold"))
    ]
    def target_compatible_peak(record, arm):
        peak = record[arm]["result_3a"]["peak"]
        value = peak.get("peak_hz")
        return bool(
            peak.get("accepted")
            and value is not None
            and TARGET_COMPATIBLE_PEAK_BAND_HZ[0]
            <= float(value)
            <= TARGET_COMPATIBLE_PEAK_BAND_HZ[1]
        )

    ieeg_target_peak = [
        value["subject"] for value in records
        if target_compatible_peak(value, "ieeg")
    ]
    scalp_target_peak = [
        value["subject"] for value in records
        if target_compatible_peak(value, "scalp")
    ]
    result["3A_C3"]["endpoint_counts"] = {
        "requested": len(records),
        "paired_spectra": sum(
            bool(value["ieeg"]["result_3a"]["endpoint_availability"]["spectrum"])
            and bool(value["scalp"]["result_3a"]["endpoint_availability"]["spectrum"])
            for value in records),
        "paired_coherence": sum(
            value["ieeg"]["result_3a"]["coherence"] is not None
            and value["scalp"]["result_3a"]["coherence"] is not None
            for value in records),
        "paired_cross_correlation": sum(
            value["ieeg"]["result_3a"]["cross_correlation"] is not None
            and value["scalp"]["result_3a"]["cross_correlation"] is not None
            for value in records),
        "both_accepted_peaks": sum(
            bool(value["ieeg"]["result_3a"]["peak"].get("accepted"))
            and bool(value["scalp"]["result_3a"]["peak"].get("accepted"))
            for value in records),
        "ieeg_accepted_peak_subjects": ieeg_peak,
        "scalp_accepted_peak_subjects": scalp_peak,
        "accepted_peak_both_subjects": sorted(
            set(ieeg_peak) & set(scalp_peak)),
        "accepted_peak_ieeg_only_subjects": sorted(
            set(ieeg_peak) - set(scalp_peak)),
        "accepted_peak_scalp_only_subjects": sorted(
            set(scalp_peak) - set(ieeg_peak)),
        "target_compatible_peak_band_hz": list(
            TARGET_COMPATIBLE_PEAK_BAND_HZ),
        "target_compatible_peak_band_status": (
            "descriptive +/-0.005-Hz window around 0.02 Hz requested for "
            "interpretation; not a Lecci-defined acceptance threshold"),
        "ieeg_target_compatible_accepted_peak_subjects": ieeg_target_peak,
        "scalp_target_compatible_accepted_peak_subjects": scalp_target_peak,
        "target_compatible_accepted_peak_both_subjects": sorted(
            set(ieeg_target_peak) & set(scalp_target_peak)),
        "target_compatible_accepted_peak_ieeg_only_subjects": sorted(
            set(ieeg_target_peak) - set(scalp_target_peak)),
        "target_compatible_accepted_peak_scalp_only_subjects": sorted(
            set(scalp_target_peak) - set(ieeg_target_peak)),
        "ieeg_coherence_above_threshold_subjects": ieeg_coherence,
        "scalp_coherence_above_threshold_subjects": scalp_coherence,
        "coherence_above_threshold_both_subjects": sorted(
            set(ieeg_coherence) & set(scalp_coherence)),
        "coherence_above_threshold_ieeg_only_subjects": sorted(
            set(ieeg_coherence) - set(scalp_coherence)),
        "coherence_above_threshold_scalp_only_subjects": sorted(
            set(scalp_coherence) - set(ieeg_coherence)),
    }
    for role in ("f3", "fz"):
        result["3B"][role] = {}
        present = [
            record["subject"] for record in records
            if bool(
                record["scalp"]["result_3b"].get(role, {}).get("available"))
        ]
        absent = [
            record["subject"] for record in records
            if not bool(
                record["scalp"]["result_3b"].get(role, {}).get("available"))
        ]
        result["3B"][role]["sensor_inventory"] = {
            "requested_participants": len(records),
            "sensor_present_n": len(present),
            "sensor_present_participants": present,
            "sensor_absent_n": len(absent),
            "sensor_absent_participants": absent,
            "note": (
                "sensor presence is distinct from whether a stage-specific "
                "endpoint has enough events and cardiac support"),
        }
        for stage in STAGES_3B:
            stage_pairs = defaultdict(list)
            participants = defaultdict(list)
            for record in records:
                scalp_role = record["scalp"]["result_3b"].get(role)
                if not scalp_role:
                    continue
                scalp_stage = scalp_role.get("stages", {}).get(stage, {})
                scalp_estimate = scalp_stage.get("estimate")
                ieeg_estimate = (
                    record["ieeg"]["result_3b"]
                    .get("stages", {}).get(stage, {}).get("estimate")
                )
                for metric in (
                    "event_locked_local_change_pct",
                    "pct_above_stage_mean",
                    "peak_lag_s",
                ):
                    left = _metric(ieeg_estimate, (metric,))
                    right = _metric(scalp_estimate, (metric,))
                    if left is not None and right is not None:
                        stage_pairs[metric].append((left, right))
                        participants[metric].append(record["subject"])
            result["3B"][role][stage] = {
                metric: {
                    **_pair_summary(stage_pairs[metric]),
                    "participants": participants[metric],
                }
                for metric in (
                    "event_locked_local_change_pct",
                    "pct_above_stage_mean",
                    "peak_lag_s",
                )
            }
    naji_means = {"N2": 12.09, "N3": 3.35}
    result["3B"]["f3"]["naji_reported_mean_context"] = {}
    for stage, reference_mean in naji_means.items():
        scalp_values = []
        participants = []
        for record in records:
            scalp_role = record["scalp"]["result_3b"].get("f3")
            estimate = (
                None if not scalp_role
                else scalp_role.get("stages", {}).get(stage, {}).get("estimate")
            )
            value = _metric(estimate, ("pct_above_stage_mean",))
            if value is not None:
                scalp_values.append(value)
                participants.append(record["subject"])
        median = (
            float(np.median(scalp_values)) if scalp_values else None)
        result["3B"]["f3"]["naji_reported_mean_context"][stage] = {
            "naji_healthy_scalp_reported_mean_pct": reference_mean,
            "our_unilateral_f3_median_pct": median,
            "our_unilateral_f3_median_over_naji_mean": (
                None if median is None else float(median / reference_mean)),
            "n_participants": len(scalp_values),
            "participants": participants,
            "warning": (
                "Naji's value is a published group mean, not an acceptance "
                "threshold; our F3 is unilateral with an undocumented reference "
                "and comes from epilepsy inpatients"),
        }
    result["interpretation_guardrail"] = (
        "These compare downstream cortical/cardiac observables. Agreement between "
        "scalp and iEEG does not identify LC as their cause.")
    return result


def _joined_support_reasons(value):
    reasons = value.get("support_reasons", []) if isinstance(value, dict) else []
    return "; ".join(str(reason) for reason in reasons) or None


def _ieeg_3b_csv_row(record, stage, comparison_role):
    ieeg_stage = (
        record["ieeg"]["result_3b"]
        .get("stages", {}).get(stage, {})
    )
    estimate = ieeg_stage.get("estimate")
    candidate_count = sum(
        int(value) for value in
        ieeg_stage.get("candidate_so_counts_by_contact", {}).values()
    )
    return {
        "subject": record["subject"],
        "question": "3B",
        "comparison_role": comparison_role,
        "sensor": "iEEG aggregate",
        "stage": stage,
        "modality": "iEEG",
        "sensor_available": True,
        "available": estimate is not None,
        "unavailability_reason": (
            None if estimate is not None
            else _joined_support_reasons(ieeg_stage)
        ),
        "n_events_or_bouts": (
            None if estimate is None else estimate.get("n_so_total")
        ),
        "pre_window_candidate_count": candidate_count,
        "peak_hz": None,
        "coherence_at_0p02_hz": None,
        "coherence_analytic_threshold": None,
        "coherence_above_threshold": None,
        "coherence_K": None,
        "shared_eeg_support_fraction": None,
        "shared_cardiac_support_fraction": None,
        "sigma_peak_window_mean": None,
        "swa_same_window_mean": None,
        "sigma_to_swa_same_peak_window_ratio": None,
        "xcorr_lecci_direction_peak_r": None,
        "xcorr_lecci_direction_peak_lag_s": None,
        "local_hr_change_pct": _metric(
            estimate, ("event_locked_local_change_pct",)),
        "pct_above_stage_mean": _metric(
            estimate, ("pct_above_stage_mean",)),
        "so_hr_peak_lag_s": _metric(estimate, ("peak_lag_s",)),
        "exploratory": True,
    }


def _role_pair_csv_rows(records):
    """Return rows keyed to an explicit scalp comparison role.

    An iEEG estimate is intentionally repeated when the same participant-stage
    contributes to both the F3 and Fz contrasts.  This table is for reconstructing
    role-specific pairs, not for pooling observations across roles.
    """
    rows = []
    for record in records:
        subject = record["subject"]
        for modality, result in (
            ("iEEG", record["ieeg"]["result_3a"]),
            ("scalp_C3", record["scalp"]["result_3a"]),
        ):
            rows.append({
                "subject": subject,
                "question": "3A",
                "comparison_role": "c3",
                "sensor": "iEEG aggregate" if modality == "iEEG" else "C3/C03",
                "stage": "NREM first 210 min",
                "modality": modality,
                "sensor_available": True,
                "available": result["endpoint_availability"]["spectrum"],
                "unavailability_reason": (
                    None
                    if result["endpoint_availability"]["spectrum"]
                    else _joined_support_reasons(result)
                ),
                "n_events_or_bouts": result.get("n_bouts"),
                "pre_window_candidate_count": None,
                "peak_hz": _metric(result, ("peak", "peak_hz")),
                "coherence_at_0p02_hz": _metric(
                    result, ("coherence", "at_0p02_hz")),
                "coherence_analytic_threshold": _metric(
                    result, ("coherence", "analytic_threshold")),
                "coherence_above_threshold": (
                    None
                    if result.get("coherence") is None
                    else bool(
                        result["coherence"]["above_analytic_threshold"])
                ),
                "coherence_K": (
                    None
                    if result.get("coherence") is None
                    else int(result["coherence"]["K"])
                ),
                "shared_eeg_support_fraction": record[
                    "shared_inputs"]["paired_3a_endpoint_support"][
                        "eeg_common"]["fraction"],
                "shared_cardiac_support_fraction": record[
                    "shared_inputs"]["paired_3a_endpoint_support"][
                        "cardiac_common"]["fraction"],
                "sigma_peak_window_mean": _metric(
                    result, ("negative_control", "sigma_window_mean")),
                "swa_same_window_mean": _metric(
                    result, ("negative_control", "swa_same_window_mean")),
                "sigma_to_swa_same_peak_window_ratio": (
                    _negative_control_ratio(result)),
                "xcorr_lecci_direction_peak_r": _metric(
                    result, (
                        "cross_correlation", "lecci_direction_peak_r")),
                "xcorr_lecci_direction_peak_lag_s": _metric(
                    result, (
                        "cross_correlation",
                        "lecci_direction_peak_lag_s")),
                "local_hr_change_pct": None,
                "pct_above_stage_mean": None,
                "so_hr_peak_lag_s": None,
                "exploratory": True,
            })
        for role in ("f3", "fz"):
            scalp_role = record["scalp"]["result_3b"].get(role)
            if not scalp_role or not bool(scalp_role.get("available")):
                reason = (
                    "missing structured sensor record"
                    if not scalp_role
                    else scalp_role.get("reason", "sensor absent")
                )
                for stage in STAGES_3B:
                    rows.append({
                        "subject": subject,
                        "question": "3B",
                        "comparison_role": role,
                        "sensor": None,
                        "stage": stage,
                        "modality": f"scalp_{role.upper()}",
                        "sensor_available": False,
                        "available": False,
                        "unavailability_reason": reason,
                        "n_events_or_bouts": None,
                        "pre_window_candidate_count": None,
                        "peak_hz": None,
                        "coherence_at_0p02_hz": None,
                        "coherence_analytic_threshold": None,
                        "coherence_above_threshold": None,
                        "coherence_K": None,
                        "shared_eeg_support_fraction": None,
                        "shared_cardiac_support_fraction": None,
                        "sigma_peak_window_mean": None,
                        "swa_same_window_mean": None,
                        "sigma_to_swa_same_peak_window_ratio": None,
                        "xcorr_lecci_direction_peak_r": None,
                        "xcorr_lecci_direction_peak_lag_s": None,
                        "local_hr_change_pct": None,
                        "pct_above_stage_mean": None,
                        "so_hr_peak_lag_s": None,
                        "exploratory": True,
                    })
                continue
            for stage in STAGES_3B:
                scalp_stage = scalp_role.get("stages", {}).get(stage, {})
                scalp_est = scalp_stage.get("estimate")
                rows.append(_ieeg_3b_csv_row(
                    record, stage, comparison_role=role))
                rows.append({
                    "subject": subject,
                    "question": "3B",
                    "comparison_role": role,
                    "sensor": scalp_role.get("channel"),
                    "stage": stage,
                    "modality": f"scalp_{role.upper()}",
                    "sensor_available": True,
                    "available": scalp_est is not None,
                    "unavailability_reason": (
                        None if scalp_est is not None
                        else _joined_support_reasons(scalp_stage)
                    ),
                    "n_events_or_bouts": (
                        None if scalp_est is None
                        else scalp_est.get("n_so_total")
                    ),
                    "pre_window_candidate_count": scalp_stage.get(
                        "candidate_so_count"),
                    "peak_hz": None,
                    "coherence_at_0p02_hz": None,
                    "coherence_analytic_threshold": None,
                    "coherence_above_threshold": None,
                    "coherence_K": None,
                    "shared_eeg_support_fraction": None,
                    "shared_cardiac_support_fraction": None,
                    "sigma_peak_window_mean": None,
                    "swa_same_window_mean": None,
                    "sigma_to_swa_same_peak_window_ratio": None,
                    "xcorr_lecci_direction_peak_r": None,
                    "xcorr_lecci_direction_peak_lag_s": None,
                    "local_hr_change_pct": _metric(
                        scalp_est, ("event_locked_local_change_pct",)),
                    "pct_above_stage_mean": _metric(
                        scalp_est, ("pct_above_stage_mean",)),
                    "so_hr_peak_lag_s": _metric(
                        scalp_est, ("peak_lag_s",)),
                    "exploratory": True,
                })
    return rows


def _normalized_csv_rows(records, role_pair_rows=None):
    """Return one observation per subject/question/stage/modality.

    In particular, each 3B iEEG participant-stage appears exactly once and is
    not duplicated merely because both an F3 and an Fz scalp comparison exists.
    The separate role-pair table retains those explicit pair memberships.
    """
    if role_pair_rows is None:
        role_pair_rows = _role_pair_csv_rows(records)

    rows = [
        dict(row) for row in role_pair_rows
        if row["question"] == "3A"
    ]
    scalp_rows = {
        (
            row["subject"],
            row["stage"],
            row["modality"],
        ): row
        for row in role_pair_rows
        if row["question"] == "3B" and row["modality"] != "iEEG"
    }
    for record in records:
        for stage in STAGES_3B:
            rows.append(_ieeg_3b_csv_row(
                record, stage, comparison_role="subject_stage"))
            for role in ("f3", "fz"):
                modality = f"scalp_{role.upper()}"
                key = (record["subject"], stage, modality)
                if key not in scalp_rows:
                    raise RuntimeError(
                        "role-pair export is missing normalized scalp row "
                        f"{key!r}"
                    )
                rows.append(dict(scalp_rows[key]))

    keys = [
        (
            row["subject"],
            row["question"],
            row["stage"],
            row["modality"],
        )
        for row in rows
    ]
    if len(keys) != len(set(keys)):
        raise RuntimeError(
            "normalized CSV would duplicate a "
            "subject/question/stage/modality observation"
        )
    return rows


def _paired_csv_fields(path):
    name = os.path.basename(os.fspath(path))
    try:
        fields = _PAIRED_CSV_FIELDS_BY_NAME[name]
    except KeyError as error:
        raise ValueError(
            f"no canonical paired-CSV schema is registered for {name!r}"
        ) from error
    if len(fields) != len(set(fields)):
        raise RuntimeError(
            f"canonical paired-CSV schema for {name!r} has duplicate columns")
    return name, fields


def _serialized_csv_rows(rows, fields, *, artifact_name):
    serialized = []
    expected_keys = set(fields)
    for index, row in enumerate(rows, start=2):
        if not isinstance(row, dict):
            raise RuntimeError(
                f"{artifact_name} expected row {index} is not an object")
        actual_keys = set(row)
        if actual_keys != expected_keys:
            missing = [
                field for field in fields if field not in actual_keys
            ]
            extra = sorted(actual_keys - expected_keys)
            raise RuntimeError(
                f"{artifact_name} expected row {index} differs from its "
                f"canonical schema: missing={missing}, extra={extra}")
        serialized.append({
            field: "" if row[field] is None else str(row[field])
            for field in fields
        })
    return serialized


def _validate_csv_artifact(path, expected_rows):
    """Require one exact canonical header and exact regenerated CSV rows."""
    artifact_name, fields = _paired_csv_fields(path)
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        actual_fields = tuple(reader.fieldnames or ())
        if actual_fields != fields:
            duplicate_fields = sorted({
                field for field in actual_fields
                if actual_fields.count(field) > 1
            })
            missing = [
                field for field in fields if field not in actual_fields
            ]
            extra = [
                field for field in actual_fields
                if field not in fields
            ]
            order_mismatch = (
                not missing and not extra and not duplicate_fields
            )
            raise RuntimeError(
                f"{artifact_name} header differs from its exact ordered "
                f"canonical schema: missing={missing}, extra={extra}, "
                f"duplicates={duplicate_fields}, "
                f"order_mismatch={order_mismatch}")
        actual_rows = list(reader)
    expected_serialized = _serialized_csv_rows(
        expected_rows,
        fields,
        artifact_name=artifact_name,
    )
    if actual_rows != expected_serialized:
        raise RuntimeError(
            f"{artifact_name} rows do not exactly match regenerated "
            "subject-results rows")
    return actual_rows


def _write_csv(path, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    artifact_name, fields = _paired_csv_fields(path)
    _serialized_csv_rows(rows, fields, artifact_name=artifact_name)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _paired_plot(ax, records, left_path, right_path, title, ylabel,
                 accepted_only=False):
    count = 0
    for index, record in enumerate(records):
        left = _metric(record["ieeg"]["result_3a"], left_path)
        right = _metric(record["scalp"]["result_3a"], right_path)
        if accepted_only:
            left_accepted = bool(
                record["ieeg"]["result_3a"]["peak"].get("accepted"))
            right_accepted = bool(
                record["scalp"]["result_3a"]["peak"].get("accepted"))
            left = left if left_accepted else None
            right = right if right_accepted else None
        if left is None or right is None:
            color = plt.cm.tab10(index % 10)
            if left is not None:
                ax.plot(0, left, "o", color=color, alpha=0.8)
                count += 1
            if right is not None:
                ax.plot(1, right, "o", color=color, alpha=0.8)
                count += 1
            continue
        color = plt.cm.tab10(index % 10)
        ax.plot([0, 1], [left, right], "-o", color=color, alpha=0.8, lw=1)
        count += 2 if accepted_only else 1
    ax.set_xticks([0, 1], ["iEEG", "scalp C3"])
    ax.set_ylabel(ylabel)
    unit = "accepted arms" if accepted_only else "paired n"
    ax.set_title(f"{title}\n{unit}={count}")
    ax.grid(axis="y", alpha=0.25)


def _paired_coherence_ratio_plot(ax, records):
    count = 0
    for index, record in enumerate(records):
        ratios = []
        for arm in ("ieeg", "scalp"):
            coherence = record[arm]["result_3a"].get("coherence")
            if not coherence:
                ratios.append(None)
                continue
            ratios.append(
                float(coherence["at_0p02_hz"])
                / float(coherence["analytic_threshold"])
            )
        left, right = ratios
        if left is None or right is None:
            continue
        color = plt.cm.tab10(index % 10)
        ax.plot([0, 1], [left, right], "-", color=color, alpha=0.8, lw=1)
        for x_value, ratio in enumerate((left, right)):
            marker = "*" if ratio > 1 else "o"
            size = 10 if marker == "*" else 5
            ax.plot(
                x_value, ratio, marker=marker, markersize=size,
                color=color, alpha=0.9)
        count += 1
    ax.axhline(1.0, ls="--", color="0.25", lw=1)
    ax.set_xticks([0, 1], ["iEEG", "scalp C3"])
    ax.set_ylabel("coherence / participant threshold")
    ax.set_title(
        f"0.02-Hz coherence relative to threshold\n"
        f"paired n={count}; star means passes (>1)")
    ax.grid(axis="y", alpha=0.25)


def _make_3a_figure(records, output_dir):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    _paired_plot(
        axes[0, 0], records,
        ("peak", "peak_hz"), ("peak", "peak_hz"),
        "Accepted peak estimates (either arm)", "Hz", accepted_only=True)
    axes[0, 0].axhspan(
        TARGET_COMPATIBLE_PEAK_BAND_HZ[0],
        TARGET_COMPATIBLE_PEAK_BAND_HZ[1],
        color="0.75", alpha=0.25,
        label="descriptive target-compatible band",
    )
    axes[0, 0].axhline(0.02, ls="--", color="0.4", lw=1)
    axes[0, 0].legend(loc="best", frameon=False, fontsize=8)
    _paired_coherence_ratio_plot(axes[0, 1], records)
    _paired_plot(
        axes[1, 0], records,
        ("cross_correlation", "lecci_direction_peak_r"),
        ("cross_correlation", "lecci_direction_peak_r"),
        "Lecci-direction HR→sigma peak", "correlation r")

    count = 0
    for index, record in enumerate(records):
        left = record["ieeg"]["raw_power"].get(
            "median_selected_contact_sigma_power")
        right = record["scalp"]["raw_power"].get("raw_sigma_median")
        if left is None or right is None or left <= 0 or right <= 0:
            continue
        axes[1, 1].plot(
            [0, 1], [left, right], "-o",
            color=plt.cm.tab10(index % 10), alpha=0.8, lw=1)
        count += 1
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xticks([0, 1], ["iEEG", "scalp C3"])
    axes[1, 1].set_ylabel("recorded-unit sigma power (log scale)")
    axes[1, 1].set_title(
        f"Descriptive raw 10–15 Hz power\npaired n={count}; not a "
        "calibrated amplitude comparison")
    axes[1, 1].grid(axis="y", alpha=0.25)
    fig.suptitle(
        "3A simultaneous sensor comparison — same participant, interval, "
        "ECG, sleep labels, and finite sample support",
        fontsize=13,
    )
    handles = [
        plt.Line2D(
            [0], [0], color=plt.cm.tab10(index % 10), marker="o", lw=1,
            label=record["subject"].replace("_phaseII", ""))
        for index, record in enumerate(records)
    ]
    fig.legend(
        handles=handles, loc="center left", ncol=1, frameon=False,
        bbox_to_anchor=(1.01, 0.5))
    for suffix in ("png", "svg"):
        fig.savefig(
            os.path.join(output_dir, f"paired_3A_C3.{suffix}"),
            dpi=180, bbox_inches="tight",
        )
    plt.close(fig)

def _make_3b_figure(records, output_dir):
    fig, axes = plt.subplots(
        2, 3, figsize=(13, 7.5), constrained_layout=True,
        sharey="row")
    for row, role in enumerate(("f3", "fz")):
        for column, stage in enumerate(STAGES_3B):
            ax = axes[row, column]
            count = 0
            for index, record in enumerate(records):
                scalp_role = record["scalp"]["result_3b"].get(role)
                if not scalp_role:
                    continue
                scalp_est = (
                    scalp_role.get("stages", {})
                    .get(stage, {}).get("estimate")
                )
                ieeg_est = (
                    record["ieeg"]["result_3b"]
                    .get("stages", {}).get(stage, {}).get("estimate")
                )
                left = _metric(
                    ieeg_est, ("event_locked_local_change_pct",))
                right = _metric(
                    scalp_est, ("event_locked_local_change_pct",))
                if left is None or right is None:
                    continue
                ax.plot(
                    [0, 1], [left, right], "-o",
                    color=plt.cm.tab10(index % 10), alpha=0.8, lw=1)
                count += 1
            ax.set_xticks([0, 1], ["iEEG", f"scalp {role.upper()}"])
            ax.set_title(f"{role.upper()} {stage}: paired n={count}")
            ax.axhline(0, color="0.5", lw=0.8)
            ax.grid(axis="y", alpha=0.25)
            if column == 0:
                ax.set_ylabel("event-locked local HR change (%)")
    fig.suptitle(
        "3B exploratory simultaneous sensor comparison\n"
        "F3 is unilateral/reference-incomplete; Fz was not a Naji sensor; "
        "inference disabled",
        fontsize=13,
    )
    handles = [
        plt.Line2D(
            [0], [0], color=plt.cm.tab10(index % 10), marker="o", lw=1,
            label=record["subject"].replace("_phaseII", ""))
        for index, record in enumerate(records)
    ]
    fig.legend(
        handles=handles, loc="center left", ncol=1, frameon=False,
        bbox_to_anchor=(1.01, 0.5))
    for suffix in ("png", "svg"):
        fig.savefig(
            os.path.join(output_dir, f"paired_3B_F3_Fz.{suffix}"),
            dpi=180, bbox_inches="tight",
        )
    plt.close(fig)
