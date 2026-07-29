"""Reusable descriptive estimators for the Naji-aligned 3B endpoint.

The production estimator averages one SO-triggered tachogram per contact,
then finds one participant-level post-trough peak.  Its whole-stage circular
shift is retained only as a diagnostic; event-locking z and p values are
disabled because that shift does not preserve local trends or event density.
"""
from __future__ import annotations

import numpy as np


FS_RR = 4.0
HALF_WIN = 5.0


def rr_baseline_hr(rr):
    """Convert mean RR to HR in the same order used for RR event curves."""
    values = np.asarray(rr, float)
    values = values[np.isfinite(values) & (values > 0)]
    return float(60.0 / values.mean()) if len(values) else np.nan


def subject_so_triggered(
        tachogram, trough_times_by_channel, stage_mean_hr, stage_pool_idx,
        n_sur=1000, rng=None, domain="hr", minimum_channels=1,
        channel_ids=None, minimum_events_per_channel=30,
        minimum_surrogate_pool_samples=100):
    """Return one descriptive SO-to-HR estimate for one participant.

    ``domain="rr"`` is the paper-aligned path: contact curves are averaged in
    RR space and minima define HR-burst timing.  ``domain="hr"`` is retained
    for explicit sensitivity analyses and withdrawn compatibility callers.
    Both observed and shifted windows must be finite and wholly within stage.
    """
    if domain not in ("rr", "hr"):
        raise ValueError("domain must be 'rr' or 'hr'")
    if int(n_sur) < 2:
        raise ValueError("at least two diagnostic surrogates are required")
    series = np.asarray(tachogram, float)
    rng = rng or np.random.RandomState(0)
    half_window = int(HALF_WIN * FS_RR)

    finite = np.isfinite(series)
    finite_cumulative = np.concatenate(([0], np.cumsum(finite)))
    candidate_centres = np.arange(
        half_window, len(series) - half_window + 1)
    full_finite_window = np.zeros(len(series), bool)
    full_finite_window[candidate_centres] = (
        finite_cumulative[candidate_centres + half_window]
        - finite_cumulative[candidate_centres - half_window]
    ) == (2 * half_window)

    stage_samples = np.asarray(stage_pool_idx, int)
    stage_samples = stage_samples[
        (stage_samples >= 0) & (stage_samples < len(series))]
    in_stage = np.zeros(len(series), bool)
    in_stage[stage_samples] = True
    stage_cumulative = np.concatenate(([0], np.cumsum(in_stage)))
    full_stage_window = np.zeros(len(series), bool)
    full_stage_window[candidate_centres] = (
        stage_cumulative[candidate_centres + half_window]
        - stage_cumulative[candidate_centres - half_window]
    ) == (2 * half_window)
    eligible = full_finite_window & full_stage_window
    pool = stage_samples[
        (stage_samples >= half_window)
        & (stage_samples <= len(series) - half_window)
    ]
    pool = np.unique(pool[eligible[pool]])
    if len(pool) < int(minimum_surrogate_pool_samples):
        return None

    if channel_ids is None:
        channel_ids = list(range(len(trough_times_by_channel)))
    if len(channel_ids) != len(trough_times_by_channel):
        raise ValueError(
            "channel_ids must align with trough_times_by_channel")

    channel_indices = []
    retained_channel_ids = []
    for channel_id, trough_times_s in zip(
            channel_ids, trough_times_by_channel):
        indices = np.round(
            np.asarray(trough_times_s, float) * FS_RR).astype(int)
        indices = indices[
            (indices >= half_window)
            & (indices <= len(series) - half_window)
        ]
        indices = np.unique(indices[eligible[indices]])
        if len(indices) >= int(minimum_events_per_channel):
            channel_indices.append(indices)
            retained_channel_ids.append(channel_id)
    if len(channel_indices) < int(minimum_channels):
        return None

    def channel_curves(index_arrays):
        return [
            np.stack([
                series[index - half_window:index + half_window]
                for index in indices
            ]).mean(axis=0)
            for indices in index_arrays
        ]

    def post_peak_index(curve, post_mask):
        values = curve[post_mask]
        return int(
            np.argmin(values) if domain == "rr" else np.argmax(values))

    def percent_above_stage_mean(curve, post_mask):
        index = post_peak_index(curve, post_mask)
        value = float(curve[post_mask][index])
        peak_hr = 60.0 / value if domain == "rr" else value
        return 100.0 * (peak_hr - stage_mean_hr) / stage_mean_hr

    contact_curves = channel_curves(channel_indices)
    curve = np.mean(contact_curves, axis=0)
    lag = (
        np.arange(len(curve)) - half_window) / FS_RR
    post = lag >= 0
    participant_peak_index = post_peak_index(curve, post)
    percent_change = percent_above_stage_mean(curve, post)

    contact_peak_lags = [
        float(lag[post][post_peak_index(contact_curve, post)])
        for contact_curve in contact_curves
    ]
    paper_peak_lag = float(np.mean(contact_peak_lags))

    event_ranks = []
    for indices in channel_indices:
        rank = np.searchsorted(pool, indices)
        if (
            np.any(rank >= len(pool))
            or not np.array_equal(pool[rank], indices)
        ):
            raise RuntimeError(
                "observed SO index was not found in the eligible stage pool")
        event_ranks.append(rank)

    # One shift preserves cross-contact coincidences and event-train structure.
    # It remains a diagnostic because a whole-stage shift does not preserve
    # local nonstationary trends or event-density clustering.
    guard = 2 * half_window
    valid_shifts = np.arange(guard, len(pool) - guard + 1)
    if not len(valid_shifts):
        valid_shifts = np.arange(1, len(pool))
    null = np.empty(int(n_sur), float)
    for surrogate_index in range(int(n_sur)):
        shift = int(rng.choice(valid_shifts))
        shifted_indices = [
            pool[(rank + shift) % len(pool)] for rank in event_ranks
        ]
        shifted_curve = np.mean(
            channel_curves(shifted_indices), axis=0)
        null[surrogate_index] = percent_above_stage_mean(
            shifted_curve, post)

    diagnostic_p = float(
        (1 + np.sum(null >= percent_change)) / (len(null) + 1))
    hr_curve = 60.0 / curve if domain == "rr" else curve
    pre = lag < 0
    local_baseline_hr = float(np.mean(hr_curve[pre]))
    local_change_pct = float(
        100.0
        * (np.max(hr_curve[post]) - local_baseline_hr)
        / local_baseline_hr
    )
    peak_to_peak_pct = float(
        100.0
        * (np.max(hr_curve) - np.min(hr_curve))
        / np.mean(hr_curve)
    )
    return {
        "n_channels": len(channel_indices),
        "event_contact_ids": retained_channel_ids,
        "n_so_total": int(sum(len(values) for values in channel_indices)),
        "pct_above_stage_mean": percent_change,
        "peak_lag_s": paper_peak_lag,
        "channel_peak_lag_s": contact_peak_lags,
        "participant_curve_peak_lag_s": float(
            lag[post][participant_peak_index]),
        "event_locked_local_change_pct": local_change_pct,
        "event_curve_peak_to_peak_pct": peak_to_peak_pct,
        "local_pre_event_mean_hr": local_baseline_hr,
        "z": None,
        "p_upper": None,
        "stage_shift_z_diagnostic": float(
            (percent_change - null.mean())
            / (null.std(ddof=1) + 1e-12)
        ),
        "stage_shift_p_upper_diagnostic": diagnostic_p,
        "null_mean_pct": float(null.mean()),
        "null_sd_pct": float(null.std(ddof=1)),
        "n_surrogates": int(len(null)),
        "n_surrogate_pool": int(len(pool)),
        "null_method": (
            "one shared circular shift in eligible stage-time across all "
            "channels"),
        "inference_status": (
            "disabled: whole-stage shifts do not preserve local "
            "nonstationary trends or event-density clustering; raw Naji "
            "magnitude and local change are descriptive"),
        "tachogram_domain": domain,
        "curve": [
            float(60.0 / value if domain == "rr" else value)
            for value in curve
        ],
        "rr_curve": (
            [float(value) for value in curve] if domain == "rr" else None),
        "lag_s": [float(value) for value in lag],
    }


def so_triggered(
        hr, trough_times_s, stage_mean_hr, stage_pool_idx,
        n_sur=200, rng=None):
    """Single-contact descriptive compatibility wrapper.

    This API is used only by withdrawn exploratory scripts.  It delegates to
    the participant estimator and never exposes inferential z or p values.
    """
    result = subject_so_triggered(
        hr,
        [np.asarray(trough_times_s, float)],
        stage_mean_hr,
        stage_pool_idx,
        n_sur=n_sur,
        rng=rng,
        domain="hr",
        minimum_channels=1,
        minimum_events_per_channel=30,
        minimum_surrogate_pool_samples=100,
    )
    if result is None:
        return None
    return {
        **result,
        "n_so": result["n_so_total"],
        "z": None,
        "p_upper": None,
        "inference_status": (
            "descriptive only; event-locking z/p disabled because whole-stage "
            "shifts do not preserve local nonstationarity"),
    }
