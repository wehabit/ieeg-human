"""Stream each subject's night ONCE and cache the derived series to disk.

Every re-analysis of 3A/3B so far required a fresh 7-hour stream from iEEG.org per subject, because
only summary statistics were ever saved. That made method fixes expensive and is the main reason
deviations from Lecci went uncorrected. This script pulls the night once and stores everything the
downstream tests need, so subsequent analyses run offline in seconds.

Cached per subject -> data/derived/lc_infraslow/<subject>.npz

    sigma_fixed   1 Hz, per-contact IED-masked sigma power, Lecci's fixed 10-15 Hz band
    sigma_fsp     compatibility copy using 12-14 Hz; individualized FSP is disabled until
                  all-clean-NREM estimation and manual quality control are available
    swa           1 Hz, IED-masked 0.5-4 Hz power (Lecci's negative control: no 0.02 Hz peak)
    *_power_numerator_by_contact / *_clean_sample_count_by_contact
                  reversible per-band 1-s support before the historical 50%-clean rule
    hr_1          1 Hz instantaneous heart rate
    hr_4          4 Hz instantaneous heart rate (Naji 2019 resolution, PCHIP of RR intervals)
    ep_dr/ep_swa/ep_clean   per-30 s-epoch staging features
    event_3d_*    20 Hz, per-contact pre-threshold spindle RMS, SO phase,
                  sample validity, and duration-qualified SO candidates
    so_candidate_*_<ch>  clean complete SO candidates and amplitudes; downstream stage-specific
                         thresholds are applied only after staging
    beats         R-peak times (s)

NaN means "not measured here" and is preserved deliberately -- downstream code must handle gaps
rather than delete-and-splice (see spectral_gapped.py for why).

    .venv/bin/python analysis/cache_lc_series.py [--subjects 165,157] [--hours 7]
"""
import argparse, concurrent.futures, json, os, time, traceback
from dataclasses import dataclass

import numpy as np
from scipy import signal, interpolate, ndimage
import neurokit2 as nk

from infraslow_rr_sigma_coherence import (
    IEEG_CONNECT_TIMEOUT_S, IEEG_READ_TIMEOUT_S, sess, notch, ROOT,
)
from signal_qc import ied_clean_mask
from spectral_gapped import contiguous_runs
from event_3d_estimators import (
    EVENT_FS as EVENT_3D_SAMPLING_HZ,
    IED_PAD_S as EVENT_3D_IED_PAD_S,
    SO_BAND as EVENT_3D_SO_BAND,
    SO_DUR as EVENT_3D_SO_DURATION_S,
    SPINDLE_BAND as EVENT_3D_SPINDLE_BAND,
    SPINDLE_RMS_WINDOW_S as EVENT_3D_RMS_WINDOW_S,
    conservative_resampled_interval,
    retain_complete_clean_event_3d_candidates,
    so_event_candidates as event_3d_so_candidates,
    spindle_rms as event_3d_spindle_rms,
)
from hup_portal import (
    COHORT, HUP_SOURCE_PIN_SCHEMA_VERSION, NIGHT_PROBE_WORKERS, cortical_channels,
    delta_ratio, expected_portal_sample_count, find_night,
    pull_continuous_exact, validate_hup_series_geometry,
    verify_hup_source_identity,
)
from staging_helpers import band_sos, EPOCH, CHUNK_S, SWA_BAND
from pipeline_version import (
    CACHE_SCHEMA_VERSION,
    atomic_savez,
    cache_code_sha256,
    file_sha256,
    git_is_dirty,
    git_revision,
    npz_scalar_text,
    require_integer_sample_rate,
    runtime_versions,
    source_tree_sha256,
    start_run_manifest,
    utc_now,
    validated_complete_run_exists,
    write_run_manifest,
)

OUT = os.path.join(ROOT, "data", "derived", "lc_infraslow")
SIGMA_FIXED = (10.0, 15.0)      # Lecci's band
SWA_BAND_L = (0.5, 4.0)       # Lecci human SWA control
SO_BAND_NAJI = (0.15, 4.0)
# Naji cites the Dang-Vu scalp detector: negative-to-positive zero crossing
# 0.3--1.5 s, followed by a positive-to-negative crossing in <1 s.  Absolute
# scalp-voltage gates do not transfer to iEEG, so amplitude is retained for an
# explicit offline percentile sensitivity analysis.
SO_NEGATIVE_HALF_DURATION_S = (0.3, 1.5)
SO_POSITIVE_HALF_MAX_S = 1.0
AMP_PCT = 75
FS_RR = 4.0
MIN_SIGNAL_COVERAGE = 0.80
MIN_CONTACT_COVERAGE = 0.80
MIN_CONTACTS = 3
MIN_CONTACT_FRACTION_PER_BIN = 0.80
FILTER_EDGE_S = 30.0
MISSING_PAD_S = 5.0
STAGING_WELCH_WINDOW_S = 4.0
STAGING_WELCH_OVERLAP_S = 2.0
# Fourteen 4-s windows at 2-s stride exactly cover a 30-s epoch.  Requiring all fourteen is the
# historical all-samples-clean rule; neutral caches retain every window so support can be
# calibrated and sensitivity-tested offline.
STAGING_REFERENCE_MIN_VALID_WINDOWS = 14
HUP_ANATOMY_SELECTION_METHOD = (
    "UNVALIDATED contact-number heuristic; lateral-contact candidates require "
    "coordinate/tissue/SOZ QC")
SIGNAL_FLAT_EPSILON_MULTIPLIER = 64.0


def finalize_ecg_cache_qc(
        ecg_failures, hr_coverage, minimum_coverage=MIN_SIGNAL_COVERAGE):
    """Enforce fatal detector exceptions while retaining low coverage as QC metadata.

    A detector exception can be state-dependent, so successful chunks cannot safely be used as a
    partial cardiac series.  This is different from a detector that completes successfully but
    yields sparse usable support: the latter remains available for endpoint-specific support
    checks and receives a warning here.
    """
    failures = list(ecg_failures)
    if failures:
        first = failures[0]
        first_error = first.get("error", "unreported detector error")
        raise RuntimeError(
            f"{len(failures)} ECG detector chunks failed; refusing status='ok' cache "
            f"(first error: {first_error})")
    warnings = []
    if hr_coverage < minimum_coverage:
        warnings.append(
            f"HR coverage {hr_coverage:.1%} is below the historical audit80 reference")
    return warnings


def sanitize_beats(beats, min_interval_s=0.25):
    """Sort/deduplicate R peaks and enforce the refractory interval sequentially."""
    values = np.unique(np.round(np.asarray(beats, float), 4))
    if not len(values):
        return values
    keep = [values[0]]
    for value in values[1:]:
        if value - keep[-1] >= min_interval_s:
            keep.append(value)
    return np.asarray(keep)


def prepare_continuous_signal(x, sf, missing_pad_s=MISSING_PAD_S):
    """Fill only for numerical filtering and return a conservative measured-data mask.

    Filters cannot consume NaNs, but replacing missing samples by zero and then forgetting the
    original mask turns acquisition gaps into plausible low-power physiology.  Interpolation here
    is only a computational scaffold: every originally missing sample and a safety margin around it
    remains ineligible for power, staging, ECG peaks, and event detection.
    """
    raw = np.asarray(x, float)
    finite = np.isfinite(raw)
    if finite.sum() >= 2:
        sample = np.arange(len(raw))
        filled = np.interp(sample, sample[finite], raw[finite])
    elif finite.any():
        filled = np.full(len(raw), float(raw[finite][0]))
    else:
        filled = np.zeros(len(raw), float)
    if missing_pad_s > 0 and (~finite).any():
        width = max(1, int(round(2 * missing_pad_s * sf)) + 1)
        contaminated = ndimage.maximum_filter1d(
            (~finite).astype(np.uint8), size=width, mode="nearest").astype(bool)
        measured = ~contaminated
    else:
        measured = finite.copy()
    return filled, measured


def empty_channel_activity_extrema(n_channels):
    """Initialize streaming raw-signal extrema for numerical flat-line rejection."""
    n_channels = int(n_channels)
    if n_channels < 0:
        raise ValueError("n_channels must be nonnegative")
    return dict(
        minimum=np.full(n_channels, np.inf),
        maximum=np.full(n_channels, -np.inf),
        maximum_absolute=np.zeros(n_channels),
        finite_count=np.zeros(n_channels, dtype=np.int64),
    )


def update_channel_activity_extrema(state, channel_by_sample):
    """Update raw per-channel extrema from one non-overlapping core data block."""
    values = np.asarray(channel_by_sample, float)
    if values.ndim != 2:
        raise ValueError("channel activity values must be channel-by-sample")
    n_channels = values.shape[0]
    required = {"minimum", "maximum", "maximum_absolute", "finite_count"}
    if set(state) != required or any(
            np.asarray(state[key]).shape != (n_channels,) for key in required):
        raise ValueError("channel activity state does not align with the data")
    for channel in range(n_channels):
        finite = values[channel][np.isfinite(values[channel])]
        if not len(finite):
            continue
        state["minimum"][channel] = min(
            state["minimum"][channel], float(np.min(finite)))
        state["maximum"][channel] = max(
            state["maximum"][channel], float(np.max(finite)))
        state["maximum_absolute"][channel] = max(
            state["maximum_absolute"][channel], float(np.max(np.abs(finite))))
        state["finite_count"][channel] += len(finite)
    return state


def finalize_channel_activity_qc(state):
    """Reject only numerical flat lines, without imposing a physiological amplitude gate."""
    minimum = np.asarray(state["minimum"], float)
    maximum = np.asarray(state["maximum"], float)
    maximum_absolute = np.asarray(state["maximum_absolute"], float)
    finite_count = np.asarray(state["finite_count"], np.int64)
    if not (
        minimum.shape == maximum.shape == maximum_absolute.shape == finite_count.shape
    ):
        raise ValueError("channel activity extrema arrays must align")
    dynamic_range = maximum - minimum
    tolerance = (
        SIGNAL_FLAT_EPSILON_MULTIPLIER
        * np.finfo(float).eps
        * np.maximum(1.0, maximum_absolute)
    )
    nonflat = (
        (finite_count >= 2)
        & np.isfinite(dynamic_range)
        & (dynamic_range > tolerance)
    )
    return dict(
        minimum=minimum,
        maximum=maximum,
        dynamic_range=dynamic_range,
        numerical_flat_tolerance=tolerance,
        finite_count=finite_count,
        nonflat_mask=nonflat,
        method=(
            "raw full-interval range must exceed 64 float64 eps times "
            "max(1, maximum absolute raw value)"
        ),
    )


def staging_epoch_features(segment, clean_mask, measured_mask, sf, *, return_details=False):
    """Gap-aware delta-ratio/SWA from complete artifact-free 4-s Welch windows.

    A contaminated sample never enters a retained periodogram.  Unlike the former fail-closed
    implementation, a brief masked event does not erase clean, temporally separate pieces of the
    entire 30-s epoch.  With a fully clean 30-s epoch, averaging the fourteen fixed 4-s Hann
    periodograms is numerically the same construction as ``scipy.signal.welch`` with its 50%
    overlap.  The neutral cache stores per-window support; an offline QC profile decides how many
    valid windows are adequate rather than hard-coding a coverage percentage here.
    """
    segment = np.asarray(segment, float)
    clean_mask = np.asarray(clean_mask, bool)
    measured_mask = np.asarray(measured_mask, bool)
    nperseg = int(round(STAGING_WELCH_WINDOW_S * float(sf)))
    noverlap = int(round(STAGING_WELCH_OVERLAP_S * float(sf)))
    step = nperseg - noverlap
    n_windows = (
        0 if len(segment) < nperseg
        else 1 + (len(segment) - nperseg) // step
    )
    empty_details = dict(
        n_valid_windows=0,
        n_total_windows=int(n_windows),
        valid_window_mask=np.zeros(n_windows, bool),
        window_swa_power=np.full(n_windows, np.nan),
        window_total_power=np.full(n_windows, np.nan),
        measured_fraction=0.0,
        clean_fraction=0.0,
        longest_valid_run_s=0.0,
        valid_window_span_s=0.0,
    )
    if (
        len(segment) == 0
        or clean_mask.shape != segment.shape
        or measured_mask.shape != segment.shape
        or nperseg < 2
        or step < 1
        or n_windows < 1
    ):
        return (np.nan, np.nan, empty_details) if return_details else (np.nan, np.nan)

    finite = np.isfinite(segment)
    measured = measured_mask & finite
    valid_samples = clean_mask & measured
    starts = np.arange(n_windows, dtype=int) * step
    window_valid = np.asarray([
        bool(valid_samples[start:start + nperseg].all())
        for start in starts
    ])

    # Nonfinite values occur only in windows that are ineligible.  Zero is numerical scaffolding
    # for the vectorized FFT and is never allowed into a retained periodogram.
    numerical_segment = np.where(finite, segment, 0.0)
    fq, _, pxx = signal.spectrogram(
        numerical_segment,
        fs=float(sf),
        window="hann",
        nperseg=nperseg,
        noverlap=noverlap,
        detrend="constant",
        scaling="density",
        mode="psd",
    )
    if pxx.shape[1] != n_windows:
        raise RuntimeError(
            f"staging Welch geometry produced {pxx.shape[1]} windows; expected {n_windows}"
        )
    swa_band = (fq >= SWA_BAND[0]) & (fq < SWA_BAND[1])
    total_band = (fq >= 0.5) & (fq < 25.0)
    window_swa = np.trapezoid(pxx[swa_band], fq[swa_band], axis=0)
    window_total = np.trapezoid(pxx[total_band], fq[total_band], axis=0)
    window_swa = np.where(window_valid, window_swa, np.nan)
    window_total = np.where(window_valid, window_total, np.nan)

    n_valid = int(window_valid.sum())
    swa = (
        float(np.nanmean(window_swa))
        if n_valid and np.isfinite(window_swa).any()
        else np.nan
    )
    total = (
        float(np.nanmean(window_total))
        if n_valid and np.isfinite(window_total).any()
        else np.nan
    )
    dr = float(swa / total) if np.isfinite(swa) and np.isfinite(total) and total > 0 else np.nan
    if not np.isfinite(dr) or not np.isfinite(swa):
        dr = swa = np.nan

    longest = 0
    if valid_samples.any():
        longest = max((stop - start for start, stop in contiguous_runs(valid_samples)), default=0)
    valid_starts = starts[window_valid]
    span = (
        float((valid_starts[-1] - valid_starts[0] + nperseg) / sf)
        if len(valid_starts)
        else 0.0
    )
    details = dict(
        n_valid_windows=n_valid,
        n_total_windows=int(n_windows),
        valid_window_mask=window_valid,
        window_swa_power=window_swa,
        window_total_power=window_total,
        measured_fraction=float(measured.mean()),
        clean_fraction=float(valid_samples.mean()),
        longest_valid_run_s=float(longest / sf),
        valid_window_span_s=span,
    )
    return (dr, swa, details) if return_details else (dr, swa)


def aggregate_staging_features(
        dr_by_contact, swa_by_contact, clean_by_contact, candidate_contact_mask,
        min_contacts=MIN_CONTACTS,
        min_contact_feature_coverage=MIN_CONTACT_COVERAGE,
        min_contact_fraction_per_epoch=MIN_CONTACT_FRACTION_PER_BIN,
        valid_window_count_by_contact=None,
        min_valid_windows=1):
    """Aggregate staging features from one fixed, normalized full-night contact set.

    Changing which contacts contribute from epoch to epoch can turn stable contact gain differences
    into an apparent two-state SWA distribution.  Candidate contacts are fixed before aggregation;
    SWA is normalized within contact over the whole available night; and every reported epoch must
    contain the profile-defined fraction of that exact set.
    """
    dr = np.asarray(dr_by_contact, float)
    swa = np.asarray(swa_by_contact, float)
    clean = np.asarray(clean_by_contact, float)
    if dr.ndim != 2 or swa.shape != dr.shape or clean.shape != dr.shape:
        raise ValueError("staging feature arrays must be aligned contact-by-epoch matrices")
    selected = np.asarray(candidate_contact_mask, bool).ravel().copy()
    if len(selected) != dr.shape[0]:
        raise ValueError("candidate_contact_mask must align with staging contacts")
    feature_valid = np.isfinite(dr) & np.isfinite(swa)
    if valid_window_count_by_contact is not None:
        window_count = np.asarray(valid_window_count_by_contact)
        if window_count.shape != dr.shape:
            raise ValueError("valid-window support must align with staging feature matrices")
        feature_valid &= window_count >= int(min_valid_windows)
        dr = np.where(feature_valid, dr, np.nan)
        swa = np.where(feature_valid, swa, np.nan)
    feature_coverage = feature_valid.mean(axis=1)
    # A contact that is stable for the 1-s sigma series can still be nearly absent from the
    # stricter, fully-clean 30-s staging features.  Letting one finite epoch qualify such a contact
    # raises the fixed-set denominator and can make every epoch unavailable.  Prequalify once over
    # the whole recording; do not change the set epoch by epoch.
    selected &= feature_coverage >= float(min_contact_feature_coverage)
    swa_scales = np.full(dr.shape[0], np.nan)
    for contact in range(dr.shape[0]):
        positive = swa[contact][np.isfinite(swa[contact]) & (swa[contact] > 0)]
        if len(positive):
            swa_scales[contact] = float(np.median(positive))
    selected &= np.isfinite(swa_scales) & (swa_scales > 0)
    n_selected = int(selected.sum())
    required = max(
        int(min_contacts),
        int(np.ceil(float(min_contact_fraction_per_epoch) * n_selected)),
    )
    ep_dr = np.full(dr.shape[1], np.nan)
    ep_swa = np.full(dr.shape[1], np.nan)
    ep_clean = np.full(dr.shape[1], np.nan)
    support = np.zeros(dr.shape[1], int)
    if n_selected >= int(min_contacts):
        normalized_swa = swa[selected] / swa_scales[selected, None]
        selected_dr = dr[selected]
        valid = np.isfinite(selected_dr) & np.isfinite(normalized_swa)
        support = valid.sum(axis=0)
        for epoch in np.where(support >= required)[0]:
            keep = valid[:, epoch]
            ep_dr[epoch] = float(np.median(selected_dr[keep, epoch]))
            ep_swa[epoch] = float(np.median(normalized_swa[keep, epoch]))
        ep_clean = np.nanmedian(clean[selected], axis=0)
    return ep_dr, ep_swa, ep_clean, dict(
        selected_contact_mask=selected,
        n_selected_contacts=n_selected,
        required_contact_count=required,
        contact_count=support,
        per_contact_feature_coverage=feature_coverage,
        minimum_contact_feature_coverage=float(min_contact_feature_coverage),
        minimum_valid_welch_windows=int(min_valid_windows),
        swa_normalization="divide each fixed contact by its full-night median clean-epoch SWA",
    )


def interpolate_tachograms(beats, total_s, fs_rr=FS_RR, max_gap_s=5.0):
    """Return 1-Hz/``fs_rr`` RR and HR series without bridging long beat gaps.

    Long gaps are explicitly masked *including their endpoint grid samples*.  Without that final
    mask, a 5.1-s gap between valid RR knots can leave only 5.0 s of internal NaNs, which a later
    "fill gaps up to 5 s" operation can silently reopen.
    """
    total_s = int(total_s)
    rr_1 = np.full(total_s, np.nan)
    rr_f = np.full(int(total_s * fs_rr), np.nan)
    hr_1 = np.full(total_s, np.nan)
    hr_f = np.full(int(total_s * fs_rr), np.nan)
    beats = sanitize_beats(beats)
    if len(beats) <= 20:
        return rr_1, rr_f, hr_1, hr_f
    intervals = np.diff(beats)
    good = (intervals >= 0.33) & (intervals <= 1.5)
    if good.sum() <= 20:
        return rr_1, rr_f, hr_1, hr_f
    tb, rrv = beats[1:][good], intervals[good]
    g1 = np.arange(total_s, dtype=float)
    gf = np.arange(0, total_s, 1.0 / fs_rr)
    long_gap_idx = np.where(np.diff(tb) > max_gap_s)[0]
    split = long_gap_idx + 1
    for run in np.split(np.arange(len(tb)), split):
        if len(run) < 4:
            continue
        m1 = (g1 >= tb[run[0]]) & (g1 <= tb[run[-1]])
        mf = (gf >= tb[run[0]]) & (gf <= tb[run[-1]])
        interp = interpolate.PchipInterpolator(tb[run], rrv[run], extrapolate=False)
        rr_1[m1] = interp(g1[m1])
        rr_f[mf] = interp(gf[mf])
    # Make the discontinuity impossible to reinterpret as a short fillable gap downstream.
    for i in long_gap_idx:
        rr_1[(g1 >= tb[i]) & (g1 <= tb[i + 1])] = np.nan
        rr_f[(gf >= tb[i]) & (gf <= tb[i + 1])] = np.nan
    with np.errstate(divide="ignore", invalid="ignore"):
        hr_1 = 60.0 / rr_1
        hr_f = 60.0 / rr_f
    return rr_1, rr_f, hr_1, hr_f


def power_from_binned_support(
        numerator, clean_sample_count, samples_per_second,
        minimum_clean_fraction_per_second=0.5):
    """Materialize per-second power at an explicit clean-sample support threshold.

    ``numerator`` is the sum of squared-envelope samples admitted by the clean mask and
    ``clean_sample_count`` is its denominator.  Keeping both quantities in the neutral cache makes
    the clean-fraction choice reversible: changing the threshold does not require filtering the
    raw recording again.  A zero-support bin is always unavailable, including when the requested
    minimum fraction is zero.
    """
    numerator = np.asarray(numerator, float)
    denominator = np.asarray(clean_sample_count, float)
    if numerator.shape != denominator.shape:
        raise ValueError("power numerator and clean-sample denominator must align")
    samples_per_second = int(samples_per_second)
    if samples_per_second < 1:
        raise ValueError("samples_per_second must be positive")
    minimum = float(minimum_clean_fraction_per_second)
    if not 0.0 <= minimum <= 1.0:
        raise ValueError("minimum_clean_fraction_per_second must lie in [0, 1]")
    eligible = (
        np.isfinite(numerator)
        & np.isfinite(denominator)
        & (denominator > 0)
        & (denominator >= minimum * samples_per_second)
    )
    return np.divide(
        numerator, denominator, out=np.full(numerator.shape, np.nan), where=eligible)


def _binned_power_values(env2, clean, sf, n_sec, *, return_support=False):
    """IED-masked mean power plus its reversible numerator/denominator support."""
    k = require_integer_sample_rate(sf, source="power-bin input")
    e2 = env2[:n_sec * k].reshape(n_sec, k)
    cm = clean[:n_sec * k].reshape(n_sec, k).astype(float)
    num, den = (e2 * cm).sum(1), cm.sum(1)
    values = power_from_binned_support(num, den, k, 0.5)
    return (values, num, den) if return_support else values


def _write_multichannel_power(x_ctx, sos, sf, n_sec, total_s, off, dest,
                              core_start_sample=0, numerator_dest=None,
                              clean_sample_count_dest=None):
    """Compute power per contact and write it to a full-night channel-by-time array.

    Averaging raw voltages before Hilbert power can cancel spindles with different phase/polarity.
    Normalization is deliberately deferred until the complete night is present: normalizing every
    streamed chunk separately would erase genuine dynamics slower than the chunk length and create
    artificial scale steps at chunk boundaries.
    """
    for j in range(x_ctx.shape[1]):
        filled, measured = prepare_continuous_signal(x_ctx[:, j], sf)
        x = notch(signal.detrend(filled), sf)
        clean = ied_clean_mask(x, sf) & measured
        env2 = np.abs(signal.hilbert(signal.sosfiltfilt(sos, x))) ** 2
        start = int(core_start_sample)
        stop = start + int(n_sec * sf)
        vals, numerator, denominator = _binned_power_values(
            env2[start:stop], clean[start:stop], sf, n_sec, return_support=True)
        sl = slice(off, min(off + n_sec, total_s))
        length = sl.stop - sl.start
        dest[j, sl] = vals[:length]
        if numerator_dest is not None:
            numerator_dest[j, sl] = numerator[:length]
        if clean_sample_count_dest is not None:
            clean_sample_count_dest[j, sl] = denominator[:length]


def _aggregate_full_night_power(
        per_channel, *, eligible_channels=None, min_contact_coverage=0.0,
        min_contacts=1, min_contact_fraction_per_bin=0.0, return_details=False):
    """Normalize once and aggregate a coverage-qualified, stable contact set.

    Without per-contact qualification, mutually disjoint contacts can yield a superficially 100%
    covered aggregate whose anatomical/contact identity changes completely over time.
    """
    per_channel = np.asarray(per_channel, float)
    scales = np.nanmedian(per_channel, axis=1)
    contact_coverage = np.isfinite(per_channel).mean(axis=1)
    valid = (
        np.isfinite(scales) & (scales > 0)
        & (contact_coverage >= float(min_contact_coverage))
    )
    if eligible_channels is not None:
        valid &= np.asarray(eligible_channels, bool)
    selected_n = int(valid.sum())
    required_count = max(
        1, int(np.ceil(float(min_contact_fraction_per_bin) * selected_n)))
    if selected_n < int(min_contacts):
        result = np.full(per_channel.shape[1], np.nan)
        count = np.zeros(per_channel.shape[1], int)
        details = dict(
            per_contact_coverage=contact_coverage,
            selected_mask=valid,
            contact_count=count,
            n_selected=selected_n,
            required_contact_count=required_count)
        return (result, details) if return_details else result
    normalized = per_channel[valid] / scales[valid, None]
    finite = np.isfinite(normalized)
    count = finite.sum(axis=0)
    total = np.nansum(normalized, axis=0)
    result = np.divide(
        total, count, out=np.full(normalized.shape[1], np.nan),
        where=count >= required_count)
    details = dict(
        per_contact_coverage=contact_coverage,
        selected_mask=valid,
        contact_count=count,
        n_selected=selected_n,
        required_contact_count=required_count)
    return (result, details) if return_details else result


def detect_so_candidates(x, sf, edge_s=5.0):
    """Full down/up-state SO candidates before the channel-night amplitude threshold.

    Thresholding here used to happen independently in every streamed chunk even though the method
    claimed a within-channel percentile.  Return amplitudes so the caller can threshold the entire
    channel-night distribution after streaming.
    """
    x = np.asarray(x, float)
    zc = np.where(np.diff(np.signbit(x)))[0]
    if len(zc) < 3:
        return []
    candidates = []
    edge = int(round(edge_s * sf))
    for a, b, c in zip(zc[:-2], zc[1:-1], zc[2:]):
        if x[a + 1] >= 0:
            continue
        negative_start = a + 1
        negative_stop = b + 1
        positive_start = b + 1
        positive_stop = c + 1
        negative_half = x[negative_start:negative_stop]
        positive_half = x[positive_start:positive_stop]
        if (
            not len(negative_half)
            or not len(positive_half)
            or not np.isfinite(negative_half).all()
            or not np.isfinite(positive_half).all()
            or not np.all(negative_half < 0)
            or not np.all(positive_half >= 0)
        ):
            continue
        negative_half_s = (b - a) / float(sf)
        positive_half_s = (c - b) / float(sf)
        if not (
            SO_NEGATIVE_HALF_DURATION_S[0]
            <= negative_half_s
            <= SO_NEGATIVE_HALF_DURATION_S[1]
            and 0 < positive_half_s < SO_POSITIVE_HALF_MAX_S
        ):
            continue
        tr = negative_start + int(np.argmin(negative_half))
        if tr < edge or tr >= len(x) - edge:
            continue
        up = float(np.max(positive_half))
        down = float(-x[tr])
        if not (down > 0 and up > 0):
            continue
        candidates.append((tr, down, up, down + up))
    return candidates


def threshold_so_candidates(candidates, percentile=AMP_PCT):
    """Apply one amplitude threshold to a complete channel-night candidate collection."""
    if not candidates:
        return np.array([], float)
    a = np.asarray(candidates, float)
    keep = ((a[:, 2] >= np.percentile(a[:, 2], percentile))
            & (a[:, 3] >= np.percentile(a[:, 3], percentile)))
    return np.sort(a[keep, 0])


def detect_so_halfwaves(x, sf):
    """Compatibility wrapper for a single continuous signal; returns trough sample indices."""
    candidates = detect_so_candidates(x, sf, edge_s=0)
    return threshold_so_candidates(candidates).astype(int)


def validated_analysis_hours(hours):
    """Require the whole-second duration assumed by cached time axes."""
    hours = float(hours)
    seconds = hours * 3600.0
    if (
        not np.isfinite(seconds)
        or seconds <= 0
        or abs(seconds - round(seconds)) > 1e-7
    ):
        raise ValueError(
            "analysis hours must describe a positive whole number of seconds")
    return hours, int(round(seconds))


@dataclass(frozen=True)
class _HupSource:
    name: str
    fp: str
    hours: float
    total_s: int
    current_cache_digest: str
    ds: object
    lab_idx: dict
    ekg: str
    ctx: list
    source_identity_json: str
    source_selection_json: str
    geometry: dict
    sf: float
    night: float
    idx: list
    fsp: float
    fsp_real: bool


@dataclass
class _CacheBuffers:
    sig_fixed_ch: np.ndarray
    sig_fsp_ch: np.ndarray
    swa_ch: np.ndarray
    sig_fixed_num_ch: np.ndarray
    sig_fsp_num_ch: np.ndarray
    swa_num_ch: np.ndarray
    sig_fixed_den_ch: np.ndarray
    sig_fsp_den_ch: np.ndarray
    swa_den_ch: np.ndarray
    ep_dr_ch: np.ndarray
    ep_swa_ch: np.ndarray
    ep_clean_ch: np.ndarray
    ep_measured_ch: np.ndarray
    ep_valid_window_ch: np.ndarray
    ep_window_swa_ch: np.ndarray
    ep_window_total_ch: np.ndarray
    ep_longest_clean_run_ch: np.ndarray
    ep_valid_window_span_ch: np.ndarray
    event_3d_rms_ch: np.ndarray
    event_3d_phase_ch: np.ndarray
    event_3d_valid_ch: np.ndarray
    event_3d_candidates: list
    beats: list
    so_candidates: dict
    failed_chunks: list
    ecg_failures: list
    acquisition_sample_counts: list
    channel_activity_state: dict


@dataclass(frozen=True)
class _FilterBank:
    fixed: np.ndarray
    fsp: np.ndarray
    swa: np.ndarray
    so: np.ndarray
    event_3d_so: np.ndarray
    event_3d_spindle: np.ndarray


@dataclass(frozen=True)
class _ChunkWindow:
    t: float
    duration_s: float
    pull_start_s: float
    core_start_sample: int
    core_sample_count: int
    core_stop_sample: int
    output_second: int


@dataclass(frozen=True)
class _PulledChunk:
    window: _ChunkWindow
    data: np.ndarray


@dataclass(frozen=True)
class _PreparedChunk:
    window: _ChunkWindow
    cortical_raw: np.ndarray
    ecg_raw: np.ndarray
    cortical_filtered: np.ndarray
    measured: np.ndarray
    clean: np.ndarray


@dataclass(frozen=True)
class _FinalizedSeries:
    channel_activity_qc: dict
    nonflat_contacts: np.ndarray
    sig_fixed: np.ndarray
    sigma_contact_qc: dict
    ep_dr: np.ndarray
    ep_swa: np.ndarray
    ep_clean: np.ndarray
    staging_contact_qc: dict
    sig_fsp: np.ndarray
    swa_1: np.ndarray
    beats: np.ndarray
    rr_1: np.ndarray
    rr_4: np.ndarray
    hr_1: np.ndarray
    hr_4: np.ndarray
    sigma_coverage: float
    hr_coverage: float
    qc_warnings: list


def _setup_hup_source(
        name, fp, hours, total_s, current_cache_digest, session):
    """Resolve the pinned source interval or write the historical skip cache."""
    ds = session.open_dataset(name)
    labels = list(ds.get_channel_labels())
    lab_idx = {label: index for index, label in enumerate(labels)}
    ekg = next((l for l in labels if l.upper().startswith(("EKG", "ECG"))), None)
    ctx = cortical_channels(labels)
    source_identity = verify_hup_source_identity(ds, ctx, ekg)
    source_identity_json = json.dumps(source_identity, sort_keys=True)
    if ekg is None or len(ctx) == 0:
        atomic_savez(fp, subject=name, status="skip",
                     cache_schema_version=CACHE_SCHEMA_VERSION,
                     cache_code_sha256=current_cache_digest,
                     source_dataset=name, source_kind="iEEG.org API",
                     source_identity_json=source_identity_json,
                     hours=float(hours),
                     reason=f"ekg={ekg} n_cortical={len(ctx)}")
        print(f"[{name}] SKIP ekg={ekg} n_cortical={len(ctx)}", flush=True)
        return None

    geometry = validate_hup_series_geometry(source_identity, ctx, ekg)
    sf = float(require_integer_sample_rate(
        geometry["sample_rate_hz"], source=f"{name} shared portal geometry"))
    total_h = geometry["duration_us"] / 3.6e9
    night, night_search_qc = find_night(
        ds, lab_idx, ctx[0], sf, total_h, required_h=hours,
        return_diagnostics=True)
    source_selection_json = json.dumps(
        dict(
            night_s=None if night is None else float(night),
            requested_hours=float(hours),
            night_probe_workers=NIGHT_PROBE_WORKERS,
            portal_timeout_s=dict(
                connect=IEEG_CONNECT_TIMEOUT_S,
                read=IEEG_READ_TIMEOUT_S),
            night_search_qc=night_search_qc,
            night_search=(
                "highest sparse-probe delta candidate (6 s sampled every 30 min "
                "within each contiguous 3 h window across the complete recording; "
                "not validated sleep onset or a dense whole-window mean)"),
        ),
        sort_keys=True)
    if night is None:
        atomic_savez(fp, subject=name, status="skip",
                     cache_schema_version=CACHE_SCHEMA_VERSION,
                     cache_code_sha256=current_cache_digest,
                     source_dataset=name, source_kind="iEEG.org API",
                     source_identity_json=source_identity_json,
                     source_selection_json=source_selection_json,
                     hours=float(hours),
                     reason="no night")
        print(f"[{name}] SKIP no night", flush=True)
        return None

    idx = [lab_idx[c] for c in ctx] + [lab_idx[ekg]]
    # Lecci visually selected FSP from all artifact-free NREM spectra. A single arbitrary 300-s
    # sample cannot reproduce that and labels 1/f noise as a peak, so individualized FSP is disabled
    # until an all-NREM, manually QC'd selection is supplied.
    fsp, fsp_real = 13.0, False
    print(f"[{name}] {sf:.0f} Hz | {len(ctx)} lateral-contact candidates | "
          f"fixed sigma primary; individualized FSP disabled | streaming {hours} h", flush=True)
    return _HupSource(
        name=name,
        fp=fp,
        hours=hours,
        total_s=total_s,
        current_cache_digest=current_cache_digest,
        ds=ds,
        lab_idx=lab_idx,
        ekg=ekg,
        ctx=ctx,
        source_identity_json=source_identity_json,
        source_selection_json=source_selection_json,
        geometry=geometry,
        sf=sf,
        night=night,
        idx=idx,
        fsp=fsp,
        fsp_real=fsp_real,
    )


def _allocate_cache_buffers(ctx, total_s):
    """Allocate the complete neutral cache workspace before streaming."""
    n_contacts = len(ctx)
    n_ep = int(total_s // EPOCH)
    n_event_3d = int(np.ceil(total_s * EVENT_3D_SAMPLING_HZ))
    return _CacheBuffers(
        sig_fixed_ch=np.full((n_contacts, total_s), np.nan),
        sig_fsp_ch=np.full((n_contacts, total_s), np.nan),
        swa_ch=np.full((n_contacts, total_s), np.nan),
        sig_fixed_num_ch=np.full((n_contacts, total_s), np.nan),
        sig_fsp_num_ch=np.full((n_contacts, total_s), np.nan),
        swa_num_ch=np.full((n_contacts, total_s), np.nan),
        sig_fixed_den_ch=np.zeros((n_contacts, total_s), dtype=np.uint32),
        sig_fsp_den_ch=np.zeros((n_contacts, total_s), dtype=np.uint32),
        swa_den_ch=np.zeros((n_contacts, total_s), dtype=np.uint32),
        ep_dr_ch=np.full((n_contacts, n_ep), np.nan),
        ep_swa_ch=np.full((n_contacts, n_ep), np.nan),
        ep_clean_ch=np.full((n_contacts, n_ep), np.nan),
        ep_measured_ch=np.full((n_contacts, n_ep), np.nan),
        ep_valid_window_ch=np.zeros(
            (n_contacts, n_ep, STAGING_REFERENCE_MIN_VALID_WINDOWS),
            dtype=np.uint8),
        ep_window_swa_ch=np.full(
            (n_contacts, n_ep, STAGING_REFERENCE_MIN_VALID_WINDOWS),
            np.nan, dtype=np.float32),
        ep_window_total_ch=np.full(
            (n_contacts, n_ep, STAGING_REFERENCE_MIN_VALID_WINDOWS),
            np.nan, dtype=np.float32),
        ep_longest_clean_run_ch=np.zeros(
            (n_contacts, n_ep), dtype=np.float32),
        ep_valid_window_span_ch=np.zeros(
            (n_contacts, n_ep), dtype=np.float32),
        event_3d_rms_ch=np.full(
            (n_contacts, n_event_3d), np.nan, dtype=np.float32),
        event_3d_phase_ch=np.full(
            (n_contacts, n_event_3d), np.nan, dtype=np.float32),
        event_3d_valid_ch=np.zeros(
            (n_contacts, n_event_3d), dtype=np.uint8),
        event_3d_candidates=[],
        beats=[],
        so_candidates={c: [] for c in ctx},
        failed_chunks=[],
        ecg_failures=[],
        acquisition_sample_counts=[],
        channel_activity_state=empty_channel_activity_extrema(n_contacts),
    )


def _build_filter_bank(sf, fsp):
    """Construct the fixed filters once for the complete source interval."""
    return _FilterBank(
        fixed=band_sos(SIGMA_FIXED, sf),
        fsp=band_sos((fsp - 1, fsp + 1), sf),
        swa=band_sos(SWA_BAND_L, sf, 3),
        so=signal.butter(
            3, list(SO_BAND_NAJI), btype="band", fs=sf, output="sos"),
        event_3d_so=band_sos(EVENT_3D_SO_BAND, sf, 3),
        event_3d_spindle=band_sos(EVENT_3D_SPINDLE_BAND, sf),
    )


def _pull_analysis_chunk(source, buffers, t, dur):
    """Acquire one padded chunk and record the exact non-overlapping core."""
    pull_start = max(0.0, t - FILTER_EDGE_S)
    pull_stop = min(float(source.total_s), t + dur + FILTER_EDGE_S)
    pull_dur = pull_stop - pull_start
    try:
        data = pull_continuous_exact(
            source.ds,
            source.idx,
            source.night + pull_start,
            pull_dur,
            source.sf,
            records=buffers.acquisition_sample_counts,
            purpose="analysis_subrequest",
        )
        core_a = expected_portal_sample_count(t - pull_start, source.sf)
        requested_core = expected_portal_sample_count(dur, source.sf)
        returned_core = max(
            0, min(len(data), core_a + requested_core) - core_a)
        core_record = {
            "purpose": "analysis_core",
            "analysis_start_s": float(t),
            "analysis_duration_s": float(dur),
            "requested_sample_count": requested_core,
            "returned_sample_count": int(returned_core),
            "status": (
                "ok" if returned_core == requested_core
                else "sample_count_mismatch"
            ),
        }
        buffers.acquisition_sample_counts.append(core_record)
        if returned_core != requested_core:
            raise RuntimeError(
                f"portal returned {returned_core} core samples; expected "
                f"{requested_core} at analysis t={t:g} s")
    except Exception as exc:
        buffers.failed_chunks.append(
            dict(start_s=float(t), duration_s=float(dur),
                 error=f"{type(exc).__name__}: {exc}"))
        return None

    window = _ChunkWindow(
        t=t,
        duration_s=dur,
        pull_start_s=pull_start,
        core_start_sample=core_a,
        core_sample_count=requested_core,
        core_stop_sample=core_a + requested_core,
        output_second=int(t),
    )
    return _PulledChunk(window=window, data=data)


def _prepare_analysis_chunk(source, buffers, pulled):
    """Prepare numerical channel signals while retaining measured/clean masks."""
    window = pulled.window
    x_ctx = pulled.data[:, :len(source.ctx)]
    x_ekg = pulled.data[:, len(source.ctx)]
    update_channel_activity_extrema(
        buffers.channel_activity_state,
        x_ctx[window.core_start_sample:window.core_stop_sample].T)
    prepared = [
        prepare_continuous_signal(x_ctx[:, ci], source.sf)
        for ci in range(len(source.ctx))
    ]
    x_channels = np.asarray([
        notch(signal.detrend(filled), source.sf)
        for filled, _ in prepared
    ])
    measured_channels = np.asarray([measured for _, measured in prepared])
    clean_channels = np.asarray([
        ied_clean_mask(x, source.sf) & measured
        for x, measured in zip(x_channels, measured_channels)
    ])
    return _PreparedChunk(
        window=window,
        cortical_raw=x_ctx,
        ecg_raw=x_ekg,
        cortical_filtered=x_channels,
        measured=measured_channels,
        clean=clean_channels,
    )


def _detect_chunk_ecg(source, buffers, chunk):
    """Detect core ECG peaks, retaining the historical per-chunk failure record."""
    window = chunk.window
    try:
        ecg_filled, ecg_measured = prepare_continuous_signal(
            chunk.ecg_raw, source.sf)
        cl = nk.ecg_clean(
            ecg_filled, sampling_rate=int(source.sf), method="neurokit")
        _, info = nk.ecg_peaks(
            cl, sampling_rate=int(source.sf), method="neurokit",
            correct_artifacts=True)
        peaks = np.asarray(info["ECG_R_Peaks"], int)
        peaks = peaks[
            (peaks >= window.core_start_sample)
            & (peaks < window.core_stop_sample)
        ]
        peaks = peaks[ecg_measured[peaks]]
        buffers.beats.extend(
            (peaks / source.sf + window.pull_start_s).tolist())
    except Exception as exc:
        buffers.ecg_failures.append(
            dict(start_s=float(window.t),
                 duration_s=float(window.duration_s),
                 error=f"{type(exc).__name__}: {exc}"))


def _write_chunk_event_3d(source, buffers, filters, chunk, ci, xc):
    """Write outcome-neutral 3D samples and complete pre-threshold SO cycles."""
    window = chunk.window
    event_clean = (
        ied_clean_mask(xc, source.sf, pad_s=EVENT_3D_IED_PAD_S)
        & chunk.measured[ci]
    )
    event_so = signal.sosfiltfilt(filters.event_3d_so, xc)
    event_phase = np.angle(signal.hilbert(event_so))
    event_spindle = signal.sosfiltfilt(filters.event_3d_spindle, xc)
    event_rms = event_3d_spindle_rms(event_spindle, source.sf)
    candidates_3d = event_3d_so_candidates(event_so, source.sf)
    clean_candidates_3d = retain_complete_clean_event_3d_candidates(
        candidates_3d,
        event_clean,
        window.core_start_sample,
        window.core_stop_sample,
    )
    if len(clean_candidates_3d):
        for trough, amplitude, cycle_start, cycle_stop in clean_candidates_3d:
            event_sample = int(round(
                (window.pull_start_s + float(trough) / source.sf)
                * EVENT_3D_SAMPLING_HZ
            ))
            (
                cycle_start_sample,
                cycle_stop_sample_exclusive,
            ) = conservative_resampled_interval(
                cycle_start,
                cycle_stop,
                source.sf,
                EVENT_3D_SAMPLING_HZ,
                offset_s=window.pull_start_s,
            )
            if (
                0 <= cycle_start_sample
                <= event_sample
                < cycle_stop_sample_exclusive
                <= buffers.event_3d_rms_ch.shape[1]
            ):
                buffers.event_3d_candidates.append(
                    (
                        ci,
                        event_sample,
                        float(amplitude),
                        cycle_start_sample,
                        cycle_stop_sample_exclusive,
                    )
                )

    event_start = int(round(window.t * EVENT_3D_SAMPLING_HZ))
    n_event_out = min(
        int(round(window.duration_s * EVENT_3D_SAMPLING_HZ)),
        buffers.event_3d_rms_ch.shape[1] - event_start,
    )
    local_samples = (
        window.core_start_sample
        + np.round(
            np.arange(n_event_out) * source.sf / EVENT_3D_SAMPLING_HZ
        ).astype(int)
    )
    in_bounds = local_samples < len(xc)
    event_indices = event_start + np.arange(n_event_out)[in_bounds]
    local_samples = local_samples[in_bounds]
    valid_event = (
        event_clean[local_samples]
        & np.isfinite(event_rms[local_samples])
        & np.isfinite(event_phase[local_samples])
    )
    event_indices = event_indices[valid_event]
    local_samples = local_samples[valid_event]
    buffers.event_3d_rms_ch[ci, event_indices] = event_rms[local_samples]
    buffers.event_3d_phase_ch[ci, event_indices] = event_phase[local_samples]
    buffers.event_3d_valid_ch[ci, event_indices] = 1


def _write_chunk_slow_oscillations(source, buffers, filters, chunk):
    """Retain the existing 3B candidates and neutral 3D event sidecar."""
    window = chunk.window
    for ci, channel in enumerate(source.ctx):
        xc = chunk.cortical_filtered[ci]
        if np.std(xc) < 1e-9:
            continue
        candidates = detect_so_candidates(
            signal.sosfiltfilt(filters.so, xc), source.sf)
        so_clean = (
            ied_clean_mask(xc, source.sf, pad_s=5.0)
            & chunk.measured[ci]
        )
        candidates = [
            value for value in candidates
            if (
                window.core_start_sample
                <= int(value[0])
                < window.core_stop_sample
                and so_clean[int(value[0])]
            )
        ]
        buffers.so_candidates[channel].extend([
            (
                trough / source.sf + window.pull_start_s,
                down,
                up,
                peak_to_peak,
            )
            for trough, down, up, peak_to_peak in candidates
        ])

        # The 3D sidecar remains pre-threshold, pre-stage, and pre-pairing.
        _write_chunk_event_3d(
            source, buffers, filters, chunk, ci, xc)


def _write_chunk_power(source, buffers, filters, chunk):
    """Write reversible per-contact, per-second band-power support."""
    window = chunk.window
    n_sec = int(window.core_sample_count // int(source.sf))
    for sos_b, dest, numerator_dest, denominator_dest in (
            (filters.fixed, buffers.sig_fixed_ch,
             buffers.sig_fixed_num_ch, buffers.sig_fixed_den_ch),
            (filters.fsp, buffers.sig_fsp_ch,
             buffers.sig_fsp_num_ch, buffers.sig_fsp_den_ch),
            (filters.swa, buffers.swa_ch,
             buffers.swa_num_ch, buffers.swa_den_ch)):
        _write_multichannel_power(
            chunk.cortical_raw,
            sos_b,
            source.sf,
            n_sec,
            source.total_s,
            window.output_second,
            dest,
            core_start_sample=window.core_start_sample,
            numerator_dest=numerator_dest,
            clean_sample_count_dest=denominator_dest,
        )


def _write_chunk_staging(source, buffers, chunk):
    """Write neutral per-contact Welch-window staging features."""
    window = chunk.window
    n_ep = buffers.ep_dr_ch.shape[1]
    samples_per_epoch = int(EPOCH * source.sf)
    for epoch in range(int(window.core_sample_count // samples_per_epoch)):
        global_epoch = int(
            (window.output_second + epoch * EPOCH) // EPOCH)
        if global_epoch >= n_ep:
            break
        for ci, (xc, clean_c, measured_c) in enumerate(zip(
                chunk.cortical_filtered, chunk.clean, chunk.measured)):
            start = window.core_start_sample + epoch * samples_per_epoch
            stop = start + samples_per_epoch
            segment = xc[start:stop]
            dr, swa_value, details = staging_epoch_features(
                segment,
                clean_c[start:stop],
                measured_c[start:stop],
                source.sf,
                return_details=True,
            )
            buffers.ep_dr_ch[ci, global_epoch] = dr
            buffers.ep_swa_ch[ci, global_epoch] = swa_value
            buffers.ep_clean_ch[ci, global_epoch] = details["clean_fraction"]
            buffers.ep_measured_ch[ci, global_epoch] = details["measured_fraction"]
            buffers.ep_valid_window_ch[ci, global_epoch] = details["valid_window_mask"]
            buffers.ep_window_swa_ch[ci, global_epoch] = details["window_swa_power"]
            buffers.ep_window_total_ch[ci, global_epoch] = details["window_total_power"]
            buffers.ep_longest_clean_run_ch[ci, global_epoch] = (
                details["longest_valid_run_s"])
            buffers.ep_valid_window_span_ch[ci, global_epoch] = (
                details["valid_window_span_s"])


def _process_analysis_chunk(source, buffers, filters, pulled):
    """Run all estimators for one successfully acquired padded chunk."""
    chunk = _prepare_analysis_chunk(source, buffers, pulled)
    _detect_chunk_ecg(source, buffers, chunk)
    _write_chunk_slow_oscillations(source, buffers, filters, chunk)
    _write_chunk_power(source, buffers, filters, chunk)
    _write_chunk_staging(source, buffers, chunk)


def _finalize_cache_series(source, buffers):
    """Aggregate complete-night scientific arrays and enforce fatal QC."""
    channel_activity_qc = finalize_channel_activity_qc(
        buffers.channel_activity_state)
    nonflat_contacts = channel_activity_qc["nonflat_mask"]
    sig_fixed, sigma_contact_qc = _aggregate_full_night_power(
        buffers.sig_fixed_ch,
        eligible_channels=nonflat_contacts,
        min_contact_coverage=MIN_CONTACT_COVERAGE,
        min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN,
        return_details=True,
    )
    eligible_contacts = sigma_contact_qc["selected_mask"]
    ep_dr, ep_swa, ep_clean, staging_contact_qc = aggregate_staging_features(
        buffers.ep_dr_ch,
        buffers.ep_swa_ch,
        buffers.ep_clean_ch,
        eligible_contacts,
        valid_window_count_by_contact=buffers.ep_valid_window_ch.sum(axis=2),
        min_valid_windows=STAGING_REFERENCE_MIN_VALID_WINDOWS,
    )
    sig_fsp = _aggregate_full_night_power(
        buffers.sig_fsp_ch,
        eligible_channels=eligible_contacts,
        min_contact_coverage=MIN_CONTACT_COVERAGE,
        min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN,
    )
    swa_1 = _aggregate_full_night_power(
        buffers.swa_ch,
        eligible_channels=eligible_contacts,
        min_contact_coverage=MIN_CONTACT_COVERAGE,
        min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN,
    )

    # R-peaks at chunk edges can be repeated or spuriously near-duplicated.
    beats = sanitize_beats(buffers.beats)
    rr_1, rr_4, hr_1, hr_4 = interpolate_tachograms(
        beats, source.total_s)
    sigma_coverage = float(np.isfinite(sig_fixed).mean())
    hr_coverage = float(np.isfinite(hr_4).mean())
    if buffers.failed_chunks:
        raise RuntimeError(
            f"{len(buffers.failed_chunks)} acquisition chunks failed")
    qc_warnings = finalize_ecg_cache_qc(
        buffers.ecg_failures, hr_coverage)
    if len(source.ctx) < MIN_CONTACTS:
        qc_warnings.append(
            f"only {len(source.ctx)} lateral-contact candidates; historical audit80 required "
            f"{MIN_CONTACTS}")

    return _FinalizedSeries(
        channel_activity_qc=channel_activity_qc,
        nonflat_contacts=nonflat_contacts,
        sig_fixed=sig_fixed,
        sigma_contact_qc=sigma_contact_qc,
        ep_dr=ep_dr,
        ep_swa=ep_swa,
        ep_clean=ep_clean,
        staging_contact_qc=staging_contact_qc,
        sig_fsp=sig_fsp,
        swa_1=swa_1,
        beats=beats,
        rr_1=rr_1,
        rr_4=rr_4,
        hr_1=hr_1,
        hr_4=hr_4,
        sigma_coverage=sigma_coverage,
        hr_coverage=hr_coverage,
        qc_warnings=qc_warnings,
    )


def _cache_provenance_fields(source, buffers, final):
    """Assemble source, runtime, acquisition, and QC provenance fields."""
    sigma_qc = final.sigma_contact_qc
    activity_qc = final.channel_activity_qc
    staging_qc = final.staging_contact_qc
    return dict(
        status="ok",
        cache_schema_version=CACHE_SCHEMA_VERSION,
        cache_code_sha256=source.current_cache_digest,
        generated_at_utc=utc_now(),
        code_revision=git_revision(ROOT),
        code_dirty=git_is_dirty(ROOT),
        source_tree_sha256=source_tree_sha256(ROOT),
        runtime_versions_json=json.dumps(runtime_versions(), sort_keys=True),
        source_dataset=source.name,
        source_kind="iEEG.org API",
        source_identity_json=source.source_identity_json,
        source_selection_json=source.source_selection_json,
        source_geometry_reference_channel=source.geometry["reference_channel"],
        failed_chunks_json=json.dumps(buffers.failed_chunks, sort_keys=True),
        ecg_failures_json=json.dumps(buffers.ecg_failures, sort_keys=True),
        acquisition_sample_counts_json=json.dumps(
            buffers.acquisition_sample_counts, sort_keys=True),
        acquisition_sample_count_semantics=(
            "every bounded portal subrequest and every non-overlapping "
            "analysis core records exact requested/returned sample "
            "counts; any mismatch is a fatal acquisition failure"),
        qc_warnings_json=json.dumps(final.qc_warnings, sort_keys=True),
        ecg_processing_method=(
            "NeuroKit2 ecg_clean/ecg_peaks method=neurokit with artifact correction; "
            "not Naji Pan-Tompkins 0.5-100 Hz; requires blinded R-peak validation"),
        ecg_visual_validation=False,
        sigma_coverage=final.sigma_coverage,
        hr_coverage=final.hr_coverage,
        sigma_meets_global_coverage_gate=bool(
            final.sigma_coverage >= MIN_SIGNAL_COVERAGE),
        sigma_per_contact_coverage=sigma_qc["per_contact_coverage"],
        sigma_selected_contact_mask=sigma_qc["selected_mask"],
        sigma_contact_count=sigma_qc["contact_count"],
        sigma_n_selected_contacts=sigma_qc["n_selected"],
        sigma_required_contact_count=sigma_qc["required_contact_count"],
        cortical_signal_nonflat_mask=final.nonflat_contacts,
        cortical_signal_raw_minimum=activity_qc["minimum"],
        cortical_signal_raw_maximum=activity_qc["maximum"],
        cortical_signal_raw_dynamic_range=activity_qc["dynamic_range"],
        cortical_signal_numerical_flat_tolerance=(
            activity_qc["numerical_flat_tolerance"]),
        cortical_signal_finite_sample_count=activity_qc["finite_count"],
        cortical_signal_activity_qc_method=activity_qc["method"],
        staging_selected_contact_mask=staging_qc["selected_contact_mask"],
        staging_contact_count=staging_qc["contact_count"],
        staging_n_selected_contacts=staging_qc["n_selected_contacts"],
        staging_required_contact_count=staging_qc["required_contact_count"],
        staging_per_contact_feature_coverage=(
            staging_qc["per_contact_feature_coverage"]),
        staging_minimum_contact_feature_coverage=(
            staging_qc["minimum_contact_feature_coverage"]),
        staging_minimum_valid_welch_windows=(
            staging_qc["minimum_valid_welch_windows"]),
        staging_swa_normalization=staging_qc["swa_normalization"],
        subject=source.name,
        sf=source.sf,
        night_s=source.night,
        hours=source.hours,
        cortical_chans=np.array(source.ctx),
        anatomy_selection_method=HUP_ANATOMY_SELECTION_METHOD,
        ekg=source.ekg,
        fsp=source.fsp,
        fsp_is_real_peak=source.fsp_real,
    )


def _cache_scientific_fields(source, buffers, final):
    """Assemble unchanged scientific series and estimator contract fields."""
    return dict(
        sigma_fixed=final.sig_fixed,
        sigma_fsp=final.sig_fsp,
        swa=final.swa_1,
        sigma_fixed_by_contact=buffers.sig_fixed_ch,
        sigma_fsp_by_contact=buffers.sig_fsp_ch,
        swa_by_contact=buffers.swa_ch,
        sigma_fixed_power_numerator_by_contact=buffers.sig_fixed_num_ch,
        sigma_fixed_clean_sample_count_by_contact=buffers.sig_fixed_den_ch,
        sigma_fsp_power_numerator_by_contact=buffers.sig_fsp_num_ch,
        sigma_fsp_clean_sample_count_by_contact=buffers.sig_fsp_den_ch,
        swa_power_numerator_by_contact=buffers.swa_num_ch,
        swa_clean_sample_count_by_contact=buffers.swa_den_ch,
        power_samples_per_second=int(source.sf),
        power_historical_minimum_clean_fraction_per_second=0.5,
        power_support_semantics=(
            "numerator=sum of clean squared-envelope samples in each nominal "
            "1-s bin; denominator=count of those clean samples; historical "
            "*_by_contact arrays require denominator >= 0.5*samples_per_second"),
        hr_1=final.hr_1,
        hr_4=final.hr_4,
        rr_1=final.rr_1,
        rr_4=final.rr_4,
        fs_rr=FS_RR,
        ep_dr=final.ep_dr,
        ep_swa=final.ep_swa,
        ep_clean=final.ep_clean,
        epoch_s=EPOCH,
        ep_dr_by_contact=buffers.ep_dr_ch,
        ep_swa_by_contact=buffers.ep_swa_ch,
        ep_clean_fraction_by_contact=buffers.ep_clean_ch,
        ep_measured_fraction_by_contact=buffers.ep_measured_ch,
        ep_valid_welch_window_mask_by_contact=buffers.ep_valid_window_ch,
        ep_window_swa_power_by_contact=buffers.ep_window_swa_ch,
        ep_window_total_power_by_contact=buffers.ep_window_total_ch,
        ep_longest_clean_run_s_by_contact=buffers.ep_longest_clean_run_ch,
        ep_valid_window_span_s_by_contact=buffers.ep_valid_window_span_ch,
        event_3b_so_band_hz=np.asarray(SO_BAND_NAJI),
        event_3b_negative_half_duration_s=np.asarray(
            SO_NEGATIVE_HALF_DURATION_S),
        event_3b_positive_half_max_s=SO_POSITIVE_HALF_MAX_S,
        event_3b_amplitude_semantics=(
            "Naji cites fixed Dang-Vu scalp-voltage gates, which do not "
            "transfer to iEEG; cache retains duration-qualified down/up "
            "amplitudes before an offline within-contact/stage percentile "
            "sensitivity rule"),
        event_3d_rms_12_16_by_contact=buffers.event_3d_rms_ch,
        event_3d_so_phase_0p16_1p25_by_contact=buffers.event_3d_phase_ch,
        event_3d_valid_sample_mask_by_contact=buffers.event_3d_valid_ch,
        event_3d_sampling_hz=EVENT_3D_SAMPLING_HZ,
        event_3d_so_band_hz=np.asarray(EVENT_3D_SO_BAND),
        event_3d_spindle_band_hz=np.asarray(EVENT_3D_SPINDLE_BAND),
        event_3d_so_duration_s=np.asarray(EVENT_3D_SO_DURATION_S),
        event_3d_rms_window_s=EVENT_3D_RMS_WINDOW_S,
        event_3d_ied_padding_s=EVENT_3D_IED_PAD_S,
        event_3d_support_semantics=(
            "20-Hz RMS and SO phase are retained only where the "
            "source sample is measured and passes the +/-2.5-s IED "
            "mask; every retained SO candidate's complete half-open "
            "cycle also passes that mask and has explicit start/stop "
            "indices for stage-containment checks; candidates remain "
            "pre-stage and pre-75th-percentile threshold"),
        beats=final.beats,
    )


def _attach_candidate_fields(payload, source, buffers):
    """Attach sorted 3D and per-contact 3B candidate arrays."""
    candidate_values = np.asarray(
        sorted(
            buffers.event_3d_candidates,
            key=lambda value: (value[0], value[1]),
        ),
        float,
    ).reshape(-1, 5)
    payload["event_3d_so_candidate_contact_index"] = (
        candidate_values[:, 0].astype(np.int16)
    )
    payload["event_3d_so_candidate_sample"] = (
        candidate_values[:, 1].astype(np.int32)
    )
    payload["event_3d_so_candidate_amplitude"] = (
        candidate_values[:, 2].astype(np.float32)
    )
    payload["event_3d_so_candidate_cycle_start_sample"] = (
        candidate_values[:, 3].astype(np.int32)
    )
    payload["event_3d_so_candidate_cycle_stop_sample_exclusive"] = (
        candidate_values[:, 4].astype(np.int32)
    )
    for channel in source.ctx:
        values = np.asarray(
            sorted(buffers.so_candidates[channel]), float).reshape(-1, 4)
        payload[f"so_candidate_t_{channel}"] = values[:, 0]
        payload[f"so_candidate_down_{channel}"] = values[:, 1]
        payload[f"so_candidate_up_{channel}"] = values[:, 2]
        payload[f"so_candidate_p2p_{channel}"] = values[:, 3]


def _assemble_cache_payload(source, buffers, final):
    """Combine provenance, scientific series, and candidate sidecars."""
    payload = _cache_provenance_fields(source, buffers, final)
    payload.update(_cache_scientific_fields(source, buffers, final))
    _attach_candidate_fields(payload, source, buffers)
    return payload


def _run_with_session(n, hours, force, s):
    """Stream and cache one participant using a caller-owned portal session."""
    hours, total_s = validated_analysis_hours(hours)
    name = f"HUP{n}_phaseII"
    fp = os.path.join(OUT, f"{name}.npz")
    current_cache_digest = cache_code_sha256(ROOT)
    if os.path.exists(fp) and not force:
        raise RuntimeError(
            f"{fp} cannot be reused outside a validated complete run; rerun with --force")

    t_start = time.time()
    source = _setup_hup_source(
        name, fp, hours, total_s, current_cache_digest, s)
    if source is None:
        return "skip"

    buffers = _allocate_cache_buffers(source.ctx, source.total_s)
    filters = _build_filter_bank(source.sf, source.fsp)
    t = 0.0
    while t < source.total_s:
        dur = min(CHUNK_S, source.total_s - t)
        pulled = _pull_analysis_chunk(source, buffers, t, dur)
        if pulled is not None:
            _process_analysis_chunk(source, buffers, filters, pulled)
        t += dur

    final = _finalize_cache_series(source, buffers)
    os.makedirs(OUT, exist_ok=True)
    payload = _assemble_cache_payload(source, buffers, final)
    atomic_savez(source.fp, **payload)
    frac = float(np.isfinite(final.sig_fixed).mean())
    print(
        f"[{source.name}] cached in {time.time()-t_start:.0f}s | "
        f"sigma coverage {frac:.1%} | {len(final.beats)} beats | "
        f"clean SO candidates "
        f"{sum(len(v) for v in buffers.so_candidates.values())} -> {source.fp}",
        flush=True,
    )
    return "ok"


def run(n, hours, force=False):
    """Regenerate one participant and always release its portal connection pool."""
    session = sess()
    try:
        return _run_with_session(n, hours, force, session)
    finally:
        session.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default=",".join(map(str, COHORT)))
    ap.add_argument("--hours", type=float, default=7.0)
    ap.add_argument("--force", action="store_true")
    ap.add_argument(
        "--jobs", type=int, default=1,
        help="number of subjects to regenerate concurrently (does not alter estimators)")
    a = ap.parse_args()
    if a.jobs < 1:
        raise ValueError("--jobs must be a positive integer")
    a.hours, _ = validated_analysis_hours(a.hours)
    os.makedirs(OUT, exist_ok=True)
    requested = [int(x) for x in a.subjects.split(",") if x.strip()]
    requested_names = [f"HUP{n}_phaseII" for n in requested]
    config = dict(hours=a.hours, cache_schema_version=CACHE_SCHEMA_VERSION,
                  cache_code_sha256=cache_code_sha256(ROOT),
                  source_pin_schema_version=HUP_SOURCE_PIN_SCHEMA_VERSION,
                  night_probe_workers=NIGHT_PROBE_WORKERS,
                  portal_timeout_s=dict(
                      connect=IEEG_CONNECT_TIMEOUT_S,
                      read=IEEG_READ_TIMEOUT_S))
    if not a.force and validated_complete_run_exists(
            OUT, pipeline="cache_lc_series", requested=requested_names, config=config,
            suffix=".npz", require_current_source_tree=False):
        print(f"validated existing complete cache run ({len(requested_names)} subjects)", flush=True)
        return
    pre_run_outputs = {}
    for subject in requested_names:
        path = os.path.join(OUT, f"{subject}.npz")
        exists = os.path.isfile(path)
        pre_run_outputs[subject] = {
            "existed": exists,
            "sha256": file_sha256(path) if exists else None,
        }
    run_id = start_run_manifest(
        OUT, pipeline="cache_lc_series", requested=requested_names, config=config)
    completed, skipped, failed = [], [], []

    def run_one(n):
        name = f"HUP{n}_phaseII"
        try:
            status = run(n, a.hours, a.force)
            if status == "skip":
                with np.load(os.path.join(OUT, f"{name}.npz"), allow_pickle=False) as record:
                    reason = npz_scalar_text(record, "reason", "unspecified")
                return "skip", dict(subject=name, reason=reason)
            return "ok", name
        except Exception as e:
            print(f"[HUP{n}] ERROR {type(e).__name__}: {e}", flush=True)
            traceback.print_exc()
            return "failed", dict(
                subject=name,
                error=f"{type(e).__name__}: {e}",
                pre_run_output=pre_run_outputs[name],
            )

    if a.jobs == 1:
        records = map(run_one, requested)
    else:
        executor = concurrent.futures.ThreadPoolExecutor(
            max_workers=a.jobs, thread_name_prefix="ieeg-cache-subject")
        records = executor.map(run_one, requested)
    try:
        for status, record in records:
            if status == "ok":
                completed.append(record)
            elif status == "skip":
                skipped.append(record)
            else:
                failed.append(record)
    finally:
        if a.jobs != 1:
            executor.shutdown(wait=True, cancel_futures=True)
    write_run_manifest(
        OUT, pipeline="cache_lc_series", requested=requested_names,
        completed=completed, skipped=skipped, failed=failed,
        config=config, run_id=run_id,
        run_state="failed" if failed else "complete",
        result_files_sha256={
            subject: file_sha256(os.path.join(OUT, f"{subject}.npz"))
            for subject in completed + [value["subject"] for value in skipped]
        })
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
