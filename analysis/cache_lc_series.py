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
    hr_1          1 Hz instantaneous heart rate
    hr_4          4 Hz instantaneous heart rate (Naji 2019 resolution, PCHIP of RR intervals)
    ep_dr/ep_swa/ep_clean   per-30 s-epoch staging features
    so_candidate_*_<ch>  clean complete SO candidates and amplitudes; downstream stage-specific
                         thresholds are applied only after staging
    beats         R-peak times (s)

NaN means "not measured here" and is preserved deliberately -- downstream code must handle gaps
rather than delete-and-splice (see spectral_gapped.py for why).

    .venv/bin/python analysis/cache_lc_series.py [--subjects 165,157] [--hours 7]
"""
import argparse, concurrent.futures, json, os, time, traceback
import numpy as np
from scipy import signal, interpolate, ndimage
import neurokit2 as nk

from infraslow_rr_sigma_coherence import (
    IEEG_CONNECT_TIMEOUT_S, IEEG_READ_TIMEOUT_S, sess, pull_continuous, notch, ROOT,
)
from results_3A_tutorial_style import ied_clean_mask
from cohort_3A_cortical import (
    COHORT, HUP_SOURCE_PIN_SCHEMA_VERSION, NIGHT_PROBE_WORKERS, cortical_channels,
    delta_ratio, find_night, verify_hup_source_identity,
)
from cohort_stages_3ABD import band_sos, EPOCH, CHUNK_S, SWA_BAND
from pipeline_version import (CACHE_SCHEMA_VERSION, atomic_savez, cache_code_sha256, file_sha256,
                              git_is_dirty, git_revision, npz_scalar_text, runtime_versions,
                              source_tree_sha256, utc_now, start_run_manifest,
                              validated_complete_run_exists, write_run_manifest)

OUT = os.path.join(ROOT, "data", "derived", "lc_infraslow")
SIGMA_FIXED = (10.0, 15.0)      # Lecci's band
SWA_BAND_L = (0.5, 4.0)       # Lecci human SWA control
SO_BAND_NAJI = (0.15, 4.0)
SO_DUR = (0.3, 1.0)
AMP_PCT = 75
FS_RR = 4.0
MIN_SIGNAL_COVERAGE = 0.80
MIN_CONTACT_COVERAGE = 0.80
MIN_CONTACTS = 3
MIN_CONTACT_FRACTION_PER_BIN = 0.80
FILTER_EDGE_S = 30.0
MISSING_PAD_S = 5.0
HUP_ANATOMY_SELECTION_METHOD = (
    "UNVALIDATED contact-number heuristic; lateral-contact candidates require "
    "coordinate/tissue/SOZ QC")


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


def staging_epoch_features(segment, clean_mask, measured_mask, sf):
    """Artifact-safe delta-ratio and SWA features for one contact/30-s epoch.

    Welch power is not a missing-data estimator.  Computing it on the unmasked raw epoch and only
    later accepting an epoch with (for example) 90% clean samples lets a short high-amplitude IED
    dominate both staging features.  Until a validated gap-aware spectral estimator is supplied,
    fail closed for this contact-epoch whenever any sample is unmeasured or artifact-masked.
    """
    segment = np.asarray(segment, float)
    clean_mask = np.asarray(clean_mask, bool)
    measured_mask = np.asarray(measured_mask, bool)
    if (
        len(segment) == 0
        or clean_mask.shape != segment.shape
        or measured_mask.shape != segment.shape
        or not measured_mask.all()
        or not clean_mask.all()
        or not np.isfinite(segment).all()
    ):
        return np.nan, np.nan
    dr = float(delta_ratio(segment, sf))
    fq, pp = signal.welch(
        segment, sf, nperseg=int(min(4 * sf, len(segment))))
    band = (fq >= SWA_BAND[0]) & (fq < SWA_BAND[1])
    swa = float(np.trapezoid(pp[band], fq[band])) if band.any() else np.nan
    if not np.isfinite(dr) or not np.isfinite(swa):
        return np.nan, np.nan
    return dr, swa


def aggregate_staging_features(
        dr_by_contact, swa_by_contact, clean_by_contact, candidate_contact_mask,
        min_contacts=MIN_CONTACTS,
        min_contact_feature_coverage=MIN_CONTACT_COVERAGE,
        min_contact_fraction_per_epoch=MIN_CONTACT_FRACTION_PER_BIN):
    """Aggregate staging features from one fixed, normalized full-night contact set.

    Changing which contacts contribute from epoch to epoch can turn stable contact gain differences
    into an apparent two-state SWA distribution.  Candidate contacts are fixed before aggregation;
    SWA is normalized within contact over the whole available night; and every reported epoch must
    contain the prespecified fraction of that exact set.
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


def _binned_power_values(env2, clean, sf, n_sec):
    """IED-masked mean of a squared envelope in 1 s bins."""
    k = int(sf)
    e2 = env2[:n_sec * k].reshape(n_sec, k)
    cm = clean[:n_sec * k].reshape(n_sec, k).astype(float)
    num, den = (e2 * cm).sum(1), cm.sum(1)
    return np.where(den >= 0.5 * k, num / np.maximum(den, 1e-12), np.nan)


def _write_multichannel_power(x_ctx, sos, sf, n_sec, total_s, off, dest,
                              core_start_sample=0):
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
        vals = _binned_power_values(env2[start:stop], clean[start:stop], sf, n_sec)
        sl = slice(off, min(off + n_sec, total_s))
        dest[j, sl] = vals[:sl.stop - sl.start]


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
    zc = np.where(np.diff(np.signbit(x)))[0]
    if len(zc) < 4:
        return []
    candidates = []
    edge = int(round(edge_s * sf))
    for a, b, c in zip(zc[:-2], zc[1:-1], zc[2:]):
        if x[a + 1] >= 0:
            continue
        if not (SO_DUR[0] <= (b - a) / sf <= SO_DUR[1]
                and SO_DUR[0] <= (c - b) / sf <= SO_DUR[1]):
            continue
        tr = a + int(np.argmin(x[a:b]))
        if tr < edge or tr >= len(x) - edge:
            continue
        up = float(np.max(x[b:c]))
        down = float(-x[tr])
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


def _run_with_session(n, hours, force, s):
    name = f"HUP{n}_phaseII"
    fp = os.path.join(OUT, f"{name}.npz")
    current_cache_digest = cache_code_sha256(ROOT)
    if os.path.exists(fp) and not force:
        raise RuntimeError(
            f"{fp} cannot be reused outside a validated complete run; rerun with --force")
    t_start = time.time()
    ds = s.open_dataset(name)
    labels = ds.get_channel_labels(); lab_idx = {l: i for i, l in enumerate(labels)}
    d0 = ds.get_time_series_details(labels[0]); sf = d0.sample_rate
    total_h = (getattr(d0, "duration", 0) or 0) / 3.6e9
    ekg = next((l for l in labels if l.upper().startswith(("EKG", "ECG"))), None)
    ctx = cortical_channels(labels)
    source_identity = verify_hup_source_identity(ds, ctx, ekg)
    source_identity_json = json.dumps(source_identity, sort_keys=True)
    if ekg is None or len(ctx) < 3:
        atomic_savez(fp, subject=name, status="skip",
                     cache_schema_version=CACHE_SCHEMA_VERSION,
                     cache_code_sha256=current_cache_digest,
                     source_dataset=name, source_kind="iEEG.org API",
                     source_identity_json=source_identity_json,
                     hours=float(hours),
                     reason=f"ekg={ekg} n_cortical={len(ctx)}")
        print(f"[{name}] SKIP ekg={ekg} n_cortical={len(ctx)}", flush=True)
        return "skip"
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
        return "skip"

    idx = [lab_idx[c] for c in ctx] + [lab_idx[ekg]]
    # Lecci visually selected FSP from all artifact-free NREM spectra. A single arbitrary 300-s
    # sample cannot reproduce that and labels 1/f noise as a peak, so individualized FSP is disabled
    # until an all-NREM, manually QC'd selection is supplied.
    fsp, fsp_real = 13.0, False
    print(f"[{name}] {sf:.0f} Hz | {len(ctx)} lateral-contact candidates | "
          f"fixed sigma primary; individualized FSP disabled | streaming {hours} h", flush=True)

    total_s = int(hours * 3600); n_ep = int(total_s // EPOCH)
    sig_fixed_ch = np.full((len(ctx), total_s), np.nan)
    sig_fsp_ch = np.full((len(ctx), total_s), np.nan)
    swa_ch = np.full((len(ctx), total_s), np.nan)
    ep_dr_ch = np.full((len(ctx), n_ep), np.nan)
    ep_swa_ch = np.full((len(ctx), n_ep), np.nan)
    ep_clean_ch = np.full((len(ctx), n_ep), np.nan)
    beats = []
    so_candidates = {c: [] for c in ctx}
    failed_chunks = []
    ecg_failures = []

    sos_fixed = band_sos(SIGMA_FIXED, sf)
    sos_fsp = band_sos((fsp - 1, fsp + 1), sf)
    sos_swa = band_sos(SWA_BAND_L, sf, 3)
    sos_so = signal.butter(3, list(SO_BAND_NAJI), btype="band", fs=sf, output="sos")

    t = 0.0
    while t < total_s:
        dur = min(CHUNK_S, total_s - t)
        pull_start = max(0.0, t - FILTER_EDGE_S)
        pull_stop = min(float(total_s), t + dur + FILTER_EDGE_S)
        pull_dur = pull_stop - pull_start
        try:
            d = pull_continuous(ds, idx, night + pull_start, pull_dur)
        except Exception as exc:
            failed_chunks.append(dict(start_s=float(t), duration_s=float(dur),
                                      error=f"{type(exc).__name__}: {exc}"))
            t += dur
            continue
        off = int(t)
        core_a = int(round((t - pull_start) * sf))
        core_n = min(int(round(dur * sf)), max(0, len(d) - core_a))
        core_b = core_a + core_n
        x_ctx, x_ekg = d[:, :len(ctx)], d[:, len(ctx)]
        prepared = [prepare_continuous_signal(x_ctx[:, ci], sf) for ci in range(len(ctx))]
        x_channels = np.asarray([
            notch(signal.detrend(filled), sf) for filled, _ in prepared
        ])
        measured_channels = np.asarray([measured for _, measured in prepared])
        clean_channels = np.asarray([
            ied_clean_mask(x, sf) & measured
            for x, measured in zip(x_channels, measured_channels)
        ])

        try:
            ecg_filled, ecg_measured = prepare_continuous_signal(x_ekg, sf)
            cl = nk.ecg_clean(ecg_filled, sampling_rate=int(sf),
                              method="neurokit")
            _, info = nk.ecg_peaks(cl, sampling_rate=int(sf), method="neurokit",
                                   correct_artifacts=True)
            peaks = np.asarray(info["ECG_R_Peaks"], int)
            peaks = peaks[(peaks >= core_a) & (peaks < core_b)]
            peaks = peaks[ecg_measured[peaks]]
            beats.extend((peaks / sf + pull_start).tolist())
        except Exception as exc:
            ecg_failures.append(dict(start_s=float(t), duration_s=float(dur),
                                     error=f"{type(exc).__name__}: {exc}"))

        # per-channel SO troughs (3B needs per-channel events, not the channel average)
        for ci, c in enumerate(ctx):
            xc = x_channels[ci]
            if np.std(xc) < 1e-9:
                continue
            cand = detect_so_candidates(signal.sosfiltfilt(sos_so, xc), sf)
            so_clean = ied_clean_mask(xc, sf, pad_s=5.0) & measured_channels[ci]
            cand = [
                value for value in cand
                if core_a <= int(value[0]) < core_b and so_clean[int(value[0])]
            ]
            so_candidates[c].extend(
                [(tr / sf + pull_start, down, up, p2p) for tr, down, up, p2p in cand])

        n_sec = int(core_n // int(sf))
        for sos_b, dest in ((sos_fixed, sig_fixed_ch), (sos_fsp, sig_fsp_ch),
                            (sos_swa, swa_ch)):
            _write_multichannel_power(
                x_ctx, sos_b, sf, n_sec, total_s, off, dest,
                core_start_sample=core_a)

        ke = int(EPOCH * sf)
        for e in range(int(core_n // ke)):
            gi = int((off + e * EPOCH) // EPOCH)
            if gi >= n_ep:
                break
            for ci, (xc, clean_c, measured_c) in enumerate(zip(
                    x_channels, clean_channels, measured_channels)):
                a = core_a + e * ke
                b = a + ke
                seg = xc[a:b]
                dr, swa_value = staging_epoch_features(
                    seg, clean_c[a:b], measured_c[a:b], sf)
                ep_dr_ch[ci, gi] = dr
                ep_swa_ch[ci, gi] = swa_value
                ep_clean_ch[ci, gi] = float(clean_c[a:b].mean())
        t += dur

    sig_fixed, sigma_contact_qc = _aggregate_full_night_power(
        sig_fixed_ch, min_contact_coverage=MIN_CONTACT_COVERAGE,
        min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN,
        return_details=True)
    eligible_contacts = sigma_contact_qc["selected_mask"]
    ep_dr, ep_swa, ep_clean, staging_contact_qc = aggregate_staging_features(
        ep_dr_ch, ep_swa_ch, ep_clean_ch, eligible_contacts)
    sig_fsp = _aggregate_full_night_power(
        sig_fsp_ch, eligible_channels=eligible_contacts,
        min_contact_coverage=MIN_CONTACT_COVERAGE, min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN)
    swa_1 = _aggregate_full_night_power(
        swa_ch, eligible_channels=eligible_contacts,
        min_contact_coverage=MIN_CONTACT_COVERAGE, min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN)

    # R-peaks at chunk edges can be repeated or spuriously near-duplicated.
    beats = sanitize_beats(beats)
    rr_1, rr_4, hr_1, hr_4 = interpolate_tachograms(beats, total_s)

    sigma_coverage = float(np.isfinite(sig_fixed).mean())
    hr_coverage = float(np.isfinite(hr_4).mean())
    if failed_chunks:
        raise RuntimeError(f"{len(failed_chunks)} acquisition chunks failed")
    if ecg_failures:
        raise RuntimeError(
            f"ECG detection failed in {len(ecg_failures)} chunks; publication cache fails closed")
    if hr_coverage < MIN_SIGNAL_COVERAGE:
        # Insufficient participant data are an expected cohort exclusion, not an execution
        # failure.  Raising here used to put the participant in ``manifest.failed`` and thereby
        # invalidate every otherwise usable participant in the same production run.
        reason = (
            f"cardiac coverage QC failed: HR={hr_coverage:.1%}; "
            f"minimum is {MIN_SIGNAL_COVERAGE:.0%}")
        atomic_savez(
            fp, subject=name, status="skip",
            cache_schema_version=CACHE_SCHEMA_VERSION,
            cache_code_sha256=current_cache_digest,
            source_dataset=name, source_kind="iEEG.org API",
            source_identity_json=source_identity_json,
            source_selection_json=source_selection_json,
            hours=float(hours), reason=reason,
            sigma_coverage=sigma_coverage, hr_coverage=hr_coverage)
        print(f"[{name}] SKIP {reason}", flush=True)
        return "skip"

    os.makedirs(OUT, exist_ok=True)
    payload = dict(status="ok", cache_schema_version=CACHE_SCHEMA_VERSION,
                   cache_code_sha256=current_cache_digest,
                   generated_at_utc=utc_now(), code_revision=git_revision(ROOT),
                   code_dirty=git_is_dirty(ROOT), source_tree_sha256=source_tree_sha256(ROOT),
                   runtime_versions_json=json.dumps(runtime_versions(), sort_keys=True),
                   source_dataset=name, source_kind="iEEG.org API",
                   source_identity_json=source_identity_json,
                   source_selection_json=source_selection_json,
                   failed_chunks_json=json.dumps(failed_chunks, sort_keys=True),
                   ecg_failures_json=json.dumps(ecg_failures, sort_keys=True),
                   ecg_processing_method=(
                       "NeuroKit2 ecg_clean/ecg_peaks method=neurokit with artifact correction; "
                       "not Naji Pan-Tompkins 0.5-100 Hz; requires blinded R-peak validation"),
                   ecg_visual_validation=False,
                   sigma_coverage=sigma_coverage, hr_coverage=hr_coverage,
                   sigma_meets_global_coverage_gate=bool(
                       sigma_coverage >= MIN_SIGNAL_COVERAGE),
                   sigma_per_contact_coverage=sigma_contact_qc["per_contact_coverage"],
                   sigma_selected_contact_mask=sigma_contact_qc["selected_mask"],
                   sigma_contact_count=sigma_contact_qc["contact_count"],
                   sigma_n_selected_contacts=sigma_contact_qc["n_selected"],
                   sigma_required_contact_count=sigma_contact_qc["required_contact_count"],
                   staging_selected_contact_mask=staging_contact_qc["selected_contact_mask"],
                   staging_contact_count=staging_contact_qc["contact_count"],
                   staging_n_selected_contacts=staging_contact_qc["n_selected_contacts"],
                   staging_required_contact_count=staging_contact_qc["required_contact_count"],
                   staging_per_contact_feature_coverage=(
                       staging_contact_qc["per_contact_feature_coverage"]),
                   staging_minimum_contact_feature_coverage=(
                       staging_contact_qc["minimum_contact_feature_coverage"]),
                   staging_swa_normalization=staging_contact_qc["swa_normalization"],
                   subject=name, sf=sf, night_s=night, hours=hours,
                   cortical_chans=np.array(ctx),
                   anatomy_selection_method=HUP_ANATOMY_SELECTION_METHOD,
                   ekg=ekg, fsp=fsp, fsp_is_real_peak=fsp_real,
                   sigma_fixed=sig_fixed, sigma_fsp=sig_fsp, swa=swa_1,
                   hr_1=hr_1, hr_4=hr_4, rr_1=rr_1, rr_4=rr_4, fs_rr=FS_RR,
                   ep_dr=ep_dr, ep_swa=ep_swa, ep_clean=ep_clean, epoch_s=EPOCH,
                   beats=beats)
    for c in ctx:
        values = np.asarray(sorted(so_candidates[c]), float).reshape(-1, 4)
        payload[f"so_candidate_t_{c}"] = values[:, 0]
        payload[f"so_candidate_down_{c}"] = values[:, 1]
        payload[f"so_candidate_up_{c}"] = values[:, 2]
        payload[f"so_candidate_p2p_{c}"] = values[:, 3]
    atomic_savez(fp, **payload)
    frac = float(np.isfinite(sig_fixed).mean())
    print(f"[{name}] cached in {time.time()-t_start:.0f}s | sigma coverage {frac:.1%} | "
          f"{len(beats)} beats | clean SO candidates "
          f"{sum(len(v) for v in so_candidates.values())} -> {fp}", flush=True)
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
            return "failed", dict(subject=name, error=f"{type(e).__name__}: {e}")

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
        result_files_sha256={
            subject: file_sha256(os.path.join(OUT, f"{subject}.npz"))
            for subject in completed + [value["subject"] for value in skipped]
        })
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
