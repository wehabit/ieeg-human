"""
Withdrawn historical 3B attempt adapted from Naji, Krishnan, McDevitt, Bazhenov & Mednick (2019),
J Cogn Neurosci 31:1484 — "Timing between cortical slow oscillations and heart rate bursts
during sleep" — the actual methods paper behind the SO<->heartbeat test.

The earlier 3B cited Chen/Mednick 2022 PNAS, which is a REVIEW with no such analysis in it; the
method used was improvised from a one-line description. This is not a step-by-step replication:
among other deviations, cache builders use NeuroKit's default cleaner/detector rather than Naji's
0.5-100 Hz ECG preprocessing, Pan-Tompkins detector, and visual R-peak confirmation.

  ECG      0.5-100 Hz Butterworth -> R-peaks -> RR intervals
  RR       resampled at 4 Hz, PIECEWISE CUBIC SPLINE          (was: 1 Hz linear)
  EEG      zero-phase bandpass 0.15-4 Hz                      (was: 0.5-1.25 Hz)
  SO       detected PER CHANNEL by zero-crossing down/up states with duration and amplitude
           criteria (Dang-Vu 2008 style)                      (was: channel average, top-quartile peaks)
  measure  average RR time-series referenced to the SO DOWN-STATE TROUGH over a 10 s window;
           report (a) the HR peak as % above that stage's mean HR, and
                  (b) the SO-HR peak-to-peak INTERVAL = time from SO trough to the HR maximum
           -- (b) is the paper's headline statistic and the previous version never computed it.

VALIDATION TARGETS from Naji 2019 (scalp F3/F4, healthy sleepers):
    Stage 2 : HR peak  12.09 +/- 1.48 %  above mean Stage-2 HR
    SWS     : HR peak   3.35 +/- 1.01 %
    HR acceleration+deceleration duration: 5.52 +/- 2.51 s
So the paper predicts HR modulation ~3.6x LARGER in N2 than in SWS -- a quantitative,
directional prediction this script can be checked against.

DEVIATIONS (unavoidable, stated rather than hidden):
  * Region: Naji used FRONTAL SCALP (F3/F4); here lateral neocortical iEEG contacts.
  * Dang-Vu amplitude criteria are absolute microvolts for scalp EEG and do not transfer to
    intracranial recordings, so amplitude thresholds are PERCENTILE-based within channel.
  * Staging is the GMM-on-slow-wave-power approximation (no EOG/EMG in iEEG), not R&K scoring.

    .venv/bin/python analysis/event_3B_mednick.py [--subjects 165,...] [--hours 7]
"""
import argparse, json, os, traceback
import numpy as np
from scipy import signal
import neurokit2 as nk

from infraslow_rr_sigma_coherence import sess, pull_continuous, notch, ROOT
from cohort_3A_cortical import COHORT, cortical_channels, delta_ratio, find_night
from cohort_stages_3ABD import stage_epochs, EPOCH, CHUNK_S, SWA_BAND
from cache_lc_series import (detect_so_candidates, threshold_so_candidates,
                             MIN_SIGNAL_COVERAGE, FILTER_EDGE_S, prepare_continuous_signal,
                             interpolate_tachograms)
from pipeline_version import (ANALYSIS_VERSION, atomic_json_dump, source_tree_sha256,
                              start_run_manifest, write_run_manifest)
from results_3A_tutorial_style import ied_clean_mask

OUT = os.path.join(ROOT, "outputs", "event_3B_mednick")
SO_BAND_NAJI = (0.15, 4.0)      # Naji 2019
FS_RR = 4.0                     # Naji 2019: RR resampled at 4 Hz
HALF_WIN = 5.0                  # Naji 2019: +/-5 s around the SO downstate (10 s window)


def rr_baseline_hr(rr):
    """Convert an average-RR baseline to HR in the same order as RR event curves."""
    values = np.asarray(rr, float)
    values = values[np.isfinite(values) & (values > 0)]
    return float(60.0 / values.mean()) if len(values) else np.nan


def rr_to_hr_4hz(beats, total_s):
    """RR intervals -> instantaneous HR on a 4 Hz grid via shape-preserving cubic interpolation.

    Long acquisition gaps are left missing.  The previous CubicSpline bridged every missing chunk
    and could overshoot between beats, creating synthetic HR dynamics.
    """
    _, _, _, hr = interpolate_tachograms(beats, total_s, fs_rr=FS_RR)
    if np.isfinite(hr).sum() < 100:
        return None
    return hr


def so_triggered(hr, trough_times_s, stage_mean_hr, stage_pool_idx, n_sur=200, rng=None):
    """Average HR around SO troughs -> % above stage mean, peak latency, and a random-trigger null.

    `stage_pool_idx` are the HR-grid sample indices belonging to THIS stage; surrogate triggers are
    drawn only from them.

    BUG FIXED HERE. Surrogates used to be drawn from the whole night (`rng.randint(hw, len(hr)-hw)`)
    while the statistic was expressed as "% above THIS STAGE's mean HR". The null therefore measured
    the offset between the stage mean and the whole-night mean rather than SO-locking. Proof from
    the HUP165 output: inverting each stage's null baseline recovered one common absolute heart rate
    (N2 -> 93.284 bpm, N3 -> 93.251 bpm, agreeing to 0.03 bpm), and the observed effect was
    near-identical in the two stages (+0.471% vs +0.478%) while z flipped from -1.75 to +4.90 purely
    because N2's mean HR sits below and N3's above the whole-night mean.

    Drawing observed and surrogate triggers from the same pool also cancels the upward bias of the
    `max()` over the post-trough window, which otherwise inflates the observed value alone.

    NaN-SAFE. The cached HR series carries NaN across gaps > 5 s (the old rr_to_hr_4hz interpolated
    across every gap, so this never arose). A trigger whose +/-5 s window overlaps such a gap would
    turn the averaged curve into NaN. Both observed and surrogate triggers are therefore restricted
    to centres whose FULL window is finite; when HR has no gaps this changes nothing.
    """
    hw = int(HALF_WIN * FS_RR)
    finite = np.isfinite(hr)
    # cumulative count of finite samples -> a window [i-hw, i+hw) is all-finite iff it contains 2*hw
    csum = np.concatenate(([0], np.cumsum(finite)))
    full_window = np.zeros(len(hr), bool)
    c = np.arange(hw, len(hr) - hw)
    full_window[c] = (csum[c + hw] - csum[c - hw]) == (2 * hw)

    idx = np.round(np.asarray(trough_times_s) * FS_RR).astype(int)
    idx = idx[(idx >= hw) & (idx < len(hr) - hw)]
    idx = idx[full_window[idx]]
    if len(idx) < 30:
        return None
    pool = np.asarray(stage_pool_idx, int)
    pool = pool[(pool >= hw) & (pool < len(hr) - hw)]
    pool = pool[full_window[pool]]
    if len(pool) < 100:
        return None
    seg = np.stack([hr[i - hw:i + hw] for i in idx])
    curve = seg.mean(axis=0)
    lag = (np.arange(len(curve)) - hw) / FS_RR
    post = lag >= 0                                   # Naji: HR peak FOLLOWS the downstate
    pk_i = int(np.argmax(curve[post]))
    pk_val, pk_lag = float(curve[post][pk_i]), float(lag[post][pk_i])
    pct = 100.0 * (pk_val - stage_mean_hr) / stage_mean_hr
    rng = rng or np.random.RandomState(0)
    null = []
    for _ in range(n_sur):
        r = rng.choice(pool, size=len(idx), replace=True)
        c2 = np.stack([hr[i - hw:i + hw] for i in r]).mean(axis=0)
        null.append(100.0 * (c2[post].max() - stage_mean_hr) / stage_mean_hr)
    null = np.array(null)
    return dict(n_so=int(len(idx)), pct_above_stage_mean=pct, peak_lag_s=pk_lag,
                z=float((pct - null.mean()) / (null.std() + 1e-12)),
                null_mean_pct=float(null.mean()), n_surrogate_pool=int(len(pool)),
                curve=[float(v) for v in curve], lag_s=[float(v) for v in lag])


def subject_so_triggered(tachogram, trough_times_by_channel, stage_mean_hr, stage_pool_idx,
                         n_sur=1000, rng=None, domain="hr", minimum_channels=1,
                         channel_ids=None, minimum_events_per_channel=30,
                         minimum_surrogate_pool_samples=100):
    """One participant-level SO->HR estimate across channels.

    Naji first formed an SO-triggered curve for each frontal electrode.  The old implementation
    maximised every intracranial channel separately and then averaged those maxima, which inflated
    the reported magnitude and made the average latency physically ambiguous.  Here channel curves
    are averaged first and the participant has exactly one peak/effect.  Each null replicate repeats
    that complete aggregation, preserving the estimator used for the observed value.

    ``domain="rr"`` is the paper-aligned path: event curves are averaged in RR space and their
    minima define HR-burst timing. ``domain="hr"`` is retained for explicit sensitivity tests.
    """
    if domain not in ("rr", "hr"):
        raise ValueError("domain must be 'rr' or 'hr'")
    series = np.asarray(tachogram, float)
    rng = rng or np.random.RandomState(0)
    hw = int(HALF_WIN * FS_RR)
    finite = np.isfinite(series)
    csum = np.concatenate(([0], np.cumsum(finite)))
    full_window = np.zeros(len(series), bool)
    centres = np.arange(hw, len(series) - hw)
    full_window[centres] = (csum[centres + hw] - csum[centres - hw]) == (2 * hw)

    pool = np.asarray(stage_pool_idx, int)
    pool = pool[(pool >= hw) & (pool < len(series) - hw)]
    # Both observed and null windows must remain entirely inside the stage.  Checking only the
    # centre allowed wake/stage-transition HR to leak into ±5 s windows at epoch boundaries.
    in_stage = np.zeros(len(series), bool)
    in_stage[pool] = True
    stage_csum = np.concatenate(([0], np.cumsum(in_stage)))
    full_stage = np.zeros(len(series), bool)
    full_stage[centres] = (
        stage_csum[centres + hw] - stage_csum[centres - hw]) == (2 * hw)
    eligible = full_window & full_stage
    pool = np.unique(pool[eligible[pool]])
    if len(pool) < int(minimum_surrogate_pool_samples):
        return None

    if channel_ids is None:
        channel_ids = list(range(len(trough_times_by_channel)))
    if len(channel_ids) != len(trough_times_by_channel):
        raise ValueError("channel_ids must align with trough_times_by_channel")
    channel_indices = []
    retained_channel_ids = []
    for channel_id, trough_times_s in zip(channel_ids, trough_times_by_channel):
        idx = np.round(np.asarray(trough_times_s, float) * FS_RR).astype(int)
        idx = idx[(idx >= hw) & (idx < len(series) - hw)]
        idx = np.unique(idx[eligible[idx]])
        if len(idx) >= int(minimum_events_per_channel):
            channel_indices.append(idx)
            retained_channel_ids.append(channel_id)
    if len(channel_indices) < int(minimum_channels):
        return None

    def channel_curves(indices):
        return [np.stack([series[i - hw:i + hw] for i in idx]).mean(axis=0)
                for idx in indices]

    def post_peak_index(curve, post_mask):
        values = curve[post_mask]
        return int(np.argmin(values) if domain == "rr" else np.argmax(values))

    def percent_above_stage_mean(curve, post_mask):
        idx = post_peak_index(curve, post_mask)
        value = float(curve[post_mask][idx])
        peak_hr = 60.0 / value if domain == "rr" else value
        return 100.0 * (peak_hr - stage_mean_hr) / stage_mean_hr

    observed_channel_curves = channel_curves(channel_indices)
    curve = np.mean(observed_channel_curves, axis=0)
    lag = (np.arange(len(curve)) - hw) / FS_RR
    post = lag >= 0
    pk_rel = post_peak_index(curve, post)
    pct = percent_above_stage_mean(curve, post)
    # Naji's timing endpoint takes a peak time for each electrode and then averages those times.
    # Keep that paper-defined latency distinct from the peak of the participant-average curve.
    channel_peak_lags = [
        float(lag[post][post_peak_index(ch_curve, post)])
        for ch_curve in observed_channel_curves
    ]
    paper_peak_lag = float(np.mean(channel_peak_lags))

    null = np.empty(n_sur, float)
    ranks = []
    for idx in channel_indices:
        rank = np.searchsorted(pool, idx)
        if np.any(rank >= len(pool)) or not np.array_equal(pool[rank], idx):
            raise RuntimeError("observed SO index was not found in the eligible stage pool")
        ranks.append(rank)
    # Shift the complete multichannel event ensemble by one shared offset in stage-time. This
    # preserves exact cross-channel coincidences and event-train dependence. Independently drawing
    # each channel made the null average artificially precise and produced a 78% false-positive
    # rate in a null simulation with identical channel event trains.
    guard = 2 * hw
    valid_shifts = np.arange(guard, len(pool) - guard + 1)
    if not len(valid_shifts):
        valid_shifts = np.arange(1, len(pool))
    for s in range(n_sur):
        shift = int(rng.choice(valid_shifts))
        random_indices = [pool[(rank + shift) % len(pool)] for rank in ranks]
        c2 = np.mean(channel_curves(random_indices), axis=0)
        null[s] = percent_above_stage_mean(c2, post)

    # Retain the former whole-stage shift only as a diagnostic. It fixes the earlier stage-offset
    # and independent-channel defects, but does not preserve local nonstationary trends when SOs
    # cluster in time, so it is not a valid event-locking p value.
    p_upper = float((1 + np.sum(null >= pct)) / (n_sur + 1))
    hr_curve = 60.0 / curve if domain == "rr" else curve
    pre = lag < 0
    local_baseline_hr = float(np.mean(hr_curve[pre]))
    local_change_pct = float(
        100.0 * (np.max(hr_curve[post]) - local_baseline_hr) / local_baseline_hr)
    peak_to_peak_pct = float(
        100.0 * (np.max(hr_curve) - np.min(hr_curve)) / np.mean(hr_curve))
    return dict(
        n_channels=len(channel_indices),
        event_contact_ids=retained_channel_ids,
        n_so_total=int(sum(len(x) for x in channel_indices)),
        pct_above_stage_mean=pct,
        peak_lag_s=paper_peak_lag,
        channel_peak_lag_s=channel_peak_lags,
        participant_curve_peak_lag_s=float(lag[post][pk_rel]),
        event_locked_local_change_pct=local_change_pct,
        event_curve_peak_to_peak_pct=peak_to_peak_pct,
        local_pre_event_mean_hr=local_baseline_hr,
        z=None,
        p_upper=None,
        stage_shift_z_diagnostic=float(
            (pct - null.mean()) / (null.std(ddof=1) + 1e-12)),
        stage_shift_p_upper_diagnostic=p_upper,
        null_mean_pct=float(null.mean()),
        null_sd_pct=float(null.std(ddof=1)),
        n_surrogates=int(n_sur),
        n_surrogate_pool=int(len(pool)),
        null_method="one shared circular shift in eligible stage-time across all channels",
        inference_status=(
            "disabled: whole-stage shifts do not preserve local nonstationary trends or "
            "event-density clustering; raw Naji magnitude and local change are descriptive"),
        tachogram_domain=domain,
        curve=[float(60.0 / v if domain == "rr" else v) for v in curve],
        rr_curve=([float(v) for v in curve] if domain == "rr" else None),
        lag_s=[float(v) for v in lag],
    )


def run(n, hours, rng, force=False, run_id=None, tree_digest=None):
    name = f"HUP{n}_phaseII"
    fp = os.path.join(OUT, f"{name}.json")
    if os.path.exists(fp) and not force:
        try:
            existing = json.load(open(fp))
        except Exception:
            existing = {}
        if (existing.get("analysis_version") == ANALYSIS_VERSION
                and np.isclose(existing.get("hours", np.nan), hours)
                and existing.get("source_tree_sha256") == tree_digest):
            existing["run_id"] = run_id
            atomic_json_dump(existing, fp)
            print(f"[{name}] current output", flush=True)
            return "cached" if existing.get("status") == "ok" else existing.get("status", "skip")
        raise RuntimeError(f"{fp} is legacy or has a different configuration; rerun with --force")
    s = sess(); ds = s.open_dataset(name)
    labels = ds.get_channel_labels(); lab_idx = {l: i for i, l in enumerate(labels)}
    d0 = ds.get_time_series_details(labels[0]); sf = d0.sample_rate
    total_h = (getattr(d0, "duration", 0) or 0) / 3.6e9
    ekg = next((l for l in labels if l.upper().startswith(("EKG", "ECG"))), None)
    ctx = cortical_channels(labels)
    if ekg is None or len(ctx) < 3:
        atomic_json_dump(
            dict(subject=name, status="skip", analysis_version=ANALYSIS_VERSION,
                 run_id=run_id, source_tree_sha256=tree_digest,
                 reason=f"ekg={ekg} n_cortical={len(ctx)}"), fp)
        print(f"[{name}] SKIP", flush=True)
        return "skip"
    night = find_night(ds, lab_idx, ctx[0], sf, total_h, required_h=hours)
    if night is None:
        atomic_json_dump(
            dict(subject=name, status="skip", analysis_version=ANALYSIS_VERSION,
                 run_id=run_id, source_tree_sha256=tree_digest, reason="no night"), fp)
        print(f"[{name}] SKIP no night", flush=True)
        return "skip"

    idx = [lab_idx[c] for c in ctx] + [lab_idx[ekg]]
    sos_so = signal.butter(3, list(SO_BAND_NAJI), btype="band", fs=sf, output="sos")
    total_s = int(hours * 3600); n_ep = int(total_s // EPOCH)
    ep = dict(dr=np.full(n_ep, np.nan), swa=np.full(n_ep, np.nan), clean=np.zeros(n_ep))
    so_candidates = {c: [] for c in ctx}
    beats = []
    failed_chunks = []
    ecg_failures = []
    print(f"[{name}] {sf:.0f} Hz | {len(ctx)} cortical ch | streaming {hours} h "
          f"| SO {SO_BAND_NAJI[0]}-{SO_BAND_NAJI[1]} Hz (Naji)", flush=True)

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
        try:
            ecg_filled, ecg_measured = prepare_continuous_signal(d[:, len(ctx)], sf)
            cl = nk.ecg_clean(ecg_filled, sampling_rate=int(sf), method="neurokit")
            _, info = nk.ecg_peaks(cl, sampling_rate=int(sf), method="neurokit", correct_artifacts=True)
            peaks = np.asarray(info["ECG_R_Peaks"], int)
            peaks = peaks[(peaks >= core_a) & (peaks < core_b)]
            peaks = peaks[ecg_measured[peaks]]
            beats.extend((peaks / sf + pull_start).tolist())
        except Exception as exc:
            ecg_failures.append(dict(start_s=float(t), duration_s=float(dur),
                                     error=f"{type(exc).__name__}: {exc}"))
        prepared = [prepare_continuous_signal(d[:, ci], sf) for ci in range(len(ctx))]
        x_channels = [
            notch(signal.detrend(filled), sf) for filled, _ in prepared
        ]
        measured_channels = [measured for _, measured in prepared]
        clean_channels = [
            ied_clean_mask(x, sf, pad_s=HALF_WIN) & measured
            for x, measured in zip(x_channels, measured_channels)
        ]
        for ci, c in enumerate(ctx):
            x = x_channels[ci]
            if np.std(x) < 1e-9:
                continue
            cand = detect_so_candidates(signal.sosfiltfilt(sos_so, x), sf)
            clean = clean_channels[ci]
            cand = [
                value for value in cand
                if core_a <= int(value[0]) < core_b and clean[int(value[0])]
            ]
            so_candidates[c].extend(
                [(tr / sf + pull_start, down, up, p2p) for tr, down, up, p2p in cand])
        ke = int(EPOCH * sf)
        for e in range(int(core_n // ke)):
            gi = int((off + e * EPOCH) // EPOCH)
            if gi >= n_ep:
                break
            dr_values, swa_values, clean_values = [], [], []
            for x, clean, measured in zip(
                    x_channels, clean_channels, measured_channels):
                a = core_a + e * ke
                b = a + ke
                seg = x[a:b]
                if measured[a:b].all():
                    dr_values.append(delta_ratio(seg, sf))
                    fq, pp = signal.welch(seg, sf, nperseg=int(min(4 * sf, len(seg))))
                    m = (fq >= SWA_BAND[0]) & (fq < SWA_BAND[1])
                    swa_values.append(float(np.trapezoid(pp[m], fq[m])))
                else:
                    dr_values.append(np.nan)
                    swa_values.append(np.nan)
                clean_values.append(float(clean[a:b].mean()))
            ep["dr"][gi] = float(np.nanmedian(dr_values))
            ep["swa"][gi] = float(np.nanmedian(swa_values))
            ep["clean"][gi] = float(np.nanmedian(clean_values))
        t += dur

    if failed_chunks:
        raise RuntimeError(f"{len(failed_chunks)} acquisition chunks failed")
    hr = rr_to_hr_4hz(np.array(sorted(beats)), total_s)
    if hr is None:
        atomic_json_dump(
            dict(subject=name, status="skip", analysis_version=ANALYSIS_VERSION,
                 run_id=run_id, source_tree_sha256=tree_digest, reason="no usable RR"), fp)
        return "skip"
    rr_series = 60.0 / hr
    coverage = float(np.isfinite(rr_series).mean())
    if coverage < MIN_SIGNAL_COVERAGE:
        raise RuntimeError(
            f"RR coverage {coverage:.1%} is below required {MIN_SIGNAL_COVERAGE:.0%}; "
            f"ECG detector failed in {len(ecg_failures)} chunks")
    lab, nrem, sep = stage_epochs(ep)
    rec = dict(subject=name, status="ok", analysis_version=ANALYSIS_VERSION,
               run_id=run_id, source_tree_sha256=tree_digest,
               hours=float(hours),
               sf=sf, cortical_chans=ctx, night_h=night / 3600,
               rr_coverage=coverage, ecg_failures=ecg_failures,
               n_N2=int((lab == "N2").sum()), n_N3=int((lab == "N3").sum()), gmm_separation=sep,
               method="Naji adaptation (SO 0.15-4 Hz, RR PCHIP at 4 Hz, +/-5 s, "
                      "% above stage mean HR; shared-shift multichannel null)")
    endpoint_reasons = {}
    for stage in ("N2", "N3"):
        keep = set(np.where(lab == stage)[0].tolist())
        if not keep:
            rec[stage] = None
            endpoint_reasons[stage] = "no staged epochs"
            continue
        sec = np.array(sorted(keep)) * EPOCH
        m = np.zeros(len(hr), bool)
        for s0 in sec:
            m[int(s0 * FS_RR):int((s0 + EPOCH) * FS_RR)] = True
        stage_rr = rr_series[m] if m.any() else rr_series
        stage_mean = rr_baseline_hr(stage_rr)
        stage_pool = np.where(m)[0]          # surrogate triggers must come from THIS stage only
        troughs_by_channel = []
        for c in ctx:
            stage_candidates = [
                value for value in so_candidates[c]
                if int(value[0] // EPOCH) in keep
            ]
            tt = threshold_so_candidates(stage_candidates)
            if len(tt) >= 30:
                troughs_by_channel.append(tt)
        result = subject_so_triggered(
            rr_series, troughs_by_channel, stage_mean, stage_pool,
            n_sur=1000, rng=rng, domain="rr")
        if result is None:
            rec[stage] = None
            endpoint_reasons[stage] = (
                "no channel has >=30 eligible SOs with complete in-stage RR windows")
            continue
        result.pop("curve")
        result.pop("lag_s")
        result["stage_mean_hr"] = stage_mean
        result["excess_over_null_pct"] = float(
            result["pct_above_stage_mean"] - result["null_mean_pct"])
        result["aggregation"] = (
            "one participant-level magnitude from the average channel curve; "
            "Naji latency is the mean of channel peak times")
        rec[stage] = result
        r = rec[stage]
        print(f"[{name}]   {stage}: HR peak {r['pct_above_stage_mean']:+.2f}% above stage mean "
              f"(Naji: N2 +12.09%, SWS +3.35%) | SO->HR lag {r['peak_lag_s']:.2f}s | "
              f"local Δ={r['event_locked_local_change_pct']:+.2f}% "
              f"(inference disabled) | {r['n_so_total']} SOs", flush=True)
    unavailable = [stage for stage in ("N2", "N3") if rec.get(stage) is None]
    if unavailable:
        rec["status"] = "partial"
        rec["reason"] = "unavailable 3B endpoint(s): " + ", ".join(
            f"{stage} ({endpoint_reasons.get(stage, 'unspecified')})"
            for stage in unavailable)
    rec["endpoint_unavailable_reasons"] = endpoint_reasons
    rec["endpoint_availability"] = {
        stage: rec.get(stage) is not None for stage in ("N2", "N3")
    }
    os.makedirs(OUT, exist_ok=True)
    atomic_json_dump(rec, fp)
    return rec["status"]


def main():
    raise SystemExit(
        "LEGACY/WITHDRAWN direct 3B entry point: use cache_lc_series.py followed by "
        "event_3B_cached.py. Shared RR/event helper functions remain importable.")
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default=",".join(map(str, COHORT)))
    ap.add_argument("--hours", type=float, default=7.0)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.RandomState(0)
    requested_n = [int(x) for x in a.subjects.split(",") if x.strip()]
    requested = [f"HUP{n}_phaseII" for n in requested_n]
    config = dict(hours=a.hours, force=a.force,
                  analysis_version=ANALYSIS_VERSION,
                  n_surrogates=1000,
                  null_method="shared circular shift in eligible stage-time")
    run_id = start_run_manifest(
        OUT, pipeline="event_3B_mednick_direct", requested=requested, config=config)
    tree_digest = source_tree_sha256(ROOT)
    completed, skipped, failed = [], [], []
    for n in requested_n:
        name = f"HUP{n}_phaseII"
        try:
            status = run(
                n, a.hours, rng, force=a.force, run_id=run_id, tree_digest=tree_digest)
            if status not in ("ok", "cached"):
                record = json.load(open(os.path.join(OUT, f"{name}.json")))
                skipped.append(dict(subject=name, reason=record.get("reason", "unspecified")))
            else:
                completed.append(name)
        except Exception as e:
            print(f"[{name}] ERROR {type(e).__name__}: {e}", flush=True)
            traceback.print_exc()
            failed.append(dict(subject=name, error=f"{type(e).__name__}: {e}"))
    write_run_manifest(
        OUT, pipeline="event_3B_mednick_direct",
        requested=requested, completed=completed,
        skipped=skipped, failed=failed,
        config=config, run_id=run_id)
    print(f"\nJSON -> {OUT}", flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    raise SystemExit(
        "LEGACY/WITHDRAWN direct 3B entry point: use cache_lc_series.py followed by "
        "event_3B_cached.py. Shared RR/event helper functions remain importable.")
