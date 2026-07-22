"""
3B reimplemented to follow Naji, Krishnan, McDevitt, Bazhenov & Mednick (2019),
J Cogn Neurosci 31:1484 — "Timing between cortical slow oscillations and heart rate bursts
during sleep" — the actual methods paper behind the SO<->heartbeat test.

The earlier 3B cited Chen/Mednick 2022 PNAS, which is a REVIEW with no such analysis in it; the
method used was improvised from a one-line description. This version follows Naji 2019 step by step:

  ECG      0.5-100 Hz Butterworth -> R-peaks -> RR intervals
  RR       resampled at 4 Hz, PIECEWISE CUBIC SPLINE          (was: 1 Hz linear)
  EEG      zero-phase bandpass 0.15-4 Hz                      (was: 0.5-1.25 Hz)
  SO       detected PER CHANNEL by zero-crossing half-waves with duration and amplitude
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
from scipy import signal, interpolate
import neurokit2 as nk

from infraslow_rr_sigma_coherence import sess, pull_continuous, notch, ROOT
from cohort_3A_cortical import COHORT, cortical_channels, delta_ratio, find_night
from cohort_stages_3ABD import stage_epochs, EPOCH, CHUNK_S, SWA_BAND

OUT = os.path.join(ROOT, "outputs", "event_3B_mednick")
SO_BAND_NAJI = (0.15, 4.0)      # Naji 2019
FS_RR = 4.0                     # Naji 2019: RR resampled at 4 Hz
HALF_WIN = 5.0                  # Naji 2019: +/-5 s around the SO downstate (10 s window)
SO_DUR = (0.3, 1.0)             # negative half-wave duration, Dang-Vu 2008 style
AMP_PCT = 75                    # percentile amplitude criterion (iEEG substitute for uV thresholds)


def detect_so_halfwaves(x, sf):
    """SO down-state troughs via zero-crossing half-waves + duration + amplitude criteria.
    Returns trough sample indices."""
    zc = np.where(np.diff(np.signbit(x)))[0]
    if len(zc) < 3:
        return np.array([], int)
    troughs, p2p, amps = [], [], []
    for a, b in zip(zc[:-1], zc[1:]):
        if x[a + 1] >= 0:                       # want NEGATIVE half-waves
            continue
        dur = (b - a) / sf
        if not (SO_DUR[0] <= dur <= SO_DUR[1]):
            continue
        seg = x[a:b]
        tr = a + int(np.argmin(seg))
        nxt = x[b:min(b + int(SO_DUR[1] * sf), len(x))]
        troughs.append(tr); amps.append(-x[tr])
        p2p.append((-x[tr]) + (nxt.max() if len(nxt) else 0.0))
    if not troughs:
        return np.array([], int)
    troughs, amps, p2p = np.array(troughs), np.array(amps), np.array(p2p)
    keep = (amps >= np.percentile(amps, AMP_PCT)) & (p2p >= np.percentile(p2p, AMP_PCT))
    return troughs[keep]


def rr_to_hr_4hz(beats, total_s):
    """RR intervals -> instantaneous HR on a 4 Hz grid via piecewise cubic spline (Naji 2019)."""
    if len(beats) < 20:
        return None
    rr = np.diff(beats)
    good = (rr >= 0.33) & (rr <= 1.5)
    t_b, rr_g = beats[1:][good], rr[good]
    if len(t_b) < 20:
        return None
    grid = np.arange(0, total_s, 1.0 / FS_RR)
    spl = interpolate.CubicSpline(t_b, 60.0 / rr_g, extrapolate=False)
    hr = spl(grid)
    ok = np.isfinite(hr)
    if ok.sum() < 100:
        return None
    hr[~ok] = np.interp(grid[~ok], grid[ok], hr[ok])
    return hr


def so_triggered(hr, trough_times_s, stage_mean_hr, n_sur=200, rng=None):
    """Average HR around SO troughs -> % above stage mean, peak latency, and a random-trigger null."""
    hw = int(HALF_WIN * FS_RR)
    idx = np.round(np.asarray(trough_times_s) * FS_RR).astype(int)
    idx = idx[(idx >= hw) & (idx < len(hr) - hw)]
    if len(idx) < 30:
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
        r = rng.randint(hw, len(hr) - hw, size=len(idx))
        c2 = np.stack([hr[i - hw:i + hw] for i in r]).mean(axis=0)
        null.append(100.0 * (c2[post].max() - stage_mean_hr) / stage_mean_hr)
    null = np.array(null)
    return dict(n_so=int(len(idx)), pct_above_stage_mean=pct, peak_lag_s=pk_lag,
                z=float((pct - null.mean()) / (null.std() + 1e-12)),
                null_mean_pct=float(null.mean()),
                curve=[float(v) for v in curve], lag_s=[float(v) for v in lag])


def run(n, hours, rng):
    name = f"HUP{n}_phaseII"
    fp = os.path.join(OUT, f"{name}.json")
    if os.path.exists(fp):
        print(f"[{name}] cached", flush=True); return
    s = sess(); ds = s.open_dataset(name)
    labels = ds.get_channel_labels(); lab_idx = {l: i for i, l in enumerate(labels)}
    d0 = ds.get_time_series_details(labels[0]); sf = d0.sample_rate
    total_h = (getattr(d0, "duration", 0) or 0) / 3.6e9
    ekg = next((l for l in labels if l.upper().startswith(("EKG", "ECG"))), None)
    ctx = cortical_channels(labels)
    if ekg is None or len(ctx) < 3:
        json.dump(dict(subject=name, status="skip"), open(fp, "w")); print(f"[{name}] SKIP", flush=True); return
    night = find_night(ds, lab_idx, ctx[0], sf, total_h)
    if night is None:
        json.dump(dict(subject=name, status="skip"), open(fp, "w")); print(f"[{name}] SKIP no night", flush=True); return

    idx = [lab_idx[c] for c in ctx] + [lab_idx[ekg]]
    sos_so = signal.butter(3, list(SO_BAND_NAJI), btype="band", fs=sf, output="sos")
    total_s = int(hours * 3600); n_ep = int(total_s // EPOCH)
    ep = dict(dr=np.full(n_ep, np.nan), swa=np.full(n_ep, np.nan), clean=np.ones(n_ep))
    so_t = {c: [] for c in ctx}
    beats = []
    print(f"[{name}] {sf:.0f} Hz | {len(ctx)} cortical ch | streaming {hours} h "
          f"| SO {SO_BAND_NAJI[0]}-{SO_BAND_NAJI[1]} Hz (Naji)", flush=True)

    t = 0.0
    while t < total_s:
        dur = min(CHUNK_S, total_s - t)
        try:
            d = pull_continuous(ds, idx, night + t, dur)
        except Exception:
            t += dur; continue
        off = int(t)
        try:
            cl = nk.ecg_clean(np.nan_to_num(d[:, len(ctx)].astype(float)), sampling_rate=int(sf), method="neurokit")
            _, info = nk.ecg_peaks(cl, sampling_rate=int(sf), method="neurokit", correct_artifacts=True)
            beats.extend((np.asarray(info["ECG_R_Peaks"], int) / sf + off).tolist())
        except Exception:
            pass
        for ci, c in enumerate(ctx):
            x = notch(signal.detrend(np.nan_to_num(d[:, ci].astype(float))), sf)
            if np.std(x) < 1e-9:
                continue
            tr = detect_so_halfwaves(signal.sosfiltfilt(sos_so, x), sf)
            if len(tr):
                so_t[c].extend((tr / sf + off).tolist())
        xa = notch(np.nan_to_num(d[:, :len(ctx)].astype(float)).mean(axis=1), sf)
        ke = int(EPOCH * sf)
        for e in range(int(len(xa) // ke)):
            gi = int((off + e * EPOCH) // EPOCH)
            if gi >= n_ep:
                break
            seg = xa[e * ke:(e + 1) * ke]
            ep["dr"][gi] = delta_ratio(seg, sf)
            fq, pp = signal.welch(seg, sf, nperseg=int(min(4 * sf, len(seg))))
            m = (fq >= SWA_BAND[0]) & (fq < SWA_BAND[1])
            ep["swa"][gi] = float(np.trapezoid(pp[m], fq[m]))
        t += dur

    hr = rr_to_hr_4hz(np.array(sorted(beats)), total_s)
    if hr is None:
        json.dump(dict(subject=name, status="skip", reason="no usable RR"), open(fp, "w")); return
    lab, nrem, sep = stage_epochs(ep)
    rec = dict(subject=name, status="ok", sf=sf, cortical_chans=ctx, night_h=night / 3600,
               n_N2=int((lab == "N2").sum()), n_N3=int((lab == "N3").sum()), gmm_separation=sep,
               method="Naji 2019 (SO 0.15-4 Hz, RR 4 Hz cubic spline, +/-5 s, % above stage mean HR)")
    for stage in ("N2", "N3"):
        keep = set(np.where(lab == stage)[0].tolist())
        if not keep:
            rec[stage] = None; continue
        sec = np.array(sorted(keep)) * EPOCH
        m = np.zeros(len(hr), bool)
        for s0 in sec:
            m[int(s0 * FS_RR):int((s0 + EPOCH) * FS_RR)] = True
        stage_mean = float(np.nanmean(hr[m])) if m.any() else float(np.nanmean(hr))
        per_ch = []
        for c in ctx:
            tt = [x for x in so_t[c] if int(x // EPOCH) in keep]
            r = so_triggered(hr, tt, stage_mean, rng=rng)
            if r:
                r.pop("curve"); r.pop("lag_s")
                r["ch"] = c; per_ch.append(r)
        if not per_ch:
            rec[stage] = None; continue
        rec[stage] = dict(stage_mean_hr=stage_mean, n_channels=len(per_ch),
                          n_so_total=int(sum(d["n_so"] for d in per_ch)),
                          pct_above_stage_mean=float(np.mean([d["pct_above_stage_mean"] for d in per_ch])),
                          peak_lag_s=float(np.mean([d["peak_lag_s"] for d in per_ch])),
                          z=float(np.mean([d["z"] for d in per_ch])), per_channel=per_ch)
        r = rec[stage]
        print(f"[{name}]   {stage}: HR peak {r['pct_above_stage_mean']:+.2f}% above stage mean "
              f"(Naji: N2 +12.09%, SWS +3.35%) | SO->HR lag {r['peak_lag_s']:.2f}s | "
              f"z={r['z']:.1f} | {r['n_so_total']} SOs", flush=True)
    os.makedirs(OUT, exist_ok=True)
    json.dump(rec, open(fp, "w"), indent=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default=",".join(map(str, COHORT)))
    ap.add_argument("--hours", type=float, default=7.0)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.RandomState(0)
    for n in [int(x) for x in a.subjects.split(",") if x.strip()]:
        try:
            run(n, a.hours, rng)
        except Exception as e:
            print(f"[HUP{n}] ERROR {type(e).__name__}: {e}", flush=True); traceback.print_exc()
    print(f"\nJSON -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
