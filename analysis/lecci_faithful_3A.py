"""3A re-implemented to follow Lecci et al. 2017 Sci Adv, not a fixed-frequency proxy for it.

WHAT THE PREVIOUS 3A ACTUALLY TESTED, AND WHY IT IS NOT LECCI'S TEST
-------------------------------------------------------------------
The cohort 3A asked: "is magnitude-squared coherence between sigma power and heart rate significant
in the single bin at 0.0195 Hz?" Lecci never ran that test. Their human analysis is two steps:

  Step 1 (Fig 1C/G) -- does sigma power oscillate infraslow AT ALL, and at what rate?
      Morlet wavelet power time course -> a SECOND wavelet transform of that time course at 0.001 Hz
      resolution over 0.001-0.12 Hz -> one spectrum per NREM bout >=120 s -> averaged across ALL
      bouts WEIGHTED BY DURATION -> normalised to its own mean -> Gaussian fit to locate THAT
      SUBJECT'S OWN PEAK. The reported 0.019 +/- 0.001 Hz (n=27) is the mean of those fitted peaks.
  Step 2 (Fig 6) -- does heart rate track it?
      Sigma and HR resampled to a common grid, 120 s intervals z-transformed, CROSS-CORRELATION with
      heart rate as source wave, averaged within subject then across subjects. The result is a lag,
      not a coherence value.

Testing one hard-coded bin is a materially different and stricter hypothesis. Lecci report a
between-subject SD of the peak frequency of ~0.0052 Hz (0.019 +/- 0.001 SEM, n=27) against a bin
width here of 0.0039 Hz: if EVERY subject had a textbook Lecci rhythm, only ~29% would have their
personal peak inside the single tested bin. The cohort observed 7/23 = 30.4%.

This script implements both steps. Deviations from Lecci that remain are forced by the data and are
listed in DEVIATIONS below.

DEVIATIONS (unavoidable, stated rather than hidden)
  * Staging: iEEG has no EOG/EMG, so S2 vs SWS cannot be scored. NREM is the GMM-on-delta-ratio
    approximation and the N2/N3 split is GMM-on-slow-wave-power. Lecci's S2 > SWS claim is therefore
    NOT directly testable here; pooled NREM is the primary analysis.
  * Signal: sigma power comes from a Hilbert envelope of a Butterworth band, binned to 1 s, rather
    than a 4-cycle Morlet sampled at 0.1 s then smoothed with a 4 s moving average. After Lecci's
    4 s smoothing both are band-limited well below the 0.5 Hz Nyquist of a 1 s grid, so the
    infraslow content is preserved; the optional --smooth-4s flag applies their moving average.
  * Region: lateral neocortical depth contacts, not parietal scalp. Lecci's 0.02 Hz oscillation is
    parietal-maximal and declines frontally.
  * Population: epilepsy patients on anti-seizure medication.

    .venv/bin/python analysis/lecci_faithful_3A.py [--band fsp|fixed] [--smooth-4s]
"""
import argparse, json, os
import numpy as np
from scipy import signal, optimize

from cohort_stages_3ABD import stage_epochs, EPOCH
from cohort_3A_cortical import COHORT
from spectral_gapped import (coherence_gapped, analytic_msc_threshold, fill_short_gaps,
                             contiguous_runs)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "data", "derived", "lc_infraslow")
OUT = os.path.join(ROOT, "outputs", "lecci_faithful_3A")

FS = 1.0                      # cached derived-series rate
MIN_BOUT_S = 120.0            # Lecci: bouts >= 120 s (>= 4 epochs)
F_LO, F_HI, F_STEP = 0.004, 0.120, 0.001      # Lecci: 0.001-0.12 Hz at 0.001 Hz resolution
N_CYCLES = 3.0                # cycles in the second-stage Morlet
PEAK_SEARCH = (0.008, 0.060)  # where an LC-type infraslow peak is allowed to be
F_LECCI = 0.02
NPERSEG = 256
ALPHA = 0.05
XCORR_WIN_S = 120.0           # Lecci: 120 s intervals, z-transformed
XCORR_MAX_LAG_S = 60.0


# ------------------------------------------------------------------ cache access
_CACHE_DIR = [CACHE]          # mutable so callers can point at an alternate cache (e.g. ds003848)


def set_cache(path):
    _CACHE_DIR[0] = path


def stages_for(d):
    """Return (lab, nrem, sep). Uses REAL scored stages when the cache carries a `stage_lab` array
    (ds003848, staged from EMG/EOG); otherwise falls back to the GMM-on-slow-wave proxy (HUP, no
    EOG/EMG). NREM = N2 or N3; Wake/REM/N1 are excluded from NREM by construction with real stages."""
    if "stage_lab" in getattr(d, "files", []):
        raw = np.asarray(d["stage_lab"]).astype(str)
        lab = np.full(len(raw), "", dtype=object)
        lab[raw == "N2"] = "N2"
        lab[raw == "N3"] = "N3"
        nrem = (lab == "N2") | (lab == "N3")
        return lab, nrem, None
    ep = dict(dr=d["ep_dr"], swa=d["ep_swa"], clean=d["ep_clean"])
    return stage_epochs(ep)


def load(subject):
    fp = os.path.join(_CACHE_DIR[0], f"{subject}.npz")
    if not os.path.exists(fp):
        return None
    d = np.load(fp, allow_pickle=True)
    if str(d["status"]) != "ok":
        return None
    return d


def nrem_bouts(nrem, min_s=MIN_BOUT_S):
    """Contiguous NREM runs, returned as (start_second, stop_second) and >= min_s long."""
    out = []
    for i0, i1 in contiguous_runs(nrem):
        s0, s1 = int(i0 * EPOCH), int(i1 * EPOCH)
        if s1 - s0 >= min_s:
            out.append((s0, s1))
    return out


# ------------------------------------------------------------------ Lecci Step 1
def morlet_kernel(f, fs, n_cycles):
    sigma_t = n_cycles / (2 * np.pi * f)
    t = np.arange(-3.5 * sigma_t, 3.5 * sigma_t + 1.0 / fs, 1.0 / fs)
    w = np.exp(2j * np.pi * f * t) * np.exp(-t ** 2 / (2 * sigma_t ** 2))
    return w / np.sqrt(np.sum(np.abs(w) ** 2))


def morlet_spectrum(x, fs, freqs, n_cycles=N_CYCLES):
    """Mean Morlet power across the time steps of one bout.

    Frequencies whose wavelet is longer than the bout are returned as NaN rather than silently
    estimated from an edge-dominated convolution.
    """
    x = np.asarray(x, float)
    x = x - np.nanmean(x)
    if not np.isfinite(x).all():
        return np.full(len(freqs), np.nan)
    out = np.full(len(freqs), np.nan)
    for i, f in enumerate(freqs):
        k = morlet_kernel(f, fs, n_cycles)
        if len(k) > len(x):
            continue
        w = signal.fftconvolve(x, k, mode="valid")
        out[i] = float(np.mean(np.abs(w) ** 2))
    return out


def subject_spectrum(sig, nrem, fs=FS, smooth_4s=False):
    """Duration-weighted mean infraslow spectrum over all NREM bouts >= 120 s (Lecci Fig 1G)."""
    freqs = np.arange(F_LO, F_HI + 1e-12, F_STEP)
    x, _, _ = fill_short_gaps(sig, fs, max_gap_s=5.0)
    if smooth_4s:
        k = int(round(4 * fs))
        if k > 1:
            valid = np.isfinite(x)
            xf = np.where(valid, x, 0.0)
            num = np.convolve(xf, np.ones(k) / k, mode="same")
            den = np.convolve(valid.astype(float), np.ones(k) / k, mode="same")
            x = np.where(den > 0.5, num / np.maximum(den, 1e-12), np.nan)
    # Split each NREM bout at any gap that survived interpolation and treat every artifact-free
    # sub-run >= MIN_BOUT_S as its own bout. Rejecting a whole bout for one long gap would keep only
    # the subset of bouts that happen to be gap-free -- on HUP165 that discarded 5490 s of 20191 s
    # and biased the spectrum toward short, unusually clean stretches.
    spectra, weights = [], []
    for s0, s1 in nrem_bouts(nrem):
        finite = np.isfinite(x[s0:s1])
        for a, b in contiguous_runs(finite):
            seg = x[s0 + a:s0 + b]
            if len(seg) < MIN_BOUT_S * fs:
                continue
            sp = morlet_spectrum(seg, fs, freqs)
            if np.isfinite(sp).sum() < 5:
                continue
            spectra.append(sp); weights.append(len(seg))
    if not spectra:
        return freqs, None, 0, 0.0
    S = np.array(spectra); W = np.array(weights, float)
    with np.errstate(invalid="ignore"):
        num = np.nansum(S * W[:, None], axis=0)
        den = np.nansum(np.where(np.isfinite(S), 1.0, 0.0) * W[:, None], axis=0)
    mean_spec = np.where(den > 0, num / np.maximum(den, 1e-300), np.nan)
    ok = np.isfinite(mean_spec)
    if ok.sum() < 10:
        return freqs, None, len(spectra), float(W.sum())
    mean_spec = mean_spec / np.nanmean(mean_spec[ok])       # Lecci: normalise to its own mean
    return freqs, mean_spec, len(spectra), float(W.sum())


def _scale_free(n, slope, rng):
    """Time series whose power spectrum follows f**slope (Lecci's scale-free control profile)."""
    x = rng.randn(n)
    f = np.fft.rfftfreq(n)
    S = np.fft.rfft(x)
    amp = np.zeros_like(f)
    amp[1:] = f[1:] ** (slope / 2.0)
    y = np.fft.irfft(S * amp, n)
    return y / (y.std() + 1e-12)


def peak_significance(sig, nrem, slope, observed_prominence, n_sur=200, seed=0, smooth_4s=False):
    """Lecci fig S3 G-J, per subject: is the observed peak more prominent than one produced by a
    SCALE-FREE profile with the same 1/f slope and the SAME bout-length distribution?

    HONEST LIMITATION: this per-subject control is WEAK. On synthetic data a planted 0.02 Hz rhythm
    yields p ~ 0.21 (never < 0.05), because a single ~2-hour night gives a noisy prominence estimate
    and matched 1/f surrogates occasionally produce comparable bumps. It is reported for
    completeness and to bound each subject, but a per-subject null must NOT be read as "no rhythm".
    The powerful test is at the COHORT level: real fitted peaks CLUSTER (simulated SD ~0.0008 Hz,
    100% in 0.015-0.025 Hz) whereas 1/f noise scatters (SD ~0.013 Hz, ~48% in band). `surrogate_peaks`
    is returned so summarize_corrected_3AB.py can compare real clustering against this cohort's own
    scale-free scatter. Lecci's own evidence is exactly that clustering: 0.019 +/- 0.001 across n=27.

    The surrogate inherits the real series' NaN pattern, so bout count, bout lengths and gap
    structure are identical by construction; only the spectral content is replaced.
    """
    rng = np.random.RandomState(seed)
    n = len(sig)
    finite = np.isfinite(sig)
    proms, peaks = [], []
    for _ in range(n_sur):
        surr = np.where(finite, _scale_free(n, slope, rng), np.nan)
        fr, sp, _, _ = subject_spectrum(surr, nrem, smooth_4s=smooth_4s)
        if sp is None:
            continue
        pk = fit_peak(fr, sp)
        p = pk.get("prominence_over_background", np.nan)
        if np.isfinite(p):
            proms.append(float(p))
            if pk.get("peak_hz"):
                peaks.append(float(pk["peak_hz"]))
    if len(proms) < 20:
        return None
    proms = np.array(proms)
    p_value = float((1 + np.sum(proms >= observed_prominence)) / (1 + len(proms)))
    return dict(n_surrogates=len(proms), observed_prominence=float(observed_prominence),
                null_median=float(np.median(proms)), null_p95=float(np.percentile(proms, 95)),
                p_value=p_value,
                # full surrogate peak-frequency list: the per-subject prominence test is weak, so
                # the powerful cohort-level test is whether REAL peaks cluster more tightly than
                # these do. Lecci's own evidence is exactly that clustering (0.019 +/- 0.001, n=27).
                surrogate_peaks=[float(p) for p in peaks],
                surrogate_peak_hz_median=(float(np.median(peaks)) if peaks else None),
                surrogate_peak_in_infraslow_frac=(
                    float(np.mean((np.array(peaks) >= 0.015) & (np.array(peaks) <= 0.025)))
                    if peaks else None))


def _gauss3(f, *p):
    return sum(p[i] * np.exp(-((f - p[i + 1]) / p[i + 2]) ** 2) for i in range(0, 9, 3))


def fit_peak(freqs, spec, search=PEAK_SEARCH):
    """Lecci's Gaussian fit, plus a 1/f-corrected peak as a robustness check.

    Returns dict with the fitted peak frequency, its prominence over the fitted 1/f background, and
    which method produced it.
    """
    ok = np.isfinite(spec) & (freqs >= F_LO) & (freqs <= F_HI)
    f, y = freqs[ok], spec[ok]
    if len(f) < 15:
        return dict(peak_hz=None, method="insufficient")

    # 1/f background in log-log, so a "peak" means a genuine local excess, not just low-frequency power
    pos = y > 0
    co = np.polyfit(np.log(f[pos]), np.log(y[pos]), 1)
    resid = y - np.exp(np.polyval(co, np.log(f)))
    band = (f >= search[0]) & (f <= search[1])
    pk_idx, _ = signal.find_peaks(resid[band])
    detrended_peak = float(f[band][pk_idx[np.argmax(resid[band][pk_idx])]]) if len(pk_idx) else None
    prominence = (float(np.max(resid[band]) / (np.std(resid[~band]) + 1e-12))
                  if (~band).sum() > 5 else float("nan"))

    peak, method = None, "detrended"
    try:
        a0 = float(np.max(y[band])) if band.any() else float(np.max(y))
        p0 = [a0, detrended_peak or 0.02, 0.01, a0 / 2, 0.05, 0.02, a0 / 2, 0.01, 0.05]
        lo = [0, F_LO, 1e-4] * 3
        hi = [10 * a0 + 1e-9, F_HI, 1.0] * 3
        popt, _ = optimize.curve_fit(_gauss3, f, y, p0=p0, bounds=(lo, hi), maxfev=20000)
        cands = [(popt[i], popt[i + 1]) for i in range(0, 9, 3)
                 if search[0] <= popt[i + 1] <= search[1]]
        if cands:
            peak, method = float(max(cands)[1]), "gauss3"
    except Exception:
        pass
    if peak is None:
        peak = detrended_peak
    return dict(peak_hz=peak, detrended_peak_hz=detrended_peak, method=method,
                prominence_over_background=prominence, slope_1_over_f=float(co[0]))


# ------------------------------------------------------------------ Lecci Step 2
def cross_correlation(sig, hr, nrem, fs=FS, win_s=XCORR_WIN_S, max_lag_s=XCORR_MAX_LAG_S):
    """Lecci Fig 6: z-transform each 120 s interval, cross-correlate with HR as source wave,
    average the cross-correlograms within subject. Positive lag => sigma FOLLOWS heart rate."""
    s, _, _ = fill_short_gaps(sig, fs, 5.0)
    h, _, _ = fill_short_gaps(hr, fs, 5.0)
    w = int(win_s * fs); ml = int(max_lag_s * fs)
    acc, n = [], 0
    for s0, s1 in nrem_bouts(nrem, min_s=win_s):
        for st in range(s0, s1 - w + 1, w):
            a, b = h[st:st + w], s[st:st + w]
            if not (np.isfinite(a).all() and np.isfinite(b).all()):
                continue
            if a.std() < 1e-9 or b.std() < 1e-9:
                continue
            a = (a - a.mean()) / a.std(); b = (b - b.mean()) / b.std()
            c = signal.correlate(b, a, mode="full") / len(a)
            lags = signal.correlation_lags(len(b), len(a), mode="full")
            keep = np.abs(lags) <= ml
            acc.append(c[keep]); n += 1
    if n < 5:
        return None
    lags_s = signal.correlation_lags(w, w, mode="full")[np.abs(
        signal.correlation_lags(w, w, mode="full")) <= ml] / fs
    m = np.mean(acc, axis=0)
    i = int(np.argmax(np.abs(m)))
    return dict(lag_s=lags_s.tolist(), xcorr=m.tolist(), n_intervals=n,
                peak_r=float(m[i]), peak_lag_s=float(lags_s[i]))


# ------------------------------------------------------------------ per subject
def analyse(subject, band="fsp", smooth_4s=False, n_sur=200):
    d = load(subject)
    if d is None:
        return None
    sig = d["sigma_fsp"] if band == "fsp" else d["sigma_fixed"]
    swa, hr = d["swa"], d["hr_1"]
    lab, nrem, sep = stages_for(d)
    if nrem.sum() < 40:
        return dict(subject=subject, status="skip", reason="insufficient NREM")

    freqs, spec, n_bouts, tot_s = subject_spectrum(sig, nrem, smooth_4s=smooth_4s)
    pk = fit_peak(freqs, spec) if spec is not None else dict(peak_hz=None, method="no spectrum")
    _, spec_swa, _, _ = subject_spectrum(swa, nrem, smooth_4s=smooth_4s)   # Lecci negative control
    pk_swa = fit_peak(freqs, spec_swa) if spec_swa is not None else dict(peak_hz=None)

    # Lecci's scale-free control: is the peak more prominent than matched 1/f noise produces?
    sig_null = None
    if spec is not None and np.isfinite(pk.get("prominence_over_background", np.nan)):
        sig_null = peak_significance(sig, nrem, pk.get("slope_1_over_f", 0.0),
                                     pk["prominence_over_background"],
                                     n_sur=n_sur, smooth_4s=smooth_4s)
        if sig_null:
            pk["peak_p_value"] = sig_null["p_value"]

    # coherence over ALL NREM (non-NREM masked out), gap-aware, read at the SUBJECT'S OWN peak
    mask = np.zeros(len(sig), bool)
    for s0, s1 in nrem_bouts(nrem):
        mask[s0:s1] = True
    sg = np.where(mask, sig, np.nan)
    hg = np.where(mask, hr, np.nan)
    co = coherence_gapped(hg, sg, fs=FS, nperseg=NPERSEG, highpass=0.005)
    coh = None
    if co is not None:
        crit = analytic_msc_threshold(co["K"], ALPHA)
        f, c = co["f"], co["cxy"]
        own = pk["peak_hz"] if pk.get("peak_hz") else None
        coh = dict(K=co["K"], crit=float(crit), n_valid=co["n_valid"], filled_frac=co["filled_frac"],
                   at_lecci=float(c[np.argmin(np.abs(f - F_LECCI))]),
                   sig_at_lecci=bool(c[np.argmin(np.abs(f - F_LECCI))] > crit),
                   at_own_peak=(float(c[np.argmin(np.abs(f - own))]) if own else None),
                   sig_at_own_peak=(bool(c[np.argmin(np.abs(f - own))] > crit) if own else None),
                   f=f.tolist(), cxy=c.tolist())

    xc = cross_correlation(sig, hr, nrem)
    return dict(subject=subject, status="ok", band=band, smooth_4s=smooth_4s,
                n_nrem_epochs=int(nrem.sum()), n_bouts=n_bouts, bout_seconds=tot_s,
                gmm_separation=(None if sep is None else float(sep)),
                spectrum=dict(f=freqs.tolist(),
                              sigma=(None if spec is None else np.where(np.isfinite(spec), spec, None).tolist()),
                              swa=(None if spec_swa is None else np.where(np.isfinite(spec_swa), spec_swa, None).tolist())),
                peak=pk, peak_swa=pk_swa, peak_null=sig_null, coherence=coh, xcorr=xc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default=",".join(map(str, COHORT)))
    ap.add_argument("--band", choices=["fsp", "fixed"], default="fsp")
    ap.add_argument("--smooth-4s", action="store_true")
    ap.add_argument("--n-sur", type=int, default=200,
                    help="scale-free surrogates for the Lecci fig S3 peak-significance control")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for n in [int(x) for x in a.subjects.split(",") if x.strip()]:
        name = f"HUP{n}_phaseII"
        try:
            r = analyse(name, a.band, a.smooth_4s, n_sur=a.n_sur)
        except Exception as e:
            print(f"[{name}] ERROR {type(e).__name__}: {e}", flush=True); continue
        if r is None:
            print(f"[{name}] no cache", flush=True); continue
        if r.get("status") != "ok":
            print(f"[{name}] {r.get('reason')}", flush=True); continue
        rows.append(r)
        json.dump(r, open(os.path.join(OUT, f"{name}.json"), "w"))
        pk = r["peak"]; co = r["coherence"]; xc = r["xcorr"]
        pv = pk.get("peak_p_value")
        print(f"[{name}] bouts={r['n_bouts']:3d} peak="
              f"{'none' if not pk.get('peak_hz') else format(pk['peak_hz'], '.4f') + ' Hz'} "
              f"(prom {pk.get('prominence_over_background', float('nan')):.2f}, "
              f"p={'n/a' if pv is None else format(pv, '.3f')}) | "
              f"coh@0.02={'n/a' if not co else format(co['at_lecci'], '.3f')} "
              f"{'*' if co and co['sig_at_lecci'] else ''} | "
              f"coh@own={'n/a' if not co or co['at_own_peak'] is None else format(co['at_own_peak'], '.3f')} "
              f"{'*' if co and co.get('sig_at_own_peak') else ''} | "
              f"xcorr r={'n/a' if not xc else format(xc['peak_r'], '+.3f')} "
              f"lag={'n/a' if not xc else format(xc['peak_lag_s'], '+.0f')}s", flush=True)
    print(f"\n{len(rows)} subjects -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
