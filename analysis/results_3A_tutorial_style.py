"""
LEGACY one-participant 3A visualization, retained for shared artifact/power helpers.

The executable analysis and its figures are quarantined because they use superseded FSP,
single-window, and coherence logic. Existing figures are historical outputs, not evidence.

TEST 3A on REAL patient data, drawn in the same visual language as the synthetic teaching
figure (analysis/tutorial_signal_walkthrough.py) so the result is legible at a glance.

Fixes carried over from the method audit (docs/DEC_3A_METHOD_NOTES.md):
  * sigma power is EPOCH-BINNED (mean envelope^2 per 1 s bin) -> inherently anti-aliased.
    The earlier np.interp decimation from 1024 Hz to 4 Hz aliased broadband envelope
    fluctuation into the infraslow band and flattened the spectrum.
  * sigma band = individual FAST-SPINDLE PEAK +/- 1 Hz (Lecci: the 0.02 Hz oscillation is
    strongest in a ~2 Hz band around the FSP and falls off in adjacent bands).
  * SWA (0.75-4 Hz) power is carried as a REAL-DATA NEGATIVE CONTROL: per Lecci Fig 1C/H,
    SWA lacks the 0.02 Hz peak that sigma shows.
  * significance = analytic MSC threshold 1-alpha^(1/(K-1)); surrogates are invalid here.

Usage:
  .venv/bin/python analysis/results_3A_tutorial_style.py \
      --dataset HUP165_phaseII --night HUP165_night1 --win-start-s 41280 --win-min 140
"""
import argparse, json, os
import numpy as np
from scipy import ndimage, signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from infraslow_rr_sigma_coherence import (sess, pull_continuous, notch, r_peaks, ROOT)

SP, HR_C, SO_C, NULL, INK = "#0f9b96", "#d6415f", "#d9871f", "#9aa0b4", "#20232f"
FS_P = 1.0            # power/HR analysis grid (Hz) -- 1 s epochs, Nyquist 0.5 Hz
INFRA = (0.01, 0.04)  # infraslow band (matches the teaching figure's shading)
F_TARGET = 0.02
SWA = (0.75, 4.0)
ALPHA = 0.05
EDGE_TRIM_S = 300.0


def robust_z(v):
    med = np.median(v); mad = np.median(np.abs(v - med)) + 1e-12
    return (v - med) / (1.4826 * mad)


def dilate_boolean_mask(mask, half_width):
    """Boolean equivalent of convolution with an odd all-ones window.

    ``np.convolve`` is quadratic in the window width for this use and made a five-second dilation
    at 2048 Hz dominate the multi-contact cache runtime.  A one-dimensional maximum filter
    performs the same centred inclusive dilation in linear time.
    """
    mask = np.asarray(mask, bool)
    half_width = int(half_width)
    if half_width <= 0 or not mask.any():
        return mask.copy()
    return ndimage.maximum_filter1d(
        mask.astype(np.uint8), size=2 * half_width + 1,
        mode="constant", cval=0).astype(bool)


def ied_clean_mask(y, sf, z_hf=5.0, z_amp=8.0, pad_s=0.5):
    """True where the channel is FREE of interictal epileptiform discharges / sharp artifact.
    Detected on 20-80 Hz (spikes are broadband-sharp) so the 11-16 Hz sigma band itself is not
    used as the detector; plus a raw-amplitude guard. Flags are dilated by +/- pad_s."""
    sos = signal.butter(4, [20.0, 80.0], btype="band", fs=sf, output="sos")
    hf = np.abs(signal.hilbert(signal.sosfiltfilt(sos, y)))
    bad = (robust_z(hf) > z_hf) | (np.abs(robust_z(y)) > z_amp)
    k = int(pad_s * sf)
    bad = dilate_boolean_mask(bad, k)
    return ~bad


def band_power_series(x_mtl, sf, band, dur, mask_ied=True):
    """Mean band power per 1 s epoch, averaged across channels (anti-aliased by binning).
    IED-contaminated samples are excluded from each bin's average, so the time base stays
    continuous (dropping samples outright would break the infraslow phase)."""
    n_bins = int(dur // (1.0 / FS_P)); k = int(sf / FS_P)
    per_ch, masked_frac = [], []
    for j in range(x_mtl.shape[1]):
        y = notch(np.nan_to_num(x_mtl[:, j].astype(float)), sf)
        # SOS, not b/a: at 1024 Hz the SWA band (0.75-4 Hz) is normalised ~0.0015 and a
        # 4th-order b/a Butterworth is numerically unstable there (returns NaN).
        sos = signal.butter(4, [band[0], band[1]], btype="band", fs=sf, output="sos")
        env2 = np.abs(signal.hilbert(signal.sosfiltfilt(sos, y))) ** 2
        clean = ied_clean_mask(y, sf) if mask_ied else np.ones(len(y), bool)
        masked_frac.append(1.0 - clean.mean())
        e = env2[:n_bins * k].reshape(n_bins, k)
        c = clean[:n_bins * k].reshape(n_bins, k).astype(float)
        num, den = (e * c).sum(axis=1), c.sum(axis=1)
        binned = np.where(den >= 0.5 * k, num / np.maximum(den, 1e-12), np.nan)  # bin needs >50% clean
        nan = np.isnan(binned)
        if nan.any():                                   # interpolate the few fully-dirty bins
            idx = np.arange(n_bins)
            binned[nan] = np.interp(idx[nan], idx[~nan], binned[~nan])
        per_ch.append(binned / (np.median(binned) + 1e-12))
    return np.mean(per_ch, axis=0) * 100.0, float(np.mean(masked_frac))   # % of median power


def fast_spindle_peak(x_mtl, sf):
    """Individual fast-spindle peak: PSD peak within 11-16 Hz, averaged across channels."""
    ps = []
    for j in range(x_mtl.shape[1]):
        f, p = signal.welch(notch(np.nan_to_num(x_mtl[:, j].astype(float)), sf), sf, nperseg=int(8 * sf))
        ps.append(p)
    p = np.mean(ps, axis=0)
    m = (f >= 11) & (f <= 16)
    return float(f[m][np.argmax(p[m])])


def hr_grid(x_ekg, sf, n_bins):
    r = r_peaks(x_ekg, sf); bt = r / sf; rr = np.diff(bt)
    bad = (rr < 0.33) | (rr > 1.5)
    t = np.arange(n_bins) / FS_P
    hr = np.interp(t, bt[1:][~bad], 60.0 / rr[~bad])          # upsampling: legitimate
    return hr, float(np.median(60.0 / rr[~bad])), int(len(r)), float(bad.mean())


def hp(x):
    return signal.sosfiltfilt(signal.butter(3, 0.005, btype="high", fs=FS_P, output="sos"), x)


def bp_infra(x):
    return signal.sosfiltfilt(signal.butter(3, list(INFRA), btype="band", fs=FS_P, output="sos"), x)


def msc(a, b):
    nperseg = int(min(512, (len(a) // 4) // 2 * 2))
    f, cxy = signal.coherence(a, b, fs=FS_P, nperseg=nperseg, noverlap=nperseg // 2)
    K = int((len(a) - nperseg // 2) // (nperseg // 2))
    band = (f >= INFRA[0]) & (f <= INFRA[1]); nb = int(band.sum())
    return dict(f=f, cxy=cxy, K=K, at=float(cxy[np.argmin(np.abs(f - F_TARGET))]),
                bmax=float(cxy[band].max()), fpk=float(f[band][np.argmax(cxy[band])]),
                crit=1 - ALPHA ** (1 / (K - 1)), crit_band=1 - (ALPHA / nb) ** (1 / (K - 1)),
                nperseg=nperseg)


def main():
    raise SystemExit(
        "LEGACY/WITHDRAWN one-participant 3A entry point: use cache_lc_series.py followed "
        "by lecci_faithful_3A.py. Shared helper functions remain importable.")
    raise SystemExit(
        "LEGACY 3A FIGURE QUARANTINED: existing figures use superseded estimators. "
        "Use the versioned corrected pipeline; ied_clean_mask remains a shared helper.")
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="HUP165_phaseII")
    ap.add_argument("--night", default="HUP165_night1")
    ap.add_argument("--win-start-s", type=float, default=41280.0)
    ap.add_argument("--win-min", type=float, default=140.0)
    ap.add_argument("--mtl-chans", default="LB2,LB4,LC2,LC4,LH2,LH4")
    ap.add_argument("--ekg", default="EKG1")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    win_s = a.win_min * 60.0

    cache_dir = os.path.join(ROOT, "data", "cache", "results_3A_tutorial_style")
    os.makedirs(cache_dir, exist_ok=True)
    cache = os.path.join(
        cache_dir,
        f"{a.dataset}_{int(a.win_start_s)}_{int(a.win_min)}_{a.mtl_chans.replace(',','-')}.npz")
    if os.path.exists(cache):
        z = np.load(cache); x_mtl, x_ekg, sf = z["x_mtl"], z["x_ekg"], float(z["sf"])
        print(f"[3A-fig] cache {x_mtl.shape}", flush=True)
    else:
        s = sess(); ds = s.open_dataset(a.dataset)
        labels = ds.get_channel_labels(); sf = ds.get_time_series_details(labels[0]).sample_rate
        want = [c for c in a.mtl_chans.split(",") if c in labels]
        ekg = a.ekg if a.ekg in labels else next(l for l in labels if l.upper().startswith(("EKG", "ECG")))
        idx = [labels.index(c) for c in want] + [labels.index(ekg)]
        data = pull_continuous(ds, idx, a.win_start_s, win_s)
        x_mtl, x_ekg = data[:, :len(want)], data[:, len(want)]
        np.savez_compressed(cache, x_mtl=x_mtl, x_ekg=x_ekg, sf=sf)
        print(f"[3A-fig] pulled {x_mtl.shape}", flush=True)

    dur = len(x_ekg) / sf
    fsp = fast_spindle_peak(x_mtl, sf)
    sig_band = (fsp - 1.0, fsp + 1.0)
    print(f"[3A-fig] fast-spindle peak = {fsp:.2f} Hz -> sigma band {sig_band}", flush=True)

    sigma_p, ied_frac = band_power_series(x_mtl, sf, sig_band, dur)
    swa_p, _ = band_power_series(x_mtl, sf, SWA, dur)
    print(f"[3A-fig] IED/artifact-masked fraction of samples = {ied_frac*100:.1f}%", flush=True)
    n_bins = len(sigma_p)
    hr, med_hr, nbeats, frac_bad = hr_grid(x_ekg, sf, n_bins)
    print(f"[3A-fig] {nbeats} beats, {med_hr:.1f} bpm, {frac_bad*100:.2f}% RR dropped", flush=True)

    k = int(EDGE_TRIM_S * FS_P)
    sg, sw, h = hp(sigma_p)[k:-k], hp(swa_p)[k:-k], hp(hr)[k:-k]
    t = (np.arange(len(sg)) + k) / FS_P

    C = msc(h, sg)          # sigma x HR  (the test)
    Cs = msc(h, sw)         # SWA x HR    (control)

    res = dict(dataset=a.dataset, window_h=[a.win_start_s / 3600, (a.win_start_s + win_s) / 3600],
               dur_min=dur / 60, fast_spindle_peak_hz=fsp, sigma_band=list(sig_band),
               nbeats=nbeats, median_hr_bpm=med_hr, rr_dropped_frac=frac_bad,
               ied_masked_frac=ied_frac,
               sigma_coh_at_0p02=C["at"], sigma_coh_band_max=C["bmax"], sigma_band_peak_hz=C["fpk"],
               crit_bin=C["crit"], crit_band=C["crit_band"], K=C["K"],
               significant_at_0p02=bool(C["at"] > C["crit"]),
               significant_band_corrected=bool(C["bmax"] > C["crit_band"]),
               swa_control_coh_at_0p02=Cs["at"], swa_control_band_max=Cs["bmax"])
    outdir = os.path.join(ROOT, "outputs", "results_3A_tutorial_style"); os.makedirs(outdir, exist_ok=True)
    json.dump(res, open(os.path.join(outdir, f"{a.dataset}{a.tag}_3A_results.json"), "w"), indent=2)
    print("[3A-fig] " + json.dumps(res, indent=2), flush=True)

    # ================= figure, teaching-figure framing =================
    fig = plt.figure(figsize=(19, 12.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 0.95, 1.05], hspace=0.42, wspace=0.18)
    verdict = "SIGNIFICANT" if res["significant_at_0p02"] else "NOT significant"
    fig.suptitle(f"TEST 3A on REAL DATA · {a.dataset} · does spindling rise & fall WITH heart rate every ~50 s?",
                 fontsize=15, fontweight="bold", color=INK, x=0.5, y=0.985)

    # (1) zoom out on the spindle envelope
    ax = fig.add_subplot(gs[0, 0])
    seg = slice(0, int(300 * FS_P))
    ax.plot(t[seg] - t[0], sigma_p[k:-k][seg], color=SP, lw=0.8, alpha=0.45)
    sm = signal.savgol_filter(sigma_p[k:-k][seg], 51, 2)
    ax.plot(t[seg] - t[0], sm, color=SP, lw=2.4)
    ax.set(title="① Zoom OUT on the spindle envelope (5 min of real data)",
           xlabel="time (s)", ylabel="sigma power (% of median)")
    ax.text(0.02, 0.94, f"band = FSP {fsp:.1f} ± 1 Hz", transform=ax.transAxes, fontsize=9, color=SP)

    # (2) EKG -> heart rate
    ax = fig.add_subplot(gs[0, 1])
    t0 = int(60 * sf); n10 = int(10 * sf)
    ek = x_ekg[t0:t0 + n10]; te = np.arange(len(ek)) / sf
    ax.plot(te, ek, color=HR_C, lw=0.9)
    rp = r_peaks(x_ekg[t0:t0 + n10], sf)
    ax.plot(rp / sf, ek[rp], "v", color=INK, ms=7)
    ax.set(title="② Turn the EKG into heart rate (real EKG, 10 s)", xlabel="time (s)", ylabel="µV")
    ax.text(0.02, 0.06, f"{nbeats} beats · median {med_hr:.0f} bpm · {frac_bad*100:.2f}% RR dropped",
            transform=ax.transAxes, fontsize=9, color=HR_C)

    # (3) overlay the two slow lines
    ax = fig.add_subplot(gs[1, :])
    z = lambda v: (v - v.mean()) / v.std()
    ax.plot(t / 60, z(bp_infra(sg)), color=SP, lw=1.1, label="spindle envelope (infraslow)")
    ax.plot(t / 60, z(bp_infra(h)), color=HR_C, lw=1.1, label="heart-rate rhythm (infraslow)")
    ax.set(title="③ Overlay the two slow lines — do they march together?", xlabel="time (min)", ylabel="z")
    ax.legend(fontsize=9, loc="upper right")

    # (4) measure it
    ax = fig.add_subplot(gs[2, 0])
    ax.axhspan(0, C["crit"], color=NULL, alpha=0.30, label=f"chance (analytic crit {C['crit']:.3f})")
    ax.axvspan(INFRA[0], INFRA[1], color=SP, alpha=0.07, label="infraslow band")
    ax.plot(C["f"], C["cxy"], color=SP, lw=1.8, label="measured coherence")
    ax.axvline(F_TARGET, color=SP, ls=":", lw=1)
    ax.plot(F_TARGET, C["at"], "v", color=INK, ms=10)
    ax.annotate(f"@0.02 Hz = {C['at']:.3f}\n{verdict}", xy=(F_TARGET, C["at"]),
                xytext=(0.035, max(C["at"], C["crit"]) + 0.10), fontsize=10, color=INK,
                arrowprops=dict(arrowstyle="->", color=INK, lw=1))
    ax.set(title="④ Measure it: infraslow coherence vs chance", xlabel="frequency (Hz)",
           ylabel="coherence", xlim=[0, 0.06], ylim=[0, max(0.35, C["cxy"][:60].max() * 1.25)])
    ax.legend(fontsize=8, loc="upper right")

    # (5) Lecci Fig-1C check: is there a ~0.02 Hz rhythm in the power series at all?
    ax = fig.add_subplot(gs[2, 1])
    for sgl, lab, col in ((sg, f"sigma power (FSP±1 Hz)", SP), (sw, "SWA 0.75–4 Hz (control)", SO_C)):
        fp, pp = signal.welch(sgl, fs=FS_P, nperseg=C["nperseg"])
        ax.plot(fp, pp / pp[(fp > 0.005) & (fp < 0.06)].max(), color=col, lw=1.8, label=lab)
    ax.axvspan(INFRA[0], INFRA[1], color=SP, alpha=0.07)
    ax.axvline(F_TARGET, color=SP, ls=":", lw=1)
    ax.set(title="⑤ Is a ~50 s (0.02 Hz) rhythm even present? (Lecci Fig 1C test)",
           xlabel="frequency (Hz)", ylabel="norm. power", xlim=[0, 0.06])
    ax.legend(fontsize=8)

    swa_txt = (f"SWA control @0.02 Hz = {Cs['at']:.3f}")
    fig.text(0.5, 0.012,
             f"WHAT YOUR EYE SHOULD SEE:  ③ the two lines do NOT lock together.  "
             f"④ coherence @0.02 Hz = {C['at']:.3f} vs chance {C['crit']:.3f} → {verdict}; band-corrected crit "
             f"{C['crit_band']:.3f} → {'passes' if res['significant_band_corrected'] else 'FAILS'}.  "
             f"⑤ neither trace shows an isolated 0.02 Hz bump.   [{swa_txt}]",
             ha="center", fontsize=10.5, style="italic", color=INK)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, f"{a.dataset}{a.tag}_3A_results.{ext}"), dpi=130, bbox_inches="tight")
    print(f"[3A-fig] -> {outdir}", flush=True)


if __name__ == "__main__":
    raise SystemExit(
        "LEGACY/WITHDRAWN one-participant 3A entry point: use cache_lc_series.py followed "
        "by lecci_faithful_3A.py. Shared helper functions remain importable.")
