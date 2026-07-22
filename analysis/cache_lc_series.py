"""Stream each subject's night ONCE and cache the derived series to disk.

Every re-analysis of 3A/3B so far required a fresh 7-hour stream from iEEG.org per subject, because
only summary statistics were ever saved. That made method fixes expensive and is the main reason
deviations from Lecci went uncorrected. This script pulls the night once and stores everything the
downstream tests need, so subsequent analyses run offline in seconds.

Cached per subject -> data/derived/lc_infraslow/<subject>.npz

    sigma_fixed   1 Hz, IED-masked mean sigma power, Lecci's fixed 10-15 Hz band
    sigma_fsp     1 Hz, IED-masked mean sigma power, individual fast-spindle peak +/- 1 Hz
    swa           1 Hz, IED-masked 0.75-4 Hz power  (Lecci's negative control: no 0.02 Hz peak)
    hr_1          1 Hz instantaneous heart rate
    hr_4          4 Hz instantaneous heart rate (Naji 2019 resolution, cubic spline)
    ep_dr/ep_swa/ep_clean   per-30 s-epoch staging features
    so_t_<ch>     slow-oscillation trough times (s) per cortical channel, Naji 0.15-4 Hz half-waves
    beats         R-peak times (s)

NaN means "not measured here" and is preserved deliberately -- downstream code must handle gaps
rather than delete-and-splice (see spectral_gapped.py for why).

    .venv/bin/python analysis/cache_lc_series.py [--subjects 165,157] [--hours 7]
"""
import argparse, os, time, traceback
import numpy as np
from scipy import signal, interpolate
import neurokit2 as nk

from infraslow_rr_sigma_coherence import sess, pull_continuous, notch, ROOT
from results_3A_tutorial_style import ied_clean_mask
from cohort_3A_cortical import COHORT, cortical_channels, delta_ratio, find_night
from cohort_stages_3ABD import fsp_from, band_sos, EPOCH, CHUNK_S, SWA_BAND

OUT = os.path.join(ROOT, "data", "derived", "lc_infraslow")
SIGMA_FIXED = (10.0, 15.0)      # Lecci's band
SWA_BAND_L = (0.75, 4.0)
SO_BAND_NAJI = (0.15, 4.0)
SO_DUR = (0.3, 1.0)
AMP_PCT = 75
FS_RR = 4.0


def _binned_power(env2, clean, sf, n_sec, total_s, off, dest):
    """IED-masked mean of a squared envelope in 1 s bins, written into `dest` at second `off`."""
    k = int(sf)
    e2 = env2[:n_sec * k].reshape(n_sec, k)
    cm = clean[:n_sec * k].reshape(n_sec, k).astype(float)
    num, den = (e2 * cm).sum(1), cm.sum(1)
    vals = np.where(den >= 0.5 * k, num / np.maximum(den, 1e-12), np.nan)
    sl = slice(off, min(off + n_sec, total_s))
    dest[sl] = vals[:sl.stop - sl.start]


def detect_so_halfwaves(x, sf):
    """Naji 2019 / Dang-Vu style negative half-waves; returns trough sample indices."""
    zc = np.where(np.diff(np.signbit(x)))[0]
    if len(zc) < 3:
        return np.array([], int)
    troughs, p2p, amps = [], [], []
    for a, b in zip(zc[:-1], zc[1:]):
        if x[a + 1] >= 0:
            continue
        if not (SO_DUR[0] <= (b - a) / sf <= SO_DUR[1]):
            continue
        tr = a + int(np.argmin(x[a:b]))
        nxt = x[b:min(b + int(SO_DUR[1] * sf), len(x))]
        troughs.append(tr); amps.append(-x[tr])
        p2p.append((-x[tr]) + (nxt.max() if len(nxt) else 0.0))
    if not troughs:
        return np.array([], int)
    troughs, amps, p2p = np.array(troughs), np.array(amps), np.array(p2p)
    keep = (amps >= np.percentile(amps, AMP_PCT)) & (p2p >= np.percentile(p2p, AMP_PCT))
    return troughs[keep]


def run(n, hours, force=False):
    name = f"HUP{n}_phaseII"
    fp = os.path.join(OUT, f"{name}.npz")
    if os.path.exists(fp) and not force:
        print(f"[{name}] cached", flush=True); return
    t_start = time.time()
    s = sess(); ds = s.open_dataset(name)
    labels = ds.get_channel_labels(); lab_idx = {l: i for i, l in enumerate(labels)}
    d0 = ds.get_time_series_details(labels[0]); sf = d0.sample_rate
    total_h = (getattr(d0, "duration", 0) or 0) / 3.6e9
    ekg = next((l for l in labels if l.upper().startswith(("EKG", "ECG"))), None)
    ctx = cortical_channels(labels)
    if ekg is None or len(ctx) < 3:
        np.savez_compressed(fp, status="skip", reason=f"ekg={ekg} n_cortical={len(ctx)}")
        print(f"[{name}] SKIP ekg={ekg} n_cortical={len(ctx)}", flush=True); return
    night = find_night(ds, lab_idx, ctx[0], sf, total_h)
    if night is None:
        np.savez_compressed(fp, status="skip", reason="no night")
        print(f"[{name}] SKIP no night", flush=True); return

    idx = [lab_idx[c] for c in ctx] + [lab_idx[ekg]]
    fsp, fsp_real = fsp_from(pull_continuous(ds, idx, night, 300.0)[:, :len(ctx)], sf)
    print(f"[{name}] {sf:.0f} Hz | {len(ctx)} cortical | FSP {fsp:.2f} Hz | streaming {hours} h",
          flush=True)

    total_s = int(hours * 3600); n_ep = int(total_s // EPOCH)
    sig_fixed = np.full(total_s, np.nan)
    sig_fsp = np.full(total_s, np.nan)
    swa_1 = np.full(total_s, np.nan)
    ep_dr = np.full(n_ep, np.nan); ep_swa = np.full(n_ep, np.nan); ep_clean = np.zeros(n_ep)
    beats = []
    so_t = {c: [] for c in ctx}

    sos_fixed = band_sos(SIGMA_FIXED, sf)
    sos_fsp = band_sos((fsp - 1, fsp + 1), sf)
    sos_swa = band_sos(SWA_BAND_L, sf, 3)
    sos_so = signal.butter(3, list(SO_BAND_NAJI), btype="band", fs=sf, output="sos")

    t = 0.0
    while t < total_s:
        dur = min(CHUNK_S, total_s - t)
        try:
            d = pull_continuous(ds, idx, night + t, dur)
        except Exception:
            t += dur; continue
        off = int(t)
        x_ctx, x_ekg = d[:, :len(ctx)], d[:, len(ctx)]
        xa = notch(np.nan_to_num(x_ctx.astype(float)).mean(axis=1), sf)

        try:
            cl = nk.ecg_clean(np.nan_to_num(x_ekg.astype(float)), sampling_rate=int(sf),
                              method="neurokit")
            _, info = nk.ecg_peaks(cl, sampling_rate=int(sf), method="neurokit",
                                   correct_artifacts=True)
            beats.extend((np.asarray(info["ECG_R_Peaks"], int) / sf + off).tolist())
        except Exception:
            pass

        # per-channel SO troughs (3B needs per-channel events, not the channel average)
        for ci, c in enumerate(ctx):
            xc = notch(signal.detrend(np.nan_to_num(x_ctx[:, ci].astype(float))), sf)
            if np.std(xc) < 1e-9:
                continue
            tr = detect_so_halfwaves(signal.sosfiltfilt(sos_so, xc), sf)
            if len(tr):
                so_t[c].extend((tr / sf + off).tolist())

        clean = ied_clean_mask(xa, sf)
        n_sec = int(len(xa) // int(sf))
        for sos_b, dest in ((sos_fixed, sig_fixed), (sos_fsp, sig_fsp), (sos_swa, swa_1)):
            env2 = np.abs(signal.hilbert(signal.sosfiltfilt(sos_b, xa))) ** 2
            _binned_power(env2, clean, sf, n_sec, total_s, off, dest)

        ke = int(EPOCH * sf)
        for e in range(int(len(xa) // ke)):
            gi = int((off + e * EPOCH) // EPOCH)
            if gi >= n_ep:
                break
            seg = xa[e * ke:(e + 1) * ke]
            ep_dr[gi] = delta_ratio(seg, sf)
            fq, pp = signal.welch(seg, sf, nperseg=int(min(4 * sf, len(seg))))
            m = (fq >= SWA_BAND[0]) & (fq < SWA_BAND[1])
            ep_swa[gi] = float(np.trapezoid(pp[m], fq[m]))
            ep_clean[gi] = float(clean[e * ke:(e + 1) * ke].mean())
        t += dur

    beats = np.array(sorted(beats))
    hr_1 = np.full(total_s, np.nan)
    hr_4 = np.full(int(total_s * FS_RR), np.nan)
    if len(beats) > 20:
        rr = np.diff(beats); good = (rr >= 0.33) & (rr <= 1.5)
        if good.sum() > 20:
            tb, hrv = beats[1:][good], 60.0 / rr[good]
            # 1 Hz grid: interpolate only INSIDE the beat range, leave gaps as NaN
            g1 = np.arange(total_s)
            hr_1 = np.where((g1 >= tb[0]) & (g1 <= tb[-1]), np.interp(g1, tb, hrv), np.nan)
            g4 = np.arange(0, total_s, 1.0 / FS_RR)
            spl = interpolate.CubicSpline(tb, hrv, extrapolate=False)
            hr_4 = spl(g4)

    os.makedirs(OUT, exist_ok=True)
    payload = dict(status="ok", subject=name, sf=sf, night_s=night, hours=hours,
                   cortical_chans=np.array(ctx), ekg=ekg, fsp=fsp, fsp_is_real_peak=fsp_real,
                   sigma_fixed=sig_fixed, sigma_fsp=sig_fsp, swa=swa_1,
                   hr_1=hr_1, hr_4=hr_4, fs_rr=FS_RR,
                   ep_dr=ep_dr, ep_swa=ep_swa, ep_clean=ep_clean, epoch_s=EPOCH,
                   beats=beats)
    for c in ctx:
        payload[f"so_t_{c}"] = np.array(sorted(so_t[c]))
    np.savez_compressed(fp, **payload)
    frac = float(np.isfinite(sig_fsp).mean())
    print(f"[{name}] cached in {time.time()-t_start:.0f}s | sigma coverage {frac:.1%} | "
          f"{len(beats)} beats | SOs {sum(len(v) for v in so_t.values())} -> {fp}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default=",".join(map(str, COHORT)))
    ap.add_argument("--hours", type=float, default=7.0)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    for n in [int(x) for x in a.subjects.split(",") if x.strip()]:
        try:
            run(n, a.hours, a.force)
        except Exception as e:
            print(f"[HUP{n}] ERROR {type(e).__name__}: {e}", flush=True)
            traceback.print_exc()


if __name__ == "__main__":
    main()
