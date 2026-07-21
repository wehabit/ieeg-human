"""Is the pooled-3A excess SPECIFIC to ~0.02 Hz, or does it appear at any frequency?

The cohort found 7/23 subjects exceeding their own alpha=0.05 coherence threshold at the
pre-specified 0.02 Hz point (binomial p=4.7e-04 against a calibrated 6% false-positive rate).
That establishes coupling exists -- but not that it sits at the LC infraslow rhythm. If the same
excess appears at 0.05 Hz, 0.1 Hz, 0.2 Hz, then it is generic shared structure between sigma power
and heart rate, not the ~50 s fingerprint of Lecci / Osorio-Forero.

This re-streams each night and stores the FULL coherence spectrum for the pooled-NREM block, so the
exceedance rate can be counted at every frequency bin and 0.02 Hz compared against the rest.
Phase-amplitude coupling is skipped entirely (3D raw MI was ~6e-5, negligible), which makes this
pass much faster than the full cohort run.

    .venv/bin/python analysis/frequency_specificity_3A.py [--hours 7]
"""
import argparse, json, os, traceback
import numpy as np
from scipy import signal
import neurokit2 as nk

from infraslow_rr_sigma_coherence import sess, pull_continuous, notch, ROOT
from results_3A_tutorial_style import ied_clean_mask
from cohort_3A_cortical import COHORT, cortical_channels, delta_ratio, find_night
from cohort_stages_3ABD import (fsp_from, band_sos, stage_epochs, dominant_block,
                                EPOCH, FS_P, CHUNK_S, SWA_BAND, ALPHA,
                                MIN_3A_MIN, POOLED_3A_CAP_MIN)

OUT = os.path.join(ROOT, "outputs", "freq_specificity_3A")


def stream_light(ds, idx, n_ctx, sf, start_s, hours, fsp):
    """Streaming pass keeping ONLY what 3A needs: per-epoch staging features, 1 Hz sigma power,
    and heart beats. No phase-amplitude histograms, no surrogates."""
    total_s = int(hours * 3600)
    n_ep = int(total_s // EPOCH)
    ep = dict(dr=np.full(n_ep, np.nan), swa=np.full(n_ep, np.nan), clean=np.zeros(n_ep))
    sig_1 = np.full(total_s, np.nan)
    beats = []
    sos_sig = band_sos((fsp - 1, fsp + 1), sf)
    t = 0.0
    while t < total_s:
        dur = min(CHUNK_S, total_s - t)
        try:
            d = pull_continuous(ds, idx, start_s + t, dur)
        except Exception:
            t += dur; continue
        x_ctx, x_ekg = d[:, :n_ctx], d[:, n_ctx]
        xa = notch(np.nan_to_num(x_ctx.astype(float)).mean(axis=1), sf)
        off = int(t)
        try:
            cl = nk.ecg_clean(np.nan_to_num(x_ekg.astype(float)), sampling_rate=int(sf), method="neurokit")
            _, info = nk.ecg_peaks(cl, sampling_rate=int(sf), method="neurokit", correct_artifacts=True)
            beats.extend((np.asarray(info["ECG_R_Peaks"], int) / sf + off).tolist())
        except Exception:
            pass
        env2 = np.abs(signal.hilbert(signal.sosfiltfilt(sos_sig, xa))) ** 2
        clean = ied_clean_mask(xa, sf)
        k = int(sf); n_sec = int(len(xa) // k)
        e2 = env2[:n_sec * k].reshape(n_sec, k); cm = clean[:n_sec * k].reshape(n_sec, k).astype(float)
        num, den = (e2 * cm).sum(1), cm.sum(1)
        vals = np.where(den >= 0.5 * k, num / np.maximum(den, 1e-12), np.nan)
        sl = slice(off, min(off + n_sec, total_s))
        sig_1[sl] = vals[:sl.stop - sl.start]
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
            ep["clean"][gi] = float(clean[e * ke:(e + 1) * ke].mean())
        t += dur
    beats = np.array(sorted(beats))
    hr_1 = np.full(total_s, np.nan)
    if len(beats) > 10:
        rr = np.diff(beats); good = (rr >= 0.33) & (rr <= 1.5)
        if good.sum() > 10:
            hr_1 = np.interp(np.arange(total_s), beats[1:][good], 60.0 / rr[good])
    return ep, sig_1, hr_1


def msc_full(hr, sig, nperseg_cap=256):
    m = np.isfinite(hr) & np.isfinite(sig)
    if m.sum() < 600:
        return None
    hr, sig = hr[m], sig[m]
    sos = signal.butter(3, 0.005, btype="high", fs=FS_P, output="sos")
    a = signal.sosfiltfilt(sos, hr - hr.mean()); b = signal.sosfiltfilt(sos, sig - sig.mean())
    nper = int(min(nperseg_cap, (len(a) // 4) // 2 * 2))
    if nper < 128:
        return None
    f, cxy = signal.coherence(a, b, fs=FS_P, nperseg=nper, noverlap=nper // 2)
    K = int((len(a) - nper // 2) // (nper // 2))
    return dict(f=f.tolist(), cxy=cxy.tolist(), K=K,
                crit=1 - ALPHA ** (1 / max(K - 1, 1)), n_sec=int(len(a)))


def run(n, hours):
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
    fsp, real = fsp_from(pull_continuous(ds, idx, night, 300.0)[:, :len(ctx)], sf)
    print(f"[{name}] {sf:.0f} Hz | FSP {fsp:.2f} | streaming {hours} h", flush=True)
    ep, sig_1, hr_1 = stream_light(ds, idx, len(ctx), sf, night, hours, fsp)
    lab, nrem, sep = stage_epochs(ep)
    nb = dominant_block(nrem, nrem)
    rec = dict(subject=name, status="ok", sf=sf, fsp=fsp, fsp_real=real,
               n_nrem=int(nrem.sum()), block_min=(nb[2] if nb else 0.0))
    if nb and nb[2] >= MIN_3A_MIN:
        cap = int(min(nb[2], POOLED_3A_CAP_MIN) * 60)
        sl = slice(int(nb[0] * EPOCH), int(nb[0] * EPOCH) + cap)
        rec["msc"] = msc_full(hr_1[sl], sig_1[sl])
        rec["used_min"] = cap / 60.0
    else:
        rec["msc"] = None
    json.dump(rec, open(fp, "w"))
    m = rec["msc"]
    if m:
        f = np.array(m["f"]); c = np.array(m["cxy"])
        at = float(c[np.argmin(np.abs(f - 0.02))])
        print(f"[{name}] block {rec['block_min']:.0f} min, K={m['K']}, crit={m['crit']:.3f}, "
              f"coh@0.02={at:.3f} {'*' if at > m['crit'] else ''}", flush=True)
    else:
        print(f"[{name}] no usable block ({rec['block_min']:.0f} min)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default=",".join(map(str, COHORT)))
    ap.add_argument("--hours", type=float, default=7.0)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    for n in [int(x) for x in a.subjects.split(",") if x.strip()]:
        try:
            run(n, a.hours)
        except Exception as e:
            print(f"[HUP{n}] ERROR {type(e).__name__}: {e}", flush=True); traceback.print_exc()
    print(f"\nJSON -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
