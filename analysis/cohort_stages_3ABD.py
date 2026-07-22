"""
Tests 3A / 3B / 3D per sleep stage (N2-like vs N3-like), across the HUP phase-II cohort.

STAGING CAVEAT, STATED UP FRONT: these are intracranial recordings with no EOG/EMG and no scalp
derivations, so true AASM N2/N3 scoring is IMPOSSIBLE here. What this does instead is the standard
data-driven substitute: score 30 s epochs as NREM by relative delta, then split NREM into two
classes by fitting a 2-component Gaussian mixture to log slow-wave (0.5-4 Hz) power within subject.
The high-SWA class is reported as N3-like, the low-SWA class as N2-like. That is the physiological
axis the AASM N2/N3 boundary tracks, but it is NOT scored staging. ds003848 (which has EOG+EMG) is
where this can be validated properly.

ARCHITECTURE -- one streaming pass over the night per subject. The data pull dominates runtime, so
the night is streamed once in 10-min chunks and only DERIVED features are kept:
    per 30 s epoch : delta ratio, slow-wave power, clean fraction, phase-amplitude histograms
    at 1 Hz        : IED-masked sigma power, heart rate
    events         : slow-oscillation trough times
Staging needs the whole night's SWA distribution, so it is done after the pass; because the
per-epoch histograms and event times are retained, all three tests are then computed offline per
stage without re-pulling anything.

WHY THIS RECOVERS SUBJECTS: 3B (SO->heartbeat) and 3D (SO->spindle) are EVENT-based and need no
temporal continuity, so they use every epoch of a stage. Only 3A (coherence) needs a contiguous
block. An earlier version demanded contiguity for all three and skipped most of the cohort.

    .venv/bin/python analysis/cohort_stages_3ABD.py [--subjects 165,157,...] [--hours 7]
"""
import argparse, json, os, re, traceback
import numpy as np
from scipy import signal
import neurokit2 as nk

from infraslow_rr_sigma_coherence import sess, get, pull_continuous, notch, ROOT
from results_3A_tutorial_style import ied_clean_mask, robust_z
from cohort_3A_cortical import COHORT, cortical_channels, delta_ratio, find_night, tort_mi

OUT = os.path.join(ROOT, "outputs", "cohort_stages_3ABD")
EPOCH = 30.0
FS_P = 1.0
CHUNK_S = 600.0
SWA_BAND = (0.5, 4.0)
SO_BAND = (0.5, 1.25)
INFRA = (0.01, 0.04)
F_TARGET = 0.02
NB = 18                 # phase bins for Tort MI
N_SUR = 20              # surrogates for 3D / 3B
NREM_DR = 0.90
MIN_3A_MIN = 20.0       # min contiguous minutes for a coherence estimate; short blocks give few
                        # Welch segments (high per-window threshold), but N2/N3 are length-matched
                        # within subject so the ~1/K bias is identical on both sides of the paired test
ALPHA = 0.05
POOLED_3A_CAP_MIN = 55.0  # common cap so K (and thus the ~1/K coherence floor) matches across subjects


# ---------------------------------------------------------------- streaming pass
def fsp_from(x_mtl, sf):
    """Individual fast-spindle peak. The raw PSD falls as 1/f, so a plain argmax over 11-16 Hz just
    returns the low edge; whiten by removing the log-log 1/f trend first and take a genuine LOCAL
    peak. Returns (freq, found_real_peak)."""
    ps = []
    for j in range(x_mtl.shape[1]):
        f, p = signal.welch(notch(np.nan_to_num(x_mtl[:, j].astype(float)), sf), sf, nperseg=int(8 * sf))
        ps.append(p)
    p = np.mean(ps, axis=0)
    fit = (f >= 2) & (f <= 30) & (p > 0)
    co = np.polyfit(np.log(f[fit]), np.log(p[fit]), 1)          # 1/f background
    resid = np.log(p[fit]) - np.polyval(co, np.log(f[fit]))
    ff = f[fit]
    band = (ff >= 10.5) & (ff <= 16.0)
    pk, _ = signal.find_peaks(resid[band])
    if len(pk):
        return float(ff[band][pk[np.argmax(resid[band][pk])]]), True
    return 13.0, False                                           # typical fast-spindle peak fallback


def band_sos(band, sf, order=4):
    return signal.butter(order, list(band), btype="band", fs=sf, output="sos")


def stream_night(ds, idx, n_ctx, sf, start_s, hours, fsp, rng):
    """One pass over the night; returns accumulated derived features."""
    total_s = int(hours * 3600)
    n_ep = int(total_s // EPOCH)
    ep = dict(dr=np.full(n_ep, np.nan), swa=np.full(n_ep, np.nan), clean=np.zeros(n_ep))
    pac_sum = np.zeros((n_ep, NB)); pac_cnt = np.zeros((n_ep, NB))
    pac_sur = np.zeros((n_ep, N_SUR, NB))
    sig_1 = np.full(total_s, np.nan)
    beats = []
    so_t = []

    sos_sig = band_sos((fsp - 1, fsp + 1), sf)
    sos_swa = band_sos(SWA_BAND, sf, 3)
    sos_so = band_sos(SO_BAND, sf, 3)

    t = 0.0
    while t < total_s:
        dur = min(CHUNK_S, total_s - t)
        try:
            d = pull_continuous(ds, idx, start_s + t, dur)
        except Exception:
            t += dur; continue
        x_ctx, x_ekg = d[:, :n_ctx], d[:, n_ctx]
        xa = np.nan_to_num(x_ctx.astype(float)).mean(axis=1)
        xa = notch(xa, sf)
        off = int(t)                                   # chunk start, seconds

        # --- heart beats -> global times
        try:
            cl = nk.ecg_clean(np.nan_to_num(x_ekg.astype(float)), sampling_rate=int(sf), method="neurokit")
            _, info = nk.ecg_peaks(cl, sampling_rate=int(sf), method="neurokit", correct_artifacts=True)
            beats.extend((np.asarray(info["ECG_R_Peaks"], int) / sf + off).tolist())
        except Exception:
            pass

        # --- SO troughs -> global times
        f_so = signal.sosfiltfilt(sos_so, xa)
        neg, _ = signal.find_peaks(-f_so, distance=int(0.5 * sf))
        if len(neg) > 10:
            thr = np.percentile(-f_so[neg], 75)
            so_t.extend((neg[-f_so[neg] >= thr] / sf + off).tolist())
        ph = np.angle(signal.hilbert(f_so))

        # --- sigma amplitude + IED mask
        amp = np.abs(signal.hilbert(signal.sosfiltfilt(sos_sig, xa)))
        clean = ied_clean_mask(xa, sf)
        env2 = amp ** 2
        swa_env = np.abs(signal.hilbert(signal.sosfiltfilt(sos_swa, xa))) ** 2

        # --- 1 Hz sigma power (IED-masked bin means)
        k = int(sf)
        n_sec = int(len(xa) // k)
        e2 = env2[:n_sec * k].reshape(n_sec, k); cm = clean[:n_sec * k].reshape(n_sec, k).astype(float)
        num, den = (e2 * cm).sum(1), cm.sum(1)
        vals = np.where(den >= 0.5 * k, num / np.maximum(den, 1e-12), np.nan)
        sl = slice(off, min(off + n_sec, total_s))
        sig_1[sl] = vals[:sl.stop - sl.start]

        # --- per-epoch features + PAC histograms
        ke = int(EPOCH * sf)
        for e in range(int(len(xa) // ke)):
            gi = int((off + e * EPOCH) // EPOCH)
            if gi >= n_ep:
                break
            s0, s1 = e * ke, (e + 1) * ke
            seg = xa[s0:s1]
            ep["dr"][gi] = delta_ratio(seg, sf)
            fq, pp = signal.welch(seg, sf, nperseg=int(min(4 * sf, len(seg))))
            m = (fq >= SWA_BAND[0]) & (fq < SWA_BAND[1])
            ep["swa"][gi] = float(np.trapezoid(pp[m], fq[m]))
            ep["clean"][gi] = float(clean[s0:s1].mean())
            pe, ae = ph[s0:s1], amp[s0:s1]
            bi = np.clip(((pe + np.pi) / (2 * np.pi) * NB).astype(int), 0, NB - 1)
            np.add.at(pac_sum[gi], bi, ae); np.add.at(pac_cnt[gi], bi, 1.0)
            for s in range(N_SUR):
                np.add.at(pac_sur[gi, s], bi, np.roll(ae, rng.randint(int(sf), max(int(sf) + 1, len(ae) - 1))))
        t += dur

    # heart rate on the 1 Hz grid
    beats = np.array(sorted(beats))
    hr_1 = np.full(total_s, np.nan)
    if len(beats) > 10:
        rr = np.diff(beats); good = (rr >= 0.33) & (rr <= 1.5)
        if good.sum() > 10:
            hr_1 = np.interp(np.arange(total_s), beats[1:][good], 60.0 / rr[good])
    return dict(ep=ep, pac_sum=pac_sum, pac_cnt=pac_cnt, pac_sur=pac_sur,
                sig_1=sig_1, hr_1=hr_1, so_t=np.array(sorted(so_t)), n_ep=n_ep)


# ---------------------------------------------------------------- staging
def nrem_mask_adaptive(dr, clean):
    """NREM epochs, thresholded WITHIN subject.

    An absolute delta-ratio cutoff does not transfer: it was tuned on one subject and discarded
    most epochs in others (101 / 48 / 354 NREM epochs where ~800 were expected), collapsing the
    staging. Fit 2 classes to this subject's own delta-ratio distribution instead."""
    ok = np.isfinite(dr)
    if ok.sum() < 60:
        return np.zeros(len(dr), bool)
    v = dr[ok].reshape(-1, 1)
    try:
        from sklearn.mixture import GaussianMixture
        g = GaussianMixture(n_components=2, random_state=0, n_init=3).fit(v)
        hi = int(np.argmax(g.means_.ravel()))
        m = np.zeros(len(dr), bool)
        m[np.where(ok)[0]] = (g.predict(v) == hi)
    except Exception:
        m = ok & (dr >= np.nanmedian(dr))
    if not (0.15 <= m.mean() <= 0.95):            # implausible split -> percentile fallback
        m = ok & (dr >= np.nanpercentile(dr[ok], 40))
    m = m & (clean >= 0.5)
    # CONSOLIDATE. Per-epoch classification flickers in and out of NREM every few epochs, which
    # leaves no contiguous block long enough for a 0.02 Hz coherence estimate (628 NREM epochs but
    # no 20-min run). Real sleep scoring smooths over time; a 5-epoch (2.5 min) majority filter
    # fills single-epoch dropouts and removes isolated epochs.
    k = 5
    return np.convolve(m.astype(float), np.ones(k) / k, mode="same") >= 0.5


def stage_epochs(ep):
    """NREM epochs -> {N2-like, N3-like} by 2-component GMM on log slow-wave power."""
    dr, swa, clean = ep["dr"], ep["swa"], ep["clean"]
    nrem = nrem_mask_adaptive(dr, clean) & np.isfinite(swa) & (swa > 0)
    lab = np.full(len(dr), "", dtype=object)
    if nrem.sum() < 40:
        return lab, nrem, None
    v = np.log(swa[nrem]).reshape(-1, 1)
    try:
        from sklearn.mixture import GaussianMixture
        g = GaussianMixture(n_components=2, random_state=0, n_init=3).fit(v)
        cl = g.predict(v)
        hi = int(np.argmax(g.means_.ravel()))
        sep = float(abs(np.diff(g.means_.ravel())[0]) / np.sqrt(g.covariances_.ravel().mean()))
    except Exception:
        thr = np.median(v); cl = (v.ravel() >= thr).astype(int); hi = 1
        sep = float("nan")
    names = np.where(cl == hi, "N3", "N2")
    lab[np.where(nrem)[0]] = names
    return lab, nrem, sep


def dominant_block(sel, nrem, win=11, frac_thr=0.60):
    """Longest STAGE-DOMINANT block, not stage-pure.

    N2 and N3 alternate on a timescale of minutes, so no 40+ min block is ever purely one stage --
    demanding purity returns nothing. Coherence needs continuity, so instead smooth the stage label
    over ~5 min and keep periods where the stage is the local majority of NREM. Returns
    (start_epoch, end_epoch, minutes, purity) with purity = actual fraction of the block in-stage."""
    s = sel.astype(float); nz = nrem.astype(float)
    k = np.ones(win)
    num = np.convolve(s, k, mode="same"); den = np.convolve(nz, k, mode="same")
    frac = np.divide(num, np.maximum(den, 1e-9))
    dom = (frac >= frac_thr) & nrem
    best = (0, 0)
    i = 0
    while i < len(dom):
        if not dom[i]:
            i += 1; continue
        j = i
        while j < len(dom) and dom[j]:
            j += 1
        if j - i > best[1] - best[0]:
            best = (i, j)
        i = j
    i0, i1 = best
    if i1 <= i0:
        return None
    seg_nrem = nrem[i0:i1]
    purity = float(sel[i0:i1][seg_nrem].mean()) if seg_nrem.any() else 0.0
    return i0, i1, (i1 - i0) * EPOCH / 60.0, purity


# ---------------------------------------------------------------- tests
def msc_block(hr, sig, nperseg_cap=256):
    m = np.isfinite(hr) & np.isfinite(sig)
    if m.sum() < 600:
        return None
    hr, sig = hr[m], sig[m]
    sos = signal.butter(3, 0.005, btype="high", fs=FS_P, output="sos")
    a, b = signal.sosfiltfilt(sos, hr - hr.mean()), signal.sosfiltfilt(sos, sig - sig.mean())
    nper = int(min(nperseg_cap, (len(a) // 4) // 2 * 2))
    if nper < 128:
        return None
    f, cxy = signal.coherence(a, b, fs=FS_P, nperseg=nper, noverlap=nper // 2)
    K = int((len(a) - nper // 2) // (nper // 2))
    band = (f >= INFRA[0]) & (f <= INFRA[1])
    return dict(at=float(cxy[np.argmin(np.abs(f - F_TARGET))]), bmax=float(cxy[band].max()),
                K=K, crit=1 - ALPHA ** (1 / max(K - 1, 1)), n_sec=int(len(a)), nperseg=nper)


def mi_from(sums, cnts):
    m = np.divide(sums, np.maximum(cnts, 1e-12))
    if (cnts <= 0).any() or m.sum() <= 0:
        return np.nan
    p = m / m.sum()
    return float((np.log(NB) + np.sum(p * np.log(p + 1e-12))) / np.log(NB))


def test_3D(F, sel):
    if sel.sum() < 20:
        return dict(mi=None, mi_z=None, n_epochs=int(sel.sum()))
    mi = mi_from(F["pac_sum"][sel].sum(0), F["pac_cnt"][sel].sum(0))
    null = [mi_from(F["pac_sur"][sel, s].sum(0), F["pac_cnt"][sel].sum(0)) for s in range(N_SUR)]
    null = np.array([v for v in null if np.isfinite(v)])
    z = float((mi - null.mean()) / (null.std() + 1e-12)) if len(null) > 5 and np.isfinite(mi) else None
    return dict(mi=None if not np.isfinite(mi) else float(mi), mi_z=z, n_epochs=int(sel.sum()))


def test_3B(F, sel, rng, half=10):
    hr = F["hr_1"]
    ep_ok = np.where(sel)[0]
    if len(ep_ok) < 20 or not np.isfinite(hr).any():
        return dict(n_so=0, modulation=None, z=None)
    keep = set(ep_ok.tolist())
    tro = np.array([t for t in F["so_t"] if int(t // EPOCH) in keep], int)
    tro = tro[(tro >= half) & (tro < len(hr) - half)]
    tro = np.array([t for t in tro if np.isfinite(hr[t - half:t + half]).all()])
    if len(tro) < 30:
        return dict(n_so=int(len(tro)), modulation=None, z=None)
    seg = np.stack([hr[t - half:t + half] for t in tro])
    curve = seg.mean(0) - seg.mean()
    mod = float(curve.max() - curve.min())
    # Surrogate triggers must be drawn from the SAME STAGE as the real ones. Drawing them from the
    # whole night (the previous behaviour) makes the null measure the stage-vs-night mean-HR offset
    # instead of SO-locking -- see the docstring of event_3B_mednick.so_triggered for the proof.
    in_stage = np.zeros(len(hr), bool)
    for e in ep_ok:
        in_stage[int(e * EPOCH):int((e + 1) * EPOCH)] = True
    ok_t = np.where(np.isfinite(hr) & in_stage)[0]
    ok_t = ok_t[(ok_t >= half) & (ok_t < len(hr) - half)]
    if len(ok_t) < 100:
        return dict(n_so=int(len(tro)), modulation=None, z=None)
    null = []
    for _ in range(200):
        r = rng.choice(ok_t, size=len(tro), replace=True)
        s2 = np.stack([hr[t - half:t + half] for t in r])
        c2 = s2.mean(0) - s2.mean(); null.append(c2.max() - c2.min())
    null = np.array(null)
    return dict(n_so=int(len(tro)), modulation=mod,
                z=float((mod - null.mean()) / (null.std() + 1e-12)),
                curve=[float(v) for v in curve])


# ---------------------------------------------------------------- per subject
def run_subject(n, hours):
    name = f"HUP{n}_phaseII"
    fp = os.path.join(OUT, f"{name}.json")
    if os.path.exists(fp):
        print(f"[{name}] cached", flush=True); return json.load(open(fp))
    rng = np.random.RandomState(0)
    rec = dict(subject=name)
    s = sess(); ds = s.open_dataset(name)
    labels = ds.get_channel_labels(); lab_idx = {l: i for i, l in enumerate(labels)}
    d0 = ds.get_time_series_details(labels[0]); sf = d0.sample_rate
    total_h = (getattr(d0, "duration", 0) or 0) / 3.6e9
    ekg = next((l for l in labels if l.upper().startswith(("EKG", "ECG"))), None)
    ctx = cortical_channels(labels)
    if ekg is None or len(ctx) < 3:
        rec.update(status="skip", reason=f"ekg={ekg} n_cortical={len(ctx)}")
        json.dump(rec, open(fp, "w"), indent=2); print(f"[{name}] SKIP {rec['reason']}", flush=True); return rec
    night = find_night(ds, lab_idx, ctx[0], sf, total_h)
    if night is None:
        rec.update(status="skip", reason="no night found")
        json.dump(rec, open(fp, "w"), indent=2); print(f"[{name}] SKIP no night", flush=True); return rec

    idx = [lab_idx[c] for c in ctx] + [lab_idx[ekg]]
    samp = pull_continuous(ds, idx, night, 300.0)
    fsp, fsp_real = fsp_from(samp[:, :len(ctx)], sf)
    print(f"[{name}] {sf:.0f} Hz | cortical={ctx} | night {night/3600:.1f} h | FSP {fsp:.2f} Hz | "
          f"streaming {hours} h", flush=True)

    F = stream_night(ds, idx, len(ctx), sf, night, hours, fsp, rng)
    lab, nrem, sep = stage_epochs(F["ep"])
    n2, n3 = (lab == "N2"), (lab == "N3")
    print(f"[{name}] NREM epochs {int(nrem.sum())} -> N2-like {int(n2.sum())}, N3-like {int(n3.sum())} "
          f"(GMM separation {sep:.2f})" if sep is not None else f"[{name}] staging failed", flush=True)
    rec.update(status="ok", sf=sf, cortical_chans=ctx, ekg=ekg, night_h=night / 3600,
               fast_spindle_peak_hz=fsp, fsp_is_real_peak=fsp_real, hours_streamed=hours,
               n_nrem_epochs=int(nrem.sum()), n_N2=int(n2.sum()), n_N3=int(n3.sum()),
               gmm_separation=sep)

    # 3A needs contiguity; length-match the two stages so the MSC bias (~1/K) is equal
    blocks = {}
    for st, sel in (("N2", n2), ("N3", n3)):
        b = dominant_block(sel, nrem)
        if b:
            blocks[st] = b
    match_min = None
    if len(blocks) == 2:
        match_min = min(blocks["N2"][2], blocks["N3"][2])
        if match_min < MIN_3A_MIN:
            match_min = None
    if blocks:
        print(f"[{name}] 3A dominant blocks: " + ", ".join(
            f"{st} {b[2]:.0f} min (purity {b[3]:.2f})" for st, b in blocks.items())
            + (f" -> matched at {match_min:.0f} min" if match_min else " -> too short for 3A"), flush=True)

    for st, sel in (("N2", n2), ("N3", n3)):
        out = dict(n_epochs=int(sel.sum()))
        out.update({f"3D_{k}": v for k, v in test_3D(F, sel).items()})
        out.update({f"3B_{k}": v for k, v in test_3B(F, sel, rng).items() if k != "curve"})
        if match_min and st in blocks:
            i0 = blocks[st][0]; nsec = int(match_min * 60)
            sl = slice(int(i0 * EPOCH), int(i0 * EPOCH) + nsec)
            m = msc_block(F["hr_1"][sl], F["sig_1"][sl])
            out["3A"] = m
            out["3A_block_min"] = match_min
        else:
            out["3A"] = None
            out["3A_block_min"] = blocks.get(st, (0, 0, 0))[2] if st in blocks else 0.0
        rec[st] = out
        a = out["3A"]
        mi_z = out["3D_mi_z"]; b_z = out["3B_z"]
        a_txt = "n/a" if a is None else f"{a['at']:.3f} (crit {a['crit']:.3f}, K={a['K']})"
        mi_txt = "n/a" if mi_z is None else f"{mi_z:.1f}"
        b_txt = "n/a" if b_z is None else f"{b_z:.1f}"
        print(f"[{name}]   {st}: 3D MI_z={mi_txt} (n_ep={out['3D_n_epochs']}) | "
              f"3B z={b_txt} (n_SO={out['3B_n_so']}) | 3A {a_txt}", flush=True)

    # 3A on POOLED NREM as well. N2/N3 alternate faster than a 0.02 Hz coherence estimate needs, so
    # stage-resolved 3A is often unavailable; the pooled block is long enough and keeps 3A reportable
    # (it just cannot speak to the stage contrast).
    nb = dominant_block(nrem, nrem)
    rec["NREM_pooled"] = dict(block_min=(nb[2] if nb else 0.0), purity=(nb[3] if nb else 0.0))
    if nb and nb[2] >= MIN_3A_MIN:
        # CAP at a common duration across subjects. Coherence's noise floor is ~1/K, so a subject
        # with a 258-min block (K=113, chance 0.026) is not comparable to one with 58 min
        # (K=23, chance 0.127) -- without capping, the cohort comparison measures block length.
        cap_s = int(min(nb[2], POOLED_3A_CAP_MIN) * 60)
        rec["NREM_pooled"]["used_min"] = cap_s / 60.0
        sl = slice(int(nb[0] * EPOCH), int(nb[0] * EPOCH) + cap_s)
        rec["NREM_pooled"]["3A"] = msc_block(F["hr_1"][sl], F["sig_1"][sl])
        rec["NREM_pooled"]["3B"] = {k: v for k, v in test_3B(F, nrem, rng).items() if k != "curve"}
        rec["NREM_pooled"]["3D"] = test_3D(F, nrem)
        a = rec["NREM_pooled"]["3A"]
        pooled_txt = "n/a" if a is None else f"{a['at']:.3f} (crit {a['crit']:.3f}, K={a['K']})"
        print(f"[{name}]   pooled NREM ({nb[2]:.0f} min, purity {nb[3]:.2f}): 3A {pooled_txt}", flush=True)
    else:
        rec["NREM_pooled"]["3A"] = None
    json.dump(rec, open(fp, "w"), indent=2)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default=",".join(map(str, COHORT)))
    ap.add_argument("--hours", type=float, default=7.0)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for n in [int(x) for x in a.subjects.split(",") if x.strip()]:
        try:
            rows.append(run_subject(n, a.hours))
        except Exception as e:
            print(f"[HUP{n}] ERROR {type(e).__name__}: {e}", flush=True); traceback.print_exc()
    ok = [r for r in rows if r.get("status") == "ok"]
    print(f"\n=== COHORT n={len(ok)} ===")
    if ok:
        import csv
        flat = []
        for r in ok:
            row = {k: r[k] for k in ("subject", "sf", "fast_spindle_peak_hz", "n_N2", "n_N3", "gmm_separation")}
            for st in ("N2", "N3"):
                d = r[st]
                row[f"{st}_3D_mi"] = d["3D_mi"]; row[f"{st}_3D_mi_z"] = d["3D_mi_z"]
                row[f"{st}_3B_z"] = d["3B_z"]; row[f"{st}_3B_n_so"] = d["3B_n_so"]
                row[f"{st}_3A_at0p02"] = (d["3A"] or {}).get("at")
                row[f"{st}_3A_crit"] = (d["3A"] or {}).get("crit")
                row[f"{st}_3A_K"] = (d["3A"] or {}).get("K")
            flat.append(row)
        keys = sorted({k for r in flat for k in r})
        with open(os.path.join(OUT, "cohort_stages_3ABD.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(flat)
        from scipy import stats
        for tag, key in (("3A coherence@0.02", "3A_at0p02"), ("3B SO->HR z", "3B_z"), ("3D SO->spindle MI_z", "3D_mi_z")):
            n2 = np.array([r[f"N2_{key}"] for r in flat], float)
            n3 = np.array([r[f"N3_{key}"] for r in flat], float)
            m = np.isfinite(n2) & np.isfinite(n3)
            if m.sum() >= 3:
                try:
                    p = stats.wilcoxon(n2[m], n3[m])[1]
                except Exception:
                    p = float("nan")
                print(f"  {tag:22s} N2 median {np.median(n2[m]):8.3f} | N3 median {np.median(n3[m]):8.3f} "
                      f"| paired n={int(m.sum())} Wilcoxon p={p:.4f} | N2>N3 in {int((n2[m]>n3[m]).sum())}")
            else:
                print(f"  {tag:22s} insufficient paired data (n={int(m.sum())})")
    print(f"\nper-subject JSON -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
