"""
De-circularized validation of the N3 SO->spindle nesting (atlas, 204 Hz).

The sorted heatmap in `nesting_by_patient` double-dips (rows sorted by the same spindle value
that is colour-plotted), so its gradient is not evidence. Here are three non-circular views:

  1. MEAN +/- CI (no sorting)            - the actual claim; bootstrap CI over patients.
  2. SPLIT-HALF patient heatmap          - rank patients by up-state spindle on ODD events,
     display the mean envelope from their EVEN events (+ Spearman odd-vs-even across patients).
     A surviving gradient = stable per-patient trait, not double-dipping.
  3. Event heatmap sorted by SO TROUGH DEPTH (an independent variable) - tests "bigger slow
     waves carry more spindle nesting" (+ Spearman trough-depth vs up-state spindle across events).

Uses only cached atlas N3 clips.
"""
import os
import numpy as np, pandas as pd
from scipy import signal
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import atlas

SF = atlas.SFREQ
HALF = int(1.5 * SF); MIN_DIST = int(0.6 * SF)
SO_BAND = (0.5, 1.25); SPINDLE = (11, 16)
OUT = os.path.join(atlas.ROOT, "outputs", "nesting_by_patient")
MAX_EV_PT = 400          # cap events kept per patient (memory)


def bp(x, lo, hi):
    return signal.sosfiltfilt(signal.butter(3, [lo/(SF/2), hi/(SF/2)], btype="band", output="sos"), x)


def notch60(x):
    b, a = signal.iirnotch(atlas.LINE_HZ, 30, SF); return signal.filtfilt(b, a, x)


def main():
    os.makedirs(OUT, exist_ok=True)
    t = np.arange(-HALF, HALF) / SF
    win = (t >= 0.0) & (t <= 0.75)          # post-trough up-state window
    meta = atlas.load_metadata(normative_only=True)
    roi = meta.dropna(subset=["roi_group"])
    roi = roi[roi.roi_group.isin(atlas.MESIOTEMPORAL)][["pt", "name", "roi_group"]]
    clips = atlas.clip_index()

    per_pt = {}
    pool_env = []; pool_depth = []          # pooled events for panel 3
    for pt in sorted(roi.pt.unique()):
        chans = list(roi[roi.pt == pt]["name"])
        cached = [j for j in sorted(clips[(clips.pt == pt) & (clips.state == "N3")].idx)
                  if os.path.exists(os.path.join(atlas.EEGDIR, f"{pt}_N3_{j}.edf"))]
        if not cached:
            continue
        envs = []; depths = []
        for j in cached:
            try:
                data, names, sf = atlas.read_clip(pt, "N3", j)
            except Exception:
                continue
            nidx = {nm: k for k, nm in enumerate(names)}
            for chan in chans:
                if chan not in nidx:
                    continue
                x = data[nidx[chan]]
                if not np.isfinite(x).all() or np.std(x) < 1e-9:
                    continue
                x = notch60(signal.detrend(x)); so = bp(x, *SO_BAND)
                spz = np.abs(signal.hilbert(bp(x, *SPINDLE)))
                spz = (spz - spz.mean()) / (spz.std() + 1e-12)
                idx, _ = signal.find_peaks(-so, height=so.std(), distance=MIN_DIST)
                idx = idx[(idx > HALF) & (idx < len(x) - HALF)]
                for e in idx:
                    envs.append(spz[e-HALF:e+HALF]); depths.append(float(-so[e]))
        if len(envs) < 20:
            continue
        envs = np.array(envs); depths = np.array(depths)
        if len(envs) > MAX_EV_PT:                       # cap per patient
            sel = np.random.default_rng(0).choice(len(envs), MAX_EV_PT, replace=False)
            envs, depths = envs[sel], depths[sel]
        odd = np.arange(len(envs)) % 2 == 1; even = ~odd
        per_pt[pt] = {
            "all": envs.mean(0), "n": len(envs),
            "env_even": envs[even].mean(0),
            "upA": envs[odd][:, win].mean(),            # rank statistic (odd)
            "upB": envs[even][:, win].mean(),           # held-out (even)
        }
        pool_env.append(envs); pool_depth.append(depths)
        print(f"  {pt}: {len(envs)} events", flush=True)

    pts = list(per_pt)
    print(f"[validation] {len(pts)} patients")

    # ---------- panel 1: mean +/- bootstrap CI over patients ----------
    allmean = np.vstack([per_pt[p]["all"] for p in pts])
    rng = np.random.default_rng(0)
    bs = allmean[rng.integers(0, len(pts), (2000, len(pts)))].mean(1)   # (2000, time)
    lo, hi = np.percentile(bs, 2.5, 0), np.percentile(bs, 97.5, 0)
    grand = allmean.mean(0)
    up_ci = np.percentile(bs[:, win].mean(1), [2.5, 97.5])

    # ---------- panel 2: split-half ----------
    upA = np.array([per_pt[p]["upA"] for p in pts])
    upB = np.array([per_pt[p]["upB"] for p in pts])
    r_split, p_split = spearmanr(upA, upB)
    order = [pts[i] for i in np.argsort(upA)[::-1]]     # rank by ODD, display EVEN
    Meven = np.vstack([per_pt[p]["env_even"] for p in order])

    # ---------- panel 3: events sorted by SO trough depth (independent) ----------
    penv = np.vstack(pool_env); pdep = np.concatenate(pool_depth)
    r_depth, p_depth = spearmanr(pdep, penv[:, win].mean(1))
    if len(penv) > 3000:
        sel = rng.choice(len(penv), 3000, replace=False); penv, pdep = penv[sel], pdep[sel]
    dorder = np.argsort(pdep)[::-1]                     # deepest troughs on top
    penv = penv[dorder]

    print(f"[validation] split-half Spearman(odd,even up-state) r={r_split:.3f} p={p_split:.3g}")
    print(f"[validation] SO-depth vs up-state spindle  Spearman r={r_depth:.3f} p={p_depth:.3g}")
    print(f"[validation] grand up-state spindle-z mean={grand[win].mean():.3f} 95%CI[{up_ci[0]:.3f},{up_ci[1]:.3f}]")

    fig = plt.figure(figsize=(15, 5))
    gs = fig.add_gridspec(1, 3, wspace=0.3)
    # 1
    a0 = fig.add_subplot(gs[0]); a0b = a0.twinx()
    so_grand = None  # SO wave not stored per patient here; show spindle mean+CI only
    a0.plot(t, grand, color="tab:orange", lw=2, label="spindle env (z), mean")
    a0.fill_between(t, lo, hi, color="tab:orange", alpha=0.25, label="95% CI (over patients)")
    a0.axvline(0, color="k", ls=":", lw=0.8); a0.axhline(0, color="0.7", lw=0.6)
    a0.axvspan(0, 0.75, color="tab:blue", alpha=0.06)
    a0.set_xlabel("time from SO trough (s)"); a0.set_ylabel("spindle env (z)")
    a0.set_title(f"1. Mean ± CI (NO sorting)\nup-state z={grand[win].mean():.3f}, "
                 f"CI[{up_ci[0]:.3f},{up_ci[1]:.3f}]"); a0.legend(fontsize=7)
    a0b.set_yticks([])
    # 2
    a1 = fig.add_subplot(gs[1])
    im1 = a1.imshow(Meven, aspect="auto", cmap="RdBu_r", vmin=-0.6, vmax=0.6,
                    extent=[t[0], t[-1], len(order), 0])
    a1.axvline(0, color="k", ls=":", lw=0.8); a1.set_xlabel("time from SO trough (s)")
    a1.set_ylabel("patient (ranked by ODD events)")
    a1.set_title(f"2. Split-half: rank on ODD, show EVEN\nSpearman(odd,even) r={r_split:.2f} "
                 f"p={p_split:.1e}"); fig.colorbar(im1, ax=a1, label="spindle z")
    # 3
    a2 = fig.add_subplot(gs[2])
    im2 = a2.imshow(penv, aspect="auto", cmap="RdBu_r", vmin=-1.2, vmax=1.2,
                    extent=[t[0], t[-1], len(penv), 0])
    a2.axvline(0, color="k", ls=":", lw=0.8); a2.set_xlabel("time from SO trough (s)")
    a2.set_ylabel("events, sorted by SO trough DEPTH (deep→shallow)")
    a2.set_title(f"3. Sorted by independent var (SO depth)\nSpearman(depth,spindle) "
                 f"r={r_depth:.2f} p={p_depth:.1e}"); fig.colorbar(im2, ax=a2, label="spindle z")
    fig.suptitle("De-circularized N3 SO→spindle nesting (atlas, 204 Hz)", y=1.02, fontsize=12)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"nesting_validation.{ext}"), dpi=150, bbox_inches="tight")
    json_path = os.path.join(OUT, "validation_stats.json")
    pd.Series({"n_patients": len(pts), "split_half_r": r_split, "split_half_p": p_split,
               "depth_spindle_r": r_depth, "depth_spindle_p": p_depth,
               "grand_upstate_z": float(grand[win].mean()),
               "upstate_ci_lo": float(up_ci[0]), "upstate_ci_hi": float(up_ci[1])}).to_json(json_path)
    print(f"[validation] figure + stats in {OUT}")


if __name__ == "__main__":
    main()
