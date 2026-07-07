"""
Slow-power - Spectral state map (the anchor result, tests H1).

Kipnis prediction: slow-wave (0.5-4 Hz) coordination in mesiotemporal cortex is
state-dependent -- it should rise wake -> N2 -> N3. This is the human analogue of the
mouse 'field amplitude by epoch' step, but now across NATURAL sleep states that the
mouse pilot never had.

For each normative channel we compute band power per 30 s clip (Welch), average clips
within a channel, then aggregate channel -> patient -> ROI group, with a patient-level
bootstrap CI. Primary contrast: relative slow-band power, wake vs N3, in entorhinal /
parahippocampal cortex (the human match to the mouse LEC arm).

Usage:
  python analysis/slow_power_by_state.py --max-pt 12 --max-clips 4
  python analysis/slow_power_by_state.py --all          # every patient/clip (~GBs)
"""
import argparse, os, json
import numpy as np, pandas as pd
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import atlas

BANDS = {"slow": (0.5, 4), "delta": (1, 4), "theta": (4, 8),
         "spindle": (11, 16), "beta": (16, 30), "gamma": (30, 80)}
DENOM = (0.5, 45)          # normalization band, below the 60 Hz line
OUT = os.path.join(atlas.ROOT, "outputs", "slow_power_by_state")


def clip_bandpowers(x, sf):
    """Return dict of relative band powers for one channel signal, or None if bad."""
    if not np.isfinite(x).all() or np.std(x) < 1e-9:
        return None
    x = signal.detrend(x)
    b, a = signal.iirnotch(atlas.LINE_HZ, 30, sf)   # 60 Hz line
    x = signal.filtfilt(b, a, x)
    f, p = signal.welch(x, sf, nperseg=int(4 * sf), noverlap=int(2 * sf))
    denom = np.trapezoid(p[(f >= DENOM[0]) & (f < DENOM[1])], f[(f >= DENOM[0]) & (f < DENOM[1])])
    if denom <= 0:
        return None
    out = {}
    for name, (lo, hi) in BANDS.items():
        m = (f >= lo) & (f < hi)
        out[name] = float(np.trapezoid(p[m], f[m]) / denom)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-pt", type=int, default=12)
    ap.add_argument("--max-clips", type=int, default=4)
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    meta = atlas.load_metadata(normative_only=True)
    roi = meta.dropna(subset=["roi_group"])[["pt", "name", "roi_group", "final_label", "hemi"]]
    clips = atlas.clip_index()

    # patients that have any mesiotemporal channel
    meso_pt = roi[roi.roi_group.isin(atlas.MESIOTEMPORAL)].pt.unique()
    pts = sorted(set(clips.pt.unique()) & set(meso_pt))
    if not args.all:
        pts = pts[:args.max_pt]
    print(f"[Slow-power] {len(pts)} patients; groups={sorted(roi.roi_group.unique())}")

    rows = []
    for i, pt in enumerate(pts):
        ch_map = roi[roi.pt == pt].set_index("name")["roi_group"].to_dict()
        for st in atlas.STATES:
            idxs = sorted(clips[(clips.pt == pt) & (clips.state == st)].idx.tolist())
            if not args.all:
                idxs = idxs[:args.max_clips]
            if not idxs:
                continue
            atlas.download([f"files/processed_eeg_final/{pt}_{st}_{j}.edf" for j in idxs])
            acc = {}   # channel -> list of band dicts
            for j in idxs:
                try:
                    data, names, sf = atlas.read_clip(pt, st, j)
                except Exception:
                    continue
                nidx = {n: k for k, n in enumerate(names)}
                for chan, grp in ch_map.items():
                    if chan not in nidx:
                        continue
                    bp = clip_bandpowers(data[nidx[chan]], sf)
                    if bp:
                        acc.setdefault((chan, grp), []).append(bp)
            for (chan, grp), lst in acc.items():
                mean = {b: float(np.mean([d[b] for d in lst])) for b in BANDS}
                rows.append({"pt": pt, "channel": chan, "roi_group": grp,
                             "state": st, "n_clips": len(lst), **mean})
        print(f"  [{i+1}/{len(pts)}] {pt} done")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "channel_state_bandpower.csv"), index=False)
    print(f"[Slow-power] wrote {len(df)} channel-state rows")

    # channel -> patient -> group means, patient-level bootstrap CI on 'slow'
    def boot_ci(vals, n=2000):
        vals = np.asarray(vals, float)
        if len(vals) < 2:
            return (np.nan, np.nan)
        idx = np.random.default_rng(0).integers(0, len(vals), (n, len(vals)))
        m = vals[idx].mean(1)
        return (float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5)))

    summ = []
    for (grp, st), g in df.groupby(["roi_group", "state"]):
        per_pt = g.groupby("pt")["slow"].mean()
        lo, hi = boot_ci(per_pt.values)
        summ.append({"roi_group": grp, "state": st, "n_pt": per_pt.size,
                     "slow_rel_mean": float(per_pt.mean()), "ci_lo": lo, "ci_hi": hi})
    sdf = pd.DataFrame(summ)
    sdf.to_csv(os.path.join(OUT, "group_state_slow.csv"), index=False)

    # H1 test: entorhinal+parahippocampal, wake vs N3, paired across patients
    from scipy.stats import wilcoxon
    ent = df[df.roi_group.isin(["entorhinal", "parahippocampal"])]
    piv = ent.groupby(["pt", "state"])["slow"].mean().unstack("state")
    print("\n[H1] entorhinal/parahippocampal relative slow (0.5-4 Hz), per patient:")
    print(piv.round(3).to_string())
    if {"W", "N3"} <= set(piv.columns):
        pair = piv[["W", "N3"]].dropna()
        if len(pair) >= 3:
            stat, p = wilcoxon(pair.W, pair.N3)
            print(f"\n[H1] N={len(pair)} patients | slow W={pair.W.mean():.3f} "
                  f"N3={pair.N3.mean():.3f} | Wilcoxon p={p:.4g} "
                  f"| N3>W in {(pair.N3>pair.W).sum()}/{len(pair)}")

    # figure: relative slow power across states, per ROI group
    order = ["W", "N2", "N3", "R"]
    groups = [g for g in ["entorhinal", "parahippocampal", "hippocampal", "amygdala",
                          "temporal_neocortex"] if g in sdf.roi_group.unique()]
    plt.figure(figsize=(8, 5))
    for grp in groups:
        s = sdf[(sdf.roi_group == grp) & (sdf.n_pt >= 3)].set_index("state").reindex(order).dropna(subset=["slow_rel_mean"])
        if s.empty:
            continue
        xs = [order.index(k) for k in s.index]
        plt.plot(xs, s.slow_rel_mean, "-o", label=grp)
        plt.fill_between(xs, s.ci_lo, s.ci_hi, alpha=0.15)
    plt.xticks(range(len(order)), order)
    plt.xlabel("sleep/wake state"); plt.ylabel("relative slow-band power (0.5-4 Hz)")
    plt.title("Slow-power: mesiotemporal slow-wave power by state\n(Kipnis H1: rises W->N2->N3)")
    plt.legend(fontsize=8); plt.tight_layout()
    for ext in ("png", "svg"):
        plt.savefig(os.path.join(OUT, f"slow_power_by_state.{ext}"), dpi=150)
    json.dump({"bands": BANDS, "denom": DENOM, "line_hz": atlas.LINE_HZ,
               "n_patients": len(pts)}, open(os.path.join(OUT, "params.json"), "w"), indent=2)
    print(f"\n[Slow-power] figure + csv in {OUT}")


if __name__ == "__main__":
    main()
