"""
M2 - Slow-oscillation -> spindle phase-amplitude coupling (tests H2, H3).

This is the human LFP analogue of the mouse spike-field PLV. In the mouse, spikes locked
to the theta-phase of the field wave. Here there are no spikes (204 Hz LFP), so we ask the
field-level version of the same question: does the phase of the slow oscillation (0.5-1.25 Hz)
organize the amplitude of sleep spindles (11-16 Hz)? That nesting IS the Kipnis 'coordinated
field wave', and the Kipnis prediction is that it is maximal in N3/N2 and minimal in wake/REM.

Metric: Tort modulation index (KL divergence of the phase-binned amplitude distribution from
uniform), with a circular-shift surrogate null -> MI_z. Also the preferred coupling phase (H3:
where on the slow cycle spindle power peaks). Aggregated channel -> patient -> ROI group -> state.

Usage:
  python analysis/m2_slow_spindle_pac.py --max-pt 60 --max-clips 8
"""
import argparse, os, json
import numpy as np, pandas as pd
from scipy import signal
from scipy.stats import wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import atlas

PHASE_BAND = (0.5, 1.25)    # slow oscillation
AMP_BAND = (11, 16)         # spindle
NBINS = 18
N_SURR = 120
OUT = os.path.join(atlas.ROOT, "outputs", "m2_slow_spindle_pac")
GROUPS = ["entorhinal", "parahippocampal", "hippocampal", "temporal_neocortex"]


def bandpass(x, sf, lo, hi):
    b, a = signal.butter(3, [lo / (sf / 2), hi / (sf / 2)], btype="band")
    return signal.filtfilt(b, a, x)


def tort_mi(phase, amp, nbins=NBINS):
    """Tort modulation index + mean amplitude per phase bin (both radians in [-pi,pi))."""
    edges = np.linspace(-np.pi, np.pi, nbins + 1)
    idx = np.digitize(phase, edges) - 1
    idx[idx == nbins] = nbins - 1
    m = np.array([amp[idx == b].mean() if np.any(idx == b) else 0.0 for b in range(nbins)])
    if m.sum() <= 0:
        return np.nan, m
    p = m / m.sum()
    p = np.clip(p, 1e-12, None)
    mi = (np.log(nbins) + (p * np.log(p)).sum()) / np.log(nbins)   # normalized KL
    return float(mi), m


def pac_for_channel(x, sf, rng):
    """Return (mi, mi_z, pref_phase) for one channel signal, or None."""
    if not np.isfinite(x).all() or np.std(x) < 1e-9:
        return None
    x = signal.detrend(x)
    ph = np.angle(signal.hilbert(bandpass(x, sf, *PHASE_BAND)))
    am = np.abs(signal.hilbert(bandpass(x, sf, *AMP_BAND)))
    mi, m = tort_mi(ph, am)
    if not np.isfinite(mi):
        return None
    # circular-shift surrogate null on the amplitude envelope
    n = len(am)
    surr = np.empty(N_SURR)
    shifts = rng.integers(int(0.1 * n), int(0.9 * n), N_SURR)
    for k, s in enumerate(shifts):
        surr[k], _ = tort_mi(ph, np.roll(am, s))
    mi_z = (mi - np.nanmean(surr)) / (np.nanstd(surr) + 1e-12)
    edges = np.linspace(-np.pi, np.pi, NBINS + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    pref = float(centers[int(np.argmax(m))])
    return mi, float(mi_z), pref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-pt", type=int, default=60)
    ap.add_argument("--max-clips", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(0)

    meta = atlas.load_metadata(normative_only=True)
    roi = meta.dropna(subset=["roi_group"])
    roi = roi[roi.roi_group.isin(GROUPS)][["pt", "name", "roi_group"]]
    clips = atlas.clip_index()
    meso_pt = roi[roi.roi_group.isin(atlas.MESIOTEMPORAL)].pt.unique()
    pts = sorted(set(clips.pt.unique()) & set(meso_pt))[:args.max_pt]
    print(f"[M2] {len(pts)} patients; phase {PHASE_BAND} Hz -> amp {AMP_BAND} Hz")

    rows = []
    for i, pt in enumerate(pts):
        ch_map = roi[roi.pt == pt].set_index("name")["roi_group"].to_dict()
        for st in atlas.STATES:
            idxs = sorted(clips[(clips.pt == pt) & (clips.state == st)].idx.tolist())[:args.max_clips]
            if not idxs:
                continue
            atlas.download([f"files/processed_eeg_final/{pt}_{st}_{j}.edf" for j in idxs])
            acc = {}
            for j in idxs:
                try:
                    data, names, sf = atlas.read_clip(pt, st, j)
                except Exception:
                    continue
                nidx = {n: k for k, n in enumerate(names)}
                for chan, grp in ch_map.items():
                    if chan not in nidx:
                        continue
                    r = pac_for_channel(data[nidx[chan]], sf, rng)
                    if r:
                        acc.setdefault((chan, grp), []).append(r)
            for (chan, grp), lst in acc.items():
                mi = np.mean([a[0] for a in lst]); miz = np.mean([a[1] for a in lst])
                pref = np.angle(np.mean([np.exp(1j * a[2]) for a in lst]))
                rows.append({"pt": pt, "channel": chan, "roi_group": grp, "state": st,
                             "n_clips": len(lst), "mi": float(mi), "mi_z": float(miz),
                             "pref_phase": float(pref)})
        print(f"  [{i+1}/{len(pts)}] {pt} done")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "channel_state_pac.csv"), index=False)
    print(f"[M2] wrote {len(df)} channel-state PAC rows")

    def boot_ci(vals, n=2000):
        vals = np.asarray(vals, float)
        if len(vals) < 2:
            return (np.nan, np.nan)
        idx = np.random.default_rng(0).integers(0, len(vals), (n, len(vals)))
        m = vals[idx].mean(1)
        return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))

    summ = []
    for (grp, st), g in df.groupby(["roi_group", "state"]):
        per_pt = g.groupby("pt")["mi"].mean()
        lo, hi = boot_ci(per_pt.values)
        summ.append({"roi_group": grp, "state": st, "n_pt": per_pt.size,
                     "mi_mean": float(per_pt.mean()), "ci_lo": lo, "ci_hi": hi,
                     "mi_z_mean": float(g.groupby("pt")["mi_z"].mean().mean())})
    sdf = pd.DataFrame(summ)
    sdf.to_csv(os.path.join(OUT, "group_state_pac.csv"), index=False)

    # H2 test: entorhinal+parahippocampal, N3 vs W, paired across patients
    ent = df[df.roi_group.isin(["entorhinal", "parahippocampal"])]
    piv = ent.groupby(["pt", "state"])["mi"].mean().unstack("state")
    print("\n[H2] entorhinal/parahippocampal slow->spindle MI, per patient:")
    print(piv.round(4).to_string())
    if {"W", "N3"} <= set(piv.columns):
        pair = piv[["W", "N3"]].dropna()
        if len(pair) >= 3:
            stat, p = wilcoxon(pair.W, pair.N3)
            print(f"\n[H2 raw MI] N={len(pair)} | W={pair.W.mean():.4f} N3={pair.N3.mean():.4f} "
                  f"| p={p:.4g} | N3>W {(pair.N3>pair.W).sum()}/{len(pair)}  "
                  f"(raw MI is wake-artifact-prone; MI_z below is primary)")
    # PRIMARY test: surrogate-corrected MI_z removes wake movement-artifact inflation
    pivz = ent.groupby(["pt", "state"])["mi_z"].mean().unstack("state")
    if {"W", "N3"} <= set(pivz.columns):
        pz = pivz[["W", "N3"]].dropna()
        if len(pz) >= 3:
            stat, p = wilcoxon(pz.W, pz.N3)
            print(f"[H2 MI_z]  N={len(pz)} | W={pz.W.mean():.4f} N3={pz.N3.mean():.4f} "
                  f"| Wilcoxon p={p:.4g} | N3>W {(pz.N3>pz.W).sum()}/{len(pz)}")

    # figure
    order = ["W", "N2", "N3", "R"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for grp in GROUPS:   # PRIMARY metric: surrogate-corrected MI_z (wake-artifact robust)
        xs, ys, los, his = [], [], [], []
        for st in order:
            pp = df[(df.roi_group == grp) & (df.state == st)].groupby("pt")["mi_z"].mean()
            if pp.size < 3:
                continue
            lo, hi = boot_ci(pp.values)
            xs.append(order.index(st)); ys.append(pp.mean()); los.append(lo); his.append(hi)
        if xs:
            ax[0].plot(xs, ys, "-o", label=grp)
            ax[0].fill_between(xs, los, his, alpha=0.15)
    ax[0].axhline(0, color="k", lw=0.7, ls=":")
    ax[0].set_xticks(range(len(order))); ax[0].set_xticklabels(order)
    ax[0].set_xlabel("state"); ax[0].set_ylabel("slow->spindle coupling (surrogate-corrected MI_z)")
    ax[0].set_title("M2: slow-oscillation -> spindle coupling by state\n(surrogate-corrected; Kipnis H2: peaks N3, absent wake/REM)")
    ax[0].legend(fontsize=8)
    # H3: preferred SO phase of spindle power, N3, mesiotemporal
    n3 = df[(df.state == "N3") & df.roi_group.isin(["entorhinal", "parahippocampal"])]
    if len(n3):
        ax[1] = plt.subplot(1, 2, 2, projection="polar")
        ax[1].hist(n3.pref_phase, bins=18, color="tab:orange", alpha=0.8)
        ax[1].set_title("M3/H3: preferred SO phase of\nspindle power (N3, mesiotemporal)")
    plt.tight_layout()
    for ext in ("png", "svg"):
        plt.savefig(os.path.join(OUT, f"m2_pac_by_state.{ext}"), dpi=150)
    json.dump({"phase_band": PHASE_BAND, "amp_band": AMP_BAND, "nbins": NBINS,
               "n_surr": N_SURR, "n_patients": len(pts)},
              open(os.path.join(OUT, "params.json"), "w"), indent=2)
    print(f"\n[M2] figure + csv in {OUT}")


if __name__ == "__main__":
    main()
