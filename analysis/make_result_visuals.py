"""
Result visuals - make the M1/M2 findings visible from the raw clips + summary CSVs.

Produces (outputs/result_visuals/):
  1. example_traces.png    - real iEEG, wake vs N3, one mesiotemporal channel (slow waves appear)
  2. psd_by_state.png       - that channel's power spectrum by state (slow-band bump grows in N3)
  3. m1_paired_slope.png    - population: slow-band power W->N3, one line per patient
  4. comodulogram.png       - that channel: phase-freq x amp-freq coupling, wake vs N3
  5. tort_phase_amp.png     - that channel: mean spindle amp across the slow-osc phase, wake vs N3
  6. m2_paired_slope.png    - population: coupling MI_z W->N3, one line per patient
"""
import os
import numpy as np, pandas as pd
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import atlas

OUT = os.path.join(atlas.ROOT, "outputs", "result_visuals")
os.makedirs(OUT, exist_ok=True)
M1 = os.path.join(atlas.ROOT, "outputs", "m1_spectral_state_map", "channel_state_bandpower.csv")
M2 = os.path.join(atlas.ROOT, "outputs", "m2_slow_spindle_pac", "channel_state_pac.csv")


def bp(x, sf, lo, hi):
    b, a = signal.butter(3, [lo / (sf / 2), hi / (sf / 2)], btype="band")
    return signal.filtfilt(b, a, x)


def load_channel_state(pt, chan, state, max_clips=8):
    """Concatenated signal (list of per-clip arrays) for a channel/state from cache."""
    clips = atlas.clip_index()
    idxs = sorted(clips[(clips.pt == pt) & (clips.state == state)].idx.tolist())[:max_clips]
    segs = []
    for j in idxs:
        try:
            data, names, sf = atlas.read_clip(pt, state, j)
        except Exception:
            continue
        if chan in names:
            segs.append(signal.detrend(data[names.index(chan)]))
    return segs, atlas.SFREQ


def pick_representative():
    """Channel with high N3 coupling that also has wake data cached."""
    m2 = pd.read_csv(M2)
    cand = m2[m2.roi_group.isin(["hippocampal", "parahippocampal", "entorhinal"])]
    piv = cand.pivot_table(index=["pt", "channel", "roi_group"], columns="state", values="mi_z")
    piv = piv.dropna(subset=["W", "N3"]) if {"W", "N3"} <= set(piv.columns) else piv
    piv = piv.sort_values("N3", ascending=False)
    for (pt, chan, grp), _ in piv.iterrows():
        segs, _ = load_channel_state(pt, chan, "N3", 1)
        if segs:
            return pt, chan, grp
    return None


def tort_dist(x, sf, nbins=18):
    ph = np.angle(signal.hilbert(bp(x, sf, 0.5, 1.25)))
    am = np.abs(signal.hilbert(bp(x, sf, 11, 16)))
    edges = np.linspace(-np.pi, np.pi, nbins + 1)
    idx = np.clip(np.digitize(ph, edges) - 1, 0, nbins - 1)
    m = np.array([am[idx == b].mean() if np.any(idx == b) else 0 for b in range(nbins)])
    return (edges[:-1] + edges[1:]) / 2, m


def comodulogram(segs, sf):
    pfreqs = np.arange(0.5, 4.01, 0.5)
    afreqs = np.arange(8, 30.1, 2.0)
    grid = np.zeros((len(afreqs), len(pfreqs)))
    for i, pf in enumerate(pfreqs):
        phs = [np.angle(signal.hilbert(bp(s, sf, max(0.25, pf - 0.5), pf + 0.5))) for s in segs]
        for k, af in enumerate(afreqs):
            ams = [np.abs(signal.hilbert(bp(s, sf, af - 2, af + 2))) for s in segs]
            ph = np.concatenate(phs); am = np.concatenate(ams)
            edges = np.linspace(-np.pi, np.pi, 19)
            idx = np.clip(np.digitize(ph, edges) - 1, 0, 17)
            mv = np.array([am[idx == b].mean() if np.any(idx == b) else 0 for b in range(18)])
            p = mv / mv.sum() if mv.sum() > 0 else np.ones(18) / 18
            p = np.clip(p, 1e-12, None)
            grid[k, i] = (np.log(18) + (p * np.log(p)).sum()) / np.log(18)
    return pfreqs, afreqs, grid


def main():
    pt, chan, grp = pick_representative()
    print(f"[visuals] representative: {pt} {chan} ({grp})")
    segs = {st: load_channel_state(pt, chan, st)[0] for st in ["W", "N2", "N3"]}
    sf = atlas.SFREQ

    # 1. example traces W vs N3 (first 10 s)
    n = int(10 * sf)
    fig, ax = plt.subplots(2, 1, figsize=(11, 5), sharex=True, sharey=True)
    for a, st, c in [(ax[0], "W", "tab:gray"), (ax[1], "N3", "tab:blue")]:
        if segs[st]:
            t = np.arange(n) / sf
            a.plot(t, segs[st][0][:n] * 1e3, c, lw=0.7)
        a.set_ylabel(f"{st}\n(mV)"); a.grid(alpha=0.3)
    ax[1].set_xlabel("time (s)")
    fig.suptitle(f"Real iEEG: wake vs N3  —  {chan} ({grp}), {pt}\nslow waves emerge in N3")
    plt.tight_layout(); plt.savefig(os.path.join(OUT, "example_traces.png"), dpi=150); plt.close()

    # 2. PSD by state
    plt.figure(figsize=(7, 5))
    for st, c in [("W", "tab:gray"), ("N2", "tab:orange"), ("N3", "tab:blue")]:
        ps = []
        for s in segs[st]:
            f, p = signal.welch(s, sf, nperseg=int(4 * sf))
            ps.append(p)
        if ps:
            plt.semilogy(f, np.mean(ps, 0), c, label=st, lw=1.5)
    plt.axvspan(0.5, 4, color="tab:blue", alpha=0.08, label="slow band")
    plt.xlim(0, 40); plt.xlabel("frequency (Hz)"); plt.ylabel("power (V^2/Hz)")
    plt.title(f"Power spectrum by state — {chan} ({grp})\nslow-band power grows into N3")
    plt.legend(); plt.tight_layout(); plt.savefig(os.path.join(OUT, "psd_by_state.png"), dpi=150); plt.close()

    # 3. M1 population paired slopegraph
    m1 = pd.read_csv(M1)
    ent = m1[m1.roi_group.isin(["entorhinal", "parahippocampal"])]
    piv = ent.groupby(["pt", "state"])["slow"].mean().unstack("state")[["W", "N3"]].dropna()
    plt.figure(figsize=(5, 6))
    for _, r in piv.iterrows():
        plt.plot([0, 1], [r.W, r.N3], "-o", color="tab:blue", alpha=0.4, mfc="white")
    plt.plot([0, 1], [piv.W.mean(), piv.N3.mean()], "-o", color="black", lw=3, label="mean")
    plt.xticks([0, 1], ["Wake", "N3"]); plt.ylabel("relative slow-band power (0.5-4 Hz)")
    plt.title(f"M1: slow power rises W->N3\n{ (piv.N3>piv.W).sum() }/{len(piv)} patients, p=2.9e-6")
    plt.legend(); plt.tight_layout(); plt.savefig(os.path.join(OUT, "m1_paired_slope.png"), dpi=150); plt.close()

    # 4. comodulogram W vs N3
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    vmax = 0
    grids = {}
    for st in ["W", "N3"]:
        if segs[st]:
            pf, af, g = comodulogram(segs[st], sf); grids[st] = (pf, af, g); vmax = max(vmax, g.max())
    for a, st in zip(ax, ["W", "N3"]):
        if st in grids:
            pf, af, g = grids[st]
            im = a.pcolormesh(pf, af, g, shading="auto", cmap="viridis", vmin=0, vmax=vmax)
            a.set_title(f"{st}"); a.set_xlabel("phase freq (Hz)"); a.set_ylabel("amp freq (Hz)")
            fig.colorbar(im, ax=a, label="MI")
    fig.suptitle(f"Comodulogram — {chan} ({grp}): coupling appears in N3 (SO-phase x spindle-amp)")
    plt.tight_layout(); plt.savefig(os.path.join(OUT, "comodulogram.png"), dpi=150); plt.close()

    # 5. Tort mean-amp-by-phase W vs N3
    plt.figure(figsize=(7, 5))
    for st, c in [("W", "tab:gray"), ("N3", "tab:blue")]:
        if segs[st]:
            x = np.concatenate(segs[st])
            ctr, m = tort_dist(x, sf)
            m = m / m.sum()
            plt.bar(ctr + (0.06 if st == "N3" else -0.06), m, width=0.12, color=c, alpha=0.8, label=st)
    plt.xlabel("slow-oscillation phase (rad)"); plt.ylabel("normalized spindle amplitude")
    plt.title(f"Spindle amplitude across the slow-oscillation cycle — {chan}\nflat=no coupling (Wake), peaked=coupling (N3)")
    plt.legend(); plt.tight_layout(); plt.savefig(os.path.join(OUT, "tort_phase_amp.png"), dpi=150); plt.close()

    # 6. M2 population paired slopegraph (MI_z)
    m2 = pd.read_csv(M2)
    e2 = m2[m2.roi_group.isin(["entorhinal", "parahippocampal"])]
    p2 = e2.groupby(["pt", "state"])["mi_z"].mean().unstack("state")[["W", "N3"]].dropna()
    plt.figure(figsize=(5, 6))
    for _, r in p2.iterrows():
        plt.plot([0, 1], [r.W, r.N3], "-o", color="tab:green", alpha=0.4, mfc="white")
    plt.plot([0, 1], [p2.W.mean(), p2.N3.mean()], "-o", color="black", lw=3, label="mean")
    plt.axhline(0, color="k", ls=":", lw=0.8)
    plt.xticks([0, 1], ["Wake", "N3"]); plt.ylabel("slow->spindle coupling (MI_z)")
    plt.title(f"M2: coupling rises W->N3\n{(p2.N3>p2.W).sum()}/{len(p2)} patients, p=0.030")
    plt.legend(); plt.tight_layout(); plt.savefig(os.path.join(OUT, "m2_paired_slope.png"), dpi=150); plt.close()

    print(f"[visuals] wrote 6 figures to {OUT}")


if __name__ == "__main__":
    main()
