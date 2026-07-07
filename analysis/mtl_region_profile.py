"""
Region-profile - Regional gradient (human analogue of the mouse `spike_50hz_interpretation` figure:
dHPC driven-up vs LEC net-suppressed = region-specific processing).

Here we ask the human version: does the Kipnis NREM coordination differ across
hippocampal vs entorhinal/parahippocampal vs temporal-neocortex sites? Built from the
already-computed Slow-power (slow-band power) and Spindle-coupling (slow->spindle MI_z) tables — no re-download.

Outputs (outputs/mtl_region_profile/):
  mtl_region_profile.png  - N3 slow power and coupling by region, per-patient + bootstrap CI
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import atlas

OUT = os.path.join(atlas.ROOT, "outputs", "mtl_region_profile")
SLOW_CSV = os.path.join(atlas.ROOT, "outputs", "slow_power_by_state", "channel_state_bandpower.csv")
PAC_CSV = os.path.join(atlas.ROOT, "outputs", "slow_spindle_coupling", "channel_state_pac.csv")
REGIONS = ["hippocampal", "parahippocampal", "entorhinal", "temporal_neocortex"]
LABELS = ["hippo", "parahip", "entorhinal", "neocortex"]


def boot_ci(v, n=2000):
    v = np.asarray(v, float)
    if len(v) < 2:
        return np.nan, np.nan
    idx = np.random.default_rng(0).integers(0, len(v), (n, len(v)))
    m = v[idx].mean(1)
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def region_stats(df, value, state="N3"):
    """Per-patient mean per region for one state; returns dict region-> array."""
    d = df[df.state == state]
    return {r: d[d.roi_group == r].groupby("pt")[value].mean().values for r in REGIONS}


def panel(ax, data, title, ylabel):
    for i, r in enumerate(REGIONS):
        v = data[r]
        if len(v) == 0:
            continue
        lo, hi = boot_ci(v)
        ax.bar(i, np.mean(v), color=f"C{i}", alpha=0.7)
        ax.errorbar(i, np.mean(v), yerr=[[np.mean(v) - lo], [hi - np.mean(v)]],
                    color="k", capsize=4)
        jitter = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.3
        ax.scatter(i + jitter, v, s=12, color="k", alpha=0.35, zorder=3)
        ax.text(i, ax.get_ylim()[1] * 0.02, f"n={len(v)}", ha="center", fontsize=8)
    ax.set_xticks(range(len(REGIONS))); ax.set_xticklabels(LABELS, rotation=15)
    ax.set_title(title); ax.set_ylabel(ylabel)


def main():
    os.makedirs(OUT, exist_ok=True)
    m1 = pd.read_csv(SLOW_CSV); m2 = pd.read_csv(PAC_CSV)
    slow = region_stats(m1, "slow", "N3")
    miz = region_stats(m2, "mi_z", "N3")

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    panel(ax[0], slow, "N3 slow-band power (0.5-4 Hz)\nby region", "relative slow power")
    ax[1].axhline(0, color="k", lw=0.7, ls=":")
    panel(ax[1], miz, "N3 slow->spindle coupling\nby region", "coupling MI_z")
    fig.suptitle("Region-profile (DESCRIPTIVE — pairwise MI_z differences are n.s.): N3 slow-power & "
                 "coupling by MTL region")
    plt.tight_layout()
    for ext in ("png", "svg"):
        plt.savefig(os.path.join(OUT, f"mtl_region_profile.{ext}"), dpi=150)

    # simple pairwise report (Mann-Whitney, unpaired across patients)
    from scipy.stats import mannwhitneyu
    rows = []
    for val, name in [(slow, "slow_power"), (miz, "mi_z")]:
        for a in REGIONS:
            for b in REGIONS:
                if a < b and len(val[a]) >= 3 and len(val[b]) >= 3:
                    try:
                        u, p = mannwhitneyu(val[a], val[b])
                    except ValueError:
                        continue
                    rows.append({"metric": name, "A": a, "B": b,
                                 "meanA": float(np.mean(val[a])), "meanB": float(np.mean(val[b])),
                                 "p": float(p)})
    rep = pd.DataFrame(rows)
    rep.to_csv(os.path.join(OUT, "region_pairwise.csv"), index=False)
    print("[Region-profile] N3 region means:")
    for r in REGIONS:
        print(f"  {r:20s} slow={np.mean(slow[r]) if len(slow[r]) else float('nan'):.3f} "
              f"(n={len(slow[r])})  MI_z={np.mean(miz[r]) if len(miz[r]) else float('nan'):.3f} "
              f"(n={len(miz[r])})")
    print(f"[Region-profile] figure + pairwise csv in {OUT}")


if __name__ == "__main__":
    main()
