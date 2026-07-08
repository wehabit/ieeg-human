"""
Headline "strength" figure: cross-validated SO-coupling modulation across all datasets, showing the
HUP165 full continuous night is the strongest — with its per-channel consistency overlaid.

Uses the committed stats:
  nesting_phase_aligned/phase_aligned_stats.json  (atlas/Falach/Zurich, cross-val modulation + CI)
  hup165_night_coupling/night_stats.json          (HUP165 full night)
  hup165_night_coupling/per_channel_modulation.csv (HUP165 per-channel dots)
"""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs", "fullnight_strength")


def main():
    os.makedirs(OUT, exist_ok=True)
    pa = json.load(open(os.path.join(ROOT, "outputs", "nesting_phase_aligned", "phase_aligned_stats.json")))
    night = json.load(open(os.path.join(ROOT, "outputs", "hup165_night_coupling", "night_stats.json")))
    perch = pd.read_csv(os.path.join(ROOT, "outputs", "hup165_night_coupling", "per_channel_modulation.csv"))

    # assemble per band: list of (label, modulation, ci_lo, ci_hi)
    def row(d, band):
        return (d[band]["centre"], d[band].get("ci", d[band].get("centre_ci"))[0],
                d[band].get("ci", d[band].get("centre_ci"))[1])
    bands = {
        "SO→spindle": [("Atlas\n204 Hz", *row(pa["Atlas 204 Hz"], "spindle")),
                       ("Falach\n1 kHz", *row(pa["Falach 1 kHz"], "spindle")),
                       ("Zurich\n2 kHz", *row(pa["Zurich 2 kHz"], "spindle")),
                       ("HUP165\nFULL NIGHT", *row(night, "spindle"))],
        "SO→ripple": [("Falach\n1 kHz", *row(pa["Falach 1 kHz"], "ripple")),
                      ("Zurich\n2 kHz", *row(pa["Zurich 2 kHz"], "ripple")),
                      ("HUP165\nFULL NIGHT", *row(night, "ripple"))],
    }
    dotcol = {"SO→spindle": "spindle", "SO→ripple": "ripple"}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    for ax, (band, rows) in zip(axes, bands.items()):
        labels = [r[0] for r in rows]; mods = [r[1] for r in rows]
        los = [r[1]-r[2] for r in rows]; his = [r[3]-r[1] for r in rows]
        xs = np.arange(len(rows))
        colors = ["0.7"]*(len(rows)-1) + ["tab:blue"]
        ax.bar(xs, mods, color=colors, yerr=[los, his], capsize=5, zorder=2)
        # overlay HUP165 per-channel dots on the last bar
        v = perch[dotcol[band]].dropna().values
        jit = (np.random.default_rng(0).random(len(v))-.5)*0.4
        ax.scatter(np.full(len(v), len(rows)-1)+jit, v, s=14, color="tab:blue",
                   edgecolor="white", linewidth=0.3, alpha=0.7, zorder=3,
                   label=f"HUP165 channels ({int((v>0).sum())}/{len(v)} > 0)")
        ax.axhline(0, color="k", lw=0.7)
        ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel("cross-validated coupling modulation")
        ax.set_title(f"{band}: strongest in the full continuous night", fontsize=11)
        ax.legend(fontsize=8, loc="upper left")
    fig.suptitle("SO→spindle / SO→ripple coupling is strongest in one full continuous NREM night "
                 "(HUP165, iEEG.org 1 kHz) — bars = cohort mean ± 95% CI; dots = HUP165 per channel",
                 y=1.0, fontsize=11)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"fullnight_strength.{ext}"), dpi=150, bbox_inches="tight")
    # console recap
    for band, rows in bands.items():
        print(f"{band}: " + " | ".join(f"{l.split(chr(10))[0]} {m:.3f}[{lo:.3f},{hi:.3f}]" for l, m, lo, hi in rows))
    print(f"[strength] figure in {OUT}")


if __name__ == "__main__":
    main()
