"""
Full hierarchical coupling (SO -> spindle -> ripple) on the ripple-capable cohorts.

The Staresina 2015 hierarchy is a 3-level cascade: SO phase gates spindles; spindle phase gates
ripples. Earlier work here tested only the two SO-> legs. This adds the missing SPINDLE->RIPPLE leg,
so we can legitimately speak of the hierarchy. Polarity-robust + cross-validated (same method as
nesting_phase_aligned): per channel, half A (interleaved 5 s blocks) gives the preferred phase;
half B is tested at it. Statistics over SUBJECTS.

Legs:
  SO->spindle      phase 0.5-1.25 Hz -> amp 11-16 Hz
  spindle->ripple  phase 11-16 Hz    -> amp 80-120 Hz   (the missing leg)
  SO->ripple       phase 0.5-1.25 Hz -> amp 80-120 Hz   (shortcut, for completeness)

The 204 Hz atlas cannot reach ripples, so the full hierarchy is Falach (1 kHz) + Zurich (2 kHz);
the atlas contributes only SO->spindle (reported elsewhere).
"""
import os
import numpy as np, pandas as pd
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nesting_validation_cohorts as nvc
from nesting_phase_aligned import bpf, NBINS, CENTRE

OUT = os.path.join(nvc.atlas.ROOT, "outputs", "hierarchical_coupling")
LEGS = [("SO→spindle", (0.5, 1.25), (11, 16)),
        ("spindle→ripple", (11, 16), (80, 120)),
        ("SO→ripple", (0.5, 1.25), (80, 120))]


def subj_aligned(signals, sf, pband, aband):
    edges = np.linspace(-np.pi, np.pi, NBINS + 1); L = int(5 * sf)
    sA = np.zeros(NBINS); cA = np.zeros(NBINS); sB = np.zeros(NBINS); cB = np.zeros(NBINS)
    for x in signals:
        if len(x) < 2 * L or not np.isfinite(x).all() or np.std(x) < 1e-9:
            continue
        phi = np.angle(signal.hilbert(bpf(x, *pband, sf)))
        amp = np.abs(signal.hilbert(bpf(x, *aband, sf)))
        amp = (amp - amp.mean()) / (amp.std() + 1e-12)
        idx = np.clip(np.digitize(phi, edges) - 1, 0, NBINS - 1)
        blk = (np.arange(len(x)) // L) % 2 == 0
        for b in range(NBINS):
            mA = (idx == b) & blk; mB = (idx == b) & ~blk
            sA[b] += amp[mA].sum(); cA[b] += mA.sum()
            sB[b] += amp[mB].sum(); cB[b] += mB.sum()
    if (cA < 5).any() or (cB < 5).any():
        return None
    a, bb = sA / cA, sB / cB
    curve = np.roll(bb, CENTRE - int(np.argmax(a)))
    return curve - curve.mean()


def cohort_leg(subs, pband, aband):
    rows = []
    for sid, (sigs, sf, _) in subs.items():
        c = subj_aligned(sigs, sf, pband, aband)
        if c is not None:
            rows.append(c)
    if len(rows) < 4:
        return None
    M = np.vstack(rows); rng = np.random.default_rng(0)
    bs = M[rng.integers(0, len(M), (2000, len(M)))].mean(1)
    ci = np.percentile(bs[:, CENTRE], [2.5, 97.5])
    return {"n": len(M), "mean": M.mean(0), "lo": np.percentile(bs, 2.5, 0),
            "hi": np.percentile(bs, 97.5, 0), "centre": float(M[:, CENTRE].mean()),
            "ci": [float(ci[0]), float(ci[1])]}


def main():
    raise SystemExit(
        "WITHDRAWN: this legacy hierarchical-coupling path has not passed the current artifact, "
        "staging, participant-unit, or provenance audit. Do not regenerate or cite its figures.")
    os.makedirs(OUT, exist_ok=True)
    cohorts = [("Falach 1 kHz", nvc.load_falach()), ("Zurich 2 kHz", nvc.load_zurich())]
    x = np.linspace(-np.pi, np.pi, NBINS)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8)); summary = {}
    for i, (cname, subs) in enumerate(cohorts):
        summary[cname] = {}
        for j, (leg, pb, ab) in enumerate(LEGS):
            R = cohort_leg(subs, pb, ab); ax = axes[i][j]
            if R is None:
                ax.text(0.5, 0.5, "insufficient", ha="center"); ax.set_xticks([]); ax.set_yticks([]); continue
            summary[cname][leg] = {"n": R["n"], "centre": R["centre"], "ci": R["ci"]}
            col = ["tab:blue", "tab:green", "tab:red"][j]
            ax.plot(x, R["mean"], color=col, lw=2, marker="o", ms=3)
            ax.fill_between(x, R["lo"], R["hi"], color=col, alpha=0.2)
            ax.axvline(0, color="k", ls=":", lw=.8); ax.axhline(0, color="0.7", lw=.6)
            ok = "✓>0" if R["ci"][0] > 0 else "n.s."
            ax.set_title(f"{cname}: {leg}\nmod={R['centre']:.3f} [{R['ci'][0]:.3f},{R['ci'][1]:.3f}] {ok} (n={R['n']})",
                         fontsize=8.5)
            ax.set_xlabel("phase vs preferred"); ax.set_ylabel("amp (z)")
            print(f"[{cname}/{leg}] n={R['n']} modulation={R['centre']:.3f} "
                  f"CI[{R['ci'][0]:.3f},{R['ci'][1]:.3f}]", flush=True)
    fig.suptitle("Hierarchical coupling SO→spindle→ripple (polarity-robust, cross-validated) — "
                 "the full Staresina 2015 cascade", y=1.0, fontsize=12)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"hierarchical_coupling.{ext}"), dpi=150, bbox_inches="tight")
    pd.Series(summary).to_json(os.path.join(OUT, "hierarchical_stats.json"))
    print(f"[hierarchical] figure + stats in {OUT}")


if __name__ == "__main__":
    main()
