"""
Polarity-robust, cross-validated SO->spindle/ripple coupling across all three cohorts.

The fixed post-trough window (nesting_validation_cohorts) is confounded by bipolar SO polarity.
Here we remove that confound: align each channel to ITS OWN preferred SO phase, cross-validated.

Per channel, split samples into interleaved 5 s blocks -> half A / half B (quasi-independent):
  - half A: amplitude-by-SO-phase curve -> preferred phase bin theta* (where amp peaks)
  - half B: amplitude-by-SO-phase curve, circularly rotated so theta* -> centre bin
Average the rotated half-B curves across subjects. A centre peak (CI above the curve mean) means
real, phase-consistent coupling that is polarity-independent. theta* from A, measured on B => not
circular. Reported stat: cross-validated modulation depth = centre value of the aligned B curve
(z-scored amplitude relative to the channel's mean), bootstrap CI over subjects.

Reuses the cohort loaders from nesting_validation_cohorts.
"""
import os
import numpy as np, pandas as pd
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nesting_validation_cohorts as nvc

OUT = os.path.join(nvc.atlas.ROOT, "outputs", "nesting_phase_aligned")
SO_BAND = (0.5, 1.25); BANDS = {"spindle": (11, 16), "ripple": (80, 120)}
NBINS = 18; CENTRE = NBINS // 2


def bpf(x, lo, hi, sf):
    return signal.filtfilt(*signal.butter(3, [lo/(sf/2), hi/(sf/2)], btype="band"), x)


def subject_aligned(signals, sf, band):
    """Rotated half-B phase-amplitude curve (centre = half-A preferred phase), or None."""
    edges = np.linspace(-np.pi, np.pi, NBINS + 1)
    L = int(5 * sf)
    sumA = np.zeros(NBINS); cntA = np.zeros(NBINS); sumB = np.zeros(NBINS); cntB = np.zeros(NBINS)
    for x in signals:
        if len(x) < 2 * L or not np.isfinite(x).all() or np.std(x) < 1e-9:
            continue
        phi = np.angle(signal.hilbert(bpf(x, *SO_BAND, sf)))
        amp = np.abs(signal.hilbert(bpf(x, *band, sf)))
        amp = (amp - amp.mean()) / (amp.std() + 1e-12)
        idx = np.clip(np.digitize(phi, edges) - 1, 0, NBINS - 1)
        blk = (np.arange(len(x)) // L) % 2 == 0        # True = half A
        for b in range(NBINS):
            mA = (idx == b) & blk; mB = (idx == b) & ~blk
            sumA[b] += amp[mA].sum(); cntA[b] += mA.sum()
            sumB[b] += amp[mB].sum(); cntB[b] += mB.sum()
    if (cntA < 5).any() or (cntB < 5).any():
        return None
    cA = sumA / cntA; cB = sumB / cntB
    theta = int(np.argmax(cA))
    aligned = np.roll(cB, CENTRE - theta)
    return aligned - aligned.mean()                    # modulation relative to channel mean


def cohort(subs, name, bands):
    res = {}
    for b in bands:
        rows = []
        for sid, (sigs, sf, _) in subs.items():
            a = subject_aligned(sigs, sf, bands[b])
            if a is not None:
                rows.append(a)
        if len(rows) < 4:
            continue
        M = np.vstack(rows)
        rng = np.random.default_rng(0)
        bs = M[rng.integers(0, len(M), (2000, len(M)))].mean(1)
        centre_ci = np.percentile(bs[:, CENTRE], [2.5, 97.5])
        res[b] = {"n": len(M), "mean": M.mean(0),
                  "lo": np.percentile(bs, 2.5, 0), "hi": np.percentile(bs, 97.5, 0),
                  "centre": float(M[:, CENTRE].mean()),
                  "centre_ci": [float(centre_ci[0]), float(centre_ci[1])]}
        print(f"[{name}/{b}] n={len(M)} centre modulation={res[b]['centre']:.3f} "
              f"CI[{centre_ci[0]:.3f},{centre_ci[1]:.3f}] "
              f"{'>0' if centre_ci[0] > 0 else 'n.s.'}", flush=True)
    return res


def main():
    os.makedirs(OUT, exist_ok=True)
    cohorts = [("Atlas 204 Hz", nvc.load_atlas(), {"spindle": BANDS["spindle"]}),
               ("Falach 1 kHz", nvc.load_falach(), BANDS),
               ("Zurich 2 kHz", nvc.load_zurich(), BANDS)]
    x = np.linspace(-np.pi, np.pi, NBINS)            # phase offset from preferred (0 = preferred)
    fig, axes = plt.subplots(3, 2, figsize=(11, 11)); summary = {}
    for i, (cname, subs, bands) in enumerate(cohorts):
        print(f"== {cname}: {len(subs)} subjects", flush=True)
        res = cohort(subs, cname, bands)
        summary[cname] = {b: {k: res[b][k] for k in ["n", "centre", "centre_ci"]} for b in res}
        for j, b in enumerate(["spindle", "ripple"]):
            ax = axes[i][j]
            if b not in res:
                ax.text(0.5, 0.5, "not resolvable\n(204 Hz Nyquist)" if b == "ripple" and i == 0
                        else "insufficient data", ha="center", va="center", color="0.5")
                ax.set_xticks([]); ax.set_yticks([]); ax.set_title(f"{cname} — SO→{b}", fontsize=9)
                continue
            R = res[b]; col = "tab:orange" if b == "spindle" else "tab:red"
            ax.plot(x, R["mean"], color=col, lw=2, marker="o", ms=3)
            ax.fill_between(x, R["lo"], R["hi"], color=col, alpha=0.22)
            ax.axvline(0, color="k", ls=":", lw=0.8); ax.axhline(0, color="0.7", lw=0.6)
            ok = "✓ centre>0" if R["centre_ci"][0] > 0 else "n.s."
            ax.set_title(f"{cname} — SO→{b}   cross-val modulation={R['centre']:.3f} "
                         f"[{R['centre_ci'][0]:.3f},{R['centre_ci'][1]:.3f}] {ok}  (n={R['n']})", fontsize=8.5)
            ax.set_xlabel("SO phase relative to preferred (rad; 0 = held-out preferred)")
            ax.set_ylabel(f"{b} amp (z, vs mean)")
    fig.suptitle("Polarity-robust, cross-validated SO coupling (align each channel to its own "
                 "preferred phase on half A, test on half B)", y=1.0, fontsize=11)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"nesting_phase_aligned.{ext}"), dpi=150, bbox_inches="tight")
    pd.Series(summary).to_json(os.path.join(OUT, "phase_aligned_stats.json"))
    print(f"[phase-aligned] figure + stats in {OUT}")


if __name__ == "__main__":
    main()
