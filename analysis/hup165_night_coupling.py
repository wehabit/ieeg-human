"""
Polarity-robust, cross-validated SO->spindle/ripple coupling on ONE FULL NIGHT of continuous
NREM (HUP165_phaseII, iEEG.org, 1024 Hz), plus how the coupling tracks ACROSS the night.

This is the thing no clip cohort could do: ~4 h of NREM spanning a single ~9 h night, from MTL
depth electrodes (LA/LB/LC/LH -> 28 bipolar). Same non-circular method as nesting_phase_aligned:
per bipolar channel, half A (interleaved 5 s blocks) gives the preferred SO phase, half B is tested
at it. Statistics are over the 28 bipolar CHANNELS (one patient), and repeated in time bins across
the night.

Memory-safe: accumulates phase-amplitude bin sums per channel (and per time bin) block by block.
"""
import os, json, re
import numpy as np, pandas as pd
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nesting_phase_aligned import bpf, SO_BAND, NBINS, CENTRE

BANDS = {"spindle": (11, 16), "ripple": (80, 120)}
LINE = 60.0                      # HUP = US 60 Hz line
NBIN_TIME = 6
DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "data", "ieeg_portal", "HUP165_night1")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "outputs", "hup165_night_coupling")


def notch(x, sf):
    for f0 in (LINE, 2 * LINE):                     # 60 + 120 (120 sits in ripple band)
        if f0 < sf / 2:
            x = signal.filtfilt(*signal.iirnotch(f0, 30, sf), x)
    return x


def bipolar(data, names):
    col = {n: i for i, n in enumerate(names)}
    out = {}
    for sh in ("LA", "LB", "LC", "LH"):
        for k in range(1, 8):
            a, b = f"{sh}{k}", f"{sh}{k+1}"
            if a in col and b in col:
                out[f"{sh}{k}-{k+1}"] = data[:, col[a]] - data[:, col[b]]
    return out


def blank():
    return {b: {"sA": np.zeros(NBINS), "cA": np.zeros(NBINS),
                "sB": np.zeros(NBINS), "cB": np.zeros(NBINS)} for b in BANDS}


def accumulate(acc, sig, sf):
    edges = np.linspace(-np.pi, np.pi, NBINS + 1)
    L = int(5 * sf)
    if len(sig) < 2 * L or not np.isfinite(sig).all() or np.std(sig) < 1e-9:
        return
    phi = np.angle(signal.hilbert(bpf(sig, *SO_BAND, sf)))
    idx = np.clip(np.digitize(phi, edges) - 1, 0, NBINS - 1)
    half = (np.arange(len(sig)) // L) % 2 == 0
    for band, (lo, hi) in BANDS.items():
        amp = np.abs(signal.hilbert(bpf(sig, lo, hi, sf)))
        amp = (amp - amp.mean()) / (amp.std() + 1e-12)
        for bn in range(NBINS):
            mA = (idx == bn) & half; mB = (idx == bn) & ~half
            acc[band]["sA"][bn] += amp[mA].sum(); acc[band]["cA"][bn] += mA.sum()
            acc[band]["sB"][bn] += amp[mB].sum(); acc[band]["cB"][bn] += mB.sum()


def aligned_curve(a):
    if (a["cA"] < 5).any() or (a["cB"] < 5).any():
        return None
    cA = a["sA"] / a["cA"]; cB = a["sB"] / a["cB"]
    curve = np.roll(cB, CENTRE - int(np.argmax(cA)))
    return curve - curve.mean()


def main():
    os.makedirs(OUT, exist_ok=True)
    idx = pd.read_csv(os.path.join(DIR, "index.csv"))
    t0, t1 = idx.start_sec.min(), idx.start_sec.max() + idx.dur_sec.iloc[-1]
    edges = np.linspace(t0, t1, NBIN_TIME + 1)
    whole = {}; bytime = {}                 # channel -> acc ; (channel,tbin) -> acc
    for _, row in idx.iterrows():
        b = np.load(os.path.join(DIR, row["file"]), allow_pickle=True)
        data = b["data"]; names = [str(n) for n in b["ch_names"]]; sf = float(b["sfreq"])
        tb = min(int(np.searchsorted(edges, row["start_sec"], "right") - 1), NBIN_TIME - 1)
        for lab, sig in bipolar(data, names).items():
            sig = notch(signal.detrend(sig), sf)
            whole.setdefault(lab, blank()); accumulate(whole[lab], sig, sf)
            bytime.setdefault((lab, tb), blank()); accumulate(bytime[(lab, tb)], sig, sf)
    print(f"[hup165] {len(idx)} blocks; {len(whole)} bipolar channels", flush=True)

    # ---- whole-night coupling (over channels) ----
    x = np.linspace(-np.pi, np.pi, NBINS); rng = np.random.default_rng(0)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    stats = {}
    for j, band in enumerate(BANDS):
        curves = [aligned_curve(whole[c][band]) for c in whole]
        curves = np.array([c for c in curves if c is not None])
        M = curves; bs = M[rng.integers(0, len(M), (2000, len(M)))].mean(1)
        cci = np.percentile(bs[:, CENTRE], [2.5, 97.5])
        stats[band] = {"n_ch": len(M), "centre": float(M[:, CENTRE].mean()),
                       "ci": [float(cci[0]), float(cci[1])]}
        ax[j].plot(x, M.mean(0), color="tab:orange" if band == "spindle" else "tab:red", lw=2, marker="o", ms=3)
        ax[j].fill_between(x, np.percentile(bs, 2.5, 0), np.percentile(bs, 97.5, 0),
                           color="tab:orange" if band == "spindle" else "tab:red", alpha=0.22)
        ax[j].axvline(0, color="k", ls=":", lw=.8); ax[j].axhline(0, color="0.7", lw=.6)
        ok = "✓>0" if cci[0] > 0 else "n.s."
        ax[j].set_title(f"Whole night — SO→{band}\ncross-val modulation={M[:,CENTRE].mean():.3f} "
                        f"[{cci[0]:.3f},{cci[1]:.3f}] {ok} (n={len(M)} ch)", fontsize=9)
        ax[j].set_xlabel("SO phase vs preferred (rad)"); ax[j].set_ylabel(f"{band} amp (z)")

    # ---- across-night time course ----
    mids = (edges[:-1] + edges[1:]) / 2 / 3600.0
    for band, col in [("spindle", "tab:orange"), ("ripple", "tab:red")]:
        ys, es = [], []
        for tb in range(NBIN_TIME):
            vals = [aligned_curve(bytime[(c, tb)][band])[CENTRE]
                    for c in whole if (c, tb) in bytime and aligned_curve(bytime[(c, tb)][band]) is not None]
            ys.append(np.mean(vals) if vals else np.nan)
            es.append(np.std(vals)/np.sqrt(len(vals)) if vals else np.nan)
        ax[2].errorbar(mids, ys, yerr=es, marker="o", color=col, label=f"SO→{band}", capsize=3)
    ax[2].axhline(0, color="0.7", lw=.6)
    ax[2].set_title("Coupling across the night"); ax[2].set_xlabel("hours into recording")
    ax[2].set_ylabel("cross-val modulation (centre)"); ax[2].legend(fontsize=8)
    fig.suptitle("HUP165 — one full NREM night (iEEG.org, 1024 Hz), polarity-robust cross-validated "
                 "SO coupling", y=1.02, fontsize=12)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"hup165_night_coupling.{ext}"), dpi=150, bbox_inches="tight")
    json.dump(stats, open(os.path.join(OUT, "night_stats.json"), "w"), indent=2)
    for b in BANDS:
        print(f"[hup165/{b}] n={stats[b]['n_ch']} ch  modulation={stats[b]['centre']:.3f} "
              f"CI[{stats[b]['ci'][0]:.3f},{stats[b]['ci'][1]:.3f}]", flush=True)
    print(f"[hup165] figure + stats in {OUT}")


if __name__ == "__main__":
    main()
