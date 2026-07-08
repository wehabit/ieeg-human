"""
Full-night SO->spindle/ripple coupling across MULTIPLE HUP subjects (iEEG.org, 1024 Hz).

Generalises hup165_night_coupling to every pulled full night (data/ieeg_portal/HUP*_night1/),
giving CROSS-SUBJECT inference on continuous full-night data — the n>1 version.

Per subject: polarity-robust cross-validated modulation per bipolar MTL channel (half A -> preferred
SO phase, half B tested), then the subject's value = mean over its channels. Cohort statistics are
over SUBJECTS. Reports per-channel positivity per subject too (pathology-robustness: coupling driven
by many channels, not a few epileptogenic ones).
"""
import os, glob, re, json
import numpy as np, pandas as pd
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from hup165_night_coupling import notch, blank, accumulate, aligned_curve, BANDS, CENTRE

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PORTAL = os.path.join(ROOT, "data", "ieeg_portal")
OUT = os.path.join(ROOT, "outputs", "hup_night_cohort")


def bipolar_generic(data, names):
    col = {n: i for i, n in enumerate(names)}
    out = {}
    for n in names:
        m = re.match(r"([A-Za-z]+)(\d+)$", n)
        if not m:
            continue
        sh, k = m.group(1), int(m.group(2)); nxt = f"{sh}{k+1}"
        if nxt in col:
            out[f"{sh}{k}-{k+1}"] = data[:, col[n]] - data[:, col[nxt]]
    return out


def subject_night(night_dir):
    idx = pd.read_csv(os.path.join(night_dir, "index.csv"))
    whole = {}
    for _, row in idx.iterrows():
        b = np.load(os.path.join(night_dir, row["file"]), allow_pickle=True)
        data = b["data"]; names = [str(n) for n in b["ch_names"]]; sf = float(b["sfreq"])
        for lab, sig in bipolar_generic(data, names).items():
            sig = notch(signal.detrend(sig), sf)
            whole.setdefault(lab, blank()); accumulate(whole[lab], sig, sf)
    # per-channel center modulation
    res = {}
    for band in BANDS:
        vals = [aligned_curve(whole[c][band]) for c in whole]
        vals = np.array([v[CENTRE] for v in vals if v is not None])
        res[band] = vals
    return res


def main():
    os.makedirs(OUT, exist_ok=True)
    dirs = sorted(glob.glob(os.path.join(PORTAL, "HUP*_night1")))
    subs = {}
    for d in dirs:
        if not os.path.exists(os.path.join(d, "index.csv")):
            continue
        sid = os.path.basename(d).split("_")[0]
        subs[sid] = subject_night(d)
        for band in BANDS:
            v = subs[sid][band]
            print(f"  {sid}/{band}: mean {v.mean():.3f}  positive {int((v>0).sum())}/{len(v)} ch", flush=True)
    ids = list(subs)
    print(f"[cohort] {len(ids)} full-night subjects: {ids}")

    fig, ax = plt.subplots(1, len(BANDS), figsize=(5 * len(BANDS), 5))
    summary = {}
    for j, band in enumerate(BANDS):
        subj_means = np.array([subs[s][band].mean() for s in ids])
        rng = np.random.default_rng(0)
        bs = subj_means[rng.integers(0, len(ids), (2000, len(ids)))].mean(1)
        ci = np.percentile(bs, [2.5, 97.5])
        summary[band] = {"n_subj": len(ids), "subj_means": {s: float(subs[s][band].mean()) for s in ids},
                         "cohort_mean": float(subj_means.mean()), "ci": [float(ci[0]), float(ci[1])],
                         "subj_pos_channels": {s: f"{int((subs[s][band]>0).sum())}/{len(subs[s][band])}" for s in ids}}
        # plot: per-subject dots + cohort mean±CI
        for i, s in enumerate(ids):
            ax[j].scatter([i] * len(subs[s][band]), subs[s][band], s=10, color="0.6", alpha=0.5)
            ax[j].plot(i, subj_means[i], "o", color="tab:blue", ms=9)
            ax[j].text(i, ax[j].get_ylim()[1], s, rotation=90, fontsize=6, va="top")
        ax[j].axhline(0, color="0.7", lw=0.7)
        ax[j].axhspan(ci[0], ci[1], color="tab:blue", alpha=0.12)
        ax[j].axhline(subj_means.mean(), color="tab:blue", lw=2)
        ok = "✓>0" if ci[0] > 0 else "n.s."
        ax[j].set_title(f"SO→{band} (full nights)\ncohort mean={subj_means.mean():.3f} "
                        f"[{ci[0]:.3f},{ci[1]:.3f}] {ok}  n={len(ids)} subj", fontsize=10)
        ax[j].set_xticks(range(len(ids))); ax[j].set_xticklabels(ids, rotation=45, fontsize=7)
        ax[j].set_ylabel("cross-val modulation (per channel = grey, subject = blue)")
    fig.suptitle("Full-night SO coupling across HUP subjects (iEEG.org, 1024 Hz, polarity-robust)",
                 y=1.0, fontsize=12)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"hup_night_cohort.{ext}"), dpi=150, bbox_inches="tight")
    json.dump(summary, open(os.path.join(OUT, "cohort_stats.json"), "w"), indent=2)
    print(f"[cohort] figure + stats in {OUT}")


if __name__ == "__main__":
    main()
