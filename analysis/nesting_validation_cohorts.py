"""
De-circularized SO-nesting validation across ALL THREE cohorts, stacked, for spindle (all 3)
and ripple (the two >=1 kHz cohorts; the 204 Hz atlas cannot reach 80-120 Hz).

Per cohort, per band, two non-circular tests (no sorting):
  - MEAN +/- bootstrap CI over subjects  -> is the up-state (0-0.75 s post SO trough) envelope > 0?
  - SPLIT-HALF Spearman(odd,even up-state) across subjects -> stable per-subject trait?

We show cohorts SEPARATELY (replication), not pooled: different rates/montages/populations make a
single pooled statistic hard to defend, and separate rows reveal whether the effect actually
replicates. Zurich is an expert-HFO cohort with no SOZ/IED mask here -> its ripple band includes
some PATHOLOGICAL HFOs (caveat); the envelope measure is far less sensitive to this than the
event-based hfo_slow_phase analysis.
"""
import os, glob, re
import numpy as np, pandas as pd
from scipy import signal
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mne
import atlas
import slow_ripple_coupling as falach
import hfo_slow_phase as zur

OUT = os.path.join(atlas.ROOT, "outputs", "nesting_cohort_validation")
SO_BAND = (0.5, 1.25); BANDS = {"spindle": (11, 16), "ripple": (80, 120)}


def bpf(x, lo, hi, sf):
    return signal.filtfilt(*signal.butter(3, [lo/(sf/2), hi/(sf/2)], btype="band"), x)


def subject_bands(signals, sf, bands, mask_fn=None):
    """Return per-band list of SO-trough-locked z-scored envelopes (+/-1.5 s)."""
    HALF = int(1.5 * sf); MIN = int(0.6 * sf)
    out = {b: [] for b in bands}
    for x in signals:
        if len(x) < 2 * HALF + 10 or not np.isfinite(x).all() or np.std(x) < 1e-9:
            continue
        so = bpf(x, *SO_BAND, sf)
        idx, _ = signal.find_peaks(-so, height=so.std(), distance=MIN)
        idx = idx[(idx > HALF) & (idx < len(x) - HALF)]
        if mask_fn is not None:
            idx = idx[[mask_fn(e) for e in idx]]
        if len(idx) == 0:
            continue
        for b in bands:
            e_env = np.abs(signal.hilbert(bpf(x, *bands[b], sf)))
            e_env = (e_env - e_env.mean()) / (e_env.std() + 1e-12)
            for e in idx:
                out[b].append(e_env[e - HALF:e + HALF])
    return out, HALF


def analyse(subs, cohort, bands):
    """subs: sid -> (signals, sf, mask_fn). Return per-band stats + time axis (by sf)."""
    per = {}; sf0 = None; HALF0 = None
    for sid, (sigs, sf, mask_fn) in subs.items():
        sf0 = sf
        ev, HALF0 = subject_bands(sigs, sf, bands, mask_fn)
        t = np.arange(-HALF0, HALF0) / sf; win = (t >= 0.0) & (t <= 0.75)
        rec = {}
        for b in bands:
            E = np.array(ev[b])
            if len(E) < 20:
                continue
            odd = np.arange(len(E)) % 2 == 1
            rec[b] = {"mean": E.mean(0), "upA": E[odd][:, win].mean(),
                      "upB": E[~odd][:, win].mean(), "n": len(E)}
        if rec:
            per[sid] = rec
    t = np.arange(-HALF0, HALF0) / sf0; win = (t >= 0.0) & (t <= 0.75)
    rng = np.random.default_rng(0); res = {}
    for b in bands:
        ids = [s for s in per if b in per[s]]
        if len(ids) < 4:
            continue
        M = np.vstack([per[s][b]["mean"] for s in ids])
        bs = M[rng.integers(0, len(ids), (2000, len(ids)))].mean(1)
        upci = np.percentile(bs[:, win].mean(1), [2.5, 97.5])
        upA = np.array([per[s][b]["upA"] for s in ids]); upB = np.array([per[s][b]["upB"] for s in ids])
        r, p = spearmanr(upA, upB)
        res[b] = {"n_subj": len(ids), "t": t, "grand": M.mean(0),
                  "lo": np.percentile(bs, 2.5, 0), "hi": np.percentile(bs, 97.5, 0),
                  "up_z": float(M[:, win].mean()), "up_ci": [float(upci[0]), float(upci[1])],
                  "split_r": float(r), "split_p": float(p)}
        print(f"[{cohort}/{b}] n={len(ids)} up z={res[b]['up_z']:.3f} "
              f"CI[{upci[0]:.3f},{upci[1]:.3f}] split-half r={r:.3f} p={p:.3g}", flush=True)
    return res


# ---------------- cohort loaders ----------------
def load_atlas():
    meta = atlas.load_metadata(True)
    roi = meta.dropna(subset=["roi_group"]); roi = roi[roi.roi_group.isin(atlas.MESIOTEMPORAL)]
    clips = atlas.clip_index(); nb, na = signal.iirnotch(atlas.LINE_HZ, 30, atlas.SFREQ)
    subs = {}
    for pt in sorted(roi.pt.unique()):
        chans = list(roi[roi.pt == pt]["name"])
        cached = [j for j in sorted(clips[(clips.pt == pt) & (clips.state == "N3")].idx)
                  if os.path.exists(os.path.join(atlas.EEGDIR, f"{pt}_N3_{j}.edf"))]
        sigs = []
        for j in cached:
            try:
                data, names, _ = atlas.read_clip(pt, "N3", j)
            except Exception:
                continue
            nidx = {nm: k for k, nm in enumerate(names)}
            for ch in chans:
                if ch in nidx and np.isfinite(data[nidx[ch]]).all() and np.std(data[nidx[ch]]) > 1e-9:
                    sigs.append(signal.filtfilt(nb, na, signal.detrend(data[nidx[ch]])))
        if sigs:
            subs[pt] = (sigs, atlas.SFREQ, None)
    return subs


def load_falach():
    subs = {}
    for sd in sorted(glob.glob(os.path.join(falach.BIDS, "sub-*"))):
        try:
            bip, sf, ied, reg = falach.load_subject(sd)
        except Exception:
            continue
        sigs = list(bip.values())
        if not sigs:
            continue
        ied = np.array(ied, int)
        def mk(ied=ied, sf=sf):
            def f(e):
                return not np.any(np.abs(ied - e) < 0.5 * sf)
            return f
        subs[os.path.basename(sd)] = (sigs, sf, mk())
    return subs


def load_zurich():
    subs = {}
    for vhdr in sorted(glob.glob(os.path.join(zur.EEG, "sub-*_ieeg.vhdr"))):
        r = mne.io.read_raw_brainvision(vhdr, preload=True, verbose="ERROR")
        sf = r.info["sfreq"]; data = r.get_data(); names = r.ch_names
        nidx = {n: i for i, n in enumerate(names)}
        sigs = []
        for n in names:
            m = re.match(r"([A-Za-z]+)(\d+)$", n)
            if not m:
                continue
            pre, num = m.group(1), int(m.group(2)); nxt = f"{pre}{num+1}"
            if nxt in nidx and zur.region_of(f"{pre}{num}-{num+1}"):
                sig = zur.notch(signal.detrend(data[nidx[n]] - data[nidx[nxt]]), sf, zur.LINE)
                sigs.append(signal.decimate(sig, 2, ftype="fir"))  # 2 kHz -> 1 kHz
        if sigs:
            subs[os.path.basename(vhdr)[:6]] = (sigs, 1000.0, None)
    return subs


def main():
    os.makedirs(OUT, exist_ok=True)
    cohorts = [("Atlas 204 Hz (n≈49)", load_atlas(), {"spindle": BANDS["spindle"]}),
               ("Falach 1 kHz", load_falach(), BANDS),
               ("Zurich 2 kHz", load_zurich(), BANDS)]
    fig, axes = plt.subplots(3, 2, figsize=(12, 11))
    summary = {}
    for i, (cname, subs, bands) in enumerate(cohorts):
        print(f"== {cname}: {len(subs)} subjects", flush=True)
        res = analyse(subs, cname, bands)
        summary[cname] = {b: {k: res[b][k] for k in ["n_subj", "up_z", "up_ci", "split_r", "split_p"]}
                          for b in res}
        for j, b in enumerate(["spindle", "ripple"]):
            ax = axes[i][j]
            if b not in res:
                ax.text(0.5, 0.5, "not resolvable\n(204 Hz Nyquist)" if b == "ripple" and i == 0
                        else "insufficient data", ha="center", va="center", fontsize=10, color="0.5")
                ax.set_xticks([]); ax.set_yticks([])
                ax.set_title(f"{cname} — SO→{b}", fontsize=9); continue
            R = res[b]; col = "tab:orange" if b == "spindle" else "tab:red"
            ax.plot(R["t"], R["grand"], color=col, lw=2)
            ax.fill_between(R["t"], R["lo"], R["hi"], color=col, alpha=0.22)
            ax.axvline(0, color="k", ls=":", lw=0.8); ax.axhline(0, color="0.7", lw=0.6)
            ax.axvspan(0, 0.75, color="tab:blue", alpha=0.06)
            sig_ci = "✓>0" if R["up_ci"][0] > 0 else "n.s."
            ax.set_title(f"{cname} — SO→{b}   up z={R['up_z']:.3f} [{R['up_ci'][0]:.3f},{R['up_ci'][1]:.3f}] {sig_ci}\n"
                         f"split-half r={R['split_r']:.2f} p={R['split_p']:.1e}  (n={R['n_subj']})", fontsize=8.5)
            ax.set_xlabel("time from SO trough (s)"); ax.set_ylabel(f"{b} env (z)")
    fig.suptitle("SO-nesting across three independent cohorts (mean ± CI over subjects; "
                 "up-state = 0–0.75 s; no sorting — de-circularized)", y=1.0, fontsize=12)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"nesting_cohort_validation.{ext}"), dpi=150, bbox_inches="tight")
    pd.Series(summary).to_json(os.path.join(OUT, "cohort_validation_stats.json"))
    print(f"[cohorts] figure + stats in {OUT}")


if __name__ == "__main__":
    main()
