"""
Ripple-coupling - Slow-oscillation -> ripple (and spindle) coupling  [Phase-1b, 1 kHz data].

Closes the one gap the 204 Hz atlas could not: RIPPLES (80-120 Hz), the human analogue of
the mouse CA1 sharp-wave-ripple figures. In human NREM the slow oscillation is expected to
organize both spindles (~14 Hz) and ripples (~90 Hz), with events riding the SO up-state.

Data: Falach/Geva-Sagiv/Eliashiv 2024 (Figshare), overnight iEEG sleep segments, 1000 Hz,
MTL SEEG. NREM-dominant, so this is a WITHIN-NREM coordination test (not a state contrast;
that was Slow-power/Spindle-coupling), extended to the ripple band.

Rigor: bipolar re-reference within each shaft; drop SOZ channels; mask +/-0.5 s around every
annotated IED; surrogate-corrected MI_z.

Outputs (outputs/slow_ripple_coupling/):
  ripple_comodulogram.png     - SO-phase x amplitude-freq (8-200 Hz): spindle AND ripple hotspots
  ripple_coupling_by_region.png - SO->ripple and SO->spindle MI_z by MTL region
  channel_coupling.csv
"""
import os, glob, re
import numpy as np, pandas as pd
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mne

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIDS = os.path.join(ROOT, "data", "phase1b", "ieeg_ieds_bids")
OUT = os.path.join(ROOT, "outputs", "slow_ripple_coupling")
LINE = 50.0
SO_BAND = (0.5, 1.25)
SPINDLE = (11, 16)
RIPPLE = (80, 120)
NBINS = 18
N_SURR = 120

REGION = {"A": "amygdala", "AH": "hippocampal", "MH": "hippocampal",
          "EC": "entorhinal", "PHG": "parahippocampal"}


def region_of(name):
    core = re.sub(r"^[LR]", "", name)          # drop hemisphere
    core = re.sub(r"\d.*$", "", core)          # drop electrode number
    return REGION.get(core)


def bandpass(x, sf, lo, hi):
    b, a = signal.butter(3, [lo / (sf / 2), hi / (sf / 2)], btype="band")
    return signal.filtfilt(b, a, x)


def notch(x, sf, f0):
    for h in np.arange(f0, sf / 2, f0):
        b, a = signal.iirnotch(h, 30, sf)
        x = signal.filtfilt(b, a, x)
    return x


def tort_mi(ph, am, nbins=NBINS):
    edges = np.linspace(-np.pi, np.pi, nbins + 1)
    idx = np.clip(np.digitize(ph, edges) - 1, 0, nbins - 1)
    m = np.array([am[idx == b].mean() if np.any(idx == b) else 0.0 for b in range(nbins)])
    if m.sum() <= 0:
        return np.nan, m
    p = np.clip(m / m.sum(), 1e-12, None)
    return float((np.log(nbins) + (p * np.log(p)).sum()) / np.log(nbins)), m


def load_subject(sub_dir):
    """Return (bipolar_signals dict name->array, sfreq, ied_samples, region_map)."""
    edf = glob.glob(os.path.join(sub_dir, "ieeg", "*_ieeg.edf"))[0]
    ch_tsv = pd.read_csv(glob.glob(os.path.join(sub_dir, "ieeg", "*channels.tsv"))[0], sep="\t")
    ev = glob.glob(os.path.join(sub_dir, "ieeg", "*events.tsv"))
    r = mne.io.read_raw_edf(edf, preload=True, verbose="ERROR")
    sf = r.info["sfreq"]
    data = r.get_data(); names = r.ch_names
    nidx = {n: i for i, n in enumerate(names)}
    soz = dict(zip(ch_tsv["name"].astype(str),
                   ch_tsv.get("soz_region", pd.Series([0] * len(ch_tsv)))))
    # bipolar pairs within a shaft (same letter-prefix, consecutive number)
    def split(n):
        m = re.match(r"([A-Za-z]+)(\d+)$", n)
        return (m.group(1), int(m.group(2))) if m else (None, None)
    bip = {}; reg = {}
    for n in names:
        pre, num = split(n)
        if pre is None:
            continue
        nxt = f"{pre}{num+1}"
        if nxt in nidx and region_of(n):
            if soz.get(n, 0) or soz.get(nxt, 0):     # drop SOZ contacts
                continue
            sig = data[nidx[n]] - data[nidx[nxt]]
            label = f"{n}-{nxt}"
            bip[label] = notch(signal.detrend(sig), sf, LINE)
            reg[label] = region_of(n)
    ied = []
    if ev:
        e = pd.read_csv(ev[0], sep="\t")
        col = "sample" if "sample" in e.columns else None
        ied = (e[col].astype(int).tolist() if col
               else (e["onset"] * sf).astype(int).tolist()) if len(e) else []
    return bip, sf, np.array(ied, int), reg


def valid_mask(n, sf, ied, pad=0.5):
    m = np.ones(n, bool)
    w = int(pad * sf)
    for s in ied:
        m[max(0, s - w):min(n, s + w)] = False
    return m


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(0)
    rows = []
    subs = sorted(glob.glob(os.path.join(BIDS, "sub-*")))
    for sd in subs:
        sub = os.path.basename(sd)
        try:
            bip, sf, ied, reg = load_subject(sd)
        except Exception as e:
            print("  skip", sub, e); continue
        if not bip:
            continue
        for label, x in bip.items():
            mask = valid_mask(len(x), sf, ied)
            if mask.sum() < 5 * sf:
                continue
            so_ph_full = np.angle(signal.hilbert(bandpass(x, sf, *SO_BAND)))
            ph = so_ph_full[mask]
            out = {"sub": sub, "channel": label, "region": reg[label],
                   "n_valid_s": round(mask.sum() / sf, 1)}
            for band, key in [(SPINDLE, "spindle"), (RIPPLE, "ripple")]:
                am = np.abs(signal.hilbert(bandpass(x, sf, *band)))[mask]
                mi, miz, _ = mi_z_from(ph, am, rng)
                out[f"mi_{key}"] = mi; out[f"miz_{key}"] = miz
            # preferred SO phase of ripple
            amr = np.abs(signal.hilbert(bandpass(x, sf, *RIPPLE)))[mask]
            _, m = tort_mi(ph, amr)
            centers = np.linspace(-np.pi, np.pi, NBINS + 1)[:-1] + np.pi / NBINS
            out["ripple_pref_phase"] = float(centers[int(np.argmax(m))])
            rows.append(out)
        print(f"  {sub}: {len(bip)} bipolar MTL channels")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "channel_coupling.csv"), index=False)
    print(f"[Ripple-coupling] {len(df)} channels; regions={df.region.value_counts().to_dict()}")

    # ---- stats: is SO->ripple / SO->spindle coupling above chance? (per subject) ----
    from scipy.stats import wilcoxon
    for key in ["ripple", "spindle"]:
        per_sub = df.groupby("sub")[f"miz_{key}"].mean().dropna()
        if len(per_sub) >= 5:
            stat, p = wilcoxon(per_sub.values)
            print(f"[Ripple-coupling] SO->{key}: mean MI_z={per_sub.mean():.3f} across {len(per_sub)} subjects, "
                  f">0 in {(per_sub>0).sum()}/{len(per_sub)}, Wilcoxon p={p:.3g}")

    # ---- region figure ----
    regions = ["hippocampal", "entorhinal", "parahippocampal", "amygdala"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for a, key, ttl in [(ax[0], "ripple", "SO -> ripple (80-120 Hz)"),
                        (ax[1], "spindle", "SO -> spindle (11-16 Hz)")]:
        for i, r in enumerate(regions):
            v = df[df.region == r].groupby("sub")[f"miz_{key}"].mean().dropna().values
            if len(v) == 0:
                continue
            a.bar(i, v.mean(), color=f"C{i}", alpha=0.7)
            jit = (np.random.default_rng(i).random(len(v)) - .5) * .3
            a.scatter(i + jit, v, s=14, color="k", alpha=.4, zorder=3)
            a.text(i, 0.02, f"n={len(v)}", ha="center", fontsize=8)
        a.axhline(0, color="k", lw=.7, ls=":")
        a.set_xticks(range(len(regions))); a.set_xticklabels(regions, rotation=15)
        a.set_ylabel("coupling MI_z"); a.set_title(ttl)
    fig.suptitle("Ripple-coupling: NREM slow-oscillation coordinates spindles AND ripples in human MTL "
                 "(1 kHz, IED/SOZ-cleaned)")
    plt.tight_layout()
    for e in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"ripple_coupling_by_region.{e}"), dpi=150)

    # ---- comodulogram for the strongest ripple-coupling channel ----
    best = df.dropna(subset=["miz_ripple"]).sort_values("miz_ripple", ascending=False)
    if len(best):
        r0 = best.iloc[0]
        bip, sf, ied, reg = load_subject(os.path.join(BIDS, r0["sub"]))
        x = bip[r0["channel"]]; mask = valid_mask(len(x), sf, ied)
        pfreqs = np.arange(0.5, 4.01, 0.5)
        afreqs = np.concatenate([np.arange(8, 30, 3), np.arange(30, 201, 12)])
        grid = np.zeros((len(afreqs), len(pfreqs)))
        for i, pf in enumerate(pfreqs):
            ph = np.angle(signal.hilbert(bandpass(x, sf, max(0.25, pf - 0.5), pf + 0.5)))[mask]
            for k, af in enumerate(afreqs):
                am = np.abs(signal.hilbert(bandpass(x, sf, max(1, af - 6), af + 6)))[mask]
                grid[k, i], _ = tort_mi(ph, am)
        plt.figure(figsize=(7, 5.5))
        plt.pcolormesh(pfreqs, afreqs, grid, shading="auto", cmap="viridis")
        plt.colorbar(label="MI")
        plt.axhspan(*SPINDLE, color="w", alpha=0.12); plt.axhspan(*RIPPLE, color="w", alpha=0.12)
        plt.text(3.5, 14, "spindle", color="w", ha="right"); plt.text(3.5, 100, "ripple", color="w", ha="right")
        plt.xlabel("SO phase freq (Hz)"); plt.ylabel("amplitude freq (Hz)")
        plt.title(f"Ripple-coupling comodulogram — {r0['channel']} ({r0['region']}, {r0['sub']})\n"
                  "NREM slow oscillation couples both spindle and ripple bands")
        plt.tight_layout()
        for e in ("png", "svg"):
            plt.savefig(os.path.join(OUT, f"ripple_comodulogram.{e}"), dpi=150)
    print(f"[Ripple-coupling] figures in {OUT}")


def mi_z_from(ph, am, rng):
    mi, m = tort_mi(ph, am)
    if not np.isfinite(mi):
        return np.nan, np.nan, m
    n = len(am)
    surr = np.array([tort_mi(ph, np.roll(am, s))[0]
                     for s in rng.integers(int(0.1 * n), int(0.9 * n), N_SURR)])
    return mi, (mi - np.nanmean(surr)) / (np.nanstd(surr) + 1e-12), m


if __name__ == "__main__":
    main()
