"""
Ripple-triggered slow-wave average in human MTL (Zurich ds003498, expert-marked ripples).

The direct human parallel to the mouse spike-triggered LFP (Kipnis Fig 1c): trigger on each
EXPERT-MARKED ripple and average the surrounding slow-wave (0.5-4 Hz) LFP. If ripples ride a
consistent slow-wave field, the average shows a clear slow wave; if not, it is flat.

Triggering on independently-marked ripples => non-circular. Bipolar polarity varies by channel, so
we work PER CHANNEL and test significance against a shuffled-ripple-time surrogate (peak-to-peak
amplitude, which is polarity-free). Complements hfo_slow_phase (the SO-phase Rayleigh test) and the
phase-amplitude coupling: this is its time-domain, waveform form.
"""
import os, glob, re
import numpy as np, pandas as pd
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mne
import hfo_slow_phase as zur

OUT = os.path.join(zur.ROOT, "outputs", "ripple_triggered_slowwave")
SO_SF = 500.0; SO_BAND = (0.5, 4.0)
HALF = int(1.5 * SO_SF); NSURR = 200; MIN_RIPPLES = 30


def so_signal(sig2k, sf):
    ds = signal.resample_poly(sig2k, 1, int(round(sf / SO_SF)))
    sos = signal.butter(3, [SO_BAND[0]/(SO_SF/2), SO_BAND[1]/(SO_SF/2)], btype="band", output="sos")
    return signal.sosfiltfilt(sos, ds)


def trig_avg(so, onsets):
    w = [so[o-HALF:o+HALF] for o in onsets if HALF <= o < len(so)-HALF]
    return np.mean(w, 0) if len(w) >= MIN_RIPPLES else None, len(w)


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(0)
    chans = []                       # (sub, label, region, avg, z, n)
    for vhdr in sorted(glob.glob(os.path.join(zur.EEG, "sub-*_ieeg.vhdr"))):
        sub = os.path.basename(vhdr)[:6]
        r = mne.io.read_raw_brainvision(vhdr, preload=True, verbose="ERROR")
        sf = r.info["sfreq"]; data = r.get_data(); names = r.ch_names
        nidx = {n: i for i, n in enumerate(names)}
        ev = pd.read_csv(os.path.join(zur.EEG, f"{sub}_run-01_events.tsv"), sep="\t")
        ev["etype"] = ev["trial_type"].astype(str).str.split("_").str[0]
        ev["chan"] = ev["trial_type"].astype(str).str.split("_", n=1).str[1]
        rip = ev[ev.etype == "ripple"]
        def split(n):
            m = re.match(r"([A-Za-z]+)(\d+)$", n); return (m.group(1), int(m.group(2))) if m else (None, None)
        for n in names:
            pre, num = split(n)
            if pre is None:
                continue
            nxt = f"{pre}{num+1}"; lab = f"{pre}{num}-{num+1}"
            reg = zur.region_of(lab)
            if nxt not in nidx or reg is None:
                continue
            onsets = rip[rip.chan == lab]["onset"].values
            if len(onsets) < MIN_RIPPLES:
                continue
            sig = zur.notch(signal.detrend(data[nidx[n]] - data[nidx[nxt]]), sf, zur.LINE)
            so = so_signal(sig, sf)
            o500 = (onsets * SO_SF).astype(int)
            avg, n_used = trig_avg(so, o500)
            if avg is None:
                continue
            ptp = avg.max() - avg.min()
            surr = []
            for _ in range(NSURR):
                r_on = rng.integers(HALF, len(so) - HALF, len(o500))
                a, _ = trig_avg(so, r_on)
                surr.append(a.max() - a.min())
            z = (ptp - np.mean(surr)) / (np.std(surr) + 1e-12)
            chans.append((sub, lab, reg, avg, float(z), n_used))
        print(f"  {sub}: {sum(1 for c in chans if c[0]==sub)} channels", flush=True)

    df = pd.DataFrame([(c[0], c[1], c[2], c[4], c[5]) for c in chans],
                      columns=["sub", "channel", "region", "z", "n_ripples"])
    df.to_csv(os.path.join(OUT, "channel_z.csv"), index=False)
    sig = (df.z > 2).sum()
    print(f"[ripple-STA] {len(df)} channels; slow-wave modulation z>2 in {sig}/{len(df)} "
          f"({100*sig/len(df):.0f}%); median z={df.z.median():.2f}")

    # figure: example ripple-triggered slow waves (top by z) + z distribution
    top = sorted(chans, key=lambda c: -c[4])[:12]
    t = np.arange(-HALF, HALF) / SO_SF
    fig = plt.figure(figsize=(15, 7))
    gs = fig.add_gridspec(3, 5)
    for i, c in enumerate(top):
        ax = fig.add_subplot(gs[i // 5, i % 5])
        ax.plot(t, c[3] * 1e6, color="tab:blue", lw=1.5)
        ax.axvline(0, color="tab:red", ls=":", lw=1)
        ax.set_title(f"{c[0]} {c[1]} ({c[2]})\nz={c[4]:.1f}, n={c[5]}", fontsize=7)
        ax.tick_params(labelsize=6); ax.set_xticks([-1, 0, 1])
    axh = fig.add_subplot(gs[:, 4] if False else gs[2, 3:])
    axh.hist(df.z, bins=25, color="0.6"); axh.axvline(2, color="tab:red", ls="--", lw=1)
    axh.set_xlabel("ripple-triggered slow-wave z (vs shuffled)"); axh.set_ylabel("channels")
    axh.set_title(f"z>2 in {sig}/{len(df)} channels", fontsize=8)
    fig.suptitle("Ripple-triggered slow-wave average in human MTL (Zurich, expert-marked ripples) — "
                 "red line = ripple; each panel = one channel", fontsize=12)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"ripple_triggered_slowwave.{ext}"), dpi=150, bbox_inches="tight")
    print(f"[ripple-STA] figure + csv in {OUT}")


if __name__ == "__main__":
    main()
