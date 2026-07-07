"""
M8 - SO-phase of EXPERT-MARKED ripples & fast-ripples  [ds003498, Zurich, 2 kHz].

Independent 3rd cohort, and a different, stronger design than M6: instead of amplitude-envelope
coupling, we use the dataset's **expert/detector-marked** ripple and fast-ripple events and ask
what slow-oscillation phase each event lands on. If events cluster at a preferred SO phase
(Rayleigh test), the slow oscillation organizes them -- the event-based version of the Kipnis /
mouse ripple result, in human MTL.

Bonus the 204 Hz atlas and even the 1 kHz Falach set could not give: **fast ripples (200-500 Hz)**.

Data: ds003498 (Fedele/Sarnthein, OpenNeuro), interictal slow-wave-sleep 5-min runs, 2000 Hz,
BrainVision. Channels are MTL depth: A=amygdala, AH/H/PH/P=hippocampus, EC/E=entorhinal.
Events: trial_type 'ripple_<chan>' / 'fr_<chan>' / 'frandr_<chan>', chan like 'AL1-2' (bipolar).

CAVEAT: these are HFO markings made for epilepsy localization, so a fraction are PATHOLOGICAL
(not physiological sleep ripples). This tests SO-coupling of marked ripples, not exclusively
physiological ones. Complementary to M6 (SOZ-excluded, IED-masked, envelope-based).
"""
import os, glob, re
import numpy as np, pandas as pd
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mne

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EEG = os.path.join(ROOT, "data", "phase1c", "eeg")
OUT = os.path.join(ROOT, "outputs", "m8_zurich_ripple_so_phase")
LINE = 50.0
SO_BAND = (0.5, 1.25)
REGIONS = ["hippocampal", "entorhinal", "amygdala"]
REGMAP = {"A": "amygdala", "AH": "hippocampal", "H": "hippocampal", "PH": "hippocampal",
          "P": "hippocampal", "EC": "entorhinal", "E": "entorhinal"}


def region_of(label):
    m = re.match(r"[A-Za-z]+", label)
    if not m:
        return None
    letters = m.group()
    core = letters[:-1] if letters[-1] in "LR" else letters
    return REGMAP.get(core)


SO_SF = 500.0   # SO phase computed at this rate (0.5-1.25 Hz is unstable at 2 kHz)


def bandpass(x, sf, lo, hi):
    sos = signal.butter(3, [lo / (sf / 2), hi / (sf / 2)], btype="band", output="sos")
    return signal.sosfiltfilt(sos, x)


def notch(x, sf, f0):
    for h in np.arange(f0, min(sf / 2, 300), f0):
        b, a = signal.iirnotch(h, 30, sf)
        x = signal.filtfilt(b, a, x)
    return x


def rayleigh(phases):
    n = len(phases)
    if n < 5:
        return np.nan, np.nan, np.nan
    C, S = np.cos(phases).sum(), np.sin(phases).sum()
    R = np.sqrt(C**2 + S**2) / n            # mean resultant length
    z = n * R**2
    p = np.exp(-z) * (1 + (2 * z - z**2) / (4 * n))   # Zar approximation
    pref = np.arctan2(S, C)
    return R, p, pref


def main():
    os.makedirs(OUT, exist_ok=True)
    phases = {r: {"ripple": [], "fr": []} for r in REGIONS}
    per_ch = []
    subs = sorted(glob.glob(os.path.join(EEG, "sub-*_ieeg.vhdr")))
    for vhdr in subs:
        sub = os.path.basename(vhdr)[:6]
        r = mne.io.read_raw_brainvision(vhdr, preload=True, verbose="ERROR")
        sf = r.info["sfreq"]
        data = r.get_data(); names = r.ch_names
        nidx = {n: i for i, n in enumerate(names)}

        def split(n):
            m = re.match(r"([A-Za-z]+)(\d+)$", n)
            return (m.group(1), int(m.group(2))) if m else (None, None)
        # bipolar SO phase per pair, labelled like the events ('AL1-2')
        so_ph = {}
        for n in names:
            pre, num = split(n)
            if pre is None:
                continue
            nxt = f"{pre}{num+1}"
            if nxt in nidx and region_of(n):
                sig = signal.detrend(data[nidx[n]] - data[nidx[nxt]])
                ds = signal.resample_poly(sig, 1, int(round(sf / SO_SF)))   # -> SO_SF Hz
                ds = notch(ds, SO_SF, LINE)
                so_ph[f"{pre}{num}-{num+1}"] = np.angle(signal.hilbert(bandpass(ds, SO_SF, *SO_BAND)))

        ev = pd.read_csv(os.path.join(EEG, f"{sub}_run-01_events.tsv"), sep="\t")
        ev["etype"] = ev["trial_type"].astype(str).str.split("_").str[0]
        ev["chan"] = ev["trial_type"].astype(str).str.split("_", n=1).str[1]
        cnt = {}
        for _, e in ev.iterrows():
            et = "fr" if e["etype"] in ("fr", "frandr") else ("ripple" if e["etype"] == "ripple" else None)
            ch = e["chan"]
            if et is None or ch not in so_ph:
                continue
            reg = region_of(ch)
            if reg not in phases:
                continue
            s = int(round(e["onset"] * SO_SF))
            if 0 <= s < so_ph[ch].shape[0] and np.isfinite(so_ph[ch][s]):
                phases[reg][et].append(so_ph[ch][s])
                cnt[(reg, et)] = cnt.get((reg, et), 0) + 1
        print(f"  {sub}: {len(so_ph)} MTL bipolar ch, events by region: "
              f"{ {k: v for k, v in cnt.items()} }")

    # ---- stats + figure ----
    rows = []
    for reg in REGIONS:
        for et in ["ripple", "fr"]:
            ph = np.array(phases[reg][et])
            R, p, pref = rayleigh(ph)
            rows.append({"region": reg, "event": et, "n": len(ph),
                         "mrl": R, "rayleigh_p": p, "pref_phase_deg": np.degrees(pref)})
    rep = pd.DataFrame(rows)
    rep.to_csv(os.path.join(OUT, "so_phase_stats.csv"), index=False)
    print("\n[M8] SO-phase clustering (Rayleigh):")
    print(rep.to_string(index=False))

    # polar histograms: ripple (top) & fast-ripple (bottom) x region
    fig, axes = plt.subplots(2, len(REGIONS), figsize=(4 * len(REGIONS), 8),
                             subplot_kw={"projection": "polar"})
    for j, reg in enumerate(REGIONS):
        for i, (et, col) in enumerate([("ripple", "tab:blue"), ("fr", "tab:red")]):
            ax = axes[i, j]; ph = np.array(phases[reg][et])
            if len(ph) >= 5:
                ax.hist(ph, bins=24, color=col, alpha=0.8)
                R, p, pref = rayleigh(ph)
                ax.plot([pref, pref], [0, ax.get_ylim()[1]], color="k", lw=2)
                ax.set_title(f"{reg} {et}\nn={len(ph)} R={R:.3f} p={p:.1e}", fontsize=9)
            else:
                ax.set_title(f"{reg} {et}\n(n={len(ph)})", fontsize=9)
    fig.suptitle("M8 (Zurich ds003498, 2 kHz): SO phase of expert-marked ripples (blue) & "
                 "fast-ripples (red)\nblack line = preferred phase; clustering = SO organizes the events",
                 fontsize=11)
    plt.tight_layout()
    for e in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"m8_ripple_so_phase.{e}"), dpi=150)
    print(f"\n[M8] figures in {OUT}")


if __name__ == "__main__":
    main()
