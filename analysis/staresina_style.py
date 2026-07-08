"""
Staresina 2015 (Nat Neurosci) — style nesting figure, reproduced on our data (HUP165 full night,
1024 Hz continuous NREM). Event-based, like the original:

  A. SO-trough-triggered spindle amplitude   (spindles ride the SO up-state)
  B. preferred SO phase of spindle peaks       (polar; Rayleigh) — example channel + % channels sig
  C. spindle-trough-triggered ripple amplitude (ripples ride spindle troughs) — the weak leg here
  D. preferred spindle phase of ripple peaks    (polar; Rayleigh) — example channel + % channels sig

Events detected by standard amplitude thresholds. Aggregate panels (A,C) use amplitude envelopes
(polarity-free). Phase panels (B,D): per-channel Rayleigh (absolute phase is montage-arbitrary, so
we report clustering, not absolute phase, and show one example channel).
"""
import os, glob, re, argparse
import numpy as np, pandas as pd
from scipy import signal
from scipy.stats import circmean
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from hup165_night_coupling import notch, DIR
from hup_night_cohort import bipolar_generic

SF = 1024.0
SO_B = (0.5, 1.25); SP_B = (11, 16); RP_B = (80, 120)
OUT = os.path.join(os.path.dirname(DIR), "..", "..", "outputs", "staresina_style")
OUT = os.path.normpath(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                    "outputs", "staresina_style"))


def bp(x, lo, hi):
    return signal.filtfilt(*signal.butter(3, [lo/(SF/2), hi/(SF/2)], btype="band"), x)


def rayleigh(ph):
    n = len(ph)
    if n < 20:
        return np.nan, np.nan
    R = np.abs(np.mean(np.exp(1j * ph))); z = n * R**2
    return R, float(np.exp(-z) * (1 + (2*z - z**2)/(4*n)))


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--sp-thr',type=float,default=1.5); ap.add_argument('--rp-thr',type=float,default=2.0); ap.add_argument('--tag',default=''); A_=ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    idx = pd.read_csv(os.path.join(DIR, "index.csv"))
    hw_so = int(1.5 * SF); hw_sp = int(0.4 * SF)
    A = {"sum": np.zeros(2*hw_so), "n": 0}          # SO-trough -> spindle env
    C = {"sum": np.zeros(2*hw_sp), "n": 0}          # spindle-trough -> ripple env
    sp_so_phase = {}                                 # channel -> list (spindle peak's SO phase)
    rp_sp_phase = {}                                 # channel -> list (ripple peak's spindle phase)

    for _, row in idx.iterrows():
        b = np.load(os.path.join(DIR, row["file"]), allow_pickle=True)
        data = b["data"]; names = [str(n) for n in b["ch_names"]]
        for lab, sig in bipolar_generic(data, names).items():
            x = notch(signal.detrend(sig), SF)
            if len(x) < 2*hw_so or np.std(x) < 1e-9:
                continue
            so = bp(x, *SO_B); so_ph = np.angle(signal.hilbert(so))
            sp = bp(x, *SP_B); sp_env = np.abs(signal.hilbert(sp)); sp_ph = np.angle(signal.hilbert(sp))
            rp_env = np.abs(signal.hilbert(bp(x, *RP_B)))
            spz = (sp_env - sp_env.mean())/(sp_env.std()+1e-12)
            rpz = (rp_env - rp_env.mean())/(rp_env.std()+1e-12)
            # events
            so_tr, _ = signal.find_peaks(-so, height=so.std(), distance=int(0.8*SF))
            sp_pk, _ = signal.find_peaks(spz, height=A_.sp_thr, distance=int(0.3*SF))
            rp_pk, _ = signal.find_peaks(rpz, height=A_.rp_thr, distance=int(0.02*SF))
            sp_tr, _ = signal.find_peaks(-sp, height=sp.std(), distance=int(0.05*SF))
            # A: SO-trough-triggered spindle env
            for e in so_tr[(so_tr > hw_so) & (so_tr < len(x)-hw_so)]:
                A["sum"] += spz[e-hw_so:e+hw_so]; A["n"] += 1
            # C: spindle-trough-triggered ripple env (only spindle troughs inside a spindle event)
            in_sp = spz > 1.0
            for e in sp_tr[(sp_tr > hw_sp) & (sp_tr < len(x)-hw_sp)]:
                if in_sp[e]:
                    C["sum"] += rpz[e-hw_sp:e+hw_sp]; C["n"] += 1
            # B/D phases
            sp_so_phase.setdefault(lab, []).extend(so_ph[sp_pk].tolist())
            rp_sp_phase.setdefault(lab, []).extend(sp_ph[rp_pk].tolist())

    # per-channel Rayleigh
    def frac_sig(d):
        rs = [(rayleigh(np.array(v))) for v in d.values() if len(v) >= 20]
        rs = [(R, p) for R, p in rs if np.isfinite(p)]
        sig = sum(1 for _, p in rs if p < 0.05)
        return sig, len(rs), rs
    sB, nB, _ = frac_sig(sp_so_phase); sD, nD, _ = frac_sig(rp_sp_phase)
    exB = max(sp_so_phase.items(), key=lambda kv: len(kv[1]))
    exD = max(rp_sp_phase.items(), key=lambda kv: len(kv[1]))

    t_so = np.arange(-hw_so, hw_so)/SF; t_sp = np.arange(-hw_sp, hw_sp)/SF
    fig = plt.figure(figsize=(14, 8))
    # A
    ax = fig.add_subplot(2, 2, 1)
    ax.plot(t_so, A["sum"]/A["n"], color="tab:orange", lw=2)
    ax.axvline(0, color="k", ls=":", lw=.8); ax.axhline(0, color="0.7", lw=.6); ax.axvspan(0, .75, color="tab:blue", alpha=.06)
    ax.set_title(f"A. SO-trough-triggered SPINDLE amplitude\n(n={A['n']} SO troughs) — spindles ride the up-state")
    ax.set_xlabel("time from SO trough (s)"); ax.set_ylabel("spindle amp (z)")
    # B polar (example channel) + summary
    axB = fig.add_subplot(2, 2, 2, projection="polar")
    ph = np.array(exB[1]); axB.hist(ph, bins=24, color="tab:orange", alpha=.8)
    R, p = rayleigh(ph)
    axB.set_title(f"B. Preferred SO phase of SPINDLE peaks\nex {exB[0]} R={R:.2f} p={p:.1e} | sig in {sB}/{nB} ch")
    # C
    axC = fig.add_subplot(2, 2, 3)
    axC.plot(t_sp, C["sum"]/max(C["n"],1), color="tab:red", lw=2)
    axC.axvline(0, color="k", ls=":", lw=.8); axC.axhline(0, color="0.7", lw=.6)
    axC.set_title(f"C. Spindle-trough-triggered RIPPLE amplitude\n(n={C['n']} spindle troughs) — weak leg in our data")
    axC.set_xlabel("time from spindle trough (s)"); axC.set_ylabel("ripple amp (z)")
    # D polar
    axD = fig.add_subplot(2, 2, 4, projection="polar")
    ph = np.array(exD[1]); axD.hist(ph, bins=24, color="tab:red", alpha=.8)
    R, p = rayleigh(ph)
    axD.set_title(f"D. Preferred SPINDLE phase of RIPPLE peaks\nex {exD[0]} R={R:.2f} p={p:.1e} | sig in {sD}/{nD} ch")
    fig.suptitle("Staresina 2015-style nesting on HUP165 full night (event-based) — "
                 "SO→spindle strong; spindle→ripple weak (as measured)", y=1.0, fontsize=12)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"staresina_style{A_.tag}.{ext}"), dpi=150, bbox_inches="tight")
    print(f"[staresina] A n={A['n']} SO troughs; C n={C['n']} spindle troughs")
    print(f"[staresina] spindle-peak SO-phase clustered in {sB}/{nB} channels; "
          f"ripple-peak spindle-phase clustered in {sD}/{nD} channels")
    print(f"[staresina] figure in {OUT}")


if __name__ == "__main__":
    main()
