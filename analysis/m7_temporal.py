"""
M7 - Temporal (time-domain) figures of the NREM nesting, on the 1 kHz M6 data.

The coupling modules (M2/M6) are frequency-domain. This one shows the phenomenon in TIME:
the slow oscillation, and the spindle/ripple bursts that ride it, aligned to the SO trough.

Panels:
  A. Example segment  - raw LFP + SO band, with spindle-band and ripple-band envelopes below,
     so you can see bursts landing on SO up-states.
  B. SO-trough-triggered grand average (all channels) - mean SO wave, mean spindle envelope,
     mean ripple envelope in +/-1.5 s. Shows the temporal cascade (SO -> spindle -> ripple).
  C. SO-trough-triggered time-frequency map (best channel) - power across 2-150 Hz vs time
     around the SO event. The spindle and ripple bands light up near the up-state.

Reuses the M6 loader (bipolar, SOZ-dropped, IED-masked, 50 Hz notched).
"""
import os, glob
import numpy as np
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import m6_so_ripple_coupling as m6

OUT = os.path.join(m6.ROOT, "outputs", "m7_temporal")
SF = 1000.0
HALF = int(1.5 * SF)          # +/-1.5 s window
MIN_DIST = int(0.6 * SF)      # min spacing between SO troughs


def zenv(x, sf, band):
    e = np.abs(signal.hilbert(m6.bandpass(x, sf, *band)))
    return (e - e.mean()) / (e.std() + 1e-12)


def so_troughs(x, sf, mask):
    """SO troughs (deep negative peaks of the SO-filtered signal) within valid samples."""
    so = m6.bandpass(x, sf, *m6.SO_BAND)
    thr = so.std()
    idx, _ = signal.find_peaks(-so, height=thr, distance=MIN_DIST)
    idx = idx[(idx > HALF) & (idx < len(x) - HALF)]
    return np.array([i for i in idx if mask[i]]), so


def main():
    os.makedirs(OUT, exist_ok=True)
    t = np.arange(-HALF, HALF) / SF

    sum_so = np.zeros(2 * HALF); sum_sp = np.zeros(2 * HALF); sum_rp = np.zeros(2 * HALF)
    freqs = np.logspace(np.log10(2), np.log10(150), 40)
    tf_sum = np.zeros((len(freqs), 2 * HALF))
    n_ev = 0
    best = {"n": -1}
    for sd in sorted(glob.glob(os.path.join(m6.BIDS, "sub-*"))):
        try:
            bip, sf, ied, reg = m6.load_subject(sd)
        except Exception:
            continue
        for label, x in bip.items():
            mask = m6.valid_mask(len(x), sf, ied)
            if mask.sum() < 5 * sf:
                continue
            ev, so = so_troughs(x, sf, mask)
            if len(ev) < 5:
                continue
            sp = zenv(x, sf, m6.SPINDLE); rp = zenv(x, sf, m6.RIPPLE)
            for e in ev:
                sl = slice(e - HALF, e + HALF)
                sum_so += so[sl]; sum_sp += sp[sl]; sum_rp += rp[sl]
            # time-frequency accumulation (z-scored envelope per freq, per channel)
            for fi, f0 in enumerate(freqs):
                bw = max(1.0, f0 * 0.3)
                env = np.abs(signal.hilbert(m6.bandpass(x, sf, max(0.5, f0 - bw), f0 + bw)))
                env = (env - env.mean()) / (env.std() + 1e-12)
                for e in ev:
                    tf_sum[fi] += env[e - HALF:e + HALF]
            n_ev += len(ev)
            if len(ev) > best["n"]:
                best = {"n": len(ev), "x": x, "so": so, "sp": sp, "rp": rp,
                        "ev": ev, "label": label, "reg": reg[label],
                        "sub": os.path.basename(sd), "sf": sf}
    print(f"[M7] {n_ev} SO troughs pooled; best channel {best['label']} "
          f"({best['reg']}, {best['sub']}) with {best['n']} events")
    avg_so, avg_sp, avg_rp = sum_so / n_ev, sum_sp / n_ev, sum_rp / n_ev

    # ---- Panel A + B ----
    fig, ax = plt.subplots(2, 1, figsize=(11, 8))
    # A: example ~6 s from best channel around a strong event
    x = best["x"]; so = best["so"]; e0 = best["ev"][len(best["ev"]) // 2]
    seg = slice(e0 - int(3 * SF), e0 + int(3 * SF))
    tt = (np.arange(seg.stop - seg.start) - int(3 * SF)) / SF
    ax[0].plot(tt, x[seg] * 1e6, color="0.6", lw=0.5, label="raw LFP")
    ax[0].plot(tt, so[seg] * 1e6, color="tab:blue", lw=1.8, label="slow osc (0.5-1.25 Hz)")
    ax[0].plot(tt, 4 * np.abs(signal.hilbert(m6.bandpass(x, SF, *m6.SPINDLE)))[seg] * 1e6 + 400,
               color="tab:orange", lw=1, label="spindle env (x4, offset)")
    ax[0].plot(tt, 8 * np.abs(signal.hilbert(m6.bandpass(x, SF, *m6.RIPPLE)))[seg] * 1e6 + 700,
               color="tab:red", lw=1, label="ripple env (x8, offset)")
    ax[0].set_title(f"A. Example segment - {best['label']} ({best['reg']}, {best['sub']}): "
                    "spindle/ripple bursts on the slow oscillation")
    ax[0].set_xlabel("time (s)"); ax[0].set_ylabel("uV (offset)"); ax[0].legend(fontsize=8, ncol=2)

    # B: SO-trough-triggered grand average
    ax2 = ax[1]; ax2b = ax2.twinx()
    ax2.plot(t, avg_so * 1e6, color="tab:blue", lw=2, label="SO wave (uV)")
    ax2b.plot(t, avg_sp, color="tab:orange", lw=1.8, label="spindle env (z)")
    ax2b.plot(t, avg_rp, color="tab:red", lw=1.8, label="ripple env (z)")
    ax2.axvline(0, color="k", ls=":", lw=0.8)
    ax2.set_xlabel("time from SO trough (s)"); ax2.set_ylabel("SO wave (uV)", color="tab:blue")
    ax2b.set_ylabel("band envelope (z)")
    ax2.set_title(f"B. SO-trough-triggered grand average ({n_ev} events, all channels): "
                  "spindle & ripple rise on the SO up-state")
    l1, la1 = ax2.get_legend_handles_labels(); l2, la2 = ax2b.get_legend_handles_labels()
    ax2.legend(l1 + l2, la1 + la2, fontsize=8, loc="upper right")
    plt.tight_layout()
    for e in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"m7_so_triggered_nesting.{e}"), dpi=150)

    # ---- Panel C: SO-triggered time-frequency (grand average, all channels) ----
    tf = tf_sum / n_ev
    tf = (tf - tf.mean(axis=1, keepdims=True)) / (tf.std(axis=1, keepdims=True) + 1e-12)
    plt.figure(figsize=(9, 5.5))
    plt.pcolormesh(t, freqs, tf, shading="auto", cmap="RdBu_r", vmin=-1.2, vmax=1.2)
    plt.yscale("log"); plt.colorbar(label="power (z per freq)")
    plt.axhspan(*m6.SPINDLE, color="k", alpha=0.06); plt.axhspan(*m6.RIPPLE, color="k", alpha=0.06)
    plt.axvline(0, color="k", ls=":", lw=0.8)
    plt.text(1.2, 13.5, "spindle", fontsize=8); plt.text(1.2, 100, "ripple", fontsize=8)
    plt.xlabel("time from SO trough (s)"); plt.ylabel("frequency (Hz)")
    plt.title(f"C. SO-trough-triggered time-frequency (grand average, {n_ev} events, all channels)")
    plt.tight_layout()
    for e in ("png", "svg"):
        plt.savefig(os.path.join(OUT, f"m7_so_triggered_tf.{e}"), dpi=150)
    print(f"[M7] figures in {OUT}")


if __name__ == "__main__":
    main()
