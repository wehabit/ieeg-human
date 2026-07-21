"""
Teaching figures: what tests 3A / 3B / 3D actually do to a signal.

Draws idealized SYNTHETIC traces (not patient data) that step through each analysis so the
method is legible at a glance. Colours are consistent everywhere:
    amber = slow oscillation (0.5-1.25 Hz)   teal = spindle (11-16 Hz)   rose = heartbeat / RR
Outputs three PNGs to outputs/signal_tutorial/.

    python analysis/tutorial_signal_walkthrough.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal

SO   = "#d9871f"   # amber
SP   = "#0f9b96"   # teal
HR   = "#d6415f"   # rose
NULL = "#9aa0b4"   # grey
INK  = "#20232f"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(ROOT, "outputs", "signal_tutorial")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#cccccc",
                     "axes.linewidth": 0.8, "figure.dpi": 130})


def envelope(x, sf, lo, hi):
    b, a = signal.butter(3, [lo/(sf/2), hi/(sf/2)], btype="band")
    return np.abs(signal.hilbert(signal.filtfilt(b, a, x)))


# ==================================================================================
# 3D  — spindles ride the up-state  (the 1-second lens)
# ==================================================================================
def fig_3d():
    sf = 400.0
    t = np.arange(0, 3, 1/sf)
    so = np.sin(2*np.pi*0.8*t - 0.4)
    gate = np.clip(np.sin(2*np.pi*0.8*t - 0.4), 0, None)**2          # up-state gate
    spenv = 0.06 + 0.94*gate
    spindle = spenv * np.sin(2*np.pi*13*t)
    lfp = so + 0.6*spindle + 0.04*np.random.default_rng(1).standard_normal(t.size)

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.9))
    fig.suptitle("TEST 3D  ·  Do spindles sit on the slow wave's up-state?   (the ~1-second lens — already done)",
                 fontsize=13, fontweight="bold", x=0.02, ha="left")

    ax[0].plot(t, so, color=SO, lw=2, alpha=.55, label="slow wave (0.5–1.25 Hz)")
    ax[0].plot(t, lfp, color=INK, lw=1, alpha=.9, label="raw LFP")
    ax[0].set_title("① Raw trace", loc="left", fontsize=11)
    ax[0].set_xlabel("time (s)"); ax[0].legend(fontsize=8, loc="upper right")
    ax[0].annotate("fast spindle bursts\nsit on the crests", xy=(0.95, 1.2), xytext=(1.5, 1.9),
                   fontsize=8.5, color=SP, ha="center",
                   arrowprops=dict(arrowstyle="->", color=SP, lw=1.2))

    ax[1].plot(t, so, color=SO, lw=2, label="slow phase")
    ax[1].plot(t, spenv, color=SP, lw=2.2, label="spindle envelope")
    ax[1].set_title("② Split the bands", loc="left", fontsize=11)
    ax[1].set_xlabel("time (s)"); ax[1].legend(fontsize=8, loc="upper right")
    ax[1].annotate("teal humps line up\nwith amber peaks", xy=(1.65, 0.95), xytext=(0.35, 1.15),
                   fontsize=8.5, color=SP, arrowprops=dict(arrowstyle="->", color=SP, lw=1.2))

    # preferred-phase polar
    ax[2].remove()
    axp = fig.add_subplot(1, 3, 3, projection="polar")
    ph = np.angle(signal.hilbert(signal.filtfilt(*signal.butter(3, [0.4/(sf/2), 1.25/(sf/2)], btype="band"), lfp)))
    am = envelope(lfp, sf, 11, 16)
    nb = 18; edges = np.linspace(-np.pi, np.pi, nb+1)
    idx = np.clip(np.digitize(ph, edges)-1, 0, nb-1)
    m = np.array([am[idx == b].mean() if np.any(idx == b) else 0 for b in range(nb)])
    centers = (edges[:-1]+edges[1:])/2
    axp.bar(centers, m, width=2*np.pi/nb, color=SP, alpha=.8, edgecolor="w", linewidth=.5)
    axp.set_title("③ Spindle power by slow phase\n(tall bump = coupling)", fontsize=10)
    axp.set_yticklabels([])

    fig.text(0.5, -0.02, "WHAT YOUR EYE SHOULD SEE:  the fast teal bursts are not scattered — they clump on the amber crests. "
             "Panel ③ turns 'they clump' into one number (modulation index) + one preferred angle.",
             ha="center", fontsize=9.5, style="italic", color="#444")
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    p = os.path.join(OUT, "test_3D_so_spindle.png")
    fig.savefig(p, bbox_inches="tight"); plt.close(fig); return p


# ==================================================================================
# 3A  — spindling breathes with the heart  (the ~50-second lens — PRIMARY)
# ==================================================================================
def fig_3a():
    fig = plt.figure(figsize=(13.5, 10.5))
    fig.suptitle("TEST 3A  ·  Does spindling rise & fall WITH heart rate every ~50 s?   (PRIMARY — the infraslow lens)",
                 fontsize=13, fontweight="bold", x=0.02, ha="left")
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1], hspace=0.5, wspace=0.22)

    # ① spindle envelope zoomed out over 240 s
    T = 240; t = np.linspace(0, T, 1600); rng = np.random.default_rng(3)
    iso = (1 + np.sin(2*np.pi*0.02*t - 1.0)) / 2
    fast = np.clip(np.sin(2*np.pi*0.22*t), 0, None)
    disc = (0.14 + 0.86*iso) * (0.25 + 0.75*fast**3) + 0.02*rng.standard_normal(t.size)
    envs = 0.14 + 0.86*iso
    a1 = fig.add_subplot(gs[0, 0])
    a1.plot(t, disc, color=SP, lw=0.7, alpha=.45)
    a1.plot(t, envs, color=SP, lw=2.4)
    a1.set_title("① Zoom OUT on the spindle envelope (4 min)", loc="left", fontsize=11)
    a1.set_xlabel("time (s)"); a1.set_ylabel("amount of spindling")
    a1.annotate("", xy=(90, 1.06), xytext=(40, 1.06), arrowprops=dict(arrowstyle="<->", color=SP))
    a1.text(65, 1.1, "~50 s / cycle", color=SP, ha="center", fontsize=9)

    # ② EKG -> RR
    a2 = fig.add_subplot(gs[0, 1])
    dur = 10; beats = [0.3]
    while beats[-1] < dur:
        beats.append(beats[-1] + 0.78 + 0.10*np.sin(2*np.pi*0.9*beats[-1]))
    beats = np.array(beats[:-1])
    for bt in beats:
        a2.plot([bt-0.03, bt, bt+0.02, bt+0.05], [0, 1, -0.25, 0], color=HR, lw=1.4)
    a2.set_ylim(-0.6, 1.4); a2.set_xlim(0, dur)
    a2.set_title("② Turn the EKG into heart rate", loc="left", fontsize=11)
    a2.set_xlabel("time (s)"); a2.set_yticks([])
    a2.text(0.2, 1.28, "each R-peak ↑ ; gap between beats = R–R interval", color=HR, fontsize=8.5)

    # ③ overlay
    a3 = fig.add_subplot(gs[1, :])
    env_n = (1 + np.sin(2*np.pi*0.02*t - 1.0)) / 2
    rr_n  = (1 + np.sin(2*np.pi*0.02*t - 1.6)) / 2
    a3.plot(t, env_n, color=SP, lw=2.6, label="spindle envelope")
    a3.plot(t, rr_n, color=HR, lw=2.6, label="heart-rate rhythm (R–R)")
    a3.axvline(37, color=SP, ls=":", lw=1); a3.axvline(41.7, color=HR, ls=":", lw=1)
    a3.annotate("fixed\nlag", xy=(39.3, 1.02), ha="center", fontsize=8.5, color="#555")
    a3.set_title("③ Overlay the two slow lines — they march together, offset by a fixed lag  (the whole point of 3A)",
                 loc="left", fontsize=11)
    a3.set_xlim(0, T); a3.set_xlabel("time (s)"); a3.set_yticks([]); a3.legend(fontsize=9, loc="upper right")

    # ④ coherence (its own full row — no overlap)
    a4 = fig.add_subplot(gs[2, :])
    f = np.linspace(0, 0.06, 300)
    C = np.clip(0.12 + 0.66*np.exp(-((f-0.02)/0.0035)**2) + 0.04*np.random.default_rng(7).standard_normal(f.size), 0, 0.95)
    a4.fill_between(f, 0, 0.34, color=NULL, alpha=.3, label="surrogate 95% (chance)")
    a4.plot(f, C, color=SP, lw=2.3, label="measured coherence")
    a4.axvline(0.02, color=SP, ls=":", lw=1)
    a4.axvspan(0.01, 0.04, color=SP, alpha=.05)
    a4.set_xlim(0, 0.06); a4.set_ylim(0, 1)
    a4.set_title("④ Measure it: infraslow coherence vs. the surrogate null  (shaded = 0.01–0.04 Hz infraslow band)",
                 loc="left", fontsize=11)
    a4.set_xlabel("frequency (Hz)"); a4.set_ylabel("coherence"); a4.legend(fontsize=9, loc="upper right")
    a4.annotate("peak @ 0.02 Hz clears the grey ✓", xy=(0.02, 0.80), xytext=(0.031, 0.88),
                fontsize=9.5, color=SP, arrowprops=dict(arrowstyle="->", color=SP))

    fig.text(0.5, 0.02, "WHAT YOUR EYE SHOULD SEE:  teal & rose in ③ are two slow waves marching together. Panel ④ is the verdict — one bump at 0.02 Hz "
             "standing clear of the grey.   WHY NIGHTS-ONLY: one cycle is 50 s, so a 1-min clip can't hold even one.",
             ha="center", fontsize=9.5, style="italic", color="#444")
    p = os.path.join(OUT, "test_3A_rr_spindle_infraslow.png")
    fig.savefig(p, bbox_inches="tight"); plt.close(fig); return p


# ==================================================================================
# 3B  — slow waves fire at a preferred moment in the heartbeat
# ==================================================================================
def fig_3b():
    fig = plt.figure(figsize=(13.5, 4.4))
    fig.suptitle("TEST 3B  ·  Do slow waves fire at a preferred point in the cardiac cycle?   (the per-beat lens)",
                 fontsize=13, fontweight="bold", x=0.02, ha="left")
    gs = fig.add_gridspec(1, 3, wspace=0.28, width_ratios=[1, 1, 1])

    sf = 400.0; dur = 6; t = np.arange(0, dur, 1/sf)
    so = np.sin(2*np.pi*0.9*t - 0.6)
    beats = [0.25]
    while beats[-1] < dur + 1:
        beats.append(beats[-1] + 0.80 + 0.05*np.sin(2*np.pi*0.5*beats[-1]))
    beats = np.array(beats)
    troughs = t[signal.argrelmin(so, order=40)[0]]
    troughs = troughs[so[np.searchsorted(t, troughs)] < -0.6]

    # ① two clocks
    a1 = fig.add_subplot(gs[0, 0])
    a1.plot(t, so, color=SO, lw=2.2)
    a1.plot(troughs, np.interp(troughs, t, so), "v", color=SO, ms=9)
    for bt in beats[beats < dur]:
        a1.plot([bt, bt], [-1.9, -1.55], color=HR, lw=1.6)
    a1.set_ylim(-2.1, 1.4); a1.set_title("① Line up the two clocks", loc="left", fontsize=11)
    a1.set_xlabel("time (s)"); a1.set_yticks([])
    a1.text(0.05, 1.15, "slow wave  ▼ = trough", color=SO, fontsize=8.5)
    a1.text(0.05, -2.05, "EKG R-peaks", color=HR, fontsize=8.5)

    # ② cardiac phase ramp
    a2 = fig.add_subplot(gs[0, 1])
    ph = np.zeros_like(t)
    for i in range(len(beats)-1):
        seg = (t >= beats[i]) & (t < beats[i+1])
        ph[seg] = (t[seg]-beats[i])/(beats[i+1]-beats[i])
    a2.plot(t, ph, color=HR, lw=1.8)
    tr_ph = []
    for tr in troughs:
        j = np.searchsorted(beats, tr)-1
        if 0 <= j < len(beats)-1:
            pv = (tr-beats[j])/(beats[j+1]-beats[j]); tr_ph.append(pv)
            a2.plot([tr, tr], [0, pv], color=SO, ls=":", lw=1)
            a2.plot(tr, pv, "o", color=SO, ms=8, mec="w")
    a2.set_ylim(-0.05, 1.05); a2.set_title("② Read cardiac phase at each trough", loc="left", fontsize=11)
    a2.set_xlabel("time (s)"); a2.set_yticks([0, 1]); a2.set_yticklabels(["0", "2π"])
    a2.text(0.1, 0.92, "each ▼ → one angle (dots)", color=SO, fontsize=8.5)

    # ③ polar rayleigh: real vs shuffled
    a3 = fig.add_subplot(gs[0, 2], projection="polar")
    rng = np.random.default_rng(4)
    real = 1.0 + rng.vonmises(0, 3.0, 160)*0 + rng.normal(0, 0.35, 160) + 1.0  # tight cluster ~2.0 rad
    real = rng.vonmises(1.0, 4.0, 160)
    shuf = rng.uniform(-np.pi, np.pi, 160)
    for ang, col, rad in [(shuf, NULL, 0.9), (real, HR, 1.0)]:
        nb = 20; edges = np.linspace(-np.pi, np.pi, nb+1)
        h, _ = np.histogram(ang, bins=edges); h = h/h.max()*rad
        centers = (edges[:-1]+edges[1:])/2
        a3.bar(centers, h, width=2*np.pi/nb, color=col, alpha=.45, edgecolor="w", linewidth=.4)
        R = np.abs(np.mean(np.exp(1j*ang)))
        mu = np.angle(np.mean(np.exp(1j*ang)))
        a3.annotate("", xy=(mu, R*rad), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color=col, lw=3))
    a3.set_title("③ Do the angles clump? (Rayleigh)", fontsize=10.5)
    a3.set_yticklabels([])
    a3.text(0.5, -0.18, "rose = REAL (clumped, long arrow)  ·  grey = SHUFFLED (even, stub arrow)",
            transform=a3.transAxes, ha="center", fontsize=8.5, color="#555")

    fig.text(0.5, -0.05, "WHAT YOUR EYE SHOULD SEE:  the rose dots pile into one arc and the arrow is long — that's the preferred cardiac phase. "
             "TRAP: EKG can bleed into brain channels — use BIPOLAR channels & check the beat-triggered average is slow-shaped, not a QRS spike.",
             ha="center", fontsize=9.5, style="italic", color="#444")
    fig.tight_layout(rect=[0, 0.03, 1, 0.9])
    p = os.path.join(OUT, "test_3B_so_heartbeat.png")
    fig.savefig(p, bbox_inches="tight"); plt.close(fig); return p


if __name__ == "__main__":
    for fn in (fig_3d, fig_3a, fig_3b):
        print("wrote", fn())
    print("\nAll three figures in:", OUT)
