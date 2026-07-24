"""
Bridge - Mouse<->Human bridge figure (the one-pager / Westover-pitch slide).

Thesis: Kipnis Fig 2b shows spike-triggered field waves split by state -- NREM = a slow
(0.5-4 Hz) wave Kipnis LINKS to CSF clearance (not measured here); REM = theta (6-10 Hz). Mouse pilot,
lacking natural NREM, only reached the REM/theta column (theta-organized firing under 50 Hz
drive). The human iEEG hits the NREM/slow column -- the clearance-relevant one. This figure
places the two side by side and states the design implication.

Embeds existing PNGs (no recompute):
  - mouse theta+slow STA waveforms (Buzsaki repo)
  - human N3 comodulogram (result_visuals)
  - human coupling-by-state (Spindle-coupling)
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import atlas

OUT = os.path.join(atlas.ROOT, "outputs", "mouse_human_bridge")
BUZ = os.environ.get("BUZSAKLI_LAB_ROOT", os.path.dirname(atlas.ROOT))
MOUSE_STA = os.path.join(BUZ, "analysis/outputs/dec4/all_dhpc_kipnis_coordination/all_dhpc_kipnis_waveforms_amp250.png")
HUMAN_COMOD = os.path.join(atlas.ROOT, "outputs/result_visuals/comodulogram.png")
HUMAN_STATE = os.path.join(atlas.ROOT, "outputs/slow_spindle_coupling/spindle_coupling_by_state.png")


def show(ax, path, title):
    ax.axis("off")
    if os.path.exists(path):
        ax.imshow(mpimg.imread(path))
    else:
        ax.text(0.5, 0.5, f"[missing]\n{os.path.basename(path)}", ha="center", va="center")
    ax.set_title(title, fontsize=11, fontweight="bold")


def main():
    raise SystemExit(
        "LEGACY/WITHDRAWN bridge figure: it embeds superseded human outputs and must be "
        "regenerated only after a corrected, quality-controlled cohort rerun.")
    os.makedirs(OUT, exist_ok=True)
    fig = plt.figure(figsize=(15, 9))
    gs = GridSpec(3, 2, height_ratios=[0.6, 1.25, 1.0], hspace=0.28, wspace=0.12)

    # top: the Kipnis thesis column-map
    axt = fig.add_subplot(gs[0, :]); axt.axis("off")
    axt.text(0.5, 0.92, "A candidate NEURAL SUBSTRATE (NREM slow-wave coordination) that Kipnis "
             "links to clearance — and the column the human data reaches",
             ha="center", va="top", fontsize=13, fontweight="bold")
    cols = [
        ("WAKE", "incoherent field", "#888888"),
        ("NREM  •  slow 0.5–4 Hz", "large slow wave organizes\nspikes/spindles/ripples\n(Kipnis Fig 2b; clearance hypothesized)", "#1f77b4"),
        ("REM  •  theta 6–10 Hz", "rhythmic theta wave", "#d62728"),
    ]
    for i, (h, sub, c) in enumerate(cols):
        x = 0.17 + i * 0.33
        axt.text(x, 0.55, h, ha="center", fontsize=12, fontweight="bold", color=c)
        axt.text(x, 0.30, sub, ha="center", fontsize=9, color="black")
    axt.annotate("", xy=(0.5, 0.05), xytext=(0.83, 0.05),
                 arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.5))
    axt.annotate("", xy=(0.5, 0.12), xytext=(0.17, 0.12),
                 arrowprops=dict(arrowstyle="->", color="#888888", lw=1.0))
    axt.text(0.5, 0.02, "MOUSE reached the REM/theta column (50 Hz drive, no natural NREM)   |   "
             "HUMAN reaches the NREM/slow column (natural N3).  No clearance is measured here.",
             ha="center", fontsize=9, style="italic")

    show(fig.add_subplot(gs[1, 0]), MOUSE_STA,
         "MOUSE (Buzsaki Dec4): spike-triggered field wave\ntheta present under drive · slow-wave coordination ABSENT")
    show(fig.add_subplot(gs[1, 1]), HUMAN_COMOD,
         "HUMAN (this study): N3 comodulogram\nslow-phase → spindle-amp coupling PRESENT, wake flat")

    # bottom-left: human state dependence
    show(fig.add_subplot(gs[2, 0]), HUMAN_STATE,
         "HUMAN: coupling is NREM-specific (peaks N3, ~0 in REM)")

    # bottom-right: design implication text
    axd = fig.add_subplot(gs[2, 1]); axd.axis("off")
    axd.text(0.02, 0.95, "Design implication", fontsize=12, fontweight="bold", va="top")
    axd.text(0.02, 0.80,
             "• Mouse result = target ENGAGEMENT: 50 Hz vibration reorganizes\n"
             "  hippocampal firing (at theta).\n\n"
             "• Human data = target STATE: the native NREM coordination is a\n"
             "  slow-oscillation -> spindle nesting (Slow-power p=2.9e-6; Spindle-coupling MI_z p=0.030)\n"
             "  and -> ripple nesting on 1 kHz data (Ripple-coupling p=0.0015).\n\n"
             "• Working HYPOTHESIS: to engage this substrate, sleep stimulation\n"
             "  should ENHANCE native slow-wave coordination (slow / Jiang-Xie\n"
             "  frequencies), not 50 Hz. Clearance itself is NOT measured here.\n\n"
             "• All human results from open data — no patient access required.",
             fontsize=9.5, va="top", family="monospace")

    fig.suptitle("Bridge — Mouse → Human bridge: from theta engagement to the NREM slow-wave "
                 "coordination state (candidate substrate; clearance not measured)",
                 fontsize=12, y=0.995)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"mouse_human_bridge.{ext}"), dpi=150, bbox_inches="tight")
    print(f"[Bridge] wrote bridge figure to {OUT}")


if __name__ == "__main__":
    raise SystemExit(
        "LEGACY/WITHDRAWN bridge figure: it embeds superseded human outputs and must be "
        "regenerated only after a corrected, quality-controlled cohort rerun.")
