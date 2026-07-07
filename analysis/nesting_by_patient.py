"""
Per-patient N3 SO-triggered nesting — does every patient show it, or does the grand
average hide per-patient mess?

For each mesiotemporal patient, pool ALL cached N3 atlas clips (204 Hz), detect slow-oscillation
troughs, and average the SO wave and the (z-scored) spindle-band envelope around them. Then:
  A. per_patient_nesting.png  - one small-multiple panel per patient (SO wave + spindle envelope)
  B. so_trough_locked_spindle_heatmap.png - SO-trough-locked heatmaps: rows = patients (top)
     and rows = individual SO events (bottom), columns = time around the SO trough, colour =
     spindle envelope.

At 204 Hz ripples (80-120 Hz) are not resolvable, so this is the SO->spindle nesting; SO->ripple
lives in the 1 kHz cohort (slow_ripple_coupling / nesting_timecourse).

Uses only CACHED atlas clips (no downloads).
"""
import os, glob
import numpy as np, pandas as pd
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import atlas

SF = atlas.SFREQ            # 204 Hz
HALF = int(1.5 * SF)        # +/-1.5 s window
MIN_DIST = int(0.6 * SF)
SO_BAND = (0.5, 1.25); SPINDLE = (11, 16)
OUT = os.path.join(atlas.ROOT, "outputs", "nesting_by_patient")


def bp(x, lo, hi):
    sos = signal.butter(3, [lo / (SF / 2), hi / (SF / 2)], btype="band", output="sos")
    return signal.sosfiltfilt(sos, x)


def notch60(x):
    b, a = signal.iirnotch(atlas.LINE_HZ, 30, SF)
    return signal.filtfilt(b, a, x)


def main():
    os.makedirs(OUT, exist_ok=True)
    meta = atlas.load_metadata(normative_only=True)
    roi = meta.dropna(subset=["roi_group"])
    roi = roi[roi.roi_group.isin(atlas.MESIOTEMPORAL)][["pt", "name", "roi_group"]]
    clips = atlas.clip_index()

    per_pt = {}         # pt -> {so_sum, sp_sum, n, events(list of spindle-env epochs, subsampled)}
    for pt in sorted(roi.pt.unique()):
        ch_map = roi[roi.pt == pt].set_index("name")["roi_group"].to_dict()
        n3 = sorted(clips[(clips.pt == pt) & (clips.state == "N3")].idx.tolist())
        cached = [j for j in n3 if os.path.exists(os.path.join(atlas.EEGDIR, f"{pt}_N3_{j}.edf"))]
        if not cached:
            continue
        so_sum = np.zeros(2 * HALF); sp_sum = np.zeros(2 * HALF); n = 0; ev = []
        for j in cached:
            try:
                data, names, sf = atlas.read_clip(pt, "N3", j)
            except Exception:
                continue
            nidx = {nm: k for k, nm in enumerate(names)}
            for chan in ch_map:
                if chan not in nidx:
                    continue
                x = data[nidx[chan]]
                if not np.isfinite(x).all() or np.std(x) < 1e-9:
                    continue
                x = notch60(signal.detrend(x))
                so = bp(x, *SO_BAND)
                spenv = np.abs(signal.hilbert(bp(x, *SPINDLE)))
                spz = (spenv - spenv.mean()) / (spenv.std() + 1e-12)
                idx, _ = signal.find_peaks(-so, height=so.std(), distance=MIN_DIST)
                idx = idx[(idx > HALF) & (idx < len(x) - HALF)]
                for e in idx:
                    sl = slice(e - HALF, e + HALF)
                    so_sum += so[sl]; sp_sum += spz[sl]; n += 1
                    if len(ev) < 4000:
                        ev.append(spz[sl])
        if n >= 20:
            per_pt[pt] = {"so": so_sum / n, "sp": sp_sum / n, "n": n, "ev": ev}
        print(f"  {pt}: {len(cached)} N3 clips, {n} SO events", flush=True)

    pts = list(per_pt)
    t = np.arange(-HALF, HALF) / SF
    print(f"[nesting-by-patient] {len(pts)} patients with >=20 events")

    # per-patient quantification: spindle-env in up-state window (post-trough 0..0.75 s) vs full mean
    win = (t >= 0.0) & (t <= 0.75)
    rows = []
    for pt in pts:
        up = per_pt[pt]["sp"][win].mean()
        rows.append({"pt": pt, "n_events": per_pt[pt]["n"], "spindle_z_upstate": float(up)})
    q = pd.DataFrame(rows).sort_values("spindle_z_upstate", ascending=False)
    q.to_csv(os.path.join(OUT, "per_patient_upstate.csv"), index=False)
    pos = (q.spindle_z_upstate > 0).sum()
    print(f"[nesting-by-patient] up-state spindle z > 0 in {pos}/{len(q)} patients "
          f"(median {q.spindle_z_upstate.median():.3f})")

    # ---- A. per-patient small multiples ----
    order = q.pt.tolist()
    ncol = 5; nrow = int(np.ceil(len(order) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3 * ncol, 2.2 * nrow), squeeze=False)
    for i, pt in enumerate(order):
        ax = axes[i // ncol][i % ncol]; axb = ax.twinx()
        ax.plot(t, per_pt[pt]["so"], color="tab:blue", lw=1)
        axb.plot(t, per_pt[pt]["sp"], color="tab:orange", lw=1)
        ax.axvline(0, color="k", ls=":", lw=0.6); axb.axhline(0, color="0.7", lw=0.5)
        ax.set_title(f"{pt}  (n={per_pt[pt]['n']})", fontsize=7)
        ax.set_yticks([]); axb.set_yticks([]); ax.set_xticks([-1, 0, 1])
        ax.tick_params(labelsize=6)
    for i in range(len(order), nrow * ncol):
        axes[i // ncol][i % ncol].axis("off")
    fig.suptitle("Per-patient N3 SO-triggered nesting (blue = SO wave, orange = spindle envelope z). "
                 "t=0 is the SO trough.", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"per_patient_nesting.{ext}"), dpi=150)

    # ---- B. SO-trough-locked heatmaps ----
    fig, ax = plt.subplots(2, 1, figsize=(9, 9), gridspec_kw={"height_ratios": [1, 1.4]})
    # rows = patients (sorted by up-state spindle), colour = spindle env z
    M = np.vstack([per_pt[pt]["sp"] for pt in order])
    im0 = ax[0].imshow(M, aspect="auto", cmap="RdBu_r", vmin=-0.6, vmax=0.6,
                       extent=[t[0], t[-1], len(order), 0])
    ax[0].axvline(0, color="k", ls=":", lw=0.8)
    ax[0].set_yticks(np.arange(len(order)) + 0.5); ax[0].set_yticklabels(order, fontsize=5)
    ax[0].set_ylabel("patient (sorted)"); ax[0].set_title(
        "SO-trough-locked spindle envelope: rows = patients, time 0 = 0.5-1.25 Hz SO trough")
    fig.colorbar(im0, ax=ax[0], label="spindle z")
    # rows = individual SO events (pooled, subsampled), sorted by up-state spindle
    allev = np.vstack([e for pt in order for e in per_pt[pt]["ev"]])
    rng = np.random.default_rng(0)
    if len(allev) > 3000:
        allev = allev[rng.choice(len(allev), 3000, replace=False)]
    allev = allev[np.argsort(allev[:, win].mean(1))[::-1]]
    im1 = ax[1].imshow(allev, aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5,
                       extent=[t[0], t[-1], len(allev), 0])
    ax[1].axvline(0, color="k", ls=":", lw=0.8)
    ax[1].set_ylabel(f"individual SO events (n={len(allev)}, sorted)")
    ax[1].set_xlabel("time from SO trough (s)")
    ax[1].set_title("rows = individual SO events")
    fig.colorbar(im1, ax=ax[1], label="spindle z")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"so_trough_locked_spindle_heatmap.{ext}"), dpi=150)
    print(f"[nesting-by-patient] figures in {OUT}")


if __name__ == "__main__":
    main()
