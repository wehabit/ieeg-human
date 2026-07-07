"""
Patient-level N3 temporal QC for slow-fast nesting.

This is the skeptical-reader companion to the group SO-triggered nesting figure.
For every patient with normative mesiotemporal N3 clips in the Pattnaik/Litt
atlas, use all available N3 clips and all available normative MTL channels to
show:

  1. a raw example segment with the slow oscillation and spindle envelope,
  2. the patient's SO-trough-locked average, and
  3. an SO-trough-locked event heatmap of slow oscillation cycles.

This is a temporal QC visualization, not a new inferential test. The locking
rhythm is the detected N3 slow oscillation: 0.5-1.25 Hz, i.e. cycles of about
0.8-2.0 s.

Important limitation: the atlas provides 30 s N3 clips, not continuous full-night
N3. "Full duration" here means all available N3 clips for that patient. The
mixed-rate atlas clips cannot cleanly resolve ripples across patients, so these
figures show SO->spindle nesting. Ripple nesting remains the 1 kHz cohort's job.
"""
from __future__ import annotations

import argparse
import math
import os

import numpy as np
import pandas as pd
from scipy import signal

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import atlas


OUT = os.path.join(atlas.ROOT, "outputs", "patient_n3_temporal_qc")
SO_BAND = (0.5, 1.25)
SPINDLE_BAND = (11, 16)
HALF_S = 1.5
MIN_DIST_S = 0.6
MAX_HEATMAP_EVENTS = 220
TARGET_SF = 200.0


def bandpass(x: np.ndarray, sf: float, lo: float, hi: float) -> np.ndarray:
    b, a = signal.butter(3, [lo / (sf / 2), hi / (sf / 2)], btype="band")
    return signal.filtfilt(b, a, x)


def clean_trace(x: np.ndarray, sf: float) -> np.ndarray | None:
    if not np.isfinite(x).all() or np.nanstd(x) < 1e-12:
        return None
    x = signal.detrend(x)
    # N3 atlas data are HUP/Penn, i.e. 60 Hz line noise.
    b, a = signal.iirnotch(atlas.LINE_HZ, 30, sf)
    return signal.filtfilt(b, a, x)


def robust_z(x: np.ndarray) -> np.ndarray:
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    scale = 1.4826 * mad if mad > 0 else np.std(x)
    return (x - med) / (scale + 1e-12)


def so_troughs(so: np.ndarray, sf: float, half: int) -> np.ndarray:
    min_dist = int(MIN_DIST_S * sf)
    # Require a clear negative slow-wave trough. Using 0.75 SD keeps enough
    # events in short atlas clips while avoiding tiny wiggles.
    idx, _ = signal.find_peaks(-so, height=0.75 * np.std(so), distance=min_dist)
    return idx[(idx > half) & (idx < len(so) - half)]


def load_patient_n3(pt: str, ch_map: dict[str, str], clip_idxs: list[int], download_missing: bool):
    """Yield cleaned traces for all requested N3 clips/channels."""
    paths = [f"files/processed_eeg_final/{pt}_N3_{j}.edf" for j in clip_idxs]
    if download_missing:
        atlas.download(paths)
    for j in clip_idxs:
        local = os.path.join(atlas.EEGDIR, f"{pt}_N3_{j}.edf")
        if not os.path.exists(local):
            continue
        try:
            data, names, sf = atlas.read_clip(pt, "N3", j)
        except Exception:
            continue
        nidx = {n: k for k, n in enumerate(names)}
        for chan, roi in ch_map.items():
            if chan not in nidx:
                continue
            x = clean_trace(data[nidx[chan]], sf)
            if x is None:
                continue
            yield {"pt": pt, "clip": j, "channel": chan, "roi_group": roi, "sf": sf, "x": x}


def analyze_patient(pt: str, ch_map: dict[str, str], clip_idxs: list[int], download_missing: bool):
    t = np.arange(-HALF_S, HALF_S, 1 / TARGET_SF)
    so_sum = np.zeros(len(t))
    sp_sum = np.zeros(len(t))
    raw_sum = np.zeros(len(t))
    heat_rows = []
    n_events = 0
    n_traces = 0
    total_s = 0.0
    example = None
    best_example_score = -np.inf

    for rec in load_patient_n3(pt, ch_map, clip_idxs, download_missing):
        x = rec["x"]
        sf = rec["sf"]
        half = int(HALF_S * sf)
        n_traces += 1
        total_s += len(x) / sf
        so = bandpass(x, sf, *SO_BAND)
        sp = np.abs(signal.hilbert(bandpass(x, sf, *SPINDLE_BAND)))
        sp_z = (sp - sp.mean()) / (sp.std() + 1e-12)
        ev = so_troughs(so, sf, half)
        if len(ev) == 0:
            continue

        # Pick a raw segment that has a big slow oscillation and visible spindle
        # envelope. This is only for the example panel.
        score = float(np.percentile(np.abs(so), 95) + 0.2 * np.percentile(sp_z, 95))
        if score > best_example_score:
            center = ev[len(ev) // 2]
            lo = max(0, center - int(4 * sf))
            hi = min(len(x), center + int(4 * sf))
            example = {**rec, "so": so, "sp_z": sp_z, "lo": lo, "hi": hi}
            best_example_score = score

        for e in ev:
            sl = slice(e - half, e + half)
            # Atlas clips appear at multiple sampling rates (commonly 204 or
            # 256 Hz). Interpolate each event window onto one shared time axis
            # so every patient can contribute to the same heatmap/average.
            local_t = (np.arange(sl.stop - sl.start) - half) / sf
            so_seg = np.interp(t, local_t, so[sl])
            sp_seg = np.interp(t, local_t, sp_z[sl])
            raw_seg = np.interp(t, local_t, x[sl])
            so_sum += so_seg
            sp_sum += sp_seg
            raw_sum += raw_seg
            n_events += 1
            if len(heat_rows) < MAX_HEATMAP_EVENTS:
                heat_rows.append(robust_z(so_seg))
            else:
                # Reservoir sample deterministically enough for stable displays.
                k = n_events % MAX_HEATMAP_EVENTS
                if k < MAX_HEATMAP_EVENTS // 20:
                    heat_rows[k] = robust_z(so_seg)

    if n_events == 0:
        return None

    return {
        "pt": pt,
        "t": t,
        "avg_so_uv": so_sum / n_events * 1e6,
        "avg_raw_uv": raw_sum / n_events * 1e6,
        "avg_sp_z": sp_sum / n_events,
        "heat": np.asarray(heat_rows),
        "n_events": n_events,
        "n_traces": n_traces,
        "n_channels": len(ch_map),
        "n_clips": len(clip_idxs),
        "duration_s": total_s,
        "example": example,
    }


def save_patient_figure(res: dict, out_dir: str):
    pt = res["pt"]
    t = res["t"]
    fig = plt.figure(figsize=(12, 9))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.15, 1.05, 1.25], hspace=0.35)

    ax0 = fig.add_subplot(gs[0])
    ex = res["example"]
    if ex is not None:
        lo, hi = ex["lo"], ex["hi"]
        tt = (np.arange(hi - lo) - (hi - lo) / 2) / ex["sf"]
        x_uv = ex["x"][lo:hi] * 1e6
        so_uv = ex["so"][lo:hi] * 1e6
        sp = ex["sp_z"][lo:hi]
        sp_scaled = sp * (0.15 * np.nanstd(x_uv)) + np.nanmedian(x_uv) + 2.5 * np.nanstd(x_uv)
        ax0.plot(tt, x_uv, color="0.65", lw=0.7, label="raw LFP")
        ax0.plot(tt, so_uv, color="tab:blue", lw=1.8, label="SO 0.5-1.25 Hz")
        ax0.plot(tt, sp_scaled, color="tab:orange", lw=1.0, label="spindle env 11-16 Hz (scaled)")
        ax0.set_title(f"{pt}: raw N3 example ({ex['channel']}, {ex['roi_group']})")
        ax0.set_ylabel("uV")
        ax0.legend(fontsize=8, ncol=3, loc="upper right")
    else:
        ax0.text(0.5, 0.5, "No example segment", ha="center", va="center")
    ax0.set_xlabel("time (s)")

    ax1 = fig.add_subplot(gs[1])
    ax1b = ax1.twinx()
    ax1.plot(t, res["avg_so_uv"], color="tab:blue", lw=2.2, label="SO-trough-locked LFP/SO wave")
    ax1b.plot(t, res["avg_sp_z"], color="tab:orange", lw=1.8, label="spindle envelope")
    ax1.axvline(0, color="k", lw=0.8, ls=":")
    ax1.set_ylabel("SO wave (uV)", color="tab:blue")
    ax1b.set_ylabel("spindle envelope (z)", color="tab:orange")
    ax1.set_xlabel("time from SO trough (s)")
    ax1.set_title(
        f"SO-trough-locked average: {res['n_events']} events, "
        f"{res['n_traces']} channel-clips, {res['duration_s']/60:.1f} min N3 clips"
    )
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax1b.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")

    ax2 = fig.add_subplot(gs[2])
    heat = res["heat"]
    if len(heat):
        vmax = np.nanpercentile(np.abs(heat), 98)
        ax2.imshow(
            heat,
            aspect="auto",
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
            extent=[t[0], t[-1], heat.shape[0], 0],
            interpolation="nearest",
        )
        ax2.axvline(0, color="k", lw=0.8, ls=":")
        ax2.set_ylabel("SO events")
        ax2.set_xlabel("time from SO trough (s; SO band 0.5-1.25 Hz)")
        ax2.set_title("SO-trough-locked event heatmap: each row is one detected N3 SO cycle")
    fig.suptitle("Patient-level N3 SO-spindle temporal QC (atlas clips; no ripple claim)", y=0.99)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"{pt}_n3_so_triggered_nesting.png"), dpi=150)
    plt.close(fig)


def save_group_figures(results: list[dict], out_dir: str):
    if not results:
        return
    t = results[0]["t"]

    # Patient heatmap: each row is one patient's SO-trough-locked slow wave.
    mat = []
    labels = []
    for r in results:
        v = r["avg_so_uv"]
        mat.append((v - np.mean(v)) / (np.std(v) + 1e-12))
        labels.append(r["pt"].replace("sub-", ""))
    mat = np.asarray(mat)
    order = np.argsort(np.argmax(mat, axis=1))
    mat = mat[order]
    labels = [labels[i] for i in order]

    fig, ax = plt.subplots(figsize=(10, max(7, 0.15 * len(results) + 2)))
    vmax = np.nanpercentile(np.abs(mat), 98)
    im = ax.imshow(
        mat,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        extent=[t[0], t[-1], len(results), 0],
        interpolation="nearest",
    )
    ax.axvline(0, color="k", lw=0.8, ls=":")
    ax.set_xlabel("time from SO trough (s)")
    ax.set_ylabel("patients")
    ax.set_title(
        "N3 SO-trough-locked slow-oscillation heatmap\n"
        "Rows = patients; time 0 = detected 0.5-1.25 Hz SO trough (~0.8-2.0 s/cycle)"
    )
    step = max(1, math.ceil(len(labels) / 28))
    ax.set_yticks(np.arange(0.5, len(labels), step))
    ax.set_yticklabels(labels[::step], fontsize=7)
    fig.colorbar(im, ax=ax, label="patient-normalized SO-filtered LFP")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "all_patients_so_trough_locked_heatmap.png"), dpi=150)
    plt.close(fig)

    # Contact sheet of all patient averages.
    n = len(results)
    cols = 4
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(14, max(8, rows * 2.1)), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, r in zip(axes, results):
        ax.plot(t, r["avg_so_uv"], color="tab:blue", lw=1.3)
        ax2 = ax.twinx()
        ax2.plot(t, r["avg_sp_z"], color="tab:orange", lw=0.9, alpha=0.8)
        ax.axvline(0, color="k", lw=0.5, ls=":")
        ax.set_title(f"{r['pt'].replace('sub-', '')} · {r['n_events']} SO", fontsize=8)
        ax.tick_params(labelsize=7)
        ax2.tick_params(labelsize=7)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle(
        "Per-patient QC: SO-trough-locked slow wave (blue) and spindle envelope (orange)",
        y=0.995,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "all_patients_so_triggered_contact_sheet.png"), dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-patients", type=int, default=None)
    ap.add_argument("--cached-only", action="store_true",
                    help="Use only already cached EDF clips; do not download missing N3 clips.")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    meta = atlas.load_metadata(normative_only=True)
    roi = meta.dropna(subset=["roi_group"])
    roi = roi[roi.roi_group.isin(atlas.MESIOTEMPORAL)]
    clips = atlas.clip_index()
    n3 = clips[clips.state == "N3"]
    pts = sorted(set(n3.pt.unique()) & set(roi.pt.unique()))
    if args.max_patients:
        pts = pts[:args.max_patients]

    patient_dir = os.path.join(OUT, "patients")
    os.makedirs(patient_dir, exist_ok=True)
    results = []
    rows = []
    for i, pt in enumerate(pts, 1):
        ch_map = roi[roi.pt == pt].set_index("name")["roi_group"].to_dict()
        idxs = sorted(n3[n3.pt == pt].idx.tolist())
        if args.cached_only:
            idxs = [j for j in idxs if os.path.exists(os.path.join(atlas.EEGDIR, f"{pt}_N3_{j}.edf"))]
        if not idxs or not ch_map:
            continue
        res = analyze_patient(pt, ch_map, idxs, download_missing=not args.cached_only)
        if res is None:
            print(f"[{i}/{len(pts)}] {pt}: no SO events")
            continue
        save_patient_figure(res, patient_dir)
        results.append(res)
        rows.append({
            "pt": pt,
            "n_clips": res["n_clips"],
            "n_channels": res["n_channels"],
            "n_channel_clips": res["n_traces"],
            "duration_s": round(res["duration_s"], 1),
            "n_so_events": res["n_events"],
            "patient_png": os.path.join("patients", f"{pt}_n3_so_triggered_nesting.png"),
        })
        print(f"[{i}/{len(pts)}] {pt}: {res['n_events']} SO events, {res['duration_s']/60:.1f} min")

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "patient_n3_temporal_summary.csv"), index=False)
    save_group_figures(results, OUT)
    print(f"[patient_n3_temporal_qc] wrote {len(results)} patient figures to {OUT}")


if __name__ == "__main__":
    main()
