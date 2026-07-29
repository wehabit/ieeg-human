"""Create novice-friendly, data-grounded examples for questions 3A, 3B, and 3D.

These figures are recomputed from the current neutral caches under the locked
``overlap11_endpoint_local`` profile.  They are not synthetic illustrations and
they are not new inferential analyses.

Selected examples:

* 3A: RESP0699 -- spectrum, fixed-0.02-Hz coherence, and cross-correlation
  are all available.
* 3B: HUP160 -- N2-like, N3-like, and pooled-NREM descriptive curves are all
  available.
* 3D: HUP172 -- pooled, N2-like, and N3-like descriptive vectors are all
  available.

Run from the repository root:

    .venv/bin/python visualization/make_participant_result_visuals.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
sys.path.insert(0, str(ANALYSIS))

import lecci_faithful_3A as lecci  # noqa: E402
import run_qc_grid as grid  # noqa: E402
from event_3B_cached import stable_stage_epoch_indices, stage_so_times  # noqa: E402
from event_3b_estimators import (  # noqa: E402
    FS_RR,
    rr_baseline_hr,
    subject_so_triggered,
)
from event_3d_estimators import (  # noqa: E402
    channel_night_events,
    stage_event_pairs,
)
from event_3d_cache_support import (  # noqa: E402
    _validated_neutral_payload,
    analyse_3d_cache_support,
)
from materialize_qc_cache import materialize  # noqa: E402
from qc_profiles import load_qc_profile  # noqa: E402
from spectral_gapped import (  # noqa: E402
    _prepare,
    analytic_msc_threshold,
    coherence_gapped,
    welch_segments,
)
from participant_visual_provenance import (  # noqa: E402
    VISUAL_RESULT_NAMES,
    assert_visual_inputs_unchanged,
    freeze_visual_inputs,
    repository_relative,
    validate_terminal_visual_manifest,
    visual_manifest_base,
    write_failed_manifest,
    write_in_progress_manifest,
    write_preflight_manifest,
    write_terminal_manifest,
)
from participant_visual_narratives import (  # noqa: E402
    render_resp0699_check,
    render_visual_readme,
)
from pipeline_version import atomic_json_dump  # noqa: E402


OUT = ROOT / "outputs" / "participant_result_visuals"
PROFILE_ID = "overlap11_endpoint_local"
LOCKED_QC_ARTIFACT = (
    ROOT
    / "outputs"
    / "qc_grid_public"
    / "locked"
    / "overlap11_endpoint_local__hup.json"
)

INK = "#172033"
MUTED = "#687085"
GRID = "#d9dde7"
SIGMA = "#008c95"
HR = "#cc3f62"
SO = "#d7831f"
PURPLE = "#7557b7"
GREEN = "#268b62"
FAIL = "#c7cbd4"

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.edgecolor": GRID,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "axes.linewidth": 0.9,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def _final_subject(grid_path: Path, subject: str) -> dict:
    with grid_path.open() as handle:
        artifact = json.load(handle)
    profile = next(
        item for item in artifact["profiles"] if item["qc_profile_id"] == PROFILE_ID
    )
    return next(item for item in profile["subjects"] if item["subject"] == subject)


def _close(actual, expected, label: str, atol: float = 1e-10) -> None:
    if not np.isclose(float(actual), float(expected), rtol=1e-9, atol=atol):
        raise RuntimeError(f"{label}: recomputed {actual!r} != final {expected!r}")


def _save(fig: plt.Figure, stem: str) -> tuple[Path, Path]:
    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / f"{stem}.png"
    svg = OUT / f"{stem}.svg"
    fig.savefig(png, dpi=190, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    return png, svg


def _smooth_sigma_4s(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    width = int(round(4 * lecci.FS))
    finite = np.isfinite(values)
    numerator = np.convolve(np.where(finite, values, 0.0), np.ones(width), mode="same")
    denominator = np.convolve(finite.astype(float), np.ones(width), mode="same")
    return np.where(
        denominator >= width / 2,
        numerator / np.maximum(denominator, 1.0),
        np.nan,
    )


def _z(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    return (values - np.nanmean(values)) / np.nanstd(values)


def make_3a(profile: dict) -> dict:
    subject = "sub-RESP0699"
    cache_path = ROOT / "data" / "derived" / "ds003848" / f"{subject}.npz"
    artifact_path = (
        ROOT
        / "outputs"
        / "qc_grid"
        / "staging_window_support_v1"
        / "respect_qc_grid.json"
    )
    final = _final_subject(artifact_path, subject)["result_3a"]

    with np.load(cache_path, allow_pickle=False) as cache:
        mat = materialize(cache, profile)
        recomputed = grid.analyse_3a(mat, profile)
        sig = np.asarray(mat["sigma_parietal"], float)
        hr = np.asarray(mat["hr_1"], float)
        labels = np.asarray(mat["stage_lab"]).astype(str)

    if recomputed["endpoint_availability"] != final["endpoint_availability"]:
        raise RuntimeError("3A availability differs from the final artifact")
    _close(recomputed["peak"]["peak_hz"], final["peak"]["peak_hz"], "3A peak")
    _close(
        recomputed["coherence"]["at_0p02_hz"],
        final["coherence"]["at_0p02_hz"],
        "3A coherence",
    )
    _close(
        recomputed["cross_correlation"]["lecci_direction_peak_r"],
        final["cross_correlation"]["lecci_direction_peak_r"],
        "3A correlation",
    )

    nrem, _ = lecci.core_study_nrem_mask(labels)
    freqs, spectrum, n_bouts, bout_seconds = lecci.subject_spectrum(sig, nrem)
    if spectrum is None:
        raise RuntimeError("RESP0699 final 3A spectrum unexpectedly unavailable")

    second_mask = np.zeros(len(sig), bool)
    for start, stop in lecci.nrem_bouts(nrem):
        second_mask[start:stop] = True
    coherence = coherence_gapped(
        np.where(second_mask, hr, np.nan),
        np.where(second_mask, sig, np.nan),
        fs=lecci.FS,
        nperseg=lecci.NPERSEG,
        highpass=0.005,
    )
    if coherence is None:
        raise RuntimeError("RESP0699 final coherence unexpectedly unavailable")
    coherence_threshold = analytic_msc_threshold(coherence["K"], lecci.ALPHA)

    xcorr = lecci.cross_correlation(
        sig,
        hr,
        nrem,
        minimum_windows=profile["endpoint_3a"]["minimum_cross_correlation_windows"],
    )
    if xcorr is None:
        raise RuntimeError("RESP0699 final cross-correlation unexpectedly unavailable")

    valid_trace = second_mask & np.isfinite(sig) & np.isfinite(hr)
    runs = [
        (start, stop)
        for start, stop in lecci.contiguous_runs(valid_trace)
        if stop - start >= 120
    ]
    if not runs:
        raise RuntimeError("RESP0699 has no finite 120-s trace for teaching panel")
    trace_start, trace_stop = max(runs, key=lambda value: value[1] - value[0])
    trace_stop = min(trace_stop, trace_start + 300)
    trace_t = np.arange(trace_stop - trace_start, dtype=float)
    sigma_trace = _z(_smooth_sigma_4s(sig)[trace_start:trace_stop])
    hr_trace = _z(hr[trace_start:trace_stop])

    fig = plt.figure(figsize=(14, 9.2))
    gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.28)
    fig.suptitle(
        "Question 3A — one real participant with all three measurements available",
        x=0.06,
        y=0.99,
        ha="left",
        fontsize=16,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.06,
        0.955,
        "RESP0699 · current v9 locked endpoint-local profile · availability pass ≠ proof of the paper prediction",
        ha="left",
        fontsize=10.5,
        color=MUTED,
    )

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(trace_t, sigma_trace, color=SIGMA, lw=1.8, label="sigma power (10–15 Hz)")
    ax.plot(trace_t, hr_trace, color=HR, lw=1.35, alpha=0.88, label="heart rate")
    ax.axhline(0, color=GRID, lw=0.8)
    ax.set_title("A. What the two slow time series look like", loc="left", fontweight="bold")
    ax.set_xlabel("seconds in an actual NREM excerpt")
    ax.set_ylabel("standardized value (z score)")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.text(
        0.02,
        0.03,
        "Illustrative excerpt only; the tests use all qualifying NREM data.",
        transform=ax.transAxes,
        color=MUTED,
        fontsize=8.5,
    )

    ax = fig.add_subplot(gs[0, 1])
    show = (freqs >= 0.001) & (freqs <= 0.06)
    ax.plot(freqs[show], spectrum[show], color=SIGMA, lw=2.2)
    ax.axvline(0.019, color=PURPLE, ls="--", lw=1.5, label="Lecci ≈0.019 Hz")
    ax.axvline(
        final["peak"]["peak_hz"],
        color=SIGMA,
        ls=":",
        lw=2,
        label=f"participant peak {final['peak']['peak_hz']:.4f} Hz",
    )
    ax.scatter(
        [final["peak"]["peak_hz"]],
        [np.interp(final["peak"]["peak_hz"], freqs, spectrum)],
        s=45,
        color=SIGMA,
        zorder=4,
    )
    ax.set_title("B. Test 1: is there an infraslow sigma peak?", loc="left", fontweight="bold")
    ax.set_xlabel("infraslow frequency (Hz)")
    ax.set_ylabel("normalized sigma spectral power")
    ax.set_xlim(0, 0.06)
    ax.legend(frameon=False, fontsize=8.7)
    ax.text(
        0.98,
        0.04,
        f"PASS: accepted peak\n{n_bouts} bouts · {bout_seconds:.0f} s total",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        color=GREEN,
        fontsize=9,
        fontweight="bold",
    )

    ax = fig.add_subplot(gs[1, 0])
    show = (coherence["f"] >= 0.0) & (coherence["f"] <= 0.06)
    ax.plot(coherence["f"][show], coherence["cxy"][show], color=PURPLE, lw=2)
    ax.axhline(
        coherence_threshold,
        color=MUTED,
        ls="--",
        lw=1.3,
        label=f"individual analytic threshold {coherence_threshold:.3f}",
    )
    ax.axvline(0.02, color=SIGMA, ls=":", lw=1.6)
    target = final["coherence"]["at_0p02_hz"]
    ax.scatter([0.02], [target], color=SIGMA, s=55, zorder=4)
    ax.annotate(
        f"0.02 Hz = {target:.3f}",
        (0.02, target),
        xytext=(0.030, min(0.94, target + 0.16)),
        arrowprops={"arrowstyle": "->", "color": SIGMA},
        color=SIGMA,
        fontsize=9,
    )
    ax.set_title(
        "C. Test 2: do sigma and heart rate share 0.02-Hz power?",
        loc="left",
        fontweight="bold",
    )
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("magnitude-squared coherence")
    ax.set_xlim(0, 0.06)
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    ax.text(
        0.02,
        0.04,
        "PASS: above this participant's threshold",
        transform=ax.transAxes,
        color=GREEN,
        fontsize=9,
        fontweight="bold",
    )

    ax = fig.add_subplot(gs[1, 1])
    lag = np.asarray(xcorr["lag_s"], float)
    values = np.asarray(xcorr["xcorr"], float)
    ax.plot(lag, values, color=HR, lw=2.2)
    ax.axhline(0, color=GRID, lw=0.8)
    ax.axvline(0, color=MUTED, lw=0.9)
    ax.axvline(5, color=PURPLE, ls="--", lw=1.5, label="paper expectation ≈+5 s")
    best_lag = final["cross_correlation"]["lecci_direction_peak_lag_s"]
    best_r = final["cross_correlation"]["lecci_direction_peak_r"]
    ax.scatter([best_lag], [best_r], color=HR, s=55, zorder=4)
    ax.annotate(
        f"best paper-direction point\nlag {best_lag:.0f} s, r={best_r:.3f}",
        (best_lag, best_r),
        xytext=(5.8, best_r + 0.12),
        arrowprops={"arrowstyle": "->", "color": HR},
        color=HR,
        fontsize=9,
    )
    ax.set_title(
        "D. Test 3: does heart rate lead sigma by ~5 seconds?",
        loc="left",
        fontweight="bold",
    )
    ax.set_xlabel("lag (s); positive = sigma follows heart rate")
    ax.set_ylabel("correlation r")
    ax.set_xlim(-15, 15)
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")
    ax.text(
        0.98,
        0.04,
        "AVAILABLE and positive,\nbut best lag is 0 s — not ~5 s",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        color=HR,
        fontsize=9,
        fontweight="bold",
    )

    fig.text(
        0.06,
        0.012,
        "Visual verdict: RESP0699 has a peak and above-threshold coherence, but its timing does not reproduce the expected ~5-s lead. "
        "This is candidate physiology, not an LC measurement.",
        color=INK,
        fontsize=10,
        fontweight="bold",
    )
    png, svg = _save(fig, "3A_RESP0699_all_three_measurements")
    return {
        "subject": subject,
        "png": repository_relative(png),
        "svg": repository_relative(svg),
        "peak_hz": final["peak"]["peak_hz"],
        "coherence_0p02": target,
        "coherence_threshold": coherence_threshold,
        "xcorr_r": best_r,
        "xcorr_lag_s": best_lag,
    }


def make_3a_threshold_explainer(profile: dict) -> dict:
    """Show exactly how RESP0699's nominal pointwise coherence cutoff was obtained."""
    subject = "sub-RESP0699"
    cache_path = ROOT / "data" / "derived" / "ds003848" / f"{subject}.npz"
    artifact_path = (
        ROOT
        / "outputs"
        / "qc_grid"
        / "staging_window_support_v1"
        / "respect_qc_grid.json"
    )
    final = _final_subject(artifact_path, subject)["result_3a"]

    with np.load(cache_path, allow_pickle=False) as cache:
        mat = materialize(cache, profile)
        sig = np.asarray(mat["sigma_parietal"], float)
        hr = np.asarray(mat["hr_1"], float)
        labels = np.asarray(mat["stage_lab"]).astype(str)

    nrem, _ = lecci.core_study_nrem_mask(labels)
    bouts = lecci.nrem_bouts(nrem)
    second_mask = np.zeros(len(sig), bool)
    for start, stop in bouts:
        second_mask[start:stop] = True
    x = np.where(second_mask, hr, np.nan)
    y = np.where(second_mask, sig, np.nan)
    coherence = coherence_gapped(
        x,
        y,
        fs=lecci.FS,
        nperseg=lecci.NPERSEG,
        highpass=0.005,
    )
    if coherence is None:
        raise RuntimeError("RESP0699 final coherence unexpectedly unavailable")

    _, _, valid, filled = _prepare(
        x,
        y,
        lecci.FS,
        max_gap_s=5.0,
        highpass=0.005,
    )
    noverlap = lecci.NPERSEG // 2
    starts = [
        start
        for start in welch_segments(valid, lecci.NPERSEG, noverlap)
        if filled[start : start + lecci.NPERSEG].mean() <= 0.25
    ]
    if starts != [990, 1118, 1246, 1650]:
        raise RuntimeError(f"unexpected RESP0699 Welch starts: {starts}")

    k = int(coherence["K"])
    alpha = float(lecci.ALPHA)
    threshold = float(analytic_msc_threshold(k, alpha))
    observed = float(final["coherence"]["at_0p02_hz"])
    index = int(np.argmin(np.abs(coherence["f"] - lecci.F_LECCI)))
    evaluated_hz = float(coherence["f"][index])
    nominal_p = float((1.0 - observed) ** (k - 1))
    bonferroni_alpha = alpha / 3.0
    bonferroni_threshold = float(analytic_msc_threshold(k, bonferroni_alpha))

    _close(threshold, final["coherence"]["analytic_threshold"], "3A threshold")
    _close(observed, coherence["cxy"][index], "3A observed coherence")

    fig = plt.figure(figsize=(15, 7.4))
    gs = fig.add_gridspec(
        1,
        3,
        left=0.055,
        right=0.98,
        bottom=0.15,
        top=0.79,
        wspace=0.31,
        width_ratios=[1.35, 1.0, 0.92],
    )
    fig.suptitle(
        "How RESP0699's 3A coherence threshold was calculated",
        x=0.055,
        y=0.975,
        ha="left",
        fontsize=16,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.055,
        0.925,
        "The 0.632 line is a nominal pointwise 5% analytic cutoff determined by four accepted spectral windows.",
        ha="left",
        fontsize=10.5,
        color=MUTED,
    )

    ax = fig.add_subplot(gs[0, 0])
    for start, stop in bouts:
        ax.broken_barh(
            [(start, stop - start)],
            (-0.34, 0.68),
            facecolors="#d7e7e7",
            edgecolors=SIGMA,
            linewidth=0.8,
        )
    window_colors = [SIGMA, SO, PURPLE, HR]
    for row, (start, color) in enumerate(zip(starts, window_colors), start=1):
        ax.broken_barh(
            [(start, lecci.NPERSEG)],
            (row - 0.34, 0.68),
            facecolors=color,
            edgecolors=color,
            alpha=0.86,
        )
        ax.text(
            start + lecci.NPERSEG / 2,
            row,
            f"{start} s",
            ha="center",
            va="center",
            fontsize=7,
            color="white",
            fontweight="bold",
        )
    ax.set_yticks(range(5))
    ax.set_yticklabels(["all qualifying\nNREM bouts", "window 1", "window 2", "window 3", "window 4"])
    ax.set_xlim(0, 2350)
    ax.set_ylim(-0.7, 4.7)
    ax.set_xlabel("recording time (seconds)")
    ax.set_title(
        "A. Only four 256-s windows survived",
        loc="left",
        fontweight="bold",
    )
    ax.grid(axis="x", color=GRID, lw=0.7, alpha=0.7)
    ax.text(
        0.02,
        0.02,
        "Windows 1–3 come from one long bout;\n"
        "adjacent windows overlap by 128 s.",
        transform=ax.transAxes,
        color=MUTED,
        fontsize=8.7,
    )

    ax = fig.add_subplot(gs[0, 1])
    show = (coherence["f"] >= 0) & (coherence["f"] <= 0.06)
    ax.plot(coherence["f"][show], coherence["cxy"][show], color=PURPLE, lw=2.2)
    ax.axhline(
        threshold,
        color=HR,
        ls="--",
        lw=1.6,
        label=f"nominal cutoff {threshold:.3f}",
    )
    ax.axhline(
        bonferroni_threshold,
        color=MUTED,
        ls=":",
        lw=1.5,
        label=f"3-participant sensitivity {bonferroni_threshold:.3f}",
    )
    ax.scatter([evaluated_hz], [observed], color=SIGMA, s=65, zorder=4)
    ax.annotate(
        f"observed {observed:.3f}\n@ {evaluated_hz:.6f} Hz",
        (evaluated_hz, observed),
        xytext=(0.030, 0.79),
        arrowprops={"arrowstyle": "->", "color": SIGMA},
        color=SIGMA,
        fontsize=9,
        fontweight="bold",
    )
    ax.set_xlim(0, 0.06)
    ax.set_ylim(0, 1)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("magnitude-squared coherence")
    ax.set_title(
        "B. The observed point exceeds 0.632",
        loc="left",
        fontweight="bold",
    )
    ax.legend(frameon=False, fontsize=8.3, loc="upper right")

    ax = fig.add_subplot(gs[0, 2])
    ax.axis("off")
    ax.set_title(
        "C. What the number means",
        loc="left",
        fontweight="bold",
        pad=10,
    )
    ax.text(
        0.02,
        0.91,
        r"$C_{crit}=1-\alpha^{1/(K-1)}$"
        "\n\n"
        rf"$=1-0.05^{{1/(4-1)}}$"
        "\n\n"
        rf"$={threshold:.6f}$",
        transform=ax.transAxes,
        color=INK,
        fontsize=14,
        va="top",
    )
    ax.text(
        0.02,
        0.54,
        f"K = {k} accepted windows\n"
        f"α = {alpha:.2f}, fixed frequency only\n"
        f"nominal uncorrected p ≈ {nominal_p:.3f}",
        transform=ax.transAxes,
        color=INK,
        fontsize=10.2,
        va="top",
        linespacing=1.5,
    )
    ax.text(
        0.02,
        0.29,
        "Interpret carefully",
        transform=ax.transAxes,
        color=HR,
        fontsize=10.5,
        fontweight="bold",
    )
    ax.text(
        0.02,
        0.24,
        "• not learned from this participant\n"
        "• not an F3/F4 measurement\n"
        "• pointwise and uncorrected across people\n"
        "• does not establish LC activity",
        transform=ax.transAxes,
        color=MUTED,
        fontsize=9.1,
        va="top",
        linespacing=1.45,
    )

    fig.text(
        0.055,
        0.055,
        "Verdict: the arithmetic reproduces exactly and the observed point passes the nominal cutoff, "
        "but it does not pass the illustrated three-participant family-wise sensitivity threshold.",
        color=INK,
        fontsize=10,
        fontweight="bold",
    )
    png, svg = _save(fig, "3A_RESP0699_coherence_threshold_explained")
    return {
        "subject": subject,
        "png": repository_relative(png),
        "svg": repository_relative(svg),
        "welch_window_starts_s": starts,
        "welch_window_seconds": int(lecci.NPERSEG),
        "welch_overlap_seconds": int(noverlap),
        "K": k,
        "alpha": alpha,
        "evaluated_frequency_hz": evaluated_hz,
        "observed_coherence": observed,
        "analytic_threshold": threshold,
        "nominal_pointwise_p": nominal_p,
        "three_participant_bonferroni_threshold_sensitivity": bonferroni_threshold,
    }


def _compute_3b_stage(cache, mat: dict, profile: dict, stage_name: str) -> dict:
    rr = np.asarray(mat["rr_4"], float)
    labels = np.asarray(mat["stage_lab"]).astype(str)
    contacts = [str(value) for value in mat["contacts"]]
    roi_mask = (
        np.asarray(mat["frontal_contact_mask"], bool)
        if mat["cohort"] == "RESPect"
        else np.ones(len(contacts), bool)
    )
    eligible_contacts = [
        contact for contact, keep in zip(contacts, roi_mask) if keep
    ]
    if stage_name == "NREM":
        pooled = np.where(np.isin(labels, ("NREM", "N2", "N3")), "NREM", "")
        stable = stable_stage_epoch_indices(pooled, "NREM")
    else:
        stable = stable_stage_epoch_indices(labels, stage_name)
    pool = grid._stage_pool(rr, stable)
    keep_epochs = set(stable.tolist())
    troughs = []
    event_contacts = []
    for contact in eligible_contacts:
        times = stage_so_times(
            cache,
            contact,
            keep_epochs,
            percentile=float(profile["endpoint_3b"]["so_amplitude_percentile"]),
        )
        if len(times) >= profile["endpoint_3b"]["minimum_so_per_contact"]:
            troughs.append(times)
            event_contacts.append(contact)
    baseline = rr_baseline_hr(rr[pool])
    estimate = subject_so_triggered(
        rr,
        troughs,
        baseline,
        pool,
        n_sur=grid.GRID_SURROGATES_3B,
        rng=np.random.RandomState(grid._seed(mat["subject"], stage_name, "3B-grid")),
        domain="rr",
        minimum_channels=profile["endpoint_3b"]["minimum_contacts"],
        channel_ids=event_contacts,
        minimum_events_per_channel=profile["endpoint_3b"]["minimum_so_per_contact"],
        minimum_surrogate_pool_samples=profile["endpoint_3b"][
            "minimum_finite_stage_samples"
        ],
    )
    if estimate is None:
        raise RuntimeError(f"{mat['subject']} {stage_name} unexpectedly unavailable")
    estimate["stage_mean_hr"] = float(baseline)
    return estimate


def make_3b(profile: dict) -> dict:
    subject = "HUP160_phaseII"
    cache_path = ROOT / "data" / "derived" / "lc_infraslow" / f"{subject}.npz"
    final = _final_subject(LOCKED_QC_ARTIFACT, subject)["result_3b"]
    with np.load(cache_path, allow_pickle=False) as cache:
        mat = materialize(cache, profile)
        recomputed_summary = grid.analyse_3b(cache, mat, profile)
        estimates = {
            stage: _compute_3b_stage(cache, mat, profile, stage)
            for stage in ("N2", "N3", "NREM")
        }

    for stage in ("N2", "N3", "NREM"):
        if not recomputed_summary["stages"][stage]["available_under_profile"]:
            raise RuntimeError(f"3B {stage} availability differs from final result")
        expected = final["stages"][stage]["estimate"]
        _close(
            estimates[stage]["event_locked_local_change_pct"],
            expected["event_locked_local_change_pct"],
            f"3B {stage} local change",
        )
        _close(
            estimates[stage]["peak_lag_s"],
            expected["peak_lag_s"],
            f"3B {stage} lag",
        )
        if estimates[stage]["n_so_total"] != expected["n_so_total"]:
            raise RuntimeError(f"3B {stage} event count differs from final artifact")

    stage_order = ["N2", "N3", "NREM"]
    stage_labels = ["N2-like", "N3-like", "pooled NREM"]
    colors = [SO, PURPLE, HR]

    fig = plt.figure(figsize=(14, 8.8))
    gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.28)
    fig.suptitle(
        "Question 3B — a real participant passing N2-like, N3-like, and pooled-NREM support",
        x=0.06,
        y=0.99,
        ha="left",
        fontsize=16,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.06,
        0.955,
        "HUP160 · curves are averaged in RR space and displayed as heart rate · descriptive only",
        ha="left",
        fontsize=10.5,
        color=MUTED,
    )

    ax = fig.add_subplot(gs[0, 0])
    for stage, label, color in zip(stage_order, stage_labels, colors):
        estimate = estimates[stage]
        lag = np.asarray(estimate["lag_s"], float)
        curve = np.asarray(estimate["curve"], float)
        baseline = estimate["local_pre_event_mean_hr"]
        relative = 100.0 * (curve - baseline) / baseline
        ax.plot(
            lag,
            relative,
            color=color,
            lw=2.2,
            label=f"{label} ({estimate['n_so_total']} SOs)",
        )
    ax.axvline(0, color=INK, ls="--", lw=1.2)
    ax.axhline(0, color=GRID, lw=0.9)
    ax.set_title(
        "A. Align every cardiac window to the SO trough at time 0",
        loc="left",
        fontweight="bold",
    )
    ax.set_xlabel("seconds from slow-oscillation trough")
    ax.set_ylabel("heart-rate change from local pre-event mean (%)")
    ax.legend(frameon=False, fontsize=8.7)
    ax.text(
        0.52,
        0.05,
        "time 0 = cortical SO down-state trough",
        transform=ax.transAxes,
        color=INK,
        fontsize=8.7,
    )

    ax = fig.add_subplot(gs[0, 1])
    estimate = estimates["NREM"]
    lag = np.asarray(estimate["lag_s"], float)
    curve = np.asarray(estimate["curve"], float)
    ax.plot(lag, curve, color=HR, lw=2.6)
    ax.axvline(0, color=INK, ls="--", lw=1.2)
    ax.axhline(
        estimate["local_pre_event_mean_hr"],
        color=MUTED,
        ls=":",
        lw=1.2,
        label=f"pre-event mean {estimate['local_pre_event_mean_hr']:.2f} bpm",
    )
    post = lag >= 0
    peak_index = np.where(post)[0][int(np.argmax(curve[post]))]
    ax.scatter([lag[peak_index]], [curve[peak_index]], color=HR, s=55, zorder=4)
    ax.annotate(
        f"participant-average HR peak\n{lag[peak_index]:.2f} s after trough",
        (lag[peak_index], curve[peak_index]),
        xytext=(0.0, curve[peak_index] + 0.25),
        arrowprops={"arrowstyle": "->", "color": HR},
        color=HR,
        fontsize=9,
    )
    ax.set_title(
        "B. The actual pooled-NREM participant curve",
        loc="left",
        fontweight="bold",
    )
    ax.set_xlabel("seconds from slow-oscillation trough")
    ax.set_ylabel("heart rate (beats/min)")
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")

    ax = fig.add_subplot(gs[1, 0])
    changes = [
        estimates[stage]["event_locked_local_change_pct"] for stage in stage_order
    ]
    bars = ax.bar(stage_labels, changes, color=colors, width=0.62)
    ax.axhline(0, color=INK, lw=0.9)
    ax.set_title(
        "C. Descriptive local HR increase for this participant",
        loc="left",
        fontweight="bold",
    )
    ax.set_ylabel("post-peak versus pre-event mean (%)")
    for bar, value in zip(bars, changes):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.025,
            f"+{value:.3f}%",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
        )
    ax.set_ylim(0, max(changes) * 1.28)

    ax = fig.add_subplot(gs[1, 1])
    lags = [estimates[stage]["peak_lag_s"] for stage in stage_order]
    bars = ax.bar(stage_labels, lags, color=colors, width=0.62)
    ax.set_title(
        "D. Naji-style timing: mean of each contact's HR-peak time",
        loc="left",
        fontweight="bold",
    )
    ax.set_ylabel("seconds after SO trough")
    for bar, value in zip(bars, lags):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.08,
            f"{value:.2f} s",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
        )
    ax.set_ylim(0, max(lags) * 1.28)

    fig.text(
        0.06,
        0.012,
        "Visual verdict: this participant shows a small HR rise after SO troughs in all three stage summaries. "
        "The p/z fields remain disabled because the current null does not preserve local HR trends and clustered SO timing.",
        color=INK,
        fontsize=10,
        fontweight="bold",
    )
    png, svg = _save(fig, "3B_HUP160_event_locked_heart_rate")
    return {
        "subject": subject,
        "png": repository_relative(png),
        "svg": repository_relative(svg),
        "stages": {
            stage: {
                "n_so_total": int(estimates[stage]["n_so_total"]),
                "local_change_pct": float(
                    estimates[stage]["event_locked_local_change_pct"]
                ),
                "mean_channel_peak_lag_s": float(estimates[stage]["peak_lag_s"]),
            }
            for stage in stage_order
        },
    }


def make_3b_respect_naji_check(profile: dict) -> dict:
    """Show the closest Naji-style check possible for RESP0699."""
    subject = "sub-RESP0699"
    cache_path = ROOT / "data" / "derived" / "ds003848" / f"{subject}.npz"
    artifact_path = (
        ROOT
        / "outputs"
        / "qc_grid"
        / "staging_window_support_v1"
        / "respect_qc_grid.json"
    )
    metadata_path = (
        ROOT
        / "data"
        / "ds003848_raw"
        / "sub-RESP0699_ses-1_task-sleep_run-030608_ieeg.json"
    )
    channels_path = (
        ROOT
        / "data"
        / "ds003848_raw"
        / "sub-RESP0699_ses-1_task-sleep_run-030608_channels.tsv"
    )
    final = _final_subject(artifact_path, subject)["result_3b"]

    with metadata_path.open() as handle:
        metadata = json.load(handle)
    with channels_path.open(newline="") as handle:
        channel_rows = list(csv.DictReader(handle, delimiter="\t"))
    scalp_labels = {
        row["name"].upper()
        for row in channel_rows
        if row.get("name", "").upper() in {"F3", "F4"}
    }
    if int(metadata.get("EEGChannelCount", -1)) != 0 or scalp_labels:
        raise RuntimeError("RESP0699 unexpectedly contains scalp F3/F4 EEG")

    with np.load(cache_path, allow_pickle=False) as cache:
        mat = materialize(cache, profile)
        recomputed_summary = grid.analyse_3b(cache, mat, profile)
        estimate = _compute_3b_stage(cache, mat, profile, "N3")

    expected_stage = final["stages"]["N3"]
    expected = expected_stage["estimate"]
    if not recomputed_summary["stages"]["N3"]["available_under_profile"]:
        raise RuntimeError("RESP0699 final N3 endpoint unexpectedly unavailable")
    if recomputed_summary["stages"]["N2"]["available_under_profile"]:
        raise RuntimeError("RESP0699 final N2 endpoint unexpectedly available")
    for key in (
        "pct_above_stage_mean",
        "event_locked_local_change_pct",
        "peak_lag_s",
        "participant_curve_peak_lag_s",
    ):
        _close(estimate[key], expected[key], f"RESP0699 Naji check {key}")
    if estimate["n_so_total"] != expected["n_so_total"]:
        raise RuntimeError("RESP0699 Naji-check event count differs from final artifact")

    lag = np.asarray(estimate["lag_s"], float)
    curve = np.asarray(estimate["curve"], float)
    local_baseline = float(estimate["local_pre_event_mean_hr"])
    relative = 100.0 * (curve - local_baseline) / local_baseline
    paper_peak_lag = float(estimate["peak_lag_s"])
    participant_peak_lag = float(estimate["participant_curve_peak_lag_s"])
    stage_mean_magnitude = float(estimate["pct_above_stage_mean"])
    local_change = float(estimate["event_locked_local_change_pct"])
    event_contacts = [str(value) for value in estimate["event_contact_ids"]]

    naji_stage2 = 12.09
    naji_stage2_pm = 1.48
    naji_sws = 3.35
    naji_sws_pm = 1.01

    fig = plt.figure(figsize=(15, 7.4))
    gs = fig.add_gridspec(
        1,
        3,
        left=0.055,
        right=0.98,
        bottom=0.17,
        top=0.79,
        wspace=0.32,
        width_ratios=[1.15, 1.0, 1.05],
    )
    fig.suptitle(
        "RESP0699 — closest available check of Naji et al. 2019",
        x=0.055,
        y=0.975,
        ha="left",
        fontsize=16,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.055,
        0.925,
        "Right-frontal ECoG adaptation only · 5 contacts · 167 SOs · 510 stable s · ECG present but zero scalp EEG channels and no F3/F4.",
        ha="left",
        fontsize=10.5,
        color=MUTED,
    )

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(lag, relative, color=HR, lw=2.6)
    ax.axvline(0, color=INK, ls="--", lw=1.2, label="SO trough")
    ax.axhline(0, color=GRID, lw=0.9)
    ax.axvline(
        paper_peak_lag,
        color=PURPLE,
        ls=":",
        lw=1.5,
        label=f"mean contact peak time {paper_peak_lag:.2f} s",
    )
    peak_index = int(np.argmin(np.abs(lag - participant_peak_lag)))
    ax.scatter(
        [lag[peak_index]],
        [relative[peak_index]],
        color=HR,
        s=58,
        zorder=4,
    )
    ax.annotate(
        f"participant-average peak\n{participant_peak_lag:.2f} s",
        (lag[peak_index], relative[peak_index]),
        xytext=(-0.35, 0.60),
        arrowprops={"arrowstyle": "->", "color": HR},
        color=HR,
        fontsize=8.7,
    )
    ax.set_xlabel("seconds from SO down-state trough")
    ax.set_ylabel("HR change from local pre-event mean (%)")
    ax.set_title(
        "A. Actual SWS/N3-like event-locked curve",
        loc="left",
        fontweight="bold",
    )
    ax.legend(frameon=False, fontsize=8.2, loc="lower right")
    ax = fig.add_subplot(gs[0, 1])
    x = np.arange(2)
    width = 0.34
    naji_values = [naji_stage2, naji_sws]
    naji_pm = [naji_stage2_pm, naji_sws_pm]
    bars = ax.bar(
        x - width / 2,
        naji_values,
        width,
        yerr=naji_pm,
        color=SO,
        alpha=0.84,
        capsize=4,
        label="Naji published group value",
    )
    ax.bar(
        [x[1] + width / 2],
        [stage_mean_magnitude],
        width,
        color=SIGMA,
        label="RESP0699",
    )
    ax.bar(
        [x[0] + width / 2],
        [0],
        width,
        facecolor="white",
        edgecolor=FAIL,
        hatch="//",
    )
    ax.text(
        x[0] + width / 2,
        0.35,
        "N2\nunavailable",
        ha="center",
        va="bottom",
        color=MUTED,
        fontsize=8.5,
    )
    for bar, value, uncertainty in zip(bars, naji_values, naji_pm):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + uncertainty + 0.35,
            f"{value:.2f}%",
            ha="center",
            fontsize=8.7,
            fontweight="bold",
            color=SO,
        )
    ax.text(
        x[1] + width / 2,
        stage_mean_magnitude + 0.35,
        f"{stage_mean_magnitude:.3f}%",
        ha="center",
        fontsize=8.7,
        fontweight="bold",
        color=SIGMA,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(["Stage 2 / N2", "SWS / N3-like"])
    ax.set_ylabel("HR peak above whole-stage mean (%)")
    ax.set_ylim(0, 15.5)
    ax.set_title(
        "B. Direct-ish magnitude comparison",
        loc="left",
        fontweight="bold",
    )
    ax.legend(frameon=False, fontsize=7.9, loc="upper right")

    ax = fig.add_subplot(gs[0, 2])
    ax.axis("off")
    ax.set_title(
        "C. What can honestly be concluded",
        loc="left",
        fontweight="bold",
        pad=10,
    )
    ax.text(
        0.02,
        0.91,
        "Partial qualitative agreement",
        transform=ax.transAxes,
        color=GREEN,
        fontsize=10.5,
        fontweight="bold",
    )
    ax.text(
        0.02,
        0.85,
        "• HR rises after the SO trough\n"
        f"• mean contact timing is {paper_peak_lag:.2f} s\n"
        f"• local-baseline HR increase is +{local_change:.3f}%\n"
        "• direction and timescale are broadly compatible",
        transform=ax.transAxes,
        color=INK,
        fontsize=9.4,
        va="top",
        linespacing=1.5,
    )
    ax.text(
        0.02,
        0.61,
        "Not a Naji replication",
        transform=ax.transAxes,
        color=HR,
        fontsize=10.5,
        fontweight="bold",
    )
    ax.text(
        0.02,
        0.55,
        "• no F3/F4 or other scalp EEG\n"
        "• right-frontal ECoG, not bilateral scalp\n"
        "• N2 endpoint is unavailable\n"
        "• smaller than Naji's SWS group magnitude\n"
        "• no valid event-locking p/z\n"
        "• no behavioral processing-speed test",
        transform=ax.transAxes,
        color=MUTED,
        fontsize=9.1,
        va="top",
        linespacing=1.42,
    )
    ax.text(
        0.02,
        0.18,
        "Used contacts",
        transform=ax.transAxes,
        color=INK,
        fontsize=9.5,
        fontweight="bold",
    )
    ax.text(
        0.02,
        0.13,
        ", ".join(event_contacts),
        transform=ax.transAxes,
        color=PURPLE,
        fontsize=9.1,
        va="top",
    )

    fig.text(
        0.055,
        0.065,
        "Verdict: RESP0699 shows a small SO-following cardiac acceleration, but the missing F3/F4 and N2 data "
        "prevent a direct replication and the result does not validate an LC proxy.",
        color=INK,
        fontsize=10,
        fontweight="bold",
    )
    png, svg = _save(fig, "3B_RESP0699_naji_frontal_ecog_check")
    return {
        "subject": subject,
        "png": repository_relative(png),
        "svg": repository_relative(svg),
        "scalp_eeg_channel_count": int(metadata["EEGChannelCount"]),
        "has_F3": "F3" in scalp_labels,
        "has_F4": "F4" in scalp_labels,
        "stage": "author-SWS-selected N3-like",
        "stable_seconds": int(expected_stage["stable_seconds"]),
        "event_contact_ids": event_contacts,
        "n_channels": int(estimate["n_channels"]),
        "n_so_total": int(estimate["n_so_total"]),
        "pct_above_stage_mean": stage_mean_magnitude,
        "event_locked_local_change_pct": local_change,
        "mean_contact_peak_lag_s": paper_peak_lag,
        "participant_curve_peak_lag_s": participant_peak_lag,
        "n2_available": False,
        "inference_enabled": False,
        "naji_reference": {
            "stage2_hr_peak_pct": naji_stage2,
            "stage2_plus_minus": naji_stage2_pm,
            "sws_hr_peak_pct": naji_sws,
            "sws_plus_minus": naji_sws_pm,
        },
    }


def make_3d(profile: dict) -> dict:
    subject = "HUP172_phaseII"
    cache_path = ROOT / "data" / "derived" / "lc_infraslow" / f"{subject}.npz"
    final = _final_subject(LOCKED_QC_ARTIFACT, subject)["result_3d"]
    with np.load(cache_path, allow_pickle=False) as cache:
        mat = materialize(cache, profile)
        reconstructed = analyse_3d_cache_support(cache, mat, profile)
        contacts = [str(value) for value in mat["contacts"]]
        payload, missing, errors = _validated_neutral_payload(cache, len(contacts))
        if payload is None or missing or errors:
            raise RuntimeError(f"3D neutral payload invalid: {missing} {errors}")
        candidate_contact = payload["candidate_contact_index"]
        events = []
        for contact_index in range(len(contacts)):
            keep = candidate_contact == contact_index
            candidates = np.column_stack(
                (
                    payload["candidate_sample"][keep],
                    payload["candidate_amplitude"][keep],
                    payload["candidate_start"][keep],
                    payload["candidate_stop"][keep],
                )
            )
            events.append(
                channel_night_events(
                    payload["rms"][contact_index],
                    payload["phase"][contact_index],
                    candidates,
                    np.asarray(mat["stage_lab"]).astype(str),
                )
            )

    if not reconstructed["descriptive_effect_available"]:
        raise RuntimeError("HUP172 final 3D effect unexpectedly unavailable")
    expected = final["descriptive_effect"]
    observed = reconstructed["descriptive_effect"]
    _close(observed["participant_R"], expected["participant_R"], "3D R")
    _close(
        observed["participant_preferred_phase_deg"],
        expected["participant_preferred_phase_deg"],
        "3D phase",
    )
    if (
        observed["n_paired_events_across_qualified_contacts"]
        != expected["n_paired_events_across_qualified_contacts"]
    ):
        raise RuntimeError("3D event count differs from final artifact")

    pooled = reconstructed["stage_descriptive_support"]["pooled_NREM"]
    pooled_epochs = set(
        np.where(np.isin(mat["stage_lab"], ("NREM", "N2", "N3")))[0].tolist()
    )
    contact_records = pooled["per_contact"]
    phases_by_contact = {}
    for contact_index, (contact, event, record) in enumerate(
        zip(contacts, events, contact_records)
    ):
        phases = np.asarray(
            stage_event_pairs(event, pooled_epochs)["phases"], float)
        if record["qualifies_for_stage_vector"]:
            phases_by_contact[contact] = phases
            if len(phases) != record["n_paired_events"]:
                raise RuntimeError(f"3D phase count mismatch for {contact}")

    qualified = list(phases_by_contact)
    paired_event_count = sum(map(len, phases_by_contact.values()))
    if paired_event_count != observed[
        "n_paired_events_across_qualified_contacts"
    ]:
        raise RuntimeError("3D plotted phase count differs from final artifact")

    colors = [SIGMA, SO, PURPLE, HR, GREEN, "#496f9b"]
    selected_colors = {
        contact: colors[index % len(colors)]
        for index, contact in enumerate(qualified)
    }
    fractions = [record["valid_pooled_nrem_fraction"] for record in contact_records]
    event_counts = [record["n_paired_events"] for record in contact_records]
    selected = [record["qualifies_for_stage_vector"] for record in contact_records]
    threshold = profile["endpoint_3d"]["minimum_valid_nrem_fraction_per_contact"]

    fig = plt.figure(figsize=(15, 7.1))
    gs = fig.add_gridspec(
        1,
        3,
        left=0.06,
        right=0.98,
        bottom=0.24,
        top=0.78,
        wspace=0.34,
        width_ratios=[1.15, 1, 1],
    )
    fig.suptitle(
        "Question 3D — a real participant with a pooled SO–spindle phase vector",
        x=0.05,
        y=0.975,
        ha="left",
        fontsize=16,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.05,
        0.925,
        f"HUP172 · {len(qualified)} of {len(contacts)} contacts qualify · "
        f"{paired_event_count:,} paired events · descriptive only",
        ha="left",
        fontsize=10.5,
        color=MUTED,
    )

    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(contacts))
    bar_colors = [
        selected_colors.get(contact, FAIL) if keep else FAIL
        for contact, keep in zip(contacts, selected)
    ]
    bars = ax.bar(x, fractions, color=bar_colors, width=0.68)
    ax.axhline(
        threshold,
        color=HR,
        ls="--",
        lw=1.4,
        label=f"base valid-NREM fraction rule = {threshold:.0%}",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(contacts, rotation=35, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("valid pooled-NREM fraction")
    ax.set_title(
        "A. First decide which contacts have enough support",
        loc="left",
        fontweight="bold",
    )
    ax.legend(frameon=False, fontsize=8.2, loc="lower right")
    for bar, count, keep in zip(bars, event_counts, selected):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            min(1.015, bar.get_height() + 0.025),
            f"{count}\npairs",
            ha="center",
            va="bottom",
            fontsize=7.8,
            color=INK if keep else MUTED,
            fontweight="bold" if keep else "normal",
        )

    ax = fig.add_subplot(gs[0, 1], projection="polar")
    edges = np.linspace(-np.pi, np.pi, 25)
    centers = (edges[:-1] + edges[1:]) / 2
    for contact, phases in phases_by_contact.items():
        hist, _ = np.histogram(phases, bins=edges)
        fraction = hist / max(hist.sum(), 1)
        closed_x = np.r_[centers, centers[0]]
        closed_y = np.r_[fraction, fraction[0]]
        ax.plot(
            closed_x,
            closed_y,
            color=selected_colors[contact],
            lw=2,
            label=f"{contact} (n={len(phases)})",
        )
    ax.set_title(
        "B. Where did each paired spindle fall\non the SO phase circle?",
        fontsize=10.2,
        fontweight="bold",
        pad=18,
    )
    ax.set_yticklabels([])
    phase_handles, phase_labels = ax.get_legend_handles_labels()
    fig.legend(
        phase_handles,
        phase_labels,
        frameon=False,
        fontsize=7.8,
        loc="center",
        bbox_to_anchor=(0.51, 0.155),
        ncol=3,
    )

    ax = fig.add_subplot(gs[0, 2], projection="polar")
    contact_vectors = []
    for contact, phases in phases_by_contact.items():
        vector = complex(np.mean(np.exp(1j * phases)))
        contact_vectors.append(vector)
        ax.annotate(
            "",
            xy=(np.angle(vector), abs(vector)),
            xytext=(0, 0),
            arrowprops={
                "arrowstyle": "-|>",
                "color": selected_colors[contact],
                "lw": 2.2,
            },
        )
        ax.text(
            np.angle(vector),
            abs(vector) + 0.025,
            f"{contact}\nR={abs(vector):.3f}",
            color=selected_colors[contact],
            ha="center",
            va="center",
            fontsize=8,
        )
    participant_vector = complex(np.mean(contact_vectors))
    ax.annotate(
        "",
        xy=(np.angle(participant_vector), abs(participant_vector)),
        xytext=(0, 0),
        arrowprops={"arrowstyle": "-|>", "color": INK, "lw": 4},
    )
    ax.set_ylim(0, 0.38)
    ax.set_rticks([0.1, 0.2, 0.3])
    ax.set_rlabel_position(225)
    ax.set_title(
        "C. Equal-weight contact arrows → one participant arrow",
        fontsize=10.2,
        fontweight="bold",
        pad=18,
    )
    fig.text(
        0.82,
        0.155,
        f"black arrow: R={abs(participant_vector):.3f}, "
        f"phase={np.degrees(np.angle(participant_vector)):.1f}°",
        ha="center",
        color=INK,
        fontsize=8.8,
        fontweight="bold",
    )

    fig.text(
        0.06,
        0.095,
        "Panel A: colored contacts are included; grey contacts are excluded.\n"
        "The 80% valid-NREM support rule is repository-defined, not paper-derived.",
        color=MUTED,
        fontsize=8.3,
    )
    fig.text(
        0.05,
        0.035,
        "Visual verdict: HUP172 has enough descriptive support and a "
        f"participant vector with R={observed['participant_R']:.3f}. "
        "Inference remains disabled because nearest-event pairing can create apparent phase structure under independence.",
        color=INK,
        fontsize=10,
        fontweight="bold",
    )
    png, svg = _save(fig, "3D_HUP172_SO_spindle_phase")
    return {
        "subject": subject,
        "png": repository_relative(png),
        "svg": repository_relative(svg),
        "qualified_contacts": qualified,
        "paired_events_by_contact": {
            contact: int(len(phases))
            for contact, phases in phases_by_contact.items()
        },
        "n_paired_events": int(
            observed["n_paired_events_across_qualified_contacts"]
        ),
        "participant_R": float(observed["participant_R"]),
        "participant_preferred_phase_deg": float(
            observed["participant_preferred_phase_deg"]
        ),
        "inference_enabled": bool(reconstructed["inference_enabled"]),
    }


def main() -> None:
    manifest_path, preflight = write_preflight_manifest(OUT)
    try:
        if preflight.get("code_dirty_at_start") is not False:
            raise RuntimeError(
                "visual generation requires clean code/config at run start; "
                "changes under outputs/participant_result_visuals are ignored"
            )
        snapshot = freeze_visual_inputs(generator_path=Path(__file__))
    except Exception as error:
        write_failed_manifest(
            manifest_path,
            preflight,
            validation_state="input_validation_failed",
            error=error,
        )
        raise
    manifest_base = visual_manifest_base(snapshot, preflight)
    write_in_progress_manifest(OUT, manifest_base)

    try:
        profile = load_qc_profile(PROFILE_ID)
        results = {
            "analysis_version": manifest_base["analysis_version"],
            "cache_schema_version": manifest_base["cache_schema_version"],
            "profile_id": PROFILE_ID,
            "profile_sha256": manifest_base["profile_sha256"],
            "profile_file_sha256": manifest_base["profile_file_sha256"],
            "figures": {
                "3A": make_3a(profile),
                "3A_threshold_explainer": make_3a_threshold_explainer(profile),
                "3B": make_3b(profile),
                "3B_RESP0699_Naji_check": make_3b_respect_naji_check(profile),
                "3D": make_3d(profile),
            },
            "interpretation_contract": {
                "figures_are_recomputed_from_current_neutral_caches": True,
                "figures_are_synthetic": False,
                "availability_pass_is_not_hypothesis_support": True,
                "3b_inference_enabled": False,
                "3d_inference_enabled": False,
            },
        }
        OUT.mkdir(parents=True, exist_ok=True)
        values_path = OUT / "figure_values.json"
        atomic_json_dump(results, str(values_path))
        (OUT / "README.md").write_text(
            render_visual_readme(results),
            encoding="utf-8",
        )
        (OUT / "RESP0699_THRESHOLD_AND_NAJI_CHECK.md").write_text(
            render_resp0699_check(results),
            encoding="utf-8",
        )
        result_paths = [
            OUT / name for name in VISUAL_RESULT_NAMES
        ]
        assert_visual_inputs_unchanged(snapshot)
        manifest_path = write_terminal_manifest(
            OUT,
            manifest_base,
            result_paths,
        )
        assert_visual_inputs_unchanged(snapshot)
        validate_terminal_visual_manifest(
            manifest_path,
            result_paths,
        )
    except Exception as error:
        write_failed_manifest(
            manifest_path,
            manifest_base,
            validation_state="terminal_validation_failed",
            error=error,
        )
        raise
    print(json.dumps(results, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
