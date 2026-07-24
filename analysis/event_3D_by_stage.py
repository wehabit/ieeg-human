"""Paper-aligned, event-based slow-oscillation -> spindle coupling.

Primary estimator:
  * slow oscillation (SO): 0.16-1.25 Hz, complete 0.8-2.0 s cycles, top 25% amplitude
  * spindle: 12-16 Hz, 200-ms RMS, 75th-percentile threshold, 0.5-3.0 s duration
  * thresholds: one artifact-free NREM distribution per participant/channel/night
  * independence: assign each spindle to its nearest SO and retain at most the maximum-amplitude
    spindle for each SO epoch
  * reporting: contacts are reduced to one equal-contact vector per participant. Per-channel
    Rayleigh and cohort participant-rotation values are not production inference: the +/-2 s pairing
    window itself can create a common phase under independent spindle/SO trains.

This is a Staresina/Helfrich hybrid adaptation, not a faithful replication: Helfrich selected the
maximum 12-16-Hz amplitude in every +/-2-s SO epoch, whereas this path first requires a
duration-qualified spindle and then retains at most one around each SO. The previous implementation
thresholded each 600-s chunk separately, admitted IED/filter-ringing events, and counted multiple
spindle peaks around one SO as independent observations. Those choices can generate severe false
positives on null and impulse-only signals.

Staging remains an unvalidated N2-like/N3-like iEEG power proxy because HUP has no EOG/EMG/scalp
EEG. Consequently the pooled-NREM result is primary and stage contrasts are exploratory.

    .venv/bin/python analysis/event_3D_by_stage.py [--subjects 165,...] [--hours 7] --force
"""
import argparse
import json
import os
import traceback

import numpy as np
from scipy import signal

from infraslow_rr_sigma_coherence import sess, pull_continuous, notch, ROOT
from cohort_3A_cortical import COHORT, cortical_channels, find_night
from cohort_stages_3ABD import band_sos, stage_epochs, EPOCH, CHUNK_S, SWA_BAND
from cache_lc_series import (
    FILTER_EDGE_S,
    MIN_CONTACT_FRACTION_PER_BIN,
    MIN_CONTACTS,
    aggregate_staging_features,
    prepare_continuous_signal,
    staging_epoch_features,
)
from results_3A_tutorial_style import ied_clean_mask
from pipeline_version import (ANALYSIS_VERSION, atomic_json_dump, file_sha256,
                              finite_float_or_none,
                              source_tree_sha256, start_run_manifest,
                              validated_complete_run_exists, write_run_manifest)


OUT = os.path.join(ROOT, "outputs", "event_3D_by_stage")
SO_BAND = (0.16, 1.25)
SPINDLE_BAND = (12.0, 16.0)
SP_DUR = (0.5, 3.0)
SO_DUR = (0.8, 2.0)
PAIR_WINDOW_S = 2.0
IED_PAD_S = 2.5
EVENT_PERCENTILE = 75
EVENT_FS = 20.0
MIN_EVENTS = 20
MIN_POOLED_CONTACTS = 3
MIN_POOLED_EVENTS = 200
MIN_POOLED_VALID_S_PER_CONTACT = 20 * 60
MIN_POOLED_NREM_COVERAGE = 0.80
PRODUCTION_3D_INFERENCE_ENABLED = False


def production_config(hours=7.0):
    """Canonical production configuration, shared with the fail-closed summary."""
    return dict(
        hours=float(hours), analysis_version=ANALYSIS_VERSION,
        so_band_hz=list(SO_BAND), spindle_band_hz=list(SPINDLE_BAND),
        so_duration_s=list(SO_DUR), spindle_duration_s=list(SP_DUR),
        event_percentile=EVENT_PERCENTILE, event_sampling_hz=EVENT_FS,
        ied_mask_padding_s=IED_PAD_S, minimum_events_per_contact=MIN_EVENTS,
        threshold_scope="channel-night pooled NREM",
        pairing=f"duration-qualified hybrid; one spindle per SO within +/-{PAIR_WINDOW_S:g} s",
        staging_qc=dict(
            minimum_contacts=MIN_CONTACTS,
            minimum_contact_fraction_per_epoch=MIN_CONTACT_FRACTION_PER_BIN,
            swa_normalization=(
                "divide each fixed contact by its full-night median clean-epoch SWA")),
        pooled_qc=dict(
            minimum_contacts=MIN_POOLED_CONTACTS,
            minimum_paired_events=MIN_POOLED_EVENTS,
            minimum_valid_nrem_seconds_per_contact=MIN_POOLED_VALID_S_PER_CONTACT,
            minimum_valid_nrem_fraction_per_contact=MIN_POOLED_NREM_COVERAGE),
        production_inference_enabled=PRODUCTION_3D_INFERENCE_ENABLED,
        group_inference=(
            "disabled; participant vectors descriptive pending a time-shift/block null that "
            "repeats event pairing and contact aggregation"))


def rayleigh(ph):
    """Resultant length, calibrated small-sample Rayleigh p value, and mean phase."""
    ph = np.asarray(ph, float)
    n = len(ph)
    if n < MIN_EVENTS:
        return np.nan, np.nan, np.nan
    C, S = np.mean(np.cos(ph)), np.mean(np.sin(ph))
    R = float(np.hypot(C, S))
    z = n * R ** 2
    p = np.exp(-z)
    if n < 50:
        # Higher-order finite-n expansion used by established circular-statistics packages.
        p *= (1 + (2 * z - z ** 2) / (4 * n)
              - (24 * z - 132 * z ** 2 + 76 * z ** 3 - 9 * z ** 4)
              / (288 * n ** 2))
    p = float(np.clip(p, np.finfo(float).tiny, 1.0))
    return R, p, float(np.arctan2(S, C))


def bh_fdr(p_values, alpha=0.05):
    """Benjamini-Hochberg decisions, returned in original order."""
    p = np.asarray(p_values, float)
    out = np.zeros(len(p), bool)
    good = np.isfinite(p)
    if not good.any():
        return out
    idx = np.where(good)[0]
    order = idx[np.argsort(p[idx])]
    passed = p[order] <= alpha * np.arange(1, len(order) + 1) / len(order)
    if passed.any():
        out[order[:np.where(passed)[0].max() + 1]] = True
    return out


def participant_rotation_test(vectors, n_sur=49999, rng=None):
    """Diagnostic rotation calculation under an assumed uniform participant-direction null.

    Each complex input is one participant's equal-contact vector. Independent random rotations
    preserve every participant's magnitude—and therefore all within-participant dependence—while
    testing whether vector directions align consistently across participants. This is not the
    production 3D null: event selection within a finite SO-centered window can itself align
    participant directions, so pairing must be repeated inside a time-shift/block null.
    """
    vectors = np.asarray(vectors, complex)
    vectors = vectors[np.isfinite(vectors.real) & np.isfinite(vectors.imag)]
    if len(vectors) < 5:
        return None
    observed_vector = np.mean(vectors)
    observed = float(np.abs(observed_vector))
    rng = rng or np.random.RandomState(0)
    rotations = rng.uniform(-np.pi, np.pi, size=(int(n_sur), len(vectors)))
    null = np.abs(np.mean(vectors[None, :] * np.exp(1j * rotations), axis=1))
    return dict(
        n_participants=int(len(vectors)), group_R=observed,
        group_preferred_phase_deg=float(np.degrees(np.angle(observed_vector))),
        p=float((1 + np.sum(null >= observed)) / (len(null) + 1)),
        n_surrogates=int(len(null)),
        null_mean_R=float(np.mean(null)), null_sd_R=float(np.std(null, ddof=1)),
        method="independent rotation of one equal-contact complex vector per participant")


def pooled_endpoint_passes_qc(endpoint):
    """Whether one participant has prespecified support for a descriptive pooled-NREM vector."""
    endpoint = endpoint or {}
    per_channel = endpoint.get("per_channel") or []
    qualified = [
        value for value in per_channel
        if value.get("n", 0) >= MIN_EVENTS
        and value.get("valid_pooled_nrem_seconds", 0)
        >= MIN_POOLED_VALID_S_PER_CONTACT
        and value.get("valid_pooled_nrem_fraction", 0)
        >= MIN_POOLED_NREM_COVERAGE
    ]
    included_events = int(sum(value.get("n", 0) for value in qualified))
    return bool(
        len(qualified) >= MIN_POOLED_CONTACTS
        and included_events >= MIN_POOLED_EVENTS
        and endpoint.get("n_channels_tested") == len(per_channel) == len(qualified)
        and endpoint.get("n_spindle_events") == included_events
    )


def so_event_candidates(so, sf):
    """Return complete SO-cycle candidates as ``(trough_sample, peak_to_peak_amplitude)``."""
    so = np.asarray(so, float)
    crossings = np.where((so[:-1] >= 0) & (so[1:] < 0))[0] + 1
    candidates = []
    for a, b in zip(crossings[:-1], crossings[1:]):
        dur = (b - a) / sf
        if not (SO_DUR[0] <= dur <= SO_DUR[1]):
            continue
        seg = so[a:b]
        trough = a + int(np.argmin(seg))
        candidates.append((trough, float(seg.max() - seg.min())))
    return np.asarray(candidates, float).reshape(-1, 2)


def select_so_events(candidates, percentile=EVENT_PERCENTILE):
    """Apply one amplitude threshold to a candidate collection."""
    candidates = np.asarray(candidates, float).reshape(-1, 2)
    if not len(candidates):
        return np.array([], int)
    threshold = np.percentile(candidates[:, 1], percentile)
    return candidates[candidates[:, 1] >= threshold, 0].astype(int)


def detect_so_events(so, sf):
    """Compatibility helper for one continuous signal."""
    return select_so_events(so_event_candidates(so, sf))


def spindle_rms(sp, sf):
    n_rms = max(1, int(round(0.2 * sf)))
    return np.sqrt(np.convolve(np.asarray(sp, float) ** 2,
                               np.ones(n_rms) / n_rms, mode="same"))


def spindle_events_from_rms(rms, sf, threshold, valid=None, return_amplitudes=False):
    """Detect duration-qualified spindle events at one precomputed channel-night threshold."""
    rms = np.asarray(rms, float)
    valid = np.isfinite(rms) if valid is None else (np.asarray(valid, bool) & np.isfinite(rms))
    above = valid & (rms > threshold)
    changes = np.diff(np.pad(above.astype(int), (1, 1)))
    starts, stops = np.where(changes == 1)[0], np.where(changes == -1)[0]
    peaks, amplitudes = [], []
    for a, b in zip(starts, stops):
        dur = (b - a) / sf
        if SP_DUR[0] <= dur <= SP_DUR[1]:
            peak = a + int(np.argmax(rms[a:b]))
            peaks.append(peak)
            amplitudes.append(float(rms[peak]))
    peaks = np.asarray(peaks, int)
    amplitudes = np.asarray(amplitudes, float)
    return (peaks, amplitudes) if return_amplitudes else peaks


def detect_spindle_events(sp, sf):
    """Compatibility helper for one continuous signal."""
    rms = spindle_rms(sp, sf)
    threshold = np.percentile(rms[np.isfinite(rms)], EVENT_PERCENTILE)
    return spindle_events_from_rms(rms, sf, threshold)


def pair_one_spindle_per_so(so_indices, spindle_indices, spindle_amplitudes,
                            max_distance_samples):
    """Assign each spindle to its nearest SO and retain one maximum-amplitude spindle per SO."""
    so_indices = np.asarray(so_indices, int)
    spindle_indices = np.asarray(spindle_indices, int)
    spindle_amplitudes = np.asarray(spindle_amplitudes, float)
    if not len(so_indices) or not len(spindle_indices):
        return np.array([], int), np.array([], int)
    order = np.argsort(so_indices)
    so_sorted = so_indices[order]
    insert = np.searchsorted(so_sorted, spindle_indices)
    assigned = np.empty(len(spindle_indices), int)
    distances = np.empty(len(spindle_indices), int)
    for i, (sp, pos) in enumerate(zip(spindle_indices, insert)):
        options = [j for j in (pos - 1, pos) if 0 <= j < len(so_sorted)]
        nearest = min(options, key=lambda j: abs(sp - so_sorted[j]))
        assigned[i] = nearest
        distances[i] = abs(sp - so_sorted[nearest])
    keep = distances <= int(max_distance_samples)
    selected_sp, selected_so = [], []
    for so_pos in np.unique(assigned[keep]):
        members = np.where(keep & (assigned == so_pos))[0]
        best = members[int(np.argmax(spindle_amplitudes[members]))]
        selected_sp.append(int(spindle_indices[best]))
        selected_so.append(int(so_sorted[so_pos]))
    return np.asarray(selected_sp, int), np.asarray(selected_so, int)


def same_stage_pair_mask(spindle_epochs, so_epochs, keep_epochs):
    """Both members of a paired event must belong to the requested stage."""
    keep = set(keep_epochs)
    return np.asarray([
        spindle_epoch in keep and so_epoch in keep
        for spindle_epoch, so_epoch in zip(spindle_epochs, so_epochs)
    ], bool)


def _nrem_grid(lab, total_samples, fs):
    out = np.zeros(total_samples, bool)
    for epoch in np.where((lab == "N2") | (lab == "N3") | (lab == "NREM"))[0]:
        a = int(round(epoch * EPOCH * fs))
        b = min(total_samples, int(round((epoch + 1) * EPOCH * fs)))
        out[a:b] = True
    return out


def channel_night_events(rms, so_phase, candidates, lab):
    """Threshold and pair one channel's events after the complete night and stages are known."""
    rms = np.asarray(rms, float)
    so_phase = np.asarray(so_phase, float)
    nrem = _nrem_grid(lab, len(rms), EVENT_FS)
    valid_rms = nrem & np.isfinite(rms) & np.isfinite(so_phase)
    nrem_samples = int(nrem.sum())
    valid_nrem_seconds = float(valid_rms.sum() / EVENT_FS)
    valid_nrem_fraction = (
        float(valid_rms.sum() / nrem_samples) if nrem_samples else 0.0)
    pooled_nrem_qc_pass = bool(
        valid_nrem_seconds >= MIN_POOLED_VALID_S_PER_CONTACT
        and valid_nrem_fraction >= MIN_POOLED_NREM_COVERAGE)
    if valid_rms.sum() < EVENT_FS * 120:
        return dict(phases=np.array([]), indices=np.array([], int),
                    epochs=np.array([], int), so_epochs=np.array([], int),
                    n_so=0, n_spindle=0,
                    valid_nrem_seconds=valid_nrem_seconds,
                    valid_nrem_fraction=valid_nrem_fraction,
                    pooled_nrem_qc_pass=pooled_nrem_qc_pass)

    candidates = np.asarray(candidates, float).reshape(-1, 2)
    if len(candidates):
        idx = candidates[:, 0].astype(int)
        keep = ((idx >= 0) & (idx < len(nrem)))
        candidates = candidates[keep]
        idx = candidates[:, 0].astype(int)
        candidates = candidates[nrem[idx] & np.isfinite(so_phase[idx])]
    so_idx = select_so_events(candidates)

    sp_threshold = float(np.percentile(rms[valid_rms], EVENT_PERCENTILE))
    sp_idx, sp_amp = spindle_events_from_rms(
        rms, EVENT_FS, sp_threshold, valid=valid_rms, return_amplitudes=True)
    selected_sp, selected_so = pair_one_spindle_per_so(
        so_idx, sp_idx, sp_amp, max_distance_samples=int(PAIR_WINDOW_S * EVENT_FS))
    phases = so_phase[selected_sp] if len(selected_sp) else np.array([])
    epochs = (selected_sp / (EVENT_FS * EPOCH)).astype(int)
    so_epochs = (selected_so / (EVENT_FS * EPOCH)).astype(int)
    return dict(phases=phases, indices=selected_sp, epochs=epochs, so_epochs=so_epochs,
                n_so=int(len(so_idx)), n_spindle=int(len(sp_idx)),
                spindle_threshold=sp_threshold,
                valid_nrem_seconds=valid_nrem_seconds,
                valid_nrem_fraction=valid_nrem_fraction,
                pooled_nrem_qc_pass=pooled_nrem_qc_pass)


def run(n, hours, force=False, run_id=None, tree_digest=None):
    name = f"HUP{n}_phaseII"
    fp = os.path.join(OUT, f"{name}.json")
    if os.path.exists(fp) and not force:
        raise RuntimeError(
            f"{fp} cannot be reused outside a validated complete run; rerun with --force")

    s = sess()
    ds = s.open_dataset(name)
    labels = ds.get_channel_labels()
    lab_idx = {label: i for i, label in enumerate(labels)}
    details = ds.get_time_series_details(labels[0])
    sf = float(details.sample_rate)
    total_h = (getattr(details, "duration", 0) or 0) / 3.6e9
    ctx = cortical_channels(labels)
    if len(ctx) < 3:
        atomic_json_dump(
            dict(subject=name, status="skip", analysis_version=ANALYSIS_VERSION,
                 run_id=run_id, source_tree_sha256=tree_digest,
                 reason=f"only {len(ctx)} cortical channels"), fp)
        print(f"[{name}] SKIP", flush=True)
        return "skip"
    night = find_night(ds, lab_idx, ctx[0], sf, total_h, required_h=hours)
    if night is None:
        atomic_json_dump(
            dict(subject=name, status="skip", analysis_version=ANALYSIS_VERSION,
                 run_id=run_id, source_tree_sha256=tree_digest,
                 reason="no candidate night"), fp)
        print(f"[{name}] SKIP no night", flush=True)
        return "skip"

    idx = [lab_idx[c] for c in ctx]
    sos_so = band_sos(SO_BAND, sf, 3)
    sos_sp = band_sos(SPINDLE_BAND, sf)
    print(f"[{name}] {sf:.0f} Hz | {len(ctx)} cortical ch | SO {SO_BAND[0]}-"
          f"{SO_BAND[1]} Hz | spindle {SPINDLE_BAND[0]}-{SPINDLE_BAND[1]} Hz | "
          f"streaming {hours} h", flush=True)

    total_s = int(hours * 3600)
    n_ep = int(total_s // EPOCH)
    n_event = int(np.ceil(total_s * EVENT_FS))
    ep_dr_ch = np.full((len(ctx), n_ep), np.nan)
    ep_swa_ch = np.full((len(ctx), n_ep), np.nan)
    ep_clean_ch = np.full((len(ctx), n_ep), np.nan)
    rms_night = {c: np.full(n_event, np.nan) for c in ctx}
    phase_night = {c: np.full(n_event, np.nan) for c in ctx}
    candidates = {c: [] for c in ctx}
    failed_chunks = []

    t = 0.0
    while t < total_s:
        dur = min(CHUNK_S, total_s - t)
        pull_start = max(0.0, t - FILTER_EDGE_S)
        pull_stop = min(float(total_s), t + dur + FILTER_EDGE_S)
        pull_dur = pull_stop - pull_start
        try:
            d = pull_continuous(ds, idx, night + pull_start, pull_dur)
        except Exception as exc:
            failed_chunks.append(dict(start_s=float(t), duration_s=float(dur),
                                      error=f"{type(exc).__name__}: {exc}"))
            t += dur
            continue

        core_a = int(round((t - pull_start) * sf))
        core_n = min(int(round(dur * sf)), max(0, len(d) - core_a))
        core_b = core_a + core_n
        processed, clean_masks, measured_masks = [], [], []
        for ci, c in enumerate(ctx):
            filled, measured = prepare_continuous_signal(d[:, ci], sf)
            x = notch(signal.detrend(filled), sf)
            clean = ied_clean_mask(x, sf, pad_s=IED_PAD_S) & measured
            processed.append(x)
            clean_masks.append(clean)
            measured_masks.append(measured)
            if np.std(x) < 1e-9:
                continue
            so = signal.sosfiltfilt(sos_so, x)
            phase = np.angle(signal.hilbert(so))
            sp = signal.sosfiltfilt(sos_sp, x)
            rms = spindle_rms(sp, sf)

            cand = so_event_candidates(so, sf)
            if len(cand):
                tr = cand[:, 0].astype(int)
                keep = ((tr >= core_a) & (tr < core_b) & clean[tr])
                for tr_i, amp in cand[keep]:
                    global_idx = int(round(
                        (pull_start + float(tr_i) / sf) * EVENT_FS))
                    if 0 <= global_idx < n_event:
                        candidates[c].append((global_idx, float(amp)))

            g0 = int(round(t * EVENT_FS))
            n_out = min(int(round(dur * EVENT_FS)), n_event - g0)
            local = core_a + np.round(np.arange(n_out) * sf / EVENT_FS).astype(int)
            in_bounds = local < len(x)
            global_idx = g0 + np.arange(n_out)[in_bounds]
            local = local[in_bounds]
            valid = clean[local]
            rms_night[c][global_idx[valid]] = rms[local[valid]]
            phase_night[c][global_idx[valid]] = phase[local[valid]]

        # Robust, per-contact staging features; no voltage cancellation before power.
        ke = int(EPOCH * sf)
        for e in range(core_n // ke):
            gi = int((t + e * EPOCH) // EPOCH)
            if gi >= n_ep:
                break
            a = core_a + e * ke
            b = a + ke
            for ci, (x, clean, measured) in enumerate(zip(
                    processed, clean_masks, measured_masks)):
                seg = x[a:b]
                dr, swa_value = staging_epoch_features(
                    seg, clean[a:b], measured[a:b], sf)
                ep_dr_ch[ci, gi] = dr
                ep_swa_ch[ci, gi] = swa_value
                ep_clean_ch[ci, gi] = float(clean[a:b].mean())
        t += dur

    if failed_chunks:
        raise RuntimeError(f"{len(failed_chunks)} acquisition chunks failed")

    staging_observed_fraction = np.isfinite(ep_clean_ch).mean(axis=1)
    staging_clean_fraction = np.divide(
        np.nansum(ep_clean_ch, axis=1),
        np.isfinite(ep_clean_ch).sum(axis=1),
        out=np.zeros(len(ctx), float),
        where=np.isfinite(ep_clean_ch).sum(axis=1) > 0,
    )
    staging_candidates = (
        (staging_observed_fraction >= 0.80)
        & (staging_clean_fraction >= 0.80)
    )
    ep_dr, ep_swa, ep_clean, staging_contact_qc = aggregate_staging_features(
        ep_dr_ch, ep_swa_ch, ep_clean_ch, staging_candidates,
        min_contact_fraction_per_epoch=MIN_CONTACT_FRACTION_PER_BIN)
    ep = dict(dr=ep_dr, swa=ep_swa, clean=ep_clean)
    lab, nrem, sep = stage_epochs(ep)
    events = {
        c: channel_night_events(rms_night[c], phase_night[c], candidates[c], lab)
        for c in ctx
    }
    production = production_config(hours)
    rec = dict(
        subject=name, status="ok", analysis_version=ANALYSIS_VERSION,
        run_id=run_id, source_tree_sha256=tree_digest,
        sf=sf, hours=float(hours), cortical_chans=ctx, night_h=float(night / 3600),
        so_band_hz=list(SO_BAND), spindle_band_hz=list(SPINDLE_BAND),
        so_duration_s=list(SO_DUR), spindle_duration_s=list(SP_DUR),
        event_percentile=EVENT_PERCENTILE, event_sampling_hz=EVENT_FS,
        minimum_events_per_contact=MIN_EVENTS,
        threshold_scope=production["threshold_scope"],
        pairing=production["pairing"],
        pairing_details=(
            f"hybrid adaptation: duration-qualified spindles within +/-{PAIR_WINDOW_S:.1f} s; "
            "each spindle assigned to nearest SO; maximum-amplitude spindle retained per SO"),
        artifact_rejection=f"IED/high-amplitude mask padded +/-{IED_PAD_S:.1f} s",
        ied_mask_padding_s=IED_PAD_S,
        production_inference_enabled=PRODUCTION_3D_INFERENCE_ENABLED,
        filter_context_s=float(FILTER_EDGE_S),
        staging_selected_contact_mask=staging_contact_qc["selected_contact_mask"].tolist(),
        staging_contact_count=staging_contact_qc["contact_count"].tolist(),
        staging_n_selected_contacts=staging_contact_qc["n_selected_contacts"],
        staging_required_contact_count=staging_contact_qc["required_contact_count"],
        staging_swa_normalization=staging_contact_qc["swa_normalization"],
        n_nrem=int(nrem.sum()), n_N2=int((lab == "N2").sum()),
        n_N3=int((lab == "N3").sum()),
        gmm_separation=finite_float_or_none(sep),
        n_so_total=int(sum(v["n_so"] for v in events.values())))

    for stage in ("N2", "N3", "ALL"):
        keep_epochs = (set(np.where(lab == stage)[0].tolist()) if stage != "ALL"
                       else set(np.where(nrem)[0].tolist()))
        per_ch, n_ev = [], 0
        for c in ctx:
            if stage == "ALL" and not events[c]["pooled_nrem_qc_pass"]:
                continue
            e = np.asarray(events[c]["epochs"], int)
            so_e = np.asarray(events[c]["so_epochs"], int)
            p = np.asarray(events[c]["phases"], float)
            selected = same_stage_pair_mask(e, so_e, keep_epochs)
            phases = p[selected]
            R, rayleigh_p, mean_phase = rayleigh(phases)
            if not np.isfinite(R):
                continue
            n_selected = len(phases)
            n_ev += int(n_selected)
            r2_unbiased_iid = max(
                (n_selected * R ** 2 - 1.0) / max(n_selected - 1, 1), 0.0)
            per_ch.append(dict(
                ch=c, n=int(n_selected), R=R,
                valid_pooled_nrem_seconds=float(events[c]["valid_nrem_seconds"]),
                valid_pooled_nrem_fraction=float(events[c]["valid_nrem_fraction"]),
                R_bias_corrected_iid_diagnostic=float(np.sqrt(r2_unbiased_iid)),
                rayleigh_p_independence_assumption=rayleigh_p,
                preferred_phase_deg=float(np.degrees(mean_phase))))

        # The raw complex mean is the estimand. The usual algebraic "unbiased R" assumes iid event
        # phases and is therefore not used here; repeated/serially dependent SOs violate that
        # assumption. Participant rotations below condition on each participant-vector magnitude.
        channel_vectors = np.asarray([
            value["R"] * np.exp(1j * np.radians(value["preferred_phase_deg"]))
            for value in per_ch
        ], complex)
        participant_vector = (
            complex(np.mean(channel_vectors)) if len(channel_vectors)
            else complex(np.nan, np.nan))
        rec[stage] = dict(
            n_coupled_so_epochs=n_ev,
            n_spindle_events=n_ev,  # backward-compatible name; one event maximum per SO
            n_channels_tested=len(per_ch),
            n_channels_passing_coverage=int(sum(
                bool(value["pooled_nrem_qc_pass"]) for value in events.values())),
            n_channels_significant=None,
            multiplicity="none; contacts and events are not inferential units",
            channel_inference=(
                "descriptive only; Rayleigh p retained solely as an independence-assumption "
                "diagnostic and never used for significance"),
            participant_vector_real=(
                float(participant_vector.real) if len(channel_vectors) else None),
            participant_vector_imag=(
                float(participant_vector.imag) if len(channel_vectors) else None),
            participant_R=(
                float(abs(participant_vector)) if len(channel_vectors) else None),
            participant_preferred_phase_deg=(
                float(np.degrees(np.angle(participant_vector)))
                if len(channel_vectors) else None),
            median_R=(float(np.median([v["R"] for v in per_ch])) if per_ch else None),
            median_R_bias_corrected_iid_diagnostic=(
                float(np.median([
                    v["R_bias_corrected_iid_diagnostic"] for v in per_ch
                ])) if per_ch else None),
            mean_preferred_phase_deg=(
                float(np.degrees(np.angle(participant_vector)))
                if len(channel_vectors) else None),
            per_channel=per_ch)
        if stage == "ALL":
            rec[stage]["quality_gate"] = dict(
                passed=pooled_endpoint_passes_qc(rec[stage]),
                minimum_contacts=MIN_POOLED_CONTACTS,
                minimum_paired_events=MIN_POOLED_EVENTS,
                minimum_valid_nrem_seconds_per_contact=MIN_POOLED_VALID_S_PER_CONTACT,
                minimum_valid_nrem_fraction_per_contact=MIN_POOLED_NREM_COVERAGE)
            rec[stage]["inference_status"] = (
                "disabled: the finite SO-centered pairing window induces a common phase under "
                "independent event trains; a valid null must shift/block-resample complete spindle "
                "trains relative to SO and repeat pairing/contact aggregation")
        print(f"[{name}]   {stage:3s}: {len(per_ch)} descriptive channels | "
              f"{n_ev} paired SO epochs | participant R={rec[stage]['participant_R']}", flush=True)

    rec["endpoint_availability"] = {
        "N2": bool((rec.get("N2") or {}).get("n_channels_tested", 0)),
        "N3": bool((rec.get("N3") or {}).get("n_channels_tested", 0)),
        "ALL": pooled_endpoint_passes_qc(rec.get("ALL")),
    }
    if not rec["endpoint_availability"]["ALL"]:
        rec["status"] = "partial"
        all_endpoint = rec["ALL"]
        rec["reason"] = (
            "pooled-NREM QC failed: "
            f"{all_endpoint['n_channels_tested']}/{MIN_POOLED_CONTACTS} contacts with "
            f">={MIN_POOLED_VALID_S_PER_CONTACT:g} s and "
            f">={MIN_POOLED_NREM_COVERAGE:.0%} NREM coverage; "
            f"{all_endpoint['n_spindle_events']}/{MIN_POOLED_EVENTS} paired events")
    atomic_json_dump(rec, fp)
    return rec["status"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--subjects", default=",".join(map(str, COHORT)))
    parser.add_argument("--hours", type=float, default=7.0)
    parser.add_argument("--force", action="store_true",
                        help="overwrite legacy outputs generated by an older estimator")
    args = parser.parse_args()
    os.makedirs(OUT, exist_ok=True)
    requested_n = [int(value) for value in args.subjects.split(",") if value.strip()]
    requested = [f"HUP{n}_phaseII" for n in requested_n]
    config = production_config(args.hours)
    if not args.force and validated_complete_run_exists(
            OUT, pipeline="event_3D_by_stage", requested=requested, config=config,
            suffix=".json", require_current_source_tree=True):
        print(f"validated existing complete 3D run ({len(requested)} subjects)", flush=True)
        return
    run_id = start_run_manifest(
        OUT, pipeline="event_3D_by_stage", requested=requested, config=config)
    tree_digest = source_tree_sha256(ROOT)
    completed, skipped, failed = [], [], []
    for n in requested_n:
        name = f"HUP{n}_phaseII"
        try:
            status = run(
                n, args.hours, force=args.force, run_id=run_id, tree_digest=tree_digest)
            if status not in ("ok", "cached"):
                record = json.load(open(os.path.join(OUT, f"{name}.json")))
                skipped.append(dict(subject=name, reason=record.get("reason", "unspecified")))
            else:
                completed.append(name)
        except Exception as exc:
            print(f"[{name}] ERROR {type(exc).__name__}: {exc}", flush=True)
            traceback.print_exc()
            failed.append(dict(subject=name, error=f"{type(exc).__name__}: {exc}"))
    write_run_manifest(
        OUT, pipeline="event_3D_by_stage",
        requested=requested, completed=completed,
        skipped=skipped, failed=failed,
        config=config, run_id=run_id,
        result_files_sha256={
            subject: file_sha256(os.path.join(OUT, f"{subject}.json"))
            for subject in completed + [value["subject"] for value in skipped]
        })
    print(f"\n{len(completed)} subjects -> {OUT}", flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
