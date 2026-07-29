"""Pure estimators for the descriptive SO-to-spindle 3D endpoint.

This module intentionally contains no portal, cache, manifest, or plotting
code. Offline neutral-cache reconstruction and historical compatibility
imports share the same event detection and pairing functions from here.

The finite +/-2-second pairing window can align phases even for independent
event trains.  Rayleigh and participant-rotation calculations are therefore
diagnostics only; production inference remains disabled.
"""
from __future__ import annotations

import numpy as np


STAGE_EPOCH_S = 30.0
SO_BAND = (0.16, 1.25)
SPINDLE_BAND = (12.0, 16.0)
SP_DUR = (0.5, 3.0)
SO_DUR = (0.8, 2.0)
SPINDLE_RMS_WINDOW_S = 0.2
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


def conservative_resampled_interval(
        start_sample, stop_sample_exclusive, source_fs, target_fs, *,
        offset_s=0.0):
    """Map a half-open interval without shrinking either boundary.

    Event centres may be mapped to the nearest target sample, but event extents
    must use floor(start) and ceil(stop). Otherwise downsampling can make a
    cycle that crosses a stage boundary appear wholly inside one stage.
    """
    source_fs = float(source_fs)
    target_fs = float(target_fs)
    start_sample = float(start_sample)
    stop_sample_exclusive = float(stop_sample_exclusive)
    if (
        source_fs <= 0
        or target_fs <= 0
        or not np.isfinite(
            [source_fs, target_fs, start_sample,
             stop_sample_exclusive, offset_s]
        ).all()
        or stop_sample_exclusive <= start_sample
    ):
        raise ValueError(
            "resampled interval geometry must be finite and ordered")

    def snap_theoretical_integer(value):
        nearest = float(np.rint(value))
        tolerance = (
            64 * np.finfo(float).eps * max(1.0, abs(value))
        )
        return nearest if abs(value - nearest) <= tolerance else value

    start = snap_theoretical_integer(
        (float(offset_s) + start_sample / source_fs) * target_fs)
    stop = snap_theoretical_integer(
        (float(offset_s) + stop_sample_exclusive / source_fs) * target_fs)
    return int(np.floor(start)), int(np.ceil(stop))


def retain_complete_clean_event_3d_candidates(
        candidates, valid_sample_mask, core_start, core_stop):
    """Keep core-owned SOs only when their complete cycle is artifact-valid."""
    values = np.asarray(candidates, float)
    if values.size == 0:
        return np.empty((0, 4), float)
    if values.ndim != 2 or values.shape[1] != 4:
        raise ValueError("3D SO candidates must have four columns")
    valid = np.asarray(valid_sample_mask, bool)
    keep = np.zeros(len(values), bool)
    for index, (
            trough, _amplitude, cycle_start,
            cycle_stop) in enumerate(values):
        trough = int(trough)
        cycle_start = int(cycle_start)
        cycle_stop = int(cycle_stop)
        keep[index] = (
            int(core_start) <= trough < int(core_stop)
            and 0 <= cycle_start <= trough < cycle_stop <= len(valid)
            and bool(valid[cycle_start:cycle_stop].all())
        )
    return values[keep]


def rayleigh(phases):
    """Return resultant length, diagnostic Rayleigh p, and mean phase."""
    phases = np.asarray(phases, float)
    n_events = len(phases)
    if n_events < MIN_EVENTS:
        return np.nan, np.nan, np.nan
    cosine = np.mean(np.cos(phases))
    sine = np.mean(np.sin(phases))
    resultant = float(np.hypot(cosine, sine))
    z_value = n_events * resultant ** 2
    p_value = np.exp(-z_value)
    if n_events < 50:
        p_value *= (
            1
            + (2 * z_value - z_value ** 2) / (4 * n_events)
            - (
                24 * z_value
                - 132 * z_value ** 2
                + 76 * z_value ** 3
                - 9 * z_value ** 4
            )
            / (288 * n_events ** 2)
        )
    p_value = float(
        np.clip(p_value, np.finfo(float).tiny, 1.0))
    return (
        resultant,
        p_value,
        float(np.arctan2(sine, cosine)),
    )


def bh_fdr(p_values, alpha=0.05):
    """Return Benjamini-Hochberg decisions in the input order."""
    values = np.asarray(p_values, float)
    decisions = np.zeros(len(values), bool)
    finite = np.isfinite(values)
    if not finite.any():
        return decisions
    finite_indices = np.where(finite)[0]
    order = finite_indices[np.argsort(values[finite_indices])]
    passed = (
        values[order]
        <= alpha * np.arange(1, len(order) + 1) / len(order)
    )
    if passed.any():
        decisions[order[:np.where(passed)[0].max() + 1]] = True
    return decisions


def participant_rotation_test(vectors, n_sur=49999, rng=None):
    """Diagnostic test under an assumed uniform participant-direction null."""
    vectors = np.asarray(vectors, complex)
    vectors = vectors[
        np.isfinite(vectors.real) & np.isfinite(vectors.imag)]
    if len(vectors) < 5:
        return None
    observed_vector = np.mean(vectors)
    observed = float(np.abs(observed_vector))
    rng = rng or np.random.RandomState(0)
    rotations = rng.uniform(
        -np.pi, np.pi, size=(int(n_sur), len(vectors)))
    null = np.abs(
        np.mean(vectors[None, :] * np.exp(1j * rotations), axis=1))
    return {
        "n_participants": int(len(vectors)),
        "group_R": observed,
        "group_preferred_phase_deg": float(
            np.degrees(np.angle(observed_vector))),
        "p": float((1 + np.sum(null >= observed)) / (len(null) + 1)),
        "n_surrogates": int(len(null)),
        "null_mean_R": float(np.mean(null)),
        "null_sd_R": float(np.std(null, ddof=1)),
        "method": (
            "independent rotation of one equal-contact complex vector "
            "per participant"),
    }


def circular_summary(phases):
    """Return a descriptive complex mean without an iid-event p value."""
    phases = np.asarray(phases, float)
    phases = phases[np.isfinite(phases)]
    if not len(phases):
        return None
    vector = complex(np.mean(np.exp(1j * phases)))
    return {
        "n_paired_events": int(len(phases)),
        "vector_real": float(vector.real),
        "vector_imag": float(vector.imag),
        "R": float(abs(vector)),
        "preferred_phase_rad": float(np.angle(vector)),
        "preferred_phase_deg": float(
            np.degrees(np.angle(vector))),
    }


def pooled_endpoint_passes_qc(endpoint):
    """Whether a participant supports a descriptive pooled-NREM vector."""
    endpoint = endpoint or {}
    per_channel = endpoint.get("per_channel") or []
    qualified = [
        value
        for value in per_channel
        if (
            value.get("n", 0) >= MIN_EVENTS
            and value.get("valid_pooled_nrem_seconds", 0)
            >= MIN_POOLED_VALID_S_PER_CONTACT
            and value.get("valid_pooled_nrem_fraction", 0)
            >= MIN_POOLED_NREM_COVERAGE
        )
    ]
    included_events = int(sum(
        value.get("n", 0) for value in qualified))
    return bool(
        len(qualified) >= MIN_POOLED_CONTACTS
        and included_events >= MIN_POOLED_EVENTS
        and endpoint.get("n_channels_tested")
        == len(per_channel)
        == len(qualified)
        and endpoint.get("n_spindle_events") == included_events
    )


def so_event_candidates(so, sampling_hz):
    """Return complete SO cycles and their half-open sample bounds.

    Columns are ``trough, peak_to_peak_amplitude, cycle_start,
    cycle_stop_exclusive``.
    """
    values = np.asarray(so, float)
    crossings = np.where(
        (values[:-1] >= 0) & (values[1:] < 0))[0] + 1
    candidates = []
    for start, stop in zip(crossings[:-1], crossings[1:]):
        duration_s = (stop - start) / float(sampling_hz)
        if not (SO_DUR[0] <= duration_s <= SO_DUR[1]):
            continue
        cycle = values[start:stop]
        trough = start + int(np.argmin(cycle))
        candidates.append((
            trough,
            float(cycle.max() - cycle.min()),
            start,
            stop,
        ))
    return np.asarray(candidates, float).reshape(-1, 4)


def select_so_events(candidates, percentile=EVENT_PERCENTILE):
    """Apply one amplitude threshold to a candidate collection."""
    candidates = np.asarray(candidates, float)
    if candidates.size == 0:
        return np.array([], int)
    if candidates.ndim != 2 or candidates.shape[1] < 2:
        raise ValueError(
            "SO candidates must contain trough and amplitude columns")
    threshold = np.percentile(candidates[:, 1], percentile)
    return candidates[
        candidates[:, 1] >= threshold, 0].astype(int)


def spindle_rms(spindle_band_signal, sampling_hz):
    """Return the fixed 200-ms RMS envelope."""
    window = max(
        1,
        int(round(SPINDLE_RMS_WINDOW_S * float(sampling_hz))),
    )
    values = np.asarray(spindle_band_signal, float)
    return np.sqrt(np.convolve(
        values ** 2,
        np.ones(window) / window,
        mode="same",
    ))


def spindle_events_from_rms(
        rms, sampling_hz, threshold, valid=None,
        return_amplitudes=False, return_bounds=False,
        require_complete_valid_extent=False):
    """Detect duration-qualified spindle events at one night-wide threshold."""
    rms = np.asarray(rms, float)
    if valid is None:
        valid = np.isfinite(rms)
    else:
        valid = np.asarray(valid, bool) & np.isfinite(rms)
    above = valid & (rms > threshold)
    changes = np.diff(np.pad(above.astype(int), (1, 1)))
    starts = np.where(changes == 1)[0]
    stops = np.where(changes == -1)[0]
    peaks = []
    amplitudes = []
    event_starts = []
    event_stops = []
    for start, stop in zip(starts, stops):
        if require_complete_valid_extent and (
            start == 0
            or stop == len(rms)
            or not valid[start - 1]
            or not valid[stop]
        ):
            continue
        duration_s = (stop - start) / float(sampling_hz)
        if not (SP_DUR[0] <= duration_s <= SP_DUR[1]):
            continue
        peak = start + int(np.argmax(rms[start:stop]))
        peaks.append(peak)
        amplitudes.append(float(rms[peak]))
        event_starts.append(int(start))
        event_stops.append(int(stop))
    peaks = np.asarray(peaks, int)
    amplitudes = np.asarray(amplitudes, float)
    event_starts = np.asarray(event_starts, int)
    event_stops = np.asarray(event_stops, int)
    if return_amplitudes and return_bounds:
        return peaks, amplitudes, event_starts, event_stops
    if return_bounds:
        return peaks, event_starts, event_stops
    if return_amplitudes:
        return peaks, amplitudes
    return peaks


def pair_one_spindle_per_so(
        so_indices, spindle_indices, spindle_amplitudes,
        max_distance_samples):
    """Assign each spindle to its nearest SO and retain one maximum per SO."""
    so_indices = np.asarray(so_indices, int)
    spindle_indices = np.asarray(spindle_indices, int)
    spindle_amplitudes = np.asarray(spindle_amplitudes, float)
    if spindle_indices.shape != spindle_amplitudes.shape:
        raise ValueError(
            "spindle indices and amplitudes must align")
    if not len(so_indices) or not len(spindle_indices):
        return np.array([], int), np.array([], int)
    so_sorted = np.sort(so_indices)
    insertion = np.searchsorted(so_sorted, spindle_indices)
    assigned = np.empty(len(spindle_indices), int)
    distances = np.empty(len(spindle_indices), int)
    for index, (spindle, position) in enumerate(
            zip(spindle_indices, insertion)):
        options = [
            candidate
            for candidate in (position - 1, position)
            if 0 <= candidate < len(so_sorted)
        ]
        nearest = min(
            options,
            key=lambda candidate: abs(
                spindle - so_sorted[candidate]),
        )
        assigned[index] = nearest
        distances[index] = abs(spindle - so_sorted[nearest])
    within = distances <= int(max_distance_samples)
    selected_spindles = []
    selected_sos = []
    for so_position in np.unique(assigned[within]):
        members = np.where(
            within & (assigned == so_position))[0]
        best = members[
            int(np.argmax(spindle_amplitudes[members]))]
        selected_spindles.append(int(spindle_indices[best]))
        selected_sos.append(int(so_sorted[so_position]))
    return (
        np.asarray(selected_spindles, int),
        np.asarray(selected_sos, int),
    )


def intervals_wholly_inside(mask, starts, stops):
    """Whether each half-open interval is valid and entirely inside ``mask``."""
    mask = np.asarray(mask, bool)
    starts = np.asarray(starts, int)
    stops = np.asarray(stops, int)
    if starts.shape != stops.shape:
        raise ValueError("event starts and stops must align")
    in_bounds = (
        (starts >= 0)
        & (stops > starts)
        & (stops <= len(mask))
    )
    result = np.zeros(starts.shape, bool)
    if in_bounds.any():
        cumulative = np.concatenate((
            [0],
            np.cumsum(mask, dtype=np.int64),
        ))
        lengths = stops[in_bounds] - starts[in_bounds]
        result[in_bounds] = (
            cumulative[stops[in_bounds]]
            - cumulative[starts[in_bounds]]
        ) == lengths
    return result


def epoch_set_sample_mask(
        keep_epochs, total_samples, sampling_hz=EVENT_FS):
    """Expand a set of 30-second epoch indices to a sample mask."""
    mask = np.zeros(int(total_samples), bool)
    samples_per_epoch = int(
        round(STAGE_EPOCH_S * float(sampling_hz)))
    for epoch in set(keep_epochs):
        start = int(epoch) * samples_per_epoch
        stop = min(len(mask), start + samples_per_epoch)
        if 0 <= start < stop:
            mask[start:stop] = True
    return mask


def _nrem_grid(labels, total_samples, sampling_hz):
    labels = np.asarray(labels).astype(str)
    keep = set(np.where(
        np.isin(labels, ("NREM", "N2", "N3")))[0])
    return epoch_set_sample_mask(
        keep, total_samples, sampling_hz)


def _empty_event_record(
        total_samples, valid_nrem_seconds,
        valid_nrem_fraction, pooled_nrem_qc_pass):
    return {
        "phases": np.array([]),
        "indices": np.array([], int),
        "epochs": np.array([], int),
        "so_epochs": np.array([], int),
        "spindle_start_samples": np.array([], int),
        "spindle_stop_samples_exclusive": np.array([], int),
        "so_start_samples": np.array([], int),
        "so_stop_samples_exclusive": np.array([], int),
        "eligible_so_indices": np.array([], int),
        "eligible_so_start_samples": np.array([], int),
        "eligible_so_stop_samples_exclusive": np.array([], int),
        "eligible_spindle_indices": np.array([], int),
        "eligible_spindle_amplitudes": np.array([], float),
        "eligible_spindle_start_samples": np.array([], int),
        "eligible_spindle_stop_samples_exclusive": np.array([], int),
        "eligible_spindle_phases": np.array([], float),
        "total_samples": int(total_samples),
        "n_so": 0,
        "n_spindle": 0,
        "valid_nrem_seconds": valid_nrem_seconds,
        "valid_nrem_fraction": valid_nrem_fraction,
        "pooled_nrem_qc_pass": pooled_nrem_qc_pass,
    }


def channel_night_events(rms, so_phase, candidates, labels):
    """Threshold and pair complete, valid pooled-NREM events for one contact."""
    rms = np.asarray(rms, float)
    so_phase = np.asarray(so_phase, float)
    if rms.shape != so_phase.shape or rms.ndim != 1:
        raise ValueError(
            "RMS and SO phase must be aligned one-dimensional arrays")
    nrem = _nrem_grid(labels, len(rms), EVENT_FS)
    jointly_valid = np.isfinite(rms) & np.isfinite(so_phase)
    valid_nrem = nrem & jointly_valid
    nrem_samples = int(nrem.sum())
    valid_nrem_seconds = float(
        valid_nrem.sum() / EVENT_FS)
    valid_nrem_fraction = (
        float(valid_nrem.sum() / nrem_samples)
        if nrem_samples
        else 0.0
    )
    pooled_nrem_qc_pass = bool(
        valid_nrem_seconds
        >= MIN_POOLED_VALID_S_PER_CONTACT
        and valid_nrem_fraction
        >= MIN_POOLED_NREM_COVERAGE
    )
    if valid_nrem.sum() < EVENT_FS * 120:
        return _empty_event_record(
            len(rms),
            valid_nrem_seconds,
            valid_nrem_fraction,
            pooled_nrem_qc_pass,
        )

    candidates = np.asarray(candidates, float)
    if candidates.size == 0:
        candidates = np.empty((0, 4), float)
    if candidates.ndim != 2 or candidates.shape[1] != 4:
        raise ValueError(
            "3D SO candidates require trough, amplitude, "
            "cycle start, and exclusive cycle stop")
    if len(candidates):
        troughs = candidates[:, 0].astype(int)
        starts = candidates[:, 2].astype(int)
        stops = candidates[:, 3].astype(int)
        keep = (
            (troughs >= 0)
            & (troughs < len(nrem))
            & (starts <= troughs)
            & (troughs < stops)
            & intervals_wholly_inside(
                valid_nrem, starts, stops)
        )
        candidates = candidates[keep]
    if len(candidates):
        so_threshold = float(np.percentile(
            candidates[:, 1], EVENT_PERCENTILE))
        selected_candidates = candidates[
            candidates[:, 1] >= so_threshold]
    else:
        selected_candidates = candidates
    so_indices = selected_candidates[:, 0].astype(int)

    spindle_threshold = float(np.percentile(
        rms[valid_nrem], EVENT_PERCENTILE))
    (
        spindle_indices,
        spindle_amplitudes,
        spindle_starts,
        spindle_stops,
    ) = spindle_events_from_rms(
        rms,
        EVENT_FS,
        spindle_threshold,
        valid=jointly_valid,
        return_amplitudes=True,
        return_bounds=True,
        require_complete_valid_extent=True,
    )
    complete_nrem = intervals_wholly_inside(
        nrem, spindle_starts, spindle_stops)
    spindle_indices = spindle_indices[complete_nrem]
    spindle_amplitudes = spindle_amplitudes[complete_nrem]
    spindle_starts = spindle_starts[complete_nrem]
    spindle_stops = spindle_stops[complete_nrem]

    selected_spindles, selected_sos = pair_one_spindle_per_so(
        so_indices,
        spindle_indices,
        spindle_amplitudes,
        max_distance_samples=int(PAIR_WINDOW_S * EVENT_FS),
    )
    spindle_position = {
        int(peak): index
        for index, peak in enumerate(spindle_indices)
    }
    so_position = {
        int(trough): index
        for index, trough in enumerate(so_indices)
    }
    selected_spindle_positions = np.asarray([
        spindle_position[int(peak)]
        for peak in selected_spindles
    ], int)
    selected_so_positions = np.asarray([
        so_position[int(trough)]
        for trough in selected_sos
    ], int)
    phases = (
        so_phase[selected_spindles]
        if len(selected_spindles)
        else np.array([])
    )
    return {
        "phases": phases,
        "indices": selected_spindles,
        "epochs": (
            selected_spindles
            / (EVENT_FS * STAGE_EPOCH_S)
        ).astype(int),
        "so_epochs": (
            selected_sos
            / (EVENT_FS * STAGE_EPOCH_S)
        ).astype(int),
        "spindle_start_samples": spindle_starts[
            selected_spindle_positions],
        "spindle_stop_samples_exclusive": spindle_stops[
            selected_spindle_positions],
        "so_start_samples": selected_candidates[
            selected_so_positions, 2].astype(int),
        "so_stop_samples_exclusive": selected_candidates[
            selected_so_positions, 3].astype(int),
        "eligible_so_indices": so_indices,
        "eligible_so_start_samples": selected_candidates[
            :, 2].astype(int),
        "eligible_so_stop_samples_exclusive": selected_candidates[
            :, 3].astype(int),
        "eligible_spindle_indices": spindle_indices,
        "eligible_spindle_amplitudes": spindle_amplitudes,
        "eligible_spindle_start_samples": spindle_starts,
        "eligible_spindle_stop_samples_exclusive": spindle_stops,
        "eligible_spindle_phases": so_phase[spindle_indices],
        "total_samples": int(len(rms)),
        "n_so": int(len(so_indices)),
        "n_spindle": int(len(spindle_indices)),
        "spindle_threshold": spindle_threshold,
        "valid_nrem_seconds": valid_nrem_seconds,
        "valid_nrem_fraction": valid_nrem_fraction,
        "pooled_nrem_qc_pass": pooled_nrem_qc_pass,
    }


def stage_event_pairs(record, keep_epochs):
    """Restrict complete event extents to a stage, then repeat pairing."""
    total_samples = int(record["total_samples"])
    stage_mask = epoch_set_sample_mask(
        keep_epochs, total_samples, EVENT_FS)
    so_keep = intervals_wholly_inside(
        stage_mask,
        record["eligible_so_start_samples"],
        record["eligible_so_stop_samples_exclusive"],
    )
    spindle_keep = intervals_wholly_inside(
        stage_mask,
        record["eligible_spindle_start_samples"],
        record["eligible_spindle_stop_samples_exclusive"],
    )
    so_indices = np.asarray(
        record["eligible_so_indices"], int)[so_keep]
    so_starts = np.asarray(
        record["eligible_so_start_samples"], int)[so_keep]
    so_stops = np.asarray(
        record["eligible_so_stop_samples_exclusive"], int)[so_keep]
    spindle_indices = np.asarray(
        record["eligible_spindle_indices"], int)[spindle_keep]
    spindle_amplitudes = np.asarray(
        record["eligible_spindle_amplitudes"], float)[spindle_keep]
    spindle_starts = np.asarray(
        record["eligible_spindle_start_samples"], int)[spindle_keep]
    spindle_stops = np.asarray(
        record["eligible_spindle_stop_samples_exclusive"], int)[spindle_keep]
    spindle_phases = np.asarray(
        record["eligible_spindle_phases"], float)[spindle_keep]

    selected_spindles, selected_sos = pair_one_spindle_per_so(
        so_indices,
        spindle_indices,
        spindle_amplitudes,
        max_distance_samples=int(PAIR_WINDOW_S * EVENT_FS),
    )
    spindle_position = {
        int(peak): index
        for index, peak in enumerate(spindle_indices)
    }
    so_position = {
        int(trough): index
        for index, trough in enumerate(so_indices)
    }
    selected_spindle_positions = np.asarray([
        spindle_position[int(peak)]
        for peak in selected_spindles
    ], int)
    selected_so_positions = np.asarray([
        so_position[int(trough)]
        for trough in selected_sos
    ], int)
    return {
        "phases": spindle_phases[selected_spindle_positions],
        "indices": selected_spindles,
        "so_indices": selected_sos,
        "epochs": (
            selected_spindles
            / (EVENT_FS * STAGE_EPOCH_S)
        ).astype(int),
        "so_epochs": (
            selected_sos
            / (EVENT_FS * STAGE_EPOCH_S)
        ).astype(int),
        "spindle_start_samples": spindle_starts[
            selected_spindle_positions],
        "spindle_stop_samples_exclusive": spindle_stops[
            selected_spindle_positions],
        "so_start_samples": so_starts[selected_so_positions],
        "so_stop_samples_exclusive": so_stops[
            selected_so_positions],
    }
