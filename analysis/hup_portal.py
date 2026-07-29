"""Fail-closed HUP source selection and iEEG.org acquisition helpers.

This module contains no endpoint estimator.  It owns the immutable source pin,
the explicitly labeled contact heuristic, exact portal-response accounting,
shared channel geometry validation, and deterministic sparse night search used
by the current HUP cache producer.
"""
from __future__ import annotations

import concurrent.futures
import json
import os
import re
import threading

import numpy as np
from scipy import signal

from infraslow_rr_sigma_coherence import MAX_REQ_S, ROOT, get, sess


COHORT = [
    116, 130, 133, 138, 139, 141, 143, 150, 151, 157, 160, 165, 171,
    172, 173, 177, 178, 182, 185, 187, 191, 199, 205, 211, 212,
]
NIGHT_PROBE_WORKERS = 4
HUP_SOURCE_PIN_PATH = os.path.join(
    ROOT, "analysis", "hup_ieeg_source_pin.json")
with open(HUP_SOURCE_PIN_PATH) as _pin_handle:
    _HUP_SOURCE_PIN_PAYLOAD = json.load(_pin_handle)
HUP_SOURCE_PIN_SCHEMA_VERSION = _HUP_SOURCE_PIN_PAYLOAD["schema_version"]
HUP_SOURCE_PINS = _HUP_SOURCE_PIN_PAYLOAD["datasets"]

STANDARD_SCALP_EEG_LABELS = frozenset({
    "FP1", "FP2", "FPZ",
    "F7", "F3", "FZ", "F4", "F8",
    "T3", "C3", "CZ", "C4", "T4",
    "T5", "P3", "PZ", "P4", "T6",
    "O1", "OZ", "O2",
    "T7", "T8", "P7", "P8",
    "A1", "A2", "M1", "M2",
})


class NightProbeSourceMismatch(RuntimeError):
    """A worker reopened a different immutable portal snapshot."""


class PortalSampleCountMismatch(RuntimeError):
    """A portal response did not contain the exact requested sample geometry."""


def expected_portal_sample_count(duration_s, sample_rate_hz):
    """Return the integer sample count for one portal request."""
    duration_s = float(duration_s)
    sample_rate_hz = float(sample_rate_hz)
    if (
        not np.isfinite([duration_s, sample_rate_hz]).all()
        or duration_s <= 0
        or sample_rate_hz <= 0
    ):
        raise ValueError(
            "portal duration and sample rate must be finite and positive")
    return int(round(duration_s * sample_rate_hz))


def pull_continuous_exact(
        ds, channel_indices, start_s, duration_s, sample_rate_hz, *,
        records=None, purpose="analysis", max_request_s=MAX_REQ_S):
    """Pull and validate every bounded portal subrequest.

    iEEG.org can return a short array without raising.  Such a response is a
    fatal acquisition mismatch, because accepting it shifts all later time
    indexing.  ``records`` receives an audit row before any mismatch is raised.
    """
    raw_indices = list(channel_indices)
    try:
        numeric_indices = [float(value) for value in raw_indices]
        indices = [int(value) for value in numeric_indices]
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(
            "portal channel indices must be finite nonnegative integers") from exc
    if (
        not indices
        or len(indices) != len(set(indices))
        or any(
            isinstance(raw, (bool, np.bool_))
            or not np.isfinite(numeric)
            or numeric != integer
            for raw, numeric, integer in zip(
                raw_indices, numeric_indices, indices)
        )
        or any(value < 0 for value in indices)
    ):
        raise ValueError(
            "portal channel indices must be a nonempty unique list of "
            "nonnegative integers")
    start_s = float(start_s)
    duration_s = float(duration_s)
    sample_rate_hz = float(sample_rate_hz)
    max_request_s = float(max_request_s)
    if (
        not np.isfinite(
            [start_s, duration_s, sample_rate_hz, max_request_s]).all()
        or start_s < 0
        or duration_s <= 0
        or sample_rate_hz <= 0
        or max_request_s <= 0
    ):
        raise ValueError(
            "portal request geometry must be finite with a nonnegative start "
            "and positive durations/sample rate")

    audit = records if records is not None else []
    blocks = []
    offset_s = 0.0
    while offset_s < duration_s:
        request_duration_s = min(max_request_s, duration_s - offset_s)
        request_start_s = start_s + offset_s
        expected = expected_portal_sample_count(
            request_duration_s, sample_rate_hz)
        record = {
            "purpose": str(purpose),
            "request_start_s": request_start_s,
            "request_duration_s": request_duration_s,
            "requested_sample_count": expected,
            "requested_channel_count": len(indices),
        }
        try:
            block = np.asarray(
                get(ds, indices, request_start_s, request_duration_s))
        except Exception as exc:
            record.update({
                "returned_sample_count": None,
                "returned_channel_count": None,
                "status": "error",
                "error": f"{type(exc).__name__}: {exc}",
            })
            audit.append(record)
            raise

        returned_samples = int(block.shape[0]) if block.ndim >= 1 else 0
        returned_channels = int(block.shape[1]) if block.ndim == 2 else None
        record.update({
            "returned_sample_count": returned_samples,
            "returned_channel_count": returned_channels,
        })
        if (
            block.ndim != 2
            or returned_channels != len(indices)
            or returned_samples != expected
        ):
            record["status"] = "sample_count_mismatch"
            audit.append(record)
            raise PortalSampleCountMismatch(
                f"portal returned shape {block.shape} for {purpose} at "
                f"{request_start_s:g} s; expected "
                f"({expected}, {len(indices)})")
        record["status"] = "ok"
        audit.append(record)
        blocks.append(block)
        offset_s += request_duration_s

    data = np.concatenate(blocks, axis=0)
    order = np.argsort(indices)
    inverse = np.empty(len(indices), dtype=int)
    inverse[order] = np.arange(len(indices))
    return data[:, inverse]


def same_series_geometry(actual, reference):
    """Whether two portal channel identities share one exact sample time base."""
    exact = ("start_time_us", "end_time_us", "number_of_samples")
    if any(actual[key] != reference[key] for key in exact):
        return False
    return bool(
        np.isclose(
            actual["duration_us"], reference["duration_us"],
            rtol=0, atol=1)
        and np.isclose(
            actual["sample_rate_hz"], reference["sample_rate_hz"],
            rtol=0, atol=1e-9)
    )


def validate_hup_series_geometry(source_identity, cortical, ekg):
    """Require all selected cortical/ECG channels to share one exact time base."""
    cortical = list(cortical)
    if not cortical or not ekg:
        raise RuntimeError(
            "shared HUP geometry requires a cortical reference and ECG channel")
    channels = source_identity.get("channels", {})
    reference_label = cortical[0]
    if reference_label not in channels:
        raise RuntimeError(
            f"cortical reference {reference_label!r} lacks source geometry")
    reference = channels[reference_label]
    for label in cortical + [ekg]:
        identity = channels.get(label)
        if identity is None:
            raise RuntimeError(f"selected channel {label!r} lacks source geometry")
        if not same_series_geometry(identity, reference):
            raise RuntimeError(
                f"selected channel {label!r} does not share the exact time "
                f"base of cortical reference {reference_label!r}")
    return {
        "reference_channel": reference_label,
        "sample_rate_hz": float(reference["sample_rate_hz"]),
        "duration_us": float(reference["duration_us"]),
        "number_of_samples": int(reference["number_of_samples"]),
        "requirement": (
            "all selected cortical and ECG channels share start/end/sample "
            "count/duration/sample rate with the first selected cortical channel"
        ),
    }


def verify_hup_source_identity(ds, cortical, ekg, pins=HUP_SOURCE_PINS):
    """Fail closed unless the opened snapshot and used revisions match the pin."""
    name = str(getattr(ds, "name", "") or "")
    if name not in pins:
        raise RuntimeError(f"no pinned iEEG.org source identity for {name!r}")
    expected = pins[name]
    labels = list(ds.get_channel_labels())
    nodes = {
        node.findtext("channelLabel"): node
        for node in getattr(ds, "ts_array", [])
    }
    used = list(cortical) + ([ekg] if ekg else [])
    channels = {}
    compact_channels = {}
    for label in used:
        if label not in labels or label not in nodes:
            raise RuntimeError(
                f"pinned source channel {label!r} is absent from {name!r}")
        detail = ds.get_time_series_details(label)
        revision_id = str(detail.portal_id)
        data_check = nodes[label].findtext("dataCheck")
        compact_channels[label] = [revision_id, data_check]
        channels[label] = {
            "revision_id": revision_id,
            "data_check": data_check,
            "start_time_us": int(detail.start_time),
            "end_time_us": int(detail.end_time),
            "duration_us": float(detail.duration),
            "number_of_samples": int(detail.number_of_samples),
            "sample_rate_hz": float(detail.sample_rate),
        }
    actual_compact = {
        "snapshot_id": str(getattr(ds, "snap_id", "") or ""),
        "cortical_channels": list(cortical),
        "ekg": ekg,
        "channels": compact_channels,
    }
    if actual_compact != expected:
        raise RuntimeError(
            f"iEEG.org source identity mismatch for {name!r}; the pinned "
            "snapshot/channel revisions must be reviewed before analysis")
    identity = {
        "pin_schema_version": HUP_SOURCE_PIN_SCHEMA_VERSION,
        "dataset_name": name,
        "snapshot_id": actual_compact["snapshot_id"],
        "cortical_channels": list(cortical),
        "ekg": ekg,
        "channels": channels,
    }
    if cortical and ekg:
        identity["shared_series_geometry"] = validate_hup_series_geometry(
            identity, cortical, ekg)
    return identity


def delta_ratio(x, sf):
    """Welch 0.5-4-Hz power divided by 0.5-25-Hz power."""
    x = np.asarray(x, float)
    finite = np.isfinite(x)
    if finite.mean() < 0.7 or np.std(x[finite]) < 1e-9:
        return np.nan
    if not finite.all():
        sample = np.arange(len(x))
        x = x.copy()
        x[~finite] = np.interp(
            sample[~finite], sample[finite], x[finite])
    frequencies, power = signal.welch(
        x, sf, nperseg=int(min(4 * sf, len(x))))
    delta = (frequencies >= 0.5) & (frequencies < 4)
    broadband = (frequencies >= 0.5) & (frequencies < 25)
    numerator = np.trapezoid(power[delta], frequencies[delta])
    denominator = np.trapezoid(power[broadband], frequencies[broadband])
    return float(numerator / denominator) if denominator > 0 else np.nan


def is_standard_scalp_eeg_label(label):
    """Return whether a label is an unadorned conventional scalp/reference name."""
    value = str(label).strip().upper()
    match = re.fullmatch(r"([A-Z]+)0*(\d+)", value)
    if match:
        value = f"{match.group(1)}{int(match.group(2))}"
    return value in STANDARD_SCALP_EEG_LABELS


def cortical_channels(labels, n_want=6):
    """Return unvalidated heuristic lateral-contact candidates.

    Contact numbers cannot establish anatomy, gray matter, SOZ exclusion, or
    signal quality.  Publication analyses still require coordinates and
    clinical/anatomical QC.
    """
    shafts = {}
    for label in labels:
        match = re.match(r"^([A-Z]{1,3})(\d+)$", label)
        if (
            match
            and not label.upper().startswith(("EKG", "ECG"))
            and not is_standard_scalp_eeg_label(label)
        ):
            shafts.setdefault(match.group(1), []).append(int(match.group(2)))
    picks = [
        (prefix, max(numbers))
        for prefix, numbers in shafts.items()
        if max(numbers) >= 6
    ]
    picks.sort(key=lambda value: -value[1])
    available = set(labels)
    return [
        f"{prefix}{number}"
        for prefix, number in picks
        if f"{prefix}{number}" in available
    ][:n_want]


def _night_probe_record(index, time_s, worker_dataset, sf):
    requested = expected_portal_sample_count(6.0, sf)
    returned = None
    try:
        worker, channel_index = worker_dataset()
        block = np.asarray(get(worker, [channel_index], time_s, 6.0))
        returned = int(block.shape[0]) if block.ndim >= 1 else 0
        returned_channels = int(block.shape[1]) if block.ndim == 2 else None
        if (
            block.ndim != 2
            or returned_channels != 1
            or returned != requested
        ):
            raise PortalSampleCountMismatch(
                f"night probe at {time_s:g} s returned shape {block.shape}; "
                f"expected ({requested}, 1)")
        value = delta_ratio(block[:, 0], sf)
    except NightProbeSourceMismatch:
        raise
    except Exception as exc:
        value = np.nan
        failed = True
        error = f"{type(exc).__name__}: {exc}"
    else:
        failed = False
        error = None
    return index, value, failed, {
        "start_s": float(time_s),
        "duration_s": 6.0,
        "requested_sample_count": requested,
        "returned_sample_count": returned,
        "status": "ok" if not failed else "failed",
        "error": error,
    }


def find_night(
        ds, lab_idx, ch, sf, total_h, scan_h=None, step_min=30.0,
        required_h=0.0, probe_workers=NIGHT_PROBE_WORKERS,
        return_diagnostics=False):
    """Select the highest-delta, coverage-qualified contiguous three-hour window."""
    if required_h < 0:
        raise ValueError("required_h must be nonnegative")
    if int(probe_workers) != probe_workers or probe_workers < 1:
        raise ValueError("probe_workers must be a positive integer")
    if total_h < required_h:
        result = None
        diagnostics = {
            "n_probes": 0,
            "n_finite_scores": 0,
            "n_probe_failures": 0,
            "probe_failure_times_s": [],
            "probe_sample_counts": [],
            "reason": "recording shorter than requested analysis interval",
        }
        return (result, diagnostics) if return_diagnostics else result

    step_s = step_min * 60
    search_h = total_h if scan_h is None else min(scan_h, total_h)
    times = np.arange(int(search_h * 3600 / step_s), dtype=float) * step_s
    scores = np.full(len(times), np.nan)
    probe_failed = np.zeros(len(times), bool)
    sample_records = [None] * len(times)

    worker_local = threading.local()
    worker_sessions = []
    sessions_lock = threading.Lock()
    use_portal_workers = bool(
        getattr(ds, "name", None) and getattr(ds, "session", None))
    expected_snapshot_id = str(getattr(ds, "snap_id", "") or "")

    def worker_dataset():
        if hasattr(worker_local, "dataset"):
            return worker_local.dataset, worker_local.channel_index
        if not use_portal_workers:
            worker_local.dataset = ds
            worker_local.channel_index = lab_idx[ch]
            return worker_local.dataset, worker_local.channel_index
        session = sess()
        try:
            worker = session.open_dataset(ds.name)
            snapshot_id = str(getattr(worker, "snap_id", "") or "")
            if expected_snapshot_id and snapshot_id != expected_snapshot_id:
                raise NightProbeSourceMismatch(
                    f"night-probe snapshot changed for {ds.name!r}: "
                    f"expected {expected_snapshot_id!r}, got {snapshot_id!r}")
            labels = worker.get_channel_labels()
            if ch not in labels:
                raise NightProbeSourceMismatch(
                    f"night-probe channel {ch!r} is absent from "
                    f"worker dataset {ds.name!r}")
        except Exception:
            session.close()
            raise
        with sessions_lock:
            worker_sessions.append(session)
        worker_local.dataset = worker
        worker_local.channel_index = labels.index(ch)
        return worker, worker_local.channel_index

    def score(index):
        return _night_probe_record(
            index, times[index], worker_dataset, sf)

    try:
        if probe_workers == 1:
            scored = map(score, range(len(times)))
            for index, value, failed, record in scored:
                scores[index] = value
                probe_failed[index] = failed
                sample_records[index] = record
        else:
            with concurrent.futures.ThreadPoolExecutor(
                    max_workers=int(probe_workers),
                    thread_name_prefix="ieeg-night-probe") as executor:
                for index, value, failed, record in executor.map(
                        score, range(len(times))):
                    scores[index] = value
                    probe_failed[index] = failed
                    sample_records[index] = record
    finally:
        for session in worker_sessions:
            session.close()

    n_failures = int(probe_failed.sum())
    diagnostics = {
        "n_probes": int(len(times)),
        "n_finite_scores": int(np.isfinite(scores).sum()),
        "n_probe_failures": n_failures,
        "probe_failure_times_s": times[probe_failed].tolist(),
        "probe_sample_counts": sample_records,
        "reason": None,
    }
    if len(times) and n_failures > 0.20 * len(times):
        raise RuntimeError(
            f"night search failed closed: {n_failures}/{len(times)} "
            "sparse probes failed")
    if np.isfinite(scores).sum() < 6:
        diagnostics["reason"] = "fewer than six finite delta-ratio probes"
        return (None, diagnostics) if return_diagnostics else None

    window_probes = max(1, int(3 * 3600 / step_s))
    if len(scores) < window_probes:
        diagnostics["reason"] = (
            "recording shorter than the three-hour search window")
        return (None, diagnostics) if return_diagnostics else None

    best_time = None
    best_score = -np.inf
    for start in range(len(scores) - window_probes + 1):
        if times[start] + required_h * 3600 > total_h * 3600:
            continue
        values = scores[start:start + window_probes]
        if np.isfinite(values).sum() < max(
                1, int(np.ceil(0.8 * window_probes))):
            continue
        score_value = float(np.nanmean(values))
        if score_value > best_score:
            best_score = score_value
            best_time = times[start]
    if best_time is None:
        diagnostics["reason"] = "no coverage-qualified candidate interval"
    return (
        (best_time, diagnostics)
        if return_diagnostics
        else best_time
    )
