"""OpenNeuro ds003848 (Utrecht RESPect long-term iEEG) -> staged derived-series cache.

The independent replication cohort for 3A/3B. Author-provided sleep/NREM/REM/SWS, transition,
artifact, and seizure annotations define the primary coarse states and exclusions. Every subject
also carries EMG + EOG, but the rule-based multimodal labels are stored only as a sensitivity
analysis within author-unknown sleep. Neither source provides expert AASM N2/N3 scoring.

Six patients, ~1 h continuous `task-[Ss]leep` runs @ 2048 Hz, 50 Hz line. Verified against the raw
channels.tsv (2026-07): all six have iEEG (3 ECoG grid + 3 SEEG depth), ECG, EMG, EOG, and (bad)
respiration belts. Channel ROLE is taken from channels.tsv by row index, which BIDS guarantees
matches the data-column order, so naming variants (ECG+/ecg1+, emg+/EMG2, orb+/Orb+) don't matter.

Pipeline per subject: download the BrainVision triple if absent -> stream in chunks with MNE ->
derive the SAME series the HUP cache stores (so lecci_faithful_3A / event_3B_cached consume it
unchanged) PLUS author-constrained `stage_lab` and proxy sensitivity fields -> save
data/derived/ds003848/.
The raw .eeg (~3.9 GB) is optionally deleted after caching so disk stays bounded.

STAGING SENSITIVITY (not primary). True AASM scoring needs scalp EEG, which this dataset lacks.
Within author-unknown sleep, the stored sensitivity proxy uses:
  * per 30 s epoch: submental EMG RMS (notch + 10-100 Hz), EOG movement variance (0.3-6 Hz),
    and slow-wave power (0.5-4 Hz), all robust-z-scored within subject. EMG/EOG retain only samples
    covered by complete finite/measured 4-s windows at 2-s stride. Fully measured epochs exactly
    reproduce the former full-epoch RMS/variance; partial epochs intentionally omit samples in
    incomplete windows instead of losing the entire epoch.
  * Wake  = high muscle tone (EMG z > 1.0)
  * REM   = muscle atonia (EMG z < -0.3) + low SWA + phasic eye movement (EOG z > 0.5)
  * NREM  = everything else with adequate delta; split N2/N3 by a 2-component GMM on log SWA
The primary array does not use that proxy: author-unknown sleep, transitions, conflicts, and
disturbed epochs remain unclassified.

    .venv/bin/python analysis/stage_ds003848.py [--subjects sub-RESP0521,...] [--delete-raw]
"""
import argparse, hashlib, json, os, time, urllib.request
import numpy as np
import csv, io
from scipy import signal
import neurokit2 as nk
import mne

from infraslow_rr_sigma_coherence import ROOT
from cohort_stages_3ABD import (fsp_from, band_sos, EPOCH, SWA_BAND, CHUNK_S)
from cache_lc_series import (detect_so_candidates, threshold_so_candidates, SIGMA_FIXED,
                             SWA_BAND_L, SO_BAND_NAJI,
                             SO_NEGATIVE_HALF_DURATION_S,
                             SO_POSITIVE_HALF_MAX_S, FS_RR, sanitize_beats,
                             _aggregate_full_night_power, MIN_SIGNAL_COVERAGE, FILTER_EDGE_S,
                             MIN_CONTACT_COVERAGE, MIN_CONTACTS,
                             MIN_CONTACT_FRACTION_PER_BIN, prepare_continuous_signal,
                             aggregate_staging_features, interpolate_tachograms,
                             staging_epoch_features, _binned_power_values,
                             empty_channel_activity_extrema,
                             update_channel_activity_extrema,
                             finalize_channel_activity_qc,
                             STAGING_REFERENCE_MIN_VALID_WINDOWS,
                             STAGING_WELCH_WINDOW_S, STAGING_WELCH_OVERLAP_S)
from results_3A_tutorial_style import ied_clean_mask
from pipeline_version import (CACHE_SCHEMA_VERSION, atomic_savez, cache_code_sha256, git_is_dirty,
                              git_revision, npz_scalar_text, runtime_versions,
                              source_tree_sha256, utc_now, start_run_manifest,
                              validated_complete_run_exists, write_run_manifest)

mne.set_log_level("ERROR")
MIN_NREM_DELTA_RATIO = 0.20
ANNOTATION_SIGNAL_PAD_S = 5.0

BASE = "https://openneuro.org/crn/datasets/ds003848/snapshots/1.0.3/files"
RAW = os.path.join(ROOT, "data", "ds003848_raw")
OUT = os.path.join(ROOT, "data", "derived", "ds003848")

# subject -> (session, sleep-run basename without extension). One representative sleep run each.
SUBJECTS = {
    "sub-RESP0521": ("ses-1", "sub-RESP0521_ses-1_task-Sleep_run-030344"),
    "sub-RESP0699": ("ses-1", "sub-RESP0699_ses-1_task-sleep_run-030608"),
    "sub-RESP0724": ("ses-1", "sub-RESP0724_ses-1_task-sleep_run-022151"),
    "sub-RESP0749": ("ses-1", "sub-RESP0749_ses-1_task-Sleep_run-030241"),
    "sub-RESP0779": ("ses-1", "sub-RESP0779_ses-1_task-sleep_run-030241"),
    "sub-RESP0800": ("ses-1", "sub-RESP0800_ses-1_task-Sleep_run-030017"),
}

SNAPSHOT_IDENTITY_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "ds003848_snapshot_1.0.3_files.json")
with open(SNAPSHOT_IDENTITY_PATH) as _snapshot_identity_handle:
    SNAPSHOT_FILES = json.load(_snapshot_identity_handle)
_expected_snapshot_files = set()
for _subject, (_session, _base) in SUBJECTS.items():
    _expected_snapshot_files.update(
        f"{_base}{extension}"
        for extension in (
            "_channels.tsv", "_ieeg.json", "_ieeg.vhdr", "_ieeg.vmrk", "_ieeg.eeg",
            "_events.tsv",
        )
    )
    _expected_snapshot_files.add(f"{_subject}_{_session}_electrodes.tsv")
if set(SNAPSHOT_FILES) != _expected_snapshot_files:
    raise RuntimeError(
        "pinned ds003848 snapshot identity manifest does not exactly cover production inputs")


# ---------------------------------------------------------------- download
def file_sha256(path, block_size=1 << 20):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_pinned_snapshot_file(filename, path):
    """Fail closed unless a local/downloaded file is the pinned snapshot-1.0.3 object."""
    identity = SNAPSHOT_FILES.get(filename)
    if identity is None:
        raise RuntimeError(f"no pinned OpenNeuro snapshot identity for {filename}")
    actual_size = os.path.getsize(path)
    if actual_size != int(identity["size"]):
        raise RuntimeError(
            f"OpenNeuro input size mismatch for {filename}: "
            f"{actual_size} != pinned {identity['size']}")
    actual_hash = file_sha256(path)
    if actual_hash != identity["sha256"]:
        raise RuntimeError(
            f"OpenNeuro input SHA-256 mismatch for {filename}: "
            f"{actual_hash} != pinned {identity['sha256']}")
    return identity


def _dl(sub, ses, filename):
    url = f"{BASE}/{sub}:{ses}:ieeg:{filename}"
    dst = os.path.join(RAW, filename)
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        verify_pinned_snapshot_file(filename, dst)
        return dst
    os.makedirs(RAW, exist_ok=True)
    tmp = dst + ".part"
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as f:
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            f.write(b)
    verify_pinned_snapshot_file(filename, tmp)
    os.replace(tmp, dst)
    return dst


def ensure_files(sub, ses, base):
    paths = {}
    for ext in ("_channels.tsv", "_ieeg.json", "_ieeg.vhdr", "_ieeg.vmrk", "_ieeg.eeg"):
        paths[ext] = _dl(sub, ses, f"{base}{ext}")
    paths["_events.tsv"] = _dl(sub, ses, f"{base}_events.tsv")
    paths["_electrodes.tsv"] = _dl(sub, ses, f"{sub}_{ses}_electrodes.tsv")
    return paths


# ---------------------------------------------------------------- channel roles
def channel_roles(channels_tsv):
    """Return dict of role -> list of row indices (0-based, = data column order), for good channels."""
    rows = list(csv.DictReader(open(channels_tsv), delimiter="\t"))
    return channel_roles_from_rows(rows)


def channel_roles_from_rows(rows):
    """Pure row parser used by ``channel_roles`` and regression tests."""
    roles = dict(ieeg=[], ecg=[], emg=[], eog=[], names=[r["name"] for r in rows])
    for i, r in enumerate(rows):
        t = (r.get("type") or "").upper()
        bad = (r.get("status") or "").lower() == "bad"
        if t in ("SEEG", "ECOG") and not bad:
            roles["ieeg"].append(i)
        elif t in ("ECG", "EKG") and not bad:
            roles["ecg"].append(i)
        elif t == "EMG" and not bad:
            roles["emg"].append(i)
        elif t == "EOG" and not bad:
            roles["eog"].append(i)
    roles["n_rows"] = len(rows)
    return roles


def selected_channel_name_mismatches(raw_names, tsv_names, selected_indices):
    """Return order/name mismatches only for channels used by the analysis.

    MNE appends suffixes such as ``-1`` when a BrainVision header contains duplicate names.
    RESP0699 has this exact condition for an unused bad placeholder channel.  Requiring every
    unused placeholder to retain an identical display name makes a valid selected-channel mapping
    fail even though BIDS TSV row order and every requested modality agree.
    """
    if len(raw_names) != len(tsv_names):
        raise ValueError("raw and TSV channel lists must have the same length")
    selected = {int(value) for value in selected_indices}
    return [
        (i, raw_name, tsv_name)
        for i, (raw_name, tsv_name) in enumerate(zip(raw_names, tsv_names))
        if i in selected
        and raw_name.strip().casefold() != tsv_name.strip().casefold()
    ]


def electrode_eligibility_from_rows(rows, channel_names):
    """Return a conservative non-pathological cortical contact set and ROI masks.

    The RESPect electrode sidecars identify seizure-onset, resected, edge, non-gray, lesion, and
    other non-cortical contacts.  A good BIDS channel status alone does not make those contacts
    suitable for normative sleep physiology.
    """
    requested = {str(value).strip().casefold() for value in channel_names}
    by_name = {}
    for row in rows:
        name = (row.get("name") or "").strip()
        key = name.casefold()
        # RESP0699 contains duplicate unused "....." placeholders.  They are irrelevant to the
        # exact-name join, but a duplicate requested contact would be ambiguous and must fail.
        if key not in requested:
            continue
        if not name or key in by_name:
            raise ValueError(f"missing or duplicate requested electrode name: {name!r}")
        by_name[key] = row

    selected, frontal, parietal, labels, reasons = [], [], [], [], {}
    pathology_fields = (
        "soz", "resected", "edge", "silicon", "screw", "csf",
        "whitematter", "lesion", "gliosis",
    )
    for name in channel_names:
        row = by_name.get(str(name).strip().casefold())
        if row is None:
            raise ValueError(f"iEEG channel {name!r} is missing from electrodes.tsv")
        why = [
            field for field in pathology_fields
            if (row.get(field) or "").strip().casefold() == "yes"
        ]
        group = (row.get("group") or "").strip().casefold()
        if group == "depth" and (row.get("graymatter") or "").strip().casefold() != "yes":
            why.append("not_gray_matter")
        label = (row.get("Destrieux_label_text") or "").strip()
        if label.casefold() in ("", "n/a", "unknown"):
            why.append("no_cortical_atlas_label")
        keep = not why
        lower = label.casefold()
        selected.append(keep)
        frontal.append(keep and "front" in lower)
        parietal.append(
            keep and any(value in lower for value in ("pariet", "postcentral", "precuneus")))
        labels.append(label)
        if why:
            reasons[str(name)] = sorted(set(why))
    return dict(
        selected_mask=np.asarray(selected, bool),
        frontal_mask=np.asarray(frontal, bool),
        parietal_mask=np.asarray(parietal, bool),
        destrieux_labels=np.asarray(labels, dtype="<U96"),
        exclusion_reasons=reasons,
    )


def _event_category(trial_type, sub_type):
    trial = str(trial_type).strip().casefold()
    subtype = str(sub_type).strip().casefold()
    if trial in ("artefact", "artifact"):
        return "artifact"
    if trial == "seizure":
        return "seizure"
    if "stimulation" in trial:
        return "stimulation"
    if trial == "sleep":
        if subtype == "nrem":
            return "sleep_nrem"
        if subtype == "rem":
            return "sleep_rem"
        return "sleep_unknown"
    if trial == "sleep-wake transition":
        return "transition"
    if trial == "sws selection":
        return "sws_selection"
    if trial == "rem selection":
        return "rem_selection"
    return "other"


def event_annotations_from_rows(rows, total_s, sf):
    """Validate and normalize RESPect event rows to clipped half-open time intervals."""
    annotations = []
    tolerance = 2.0 / float(sf)
    for row_number, row in enumerate(rows, start=2):
        try:
            onset = float(row["onset"])
            duration = float(row["duration"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid event onset/duration at row {row_number}") from exc
        if not np.isfinite(onset) or not np.isfinite(duration) or duration <= 0:
            raise ValueError(f"non-positive or non-finite event at row {row_number}")
        stop = onset + duration
        if row.get("offset") not in (None, "", "n/a"):
            if abs(float(row["offset"]) - stop) > tolerance:
                raise ValueError(f"event offset disagrees with onset+duration at row {row_number}")
        for key, expected in (("sample_start", onset), ("sample_end", stop)):
            if row.get(key) not in (None, "", "n/a"):
                if abs(float(row[key]) / float(sf) - expected) > tolerance:
                    raise ValueError(f"{key} disagrees with event time at row {row_number}")
        clipped_start = max(0.0, onset)
        clipped_stop = min(float(total_s), stop)
        if clipped_stop <= clipped_start:
            continue
        contact_text = ",".join([
            str(row.get("electrodes_involved_onset") or ""),
            str(row.get("electrodes_involved_offset") or ""),
        ])
        contacts = sorted({
            value.strip() for value in contact_text.split(",")
            if value.strip() and value.strip().casefold() not in ("n/a",)
        })
        annotations.append(dict(
            onset_s=float(clipped_start), stop_s=float(clipped_stop),
            trial_type=str(row.get("trial_type") or ""),
            sub_type=str(row.get("sub_type") or ""),
            category=_event_category(row.get("trial_type"), row.get("sub_type")),
            contacts=contacts,
        ))
    return annotations


def event_exclusion_mask(
        annotations, start_s, n_samples, sf, channel=None,
        pad_s=0.0):
    """True for samples excluded by author artifact/seizure/stimulation annotations."""
    bad = np.zeros(int(n_samples), bool)
    channel_key = None if channel is None else str(channel).strip().casefold()
    for event in annotations:
        category = event["category"]
        contacts = {value.casefold() for value in event["contacts"]}
        applies = category in ("seizure", "stimulation")
        if category == "artifact":
            applies = "all" in contacts or (
                channel_key is not None and channel_key in contacts)
        if not applies:
            continue
        first = max(
            0, int(np.floor((event["onset_s"] - float(pad_s) - start_s) * sf)))
        last = min(
            len(bad), int(np.ceil((event["stop_s"] + float(pad_s) - start_s) * sf)))
        if last > first:
            bad[first:last] = True
    return bad


def _covered_seconds(intervals, start, stop):
    clipped = sorted(
        (max(start, value[0]), min(stop, value[1]))
        for value in intervals
        if value[1] > start and value[0] < stop
    )
    if not clipped:
        return 0.0
    total = 0.0
    left, right = clipped[0]
    for next_left, next_right in clipped[1:]:
        if next_left <= right:
            right = max(right, next_right)
        else:
            total += right - left
            left, right = next_left, next_right
    return float(total + right - left)


def annotation_epoch_context(annotations, n_epochs, epoch_s=EPOCH, boundary_tolerance_s=1.0):
    """Author-derived coarse epoch context; boundary/transition conflicts remain unclassified."""
    categories = {}
    for event in annotations:
        categories.setdefault(event["category"], []).append(
            (event["onset_s"], event["stop_s"]))
    sleep_intervals = sum(
        (categories.get(value, []) for value in
         ("sleep_nrem", "sleep_rem", "sleep_unknown",
          "sws_selection", "rem_selection")), [])
    transition_intervals = categories.get("transition", [])
    labels = np.full(int(n_epochs), "", dtype="<U5")
    sources = np.full(int(n_epochs), "boundary", dtype="<U32")
    required = float(epoch_s) - float(boundary_tolerance_s)
    for epoch in range(int(n_epochs)):
        start, stop = epoch * float(epoch_s), (epoch + 1) * float(epoch_s)
        coverage = {
            key: _covered_seconds(values, start, stop)
            for key, values in categories.items()
        }
        if any(
            coverage.get(value, 0.0) > 0.0
            for value in ("artifact", "seizure", "stimulation")
        ):
            sources[epoch] = "author_disturbance"
        elif coverage.get("sws_selection", 0.0) >= required:
            labels[epoch], sources[epoch] = "N3", "author_sws_selection"
        elif (
            coverage.get("rem_selection", 0.0) >= required
            or coverage.get("sleep_rem", 0.0) >= required
        ):
            labels[epoch], sources[epoch] = "R", "author_rem"
        elif coverage.get("sleep_nrem", 0.0) >= required:
            labels[epoch], sources[epoch] = "NREM", "author_nrem"
        elif coverage.get("sleep_unknown", 0.0) >= required:
            labels[epoch], sources[epoch] = "SLEEP", "author_sleep_unknown"
        elif (
            _covered_seconds(sleep_intervals, start, stop)
            + _covered_seconds(transition_intervals, start, stop)
            <= float(boundary_tolerance_s)
        ):
            labels[epoch], sources[epoch] = "W", "author_awake"
        elif _covered_seconds(transition_intervals, start, stop) > 0:
            sources[epoch] = "author_transition"
    return labels, sources


def annotation_second_masks(annotations, total_s):
    """One-second reporting masks for exact source annotations (no processing padding)."""
    keys = (
        "artifact", "seizure", "stimulation", "sleep_nrem", "sleep_rem",
        "sleep_unknown", "transition", "sws_selection", "rem_selection",
    )
    masks = {key: np.zeros(int(total_s), bool) for key in keys}
    for event in annotations:
        if event["category"] not in masks:
            continue
        start = max(0, int(np.floor(event["onset_s"])))
        stop = min(int(total_s), int(np.ceil(event["stop_s"])))
        masks[event["category"]][start:stop] = True
    sleep = (
        masks["sleep_nrem"] | masks["sleep_rem"] | masks["sleep_unknown"]
        | masks["sws_selection"] | masks["rem_selection"]
    )
    masks["awake"] = ~(sleep | masks["transition"])
    return masks


def constrain_proxy_to_annotations(
        proxy_labels, annotation_labels, annotation_sources,
        allow_proxy_unknown_sleep=False):
    """Use author states as primary; proxy-only unknown-sleep labels are sensitivity data."""
    proxy = np.asarray(proxy_labels).astype(str)
    annotation = np.asarray(annotation_labels).astype(str)
    sources = np.asarray(annotation_sources).astype(str).copy()
    if proxy.shape != annotation.shape or proxy.shape != sources.shape:
        raise ValueError("proxy and annotation epoch arrays must align")
    out = np.full(len(proxy), "", dtype="<U4")
    for value in ("W", "R", "NREM", "N3"):
        out[annotation == value] = value
    unknown_sleep = annotation == "SLEEP"
    usable_proxy = unknown_sleep & np.isin(proxy, ("R", "NREM", "N2", "N3"))
    if allow_proxy_unknown_sleep:
        out[usable_proxy] = proxy[usable_proxy]
        sources[usable_proxy] = "proxy_within_author_sleep"
    sources[unknown_sleep & ~usable_proxy] = "unclassified_author_sleep"
    if not allow_proxy_unknown_sleep:
        sources[unknown_sleep] = "unclassified_author_sleep"
    return out, sources


def aggregate_auxiliary_epoch_features(
        values_by_channel, min_channel_coverage=MIN_CONTACT_COVERAGE,
        min_channel_fraction_per_epoch=MIN_CONTACT_FRACTION_PER_BIN,
        valid_window_count_by_channel=None, min_valid_windows=1):
    """Scale-normalize and combine a fixed set of sufficiently supported EMG/EOG channels.

    Window support is optional so historical callers remain compatible.  Neutral RESPect caches
    provide it, allowing an offline profile to choose both its complete-window requirement and its
    channel coverage/fraction rules rather than inheriting an irreversible 80% gate.
    """
    values = np.asarray(values_by_channel, float)
    if values.ndim != 2:
        raise ValueError("auxiliary features must be channel-by-epoch")
    valid = np.isfinite(values) & (values > 0)
    if valid_window_count_by_channel is not None:
        window_count = np.asarray(valid_window_count_by_channel)
        if window_count.shape != values.shape:
            raise ValueError("auxiliary valid-window support must align with feature values")
        valid &= window_count >= int(min_valid_windows)
        values = np.where(valid, values, np.nan)
    coverage = valid.mean(axis=1)
    selected = coverage >= float(min_channel_coverage)
    out = np.full(values.shape[1], np.nan)
    scales = np.full(values.shape[0], np.nan)
    for channel in np.where(selected)[0]:
        scales[channel] = float(np.median(values[channel, valid[channel]]))
    n_selected = int(selected.sum())
    required = int(np.ceil(float(min_channel_fraction_per_epoch) * n_selected))
    support = np.zeros(values.shape[1], int)
    if n_selected:
        normalized = values[selected] / scales[selected, None]
        finite = np.isfinite(normalized) & (normalized > 0)
        support = finite.sum(axis=0)
        for epoch in np.where(support >= max(1, required))[0]:
            out[epoch] = float(np.exp(np.median(np.log(normalized[finite[:, epoch], epoch]))))
    return out, dict(
        selected_channel_mask=selected,
        per_channel_coverage=coverage,
        n_selected_channels=n_selected,
        required_channel_count=max(1, required) if n_selected else 0,
        support_count=support,
        minimum_valid_windows=int(min_valid_windows),
        minimum_channel_coverage=float(min_channel_coverage),
        minimum_channel_fraction_per_epoch=float(min_channel_fraction_per_epoch),
        normalization="divide each channel by its full-record median then geometric median",
    )


def auxiliary_epoch_window_features(
        segment, measured_mask, sf, kind, *, return_details=False):
    """Gap-aware EMG/EOG feature from complete fixed 4-s windows at 2-s stride.

    Only samples belonging to at least one fully measured, finite window enter the default epoch
    value, and each retained sample is counted once despite overlapping windows.  Thus a fully
    measured 30-s epoch is *exactly* the former full-epoch RMS/variance.  In a partial epoch, the
    result differs intentionally: incomplete-window samples are omitted instead of invalidating
    the whole epoch.  Per-window mean-square (EMG) or variance (EOG) is also returned for offline
    support thresholds and alternative equal-window aggregation.
    """
    segment = np.asarray(segment, float)
    measured_mask = np.asarray(measured_mask, bool)
    if kind not in ("emg_rms", "eog_variance"):
        raise ValueError("kind must be 'emg_rms' or 'eog_variance'")
    nperseg = int(round(STAGING_WELCH_WINDOW_S * float(sf)))
    noverlap = int(round(STAGING_WELCH_OVERLAP_S * float(sf)))
    step = nperseg - noverlap
    n_windows = (
        0 if len(segment) < nperseg
        else 1 + (len(segment) - nperseg) // step
    )
    empty_details = dict(
        n_valid_windows=0,
        n_total_windows=int(n_windows),
        valid_window_mask=np.zeros(n_windows, bool),
        window_power=np.full(n_windows, np.nan),
        covered_sample_fraction=0.0,
    )
    if (
        len(segment) == 0
        or measured_mask.shape != segment.shape
        or nperseg < 1
        or step < 1
        or n_windows < 1
    ):
        return (np.nan, empty_details) if return_details else np.nan

    finite_measured = measured_mask & np.isfinite(segment)
    starts = np.arange(n_windows, dtype=int) * step
    valid_windows = np.asarray([
        bool(finite_measured[start:start + nperseg].all())
        for start in starts
    ])
    window_power = np.full(n_windows, np.nan)
    covered = np.zeros(len(segment), bool)
    for window, start in enumerate(starts):
        if not valid_windows[window]:
            continue
        values = segment[start:start + nperseg]
        window_power[window] = (
            float(np.mean(values ** 2))
            if kind == "emg_rms"
            else float(np.var(values))
        )
        covered[start:start + nperseg] = True

    if covered.any():
        retained = segment[covered]
        value = (
            float(np.sqrt(np.mean(retained ** 2)))
            if kind == "emg_rms"
            else float(np.var(retained))
        )
    else:
        value = np.nan
    details = dict(
        n_valid_windows=int(valid_windows.sum()),
        n_total_windows=int(n_windows),
        valid_window_mask=valid_windows,
        window_power=window_power,
        covered_sample_fraction=float(covered.mean()),
    )
    return (value, details) if return_details else value


# ---------------------------------------------------------------- staging
def robust_z(x, reference=None):
    """Robust z score using only prespecified eligible reference epochs."""
    x = np.asarray(x, float)
    reference = np.isfinite(x) if reference is None else (
        np.asarray(reference, bool) & np.isfinite(x))
    if not reference.any():
        return np.full(len(x), np.nan)
    m = np.median(x[reference])
    s = np.median(np.abs(x[reference] - m)) * 1.4826 + 1e-12
    return (x - m) / s


def slow_wave_candidate_mask(ep_dr, ep_clean, eligibility=None):
    """Conservative high-delta component for the unvalidated RESPect stage proxy.

    EMG/EOG rules can exclude obvious wake/REM but cannot make every remaining epoch NREM. Require
    an actually separable high-delta component plus slow-frequency enrichment above the flat-
    spectrum bandwidth ratio. If the recording does not support that distinction, return no NREM
    rather than manufacturing labels.
    """
    dr = np.asarray(ep_dr, float)
    valid = np.isfinite(dr) & (np.asarray(ep_clean, float) >= 0.5)
    if eligibility is not None:
        valid &= np.asarray(eligibility, bool)
    mask = np.zeros(len(dr), bool)
    if valid.sum() < 40:
        return mask, dict(reason="fewer than 40 clean epochs")
    values = dr[valid]
    if np.nanpercentile(values, 90) - np.nanpercentile(values, 10) < 1e-3:
        return mask, dict(reason="delta ratio is effectively constant")
    try:
        from sklearn.mixture import GaussianMixture
        value_matrix = values.reshape(-1, 1)
        one = GaussianMixture(n_components=1, random_state=0, n_init=5).fit(value_matrix)
        model = GaussianMixture(n_components=2, random_state=0, n_init=10).fit(value_matrix)
        labels = model.predict(values.reshape(-1, 1))
        means = model.means_.ravel()
        hi = int(np.argmax(means))
        separation = float(
            abs(np.diff(means)[0]) / np.sqrt(np.mean(model.covariances_.ravel())))
        selected = labels == hi
        fraction = float(np.mean(selected))
        bic_gain = float(one.bic(value_matrix) - model.bic(value_matrix))

        # A 2-component GMM can always split a skewed unimodal distribution. These gates establish
        # only a reproducible high-tail partition, not latent sleep states or stage validity.
        cv_gains, split_means = [], []
        for train_idx, test_idx in (
                (np.arange(0, len(values), 2), np.arange(1, len(values), 2)),
                (np.arange(1, len(values), 2), np.arange(0, len(values), 2))):
            if len(train_idx) < 20 or len(test_idx) < 20:
                return mask, dict(reason="insufficient epochs for split-half validation")
            one_split = GaussianMixture(
                n_components=1, random_state=0, n_init=5).fit(value_matrix[train_idx])
            two_split = GaussianMixture(
                n_components=2, random_state=0, n_init=10).fit(value_matrix[train_idx])
            cv_gains.append(float(
                (two_split.score(value_matrix[test_idx])
                 - one_split.score(value_matrix[test_idx])) * len(test_idx)))
            split_means.append(np.sort(two_split.means_.ravel()))
        full_means = np.sort(means)
        mean_gap = float(np.diff(full_means)[0])
        stable = all(np.max(np.abs(value - full_means)) <= max(0.05, 0.35 * mean_gap)
                     for value in split_means)
        diagnostic = dict(
            separation=separation, bic_gain_2_vs_1=bic_gain,
            split_half_loglik_gains=cv_gains, stable_split_means=bool(stable),
            high_component_fraction=fraction,
            component_means=[float(value) for value in full_means])
        if (bic_gain <= 10.0 or min(cv_gains) <= 0.0 or not stable
                or separation < 0.75 or not 0.10 <= fraction <= 0.90
                or means[hi] < MIN_NREM_DELTA_RATIO):
            diagnostic["reason"] = "no stable two-Gaussian high-delta partition"
            return mask, diagnostic
        selected &= values >= MIN_NREM_DELTA_RATIO
        mask[np.where(valid)[0]] = selected
        diagnostic["reason"] = "accepted stable two-Gaussian high-delta partition"
        return mask, diagnostic
    except Exception as exc:
        return mask, dict(reason=f"delta model failed: {type(exc).__name__}")


def reliable_two_state_split(values):
    """Return a stable high-tail mixture partition, not a validated physiological state split."""
    values = np.asarray(values, float)
    if len(values) < 40 or np.nanpercentile(values, 90) - np.nanpercentile(values, 10) < 1e-3:
        return None, dict(reason="insufficient variation for two states")
    try:
        from sklearn.mixture import GaussianMixture
        matrix = values.reshape(-1, 1)
        one = GaussianMixture(n_components=1, random_state=0, n_init=5).fit(matrix)
        two = GaussianMixture(n_components=2, random_state=0, n_init=10).fit(matrix)
        labels = two.predict(matrix)
        means = two.means_.ravel()
        hi = int(np.argmax(means))
        fraction = float(np.mean(labels == hi))
        separation = float(
            abs(np.diff(means)[0]) / np.sqrt(np.mean(two.covariances_.ravel())))
        bic_gain = float(one.bic(matrix) - two.bic(matrix))
        cv_gains, split_means = [], []
        for train_idx, test_idx in (
                (np.arange(0, len(values), 2), np.arange(1, len(values), 2)),
                (np.arange(1, len(values), 2), np.arange(0, len(values), 2))):
            one_split = GaussianMixture(
                n_components=1, random_state=0, n_init=5).fit(matrix[train_idx])
            two_split = GaussianMixture(
                n_components=2, random_state=0, n_init=10).fit(matrix[train_idx])
            cv_gains.append(float(
                (two_split.score(matrix[test_idx]) - one_split.score(matrix[test_idx]))
                * len(test_idx)))
            split_means.append(np.sort(two_split.means_.ravel()))
        full_means = np.sort(means)
        gap = float(np.diff(full_means)[0])
        stable = all(np.max(np.abs(value - full_means)) <= max(0.05, 0.35 * gap)
                     for value in split_means)
        diagnostic = dict(
            bic_gain_2_vs_1=bic_gain, split_half_loglik_gains=cv_gains,
            separation=separation, high_component_fraction=fraction,
            stable_split_means=bool(stable),
            component_means=[float(value) for value in full_means])
        accepted = (bic_gain > 10.0 and min(cv_gains) > 0.0 and stable
                    and separation >= 0.75 and 0.10 <= fraction <= 0.90)
        diagnostic["reason"] = (
            "accepted stable two-Gaussian high-tail partition" if accepted
            else "no stable two-Gaussian high-tail partition")
        return (labels == hi) if accepted else None, diagnostic
    except Exception as exc:
        return None, dict(reason=f"high-tail partition failed: {type(exc).__name__}")


def score_stages(ep_swa, ep_emg, ep_eog, ep_dr, ep_clean):
    """Unvalidated rule-based stage proxies from EMG, EOG, and iEEG features.

    These are not expert-scored AASM stages and must not be described as "real staging".
    """
    n = len(ep_swa)
    # Persist labels as fixed-width Unicode.  Object arrays require pickle on load, but all
    # production cache readers deliberately use ``allow_pickle=False`` for safety.
    lab = np.full(n, "", dtype="<U4")
    log_swa = np.log(np.clip(ep_swa, 1e-12, None))
    log_emg = np.log(np.clip(ep_emg, 1e-12, None))
    log_eog = np.log(np.clip(ep_eog, 1e-12, None))
    # Dirty epochs must not influence the centre/MAD that defines clean-epoch wake/REM thresholds.
    # Applying the clean mask only after z-scoring lets extreme excluded data relabel retained data.
    valid = (
        (np.asarray(ep_clean, float) >= 0.5)
        & np.isfinite(log_swa) & np.isfinite(log_emg) & np.isfinite(log_eog)
    )
    swa_z = robust_z(log_swa, valid)
    emg_z = robust_z(log_emg, valid)
    eog_z = robust_z(log_eog, valid)

    wake = valid & (emg_z > 1.0)                                   # high muscle tone
    rem = valid & ~wake & (emg_z < -0.3) & (swa_z < 0.0) & (eog_z > 0.5)   # atonia + phasic eyes + low SWA
    # The delta model must use exactly the multimodally eligible reference set. Otherwise epochs
    # with missing EMG/EOG can move the fitted components despite never being label-eligible.
    slow_wave, delta_diagnostic = slow_wave_candidate_mask(
        ep_dr, ep_clean, eligibility=valid)
    nrem = valid & ~wake & ~rem & slow_wave
    lab[wake] = "W"
    lab[rem] = "R"
    lab[nrem] = "NREM"
    swa_diagnostic = dict(reason="fewer than 40 NREM candidates")
    if nrem.sum() >= 40:
        high_swa, swa_diagnostic = reliable_two_state_split(swa_z[nrem])
        if high_swa is not None:
            lab[np.where(nrem)[0]] = np.where(high_swa, "N3", "N2")
    return lab, dict(n_wake=int(wake.sum()), n_rem=int(rem.sum()), n_nrem=int(nrem.sum()),
                     n_nrem_unsplit=int((lab == "NREM").sum()),
                     n_valid=int(valid.sum()), delta_model=delta_diagnostic,
                     swa_stage_model=swa_diagnostic)


# ---------------------------------------------------------------- main derive
def run(subject, delete_raw=False, force=False):
    fp = os.path.join(OUT, f"{subject}.npz")
    current_cache_digest = cache_code_sha256(ROOT)
    if os.path.exists(fp) and not force:
        raise RuntimeError(
            f"{fp} cannot be reused outside a validated complete run; rerun with --force")
    ses, base = SUBJECTS[subject]
    t0 = time.time()
    print(f"[{subject}] fetching {base} ...", flush=True)
    paths = ensure_files(subject, ses, base)
    source_snapshot_identities = {
        os.path.basename(path): SNAPSHOT_FILES[os.path.basename(path)]
        for path in paths.values()
    }
    roles = channel_roles(paths["_channels.tsv"])
    if not roles["ieeg"] or not roles["ecg"]:
        reason = (f"extraction inputs missing: iEEG={bool(roles['ieeg'])}, "
                  f"ECG={bool(roles['ecg'])}; EMG/EOG are retained as optional modalities")
        atomic_savez(fp, subject=subject, status="skip",
                     cache_schema_version=CACHE_SCHEMA_VERSION,
                     cache_code_sha256=current_cache_digest, reason=reason)
        print(f"[{subject}] SKIP {reason}", flush=True)
        return "skip"

    raw = mne.io.read_raw_brainvision(paths["_ieeg.vhdr"], preload=False, verbose="ERROR")
    raw_sf = float(raw.info["sfreq"])
    with open(paths["_ieeg.json"]) as handle:
        sidecar = json.load(handle)
    sf = float(sidecar.get("SamplingFrequency", raw_sf))
    if not np.isfinite(sf) or sf <= 0 or abs(raw_sf - sf) / sf > 1e-5:
        raise RuntimeError(
            f"BrainVision/BIDS sampling-frequency mismatch: {raw_sf} vs {sf}")
    n_samp = raw.n_times
    if len(raw.ch_names) != roles["n_rows"]:
        raise RuntimeError(
            f"channel count mismatch: raw {len(raw.ch_names)} vs TSV {roles['n_rows']}")
    selected_role_indices = (
        roles["ieeg"] + roles["ecg"] + roles["emg"] + roles["eog"])
    mismatched = selected_channel_name_mismatches(
        raw.ch_names, roles["names"], selected_role_indices)
    if mismatched:
        raise RuntimeError(
            f"selected BrainVision/TSV channel-order mismatch; first mismatch={mismatched[0]}")
    names = roles["names"]
    all_good_ie = roles["ieeg"]
    with open(paths["_electrodes.tsv"]) as handle:
        electrode_rows = list(csv.DictReader(handle, delimiter="\t"))
    anatomy_qc = electrode_eligibility_from_rows(
        electrode_rows, [names[index] for index in all_good_ie])
    ie = [
        index for index, keep in zip(all_good_ie, anatomy_qc["selected_mask"])
        if keep
    ]
    if len(ie) == 0:
        reason = (
            "no non-SOZ, non-resected, non-edge cortical gray-matter contacts")
        atomic_savez(
            fp, subject=subject, status="skip",
            cache_schema_version=CACHE_SCHEMA_VERSION,
            cache_code_sha256=current_cache_digest, reason=reason)
        print(f"[{subject}] SKIP {reason}", flush=True)
        return "skip"
    retained_positions = np.where(anatomy_qc["selected_mask"])[0]
    frontal_contact_mask = anatomy_qc["frontal_mask"][retained_positions]
    parietal_contact_mask = anatomy_qc["parietal_mask"][retained_positions]
    destrieux_labels = anatomy_qc["destrieux_labels"][retained_positions]
    ecg_i = roles["ecg"][0]
    emg_indices = list(roles["emg"])
    eog_indices = list(roles["eog"])
    ctx = [names[i] for i in ie]

    total_s = int(n_samp / sf)
    n_ep = int(total_s // EPOCH)
    with open(paths["_events.tsv"]) as handle:
        event_rows = list(csv.DictReader(handle, delimiter="\t"))
    annotations = event_annotations_from_rows(event_rows, total_s, sf)

    # Individual FSP is disabled until it is estimated from all clean NREM and manually QC'd, as
    # in Lecci. The fixed 10-15 Hz analysis is primary.
    fsp, fsp_real = 13.0, False
    print(
        f"[{subject}] {sf:.0f} Hz ({raw_sf:.6f} header), {len(ie)}/{len(all_good_ie)} "
        f"anatomy-eligible iEEG, ECG@{ecg_i}, EMG@{emg_indices}, EOG@{eog_indices}, "
        f"FSP {fsp:.2f} Hz, {total_s} s", flush=True)

    sig_fixed_ch = np.full((len(ie), total_s), np.nan)
    sig_fsp_ch = np.full((len(ie), total_s), np.nan)
    swa_ch = np.full((len(ie), total_s), np.nan)
    sig_fixed_num_ch = np.full((len(ie), total_s), np.nan)
    sig_fsp_num_ch = np.full((len(ie), total_s), np.nan)
    swa_num_ch = np.full((len(ie), total_s), np.nan)
    sig_fixed_den_ch = np.zeros((len(ie), total_s), dtype=np.uint32)
    sig_fsp_den_ch = np.zeros((len(ie), total_s), dtype=np.uint32)
    swa_den_ch = np.zeros((len(ie), total_s), dtype=np.uint32)
    ep_dr_ch = np.full((len(ie), n_ep), np.nan)
    ep_swa_ch = np.full((len(ie), n_ep), np.nan)
    ep_clean_ch = np.full((len(ie), n_ep), np.nan)
    ep_measured_ch = np.full((len(ie), n_ep), np.nan)
    ep_valid_window_ch = np.zeros(
        (len(ie), n_ep, STAGING_REFERENCE_MIN_VALID_WINDOWS), dtype=np.uint8)
    ep_window_swa_ch = np.full(
        (len(ie), n_ep, STAGING_REFERENCE_MIN_VALID_WINDOWS), np.nan, dtype=np.float32)
    ep_window_total_ch = np.full(
        (len(ie), n_ep, STAGING_REFERENCE_MIN_VALID_WINDOWS), np.nan, dtype=np.float32)
    ep_longest_clean_run_ch = np.zeros((len(ie), n_ep), dtype=np.float32)
    ep_valid_window_span_ch = np.zeros((len(ie), n_ep), dtype=np.float32)
    ep_emg_ch = np.full((len(emg_indices), n_ep), np.nan)
    ep_eog_ch = np.full((len(eog_indices), n_ep), np.nan)
    ep_emg_valid_window_ch = np.zeros(
        (len(emg_indices), n_ep, STAGING_REFERENCE_MIN_VALID_WINDOWS), dtype=np.uint8)
    ep_eog_valid_window_ch = np.zeros(
        (len(eog_indices), n_ep, STAGING_REFERENCE_MIN_VALID_WINDOWS), dtype=np.uint8)
    ep_emg_window_power_ch = np.full(
        (len(emg_indices), n_ep, STAGING_REFERENCE_MIN_VALID_WINDOWS),
        np.nan, dtype=np.float32)
    ep_eog_window_power_ch = np.full(
        (len(eog_indices), n_ep, STAGING_REFERENCE_MIN_VALID_WINDOWS),
        np.nan, dtype=np.float32)
    beats = []
    so_candidates = {c: [] for c in ctx}
    failed_chunks = []
    ecg_failures = []
    channel_activity_state = empty_channel_activity_extrema(len(ctx))

    sos_fixed = band_sos(SIGMA_FIXED, sf)
    sos_fsp = band_sos((fsp - 1, fsp + 1), sf)
    sos_swa = band_sos(SWA_BAND_L, sf, 3)
    sos_so = signal.butter(3, list(SO_BAND_NAJI), btype="band", fs=sf, output="sos")
    sos_emg = signal.butter(3, [10.0, min(100.0, sf / 2 - 10)], btype="band", fs=sf, output="sos")
    sos_eog = signal.butter(3, [0.3, 6.0], btype="band", fs=sf, output="sos")

    # 50 Hz line notch (+ harmonics) as a precomputed IIR cascade -- reused per channel/chunk, far
    # cheaper than re-designing an FIR (mne.filter.notch_filter) on every call.
    notch_sos = np.vstack([signal.tf2sos(*signal.iirnotch(f0, 30.0, sf))
                           for f0 in np.arange(50.0, sf / 2 - 1, 50.0)])

    def notch50(x):
        return signal.sosfiltfilt(notch_sos, x)

    picks_all = list(ie) + [ecg_i] + emg_indices + eog_indices
    t = 0.0
    while t < total_s:
        dur = min(CHUNK_S, total_s - t)
        pull_start = max(0.0, t - FILTER_EDGE_S)
        pull_stop = min(float(total_s), t + dur + FILTER_EDGE_S)
        a = int(round(pull_start * sf))
        b = min(int(round(pull_stop * sf)), n_samp)
        try:
            d = raw.get_data(picks=picks_all, start=a, stop=b) * 1e6  # (n_pick, n), uV
        except Exception as exc:
            failed_chunks.append(dict(start_s=float(t), duration_s=float(dur),
                                      error=f"{type(exc).__name__}: {exc}"))
            t += dur
            continue
        off = int(t)
        core_a = int(round((t - pull_start) * sf))
        core_n = min(int(round(dur * sf)), max(0, d.shape[1] - core_a))
        core_b = core_a + core_n
        x_ie = d[:len(ie)]
        x_ecg = d[len(ie)]
        update_channel_activity_extrema(
            channel_activity_state, x_ie[:, core_a:core_b])
        col = len(ie) + 1
        x_emg = d[col:col + len(emg_indices)]
        col += len(emg_indices)
        x_eog = d[col:col + len(eog_indices)]
        global_event_bad = event_exclusion_mask(
            annotations, pull_start, d.shape[1], sf, channel=None)

        # Keep contacts separate through filtering, artifact masking, and power extraction.
        # Raw-voltage averaging can cancel equal power with opposite polarity/phase.
        prepared = []
        for ci, channel in enumerate(ctx):
            annotated_bad = event_exclusion_mask(
                annotations, pull_start, d.shape[1], sf, channel=channel)
            prepared.append(prepare_continuous_signal(
                np.where(annotated_bad, np.nan, x_ie[ci]), sf))
        x_channels = np.asarray([
            notch50(signal.detrend(filled)) for filled, _ in prepared
        ])
        measured_channels = np.asarray([measured for _, measured in prepared])
        clean_channels = np.asarray([
            ied_clean_mask(x, sf) & measured
            for x, measured in zip(x_channels, measured_channels)
        ])
        # ECG -> beats
        try:
            ecg_filled, ecg_measured = prepare_continuous_signal(
                np.where(global_event_bad, np.nan, x_ecg), sf)
            cl = nk.ecg_clean(ecg_filled, sampling_rate=int(sf), method="neurokit")
            _, info = nk.ecg_peaks(cl, sampling_rate=int(sf), method="neurokit", correct_artifacts=True)
            peaks = np.asarray(info["ECG_R_Peaks"], int)
            peaks = peaks[(peaks >= core_a) & (peaks < core_b)]
            peaks = peaks[ecg_measured[peaks]]
            beats.extend((peaks / sf + pull_start).tolist())
        except Exception as exc:
            ecg_failures.append(dict(start_s=float(t), duration_s=float(dur),
                                     error=f"{type(exc).__name__}: {exc}"))

        # per-channel SO troughs
        for ci, c in enumerate(ctx):
            xc = x_channels[ci]
            if np.std(xc) < 1e-9:
                continue
            cand = detect_so_candidates(signal.sosfiltfilt(sos_so, xc), sf)
            so_clean = ied_clean_mask(xc, sf, pad_s=5.0) & measured_channels[ci]
            cand = [
                value for value in cand
                if core_a <= int(value[0]) < core_b and so_clean[int(value[0])]
            ]
            so_candidates[c].extend(
                [(tr / sf + pull_start, down, up, p2p) for tr, down, up, p2p in cand])

        n_sec = int(core_n // int(sf))
        k = int(sf)
        for sos_b, dest, numerator_dest, denominator_dest in (
                (sos_fixed, sig_fixed_ch, sig_fixed_num_ch, sig_fixed_den_ch),
                (sos_fsp, sig_fsp_ch, sig_fsp_num_ch, sig_fsp_den_ch),
                (sos_swa, swa_ch, swa_num_ch, swa_den_ch)):
            for ci, (xc, clean_c) in enumerate(zip(x_channels, clean_channels)):
                env2 = np.abs(signal.hilbert(signal.sosfiltfilt(sos_b, xc))) ** 2
                start = core_a
                stop = start + n_sec * k
                vals_c, numerator, denominator = _binned_power_values(
                    env2[start:stop], clean_c[start:stop], sf, n_sec,
                    return_support=True)
                sl = slice(off, min(off + n_sec, total_s))
                length = sl.stop - sl.start
                dest[ci, sl] = vals_c[:length]
                numerator_dest[ci, sl] = numerator[:length]
                denominator_dest[ci, sl] = denominator[:length]

        # EMG / EOG epoch features
        emg_filtered, emg_measured = [], []
        for channel_signal in x_emg:
            filled, measured = prepare_continuous_signal(
                np.where(global_event_bad, np.nan, channel_signal), sf)
            emg_filtered.append(
                np.abs(signal.sosfiltfilt(sos_emg, notch50(filled))))
            emg_measured.append(measured)
        eog_filtered, eog_measured = [], []
        for channel_signal in x_eog:
            filled, measured = prepare_continuous_signal(
                np.where(global_event_bad, np.nan, channel_signal), sf)
            eog_filtered.append(signal.sosfiltfilt(sos_eog, filled))
            eog_measured.append(measured)
        emg_filtered = np.asarray(emg_filtered)
        emg_measured = np.asarray(emg_measured)
        eog_filtered = np.asarray(eog_filtered)
        eog_measured = np.asarray(eog_measured)

        ke = int(EPOCH * sf)
        for e in range(int(core_n // ke)):
            gi = int((off + e * EPOCH) // EPOCH)
            if gi >= n_ep:
                break
            epoch_a = core_a + e * ke
            epoch_b = epoch_a + ke
            for ci, (xc, clean_c, measured_c) in enumerate(zip(
                    x_channels, clean_channels, measured_channels)):
                seg = xc[epoch_a:epoch_b]
                dr, swa_value, staging_details = staging_epoch_features(
                    seg, clean_c[epoch_a:epoch_b], measured_c[epoch_a:epoch_b], sf,
                    return_details=True)
                ep_dr_ch[ci, gi] = dr
                ep_swa_ch[ci, gi] = swa_value
                ep_clean_ch[ci, gi] = staging_details["clean_fraction"]
                ep_measured_ch[ci, gi] = staging_details["measured_fraction"]
                ep_valid_window_ch[ci, gi] = staging_details["valid_window_mask"]
                ep_window_swa_ch[ci, gi] = staging_details["window_swa_power"]
                ep_window_total_ch[ci, gi] = staging_details["window_total_power"]
                ep_longest_clean_run_ch[ci, gi] = staging_details["longest_valid_run_s"]
                ep_valid_window_span_ch[ci, gi] = staging_details["valid_window_span_s"]
            for channel in range(len(emg_indices)):
                value, details = auxiliary_epoch_window_features(
                    emg_filtered[channel, epoch_a:epoch_b],
                    emg_measured[channel, epoch_a:epoch_b],
                    sf, "emg_rms", return_details=True)
                ep_emg_ch[channel, gi] = value
                ep_emg_valid_window_ch[channel, gi] = details["valid_window_mask"]
                ep_emg_window_power_ch[channel, gi] = details["window_power"]
            for channel in range(len(eog_indices)):
                value, details = auxiliary_epoch_window_features(
                    eog_filtered[channel, epoch_a:epoch_b],
                    eog_measured[channel, epoch_a:epoch_b],
                    sf, "eog_variance", return_details=True)
                ep_eog_ch[channel, gi] = value
                ep_eog_valid_window_ch[channel, gi] = details["valid_window_mask"]
                ep_eog_window_power_ch[channel, gi] = details["window_power"]
        t += dur

    channel_activity_qc = finalize_channel_activity_qc(channel_activity_state)
    nonflat_contacts = channel_activity_qc["nonflat_mask"]
    sig_fixed, sigma_contact_qc = _aggregate_full_night_power(
        sig_fixed_ch, eligible_channels=nonflat_contacts,
        min_contact_coverage=MIN_CONTACT_COVERAGE,
        min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN,
        return_details=True)
    eligible_contacts = sigma_contact_qc["selected_mask"]
    sig_fixed_parietal, parietal_contact_qc = _aggregate_full_night_power(
        sig_fixed_ch, eligible_channels=(parietal_contact_mask & nonflat_contacts),
        min_contact_coverage=MIN_CONTACT_COVERAGE,
        min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN,
        return_details=True)
    ep_dr, ep_swa, ep_clean, staging_contact_qc = aggregate_staging_features(
        ep_dr_ch, ep_swa_ch, ep_clean_ch, eligible_contacts,
        valid_window_count_by_contact=ep_valid_window_ch.sum(axis=2),
        min_valid_windows=STAGING_REFERENCE_MIN_VALID_WINDOWS)
    sig_fsp = _aggregate_full_night_power(
        sig_fsp_ch, eligible_channels=eligible_contacts,
        min_contact_coverage=MIN_CONTACT_COVERAGE, min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN)
    sig_fsp_parietal = _aggregate_full_night_power(
        sig_fsp_ch, eligible_channels=parietal_contact_qc["selected_mask"],
        min_contact_coverage=MIN_CONTACT_COVERAGE, min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN)
    swa_1 = _aggregate_full_night_power(
        swa_ch, eligible_channels=eligible_contacts,
        min_contact_coverage=MIN_CONTACT_COVERAGE, min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN)
    swa_parietal = _aggregate_full_night_power(
        swa_ch, eligible_channels=parietal_contact_qc["selected_mask"],
        min_contact_coverage=MIN_CONTACT_COVERAGE, min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN)
    # These two aggregate arrays remain the historical all-measured-epoch/audit80 reference.
    # Neutral per-channel values and complete-window support below permit less destructive offline
    # profiles without changing this compatibility output.
    ep_emg, emg_contact_qc = aggregate_auxiliary_epoch_features(
        ep_emg_ch,
        valid_window_count_by_channel=ep_emg_valid_window_ch.sum(axis=2),
        min_valid_windows=STAGING_REFERENCE_MIN_VALID_WINDOWS)
    ep_eog, eog_contact_qc = aggregate_auxiliary_epoch_features(
        ep_eog_ch,
        valid_window_count_by_channel=ep_eog_valid_window_ch.sum(axis=2),
        min_valid_windows=STAGING_REFERENCE_MIN_VALID_WINDOWS)

    # heart rate grids
    beats = sanitize_beats(beats)
    rr_1, rr_4, hr_1, hr_4 = interpolate_tachograms(beats, total_s)

    sigma_coverage = float(np.isfinite(sig_fixed).mean())
    hr_coverage = float(np.isfinite(hr_4).mean())
    if failed_chunks:
        raise RuntimeError(f"{len(failed_chunks)} acquisition chunks failed")
    qc_warnings = []
    if ecg_failures:
        rr_1[:] = np.nan
        rr_4[:] = np.nan
        hr_1[:] = np.nan
        hr_4[:] = np.nan
        hr_coverage = 0.0
        qc_warnings.append(
            f"cardiac series invalidated after {len(ecg_failures)} ECG detector exceptions")
    if hr_coverage < MIN_SIGNAL_COVERAGE:
        qc_warnings.append(
            f"HR coverage {hr_coverage:.1%} is below the historical audit80 reference")
    if len(ctx) < MIN_CONTACTS:
        qc_warnings.append(
            f"only {len(ctx)} anatomy-eligible contacts; historical audit80 required "
            f"{MIN_CONTACTS}")
    if not emg_indices or not eog_indices:
        qc_warnings.append(
            f"optional staging modalities missing: EMG={bool(emg_indices)}, "
            f"EOG={bool(eog_indices)}")

    proxy_lab, stage_counts = score_stages(ep_swa, ep_emg, ep_eog, ep_dr, ep_clean)
    annotation_lab, annotation_sources = annotation_epoch_context(
        annotations, n_ep, epoch_s=EPOCH)
    stage_lab, stage_sources = constrain_proxy_to_annotations(
        proxy_lab, annotation_lab, annotation_sources)
    stage_lab_proxy_sensitivity, _ = constrain_proxy_to_annotations(
        proxy_lab, annotation_lab, annotation_sources,
        allow_proxy_unknown_sleep=True)
    annotation_masks = annotation_second_masks(annotations, total_s)
    final_stage_counts = {
        value: int((stage_lab == value).sum())
        for value in ("W", "R", "NREM", "N2", "N3")
    }
    stage_counts["final_annotation_constrained_counts"] = final_stage_counts

    os.makedirs(OUT, exist_ok=True)
    source_hashes = {os.path.basename(path): file_sha256(path) for path in paths.values()}
    payload = dict(status="ok", cache_schema_version=CACHE_SCHEMA_VERSION,
                   cache_code_sha256=current_cache_digest,
                   generated_at_utc=utc_now(), code_revision=git_revision(ROOT),
                   code_dirty=git_is_dirty(ROOT), source_tree_sha256=source_tree_sha256(ROOT),
                   runtime_versions_json=json.dumps(runtime_versions(), sort_keys=True),
                   source_dataset="OpenNeuro ds003848 snapshot 1.0.3",
                   source_files_sha256_json=json.dumps(source_hashes, sort_keys=True),
                   source_files_snapshot_identity_json=json.dumps(
                       source_snapshot_identities, sort_keys=True),
                   failed_chunks_json=json.dumps(failed_chunks, sort_keys=True),
                   ecg_failures_json=json.dumps(ecg_failures, sort_keys=True),
                   qc_warnings_json=json.dumps(qc_warnings, sort_keys=True),
                   ecg_processing_method=(
                       "NeuroKit2 ecg_clean/ecg_peaks method=neurokit with artifact correction; "
                       "not Naji Pan-Tompkins 0.5-100 Hz; requires blinded R-peak validation"),
                   ecg_visual_validation=False,
                   raw_header_sampling_frequency=raw_sf,
                   bids_sampling_frequency=sf,
                   sample_timing_policy=(
                       "BIDS nominal SamplingFrequency after <=10-ppm agreement with BrainVision "
                       "header; global rounded sample boundaries"),
                   sigma_coverage=sigma_coverage, hr_coverage=hr_coverage,
                   sigma_meets_global_coverage_gate=bool(
                       sigma_coverage >= MIN_SIGNAL_COVERAGE),
                   sigma_per_contact_coverage=sigma_contact_qc["per_contact_coverage"],
                   sigma_selected_contact_mask=sigma_contact_qc["selected_mask"],
                   sigma_contact_count=sigma_contact_qc["contact_count"],
                   sigma_n_selected_contacts=sigma_contact_qc["n_selected"],
                   sigma_required_contact_count=sigma_contact_qc["required_contact_count"],
                   cortical_signal_nonflat_mask=nonflat_contacts,
                   cortical_signal_raw_minimum=channel_activity_qc["minimum"],
                   cortical_signal_raw_maximum=channel_activity_qc["maximum"],
                   cortical_signal_raw_dynamic_range=channel_activity_qc["dynamic_range"],
                   cortical_signal_numerical_flat_tolerance=(
                       channel_activity_qc["numerical_flat_tolerance"]),
                   cortical_signal_finite_sample_count=channel_activity_qc["finite_count"],
                   cortical_signal_activity_qc_method=channel_activity_qc["method"],
                   parietal_sigma_coverage=float(np.isfinite(sig_fixed_parietal).mean()),
                   parietal_selected_contact_mask=parietal_contact_qc["selected_mask"],
                   parietal_n_selected_contacts=parietal_contact_qc["n_selected"],
                   frontal_selected_contact_mask=(
                       frontal_contact_mask & nonflat_contacts
                       & sigma_contact_qc["selected_mask"]),
                   staging_selected_contact_mask=staging_contact_qc["selected_contact_mask"],
                   staging_contact_count=staging_contact_qc["contact_count"],
                   staging_n_selected_contacts=staging_contact_qc["n_selected_contacts"],
                   staging_required_contact_count=staging_contact_qc["required_contact_count"],
                   staging_per_contact_feature_coverage=(
                       staging_contact_qc["per_contact_feature_coverage"]),
                   staging_minimum_contact_feature_coverage=(
                       staging_contact_qc["minimum_contact_feature_coverage"]),
                   staging_minimum_valid_welch_windows=(
                       staging_contact_qc["minimum_valid_welch_windows"]),
                   staging_swa_normalization=staging_contact_qc["swa_normalization"],
                   emg_channel_names=np.asarray(
                       [names[index] for index in emg_indices], dtype="<U96"),
                   eog_channel_names=np.asarray(
                       [names[index] for index in eog_indices], dtype="<U96"),
                   emg_selected_channel_mask=emg_contact_qc["selected_channel_mask"],
                   eog_selected_channel_mask=eog_contact_qc["selected_channel_mask"],
                   emg_per_channel_coverage=emg_contact_qc["per_channel_coverage"],
                   eog_per_channel_coverage=eog_contact_qc["per_channel_coverage"],
                   auxiliary_channel_aggregation=emg_contact_qc["normalization"],
                   auxiliary_historical_minimum_valid_windows=(
                       STAGING_REFERENCE_MIN_VALID_WINDOWS),
                   auxiliary_window_duration_s=STAGING_WELCH_WINDOW_S,
                   auxiliary_window_overlap_s=STAGING_WELCH_OVERLAP_S,
                   auxiliary_epoch_value_semantics=(
                       "EMG RMS or EOG variance over the union of samples belonging to at "
                       "least one complete finite/measured 4-s window; overlapping samples "
                       "count once; fully measured epochs exactly equal the former full-epoch "
                       "feature; partial epochs intentionally omit incomplete-window samples"),
                   auxiliary_window_power_semantics=(
                       "EMG stores per-window mean square (take square root after averaging "
                       "to obtain RMS); EOG stores per-window variance; invalid windows are NaN"),
                   subject=subject, sf=sf, night_s=0.0, hours=total_s / 3600.0,
                   cortical_chans=np.asarray(ctx, dtype="<U96"),
                   all_good_ieeg_channels=np.asarray(
                       [names[index] for index in all_good_ie], dtype="<U96"),
                   destrieux_labels=destrieux_labels,
                   frontal_contact_mask=frontal_contact_mask,
                   parietal_contact_mask=parietal_contact_mask,
                   anatomy_exclusion_reasons_json=json.dumps(
                       anatomy_qc["exclusion_reasons"], sort_keys=True),
                   anatomy_selection_method=(
                       "exact channels/electrodes.tsv name join; BIDS-good iEEG; exclude SOZ, "
                       "resected, edge, silicon, screw, CSF, white matter, lesion, gliosis, "
                       "non-gray depth, and unknown cortical-atlas contacts; endpoint ROI masks "
                       "from Destrieux labels"),
                   ekg=names[ecg_i], fsp=fsp, fsp_is_real_peak=fsp_real,
                   sigma_fixed=sig_fixed, sigma_fsp=sig_fsp, swa=swa_1,
                   sigma_fixed_by_contact=sig_fixed_ch,
                   sigma_fsp_by_contact=sig_fsp_ch,
                   swa_by_contact=swa_ch,
                   sigma_fixed_power_numerator_by_contact=sig_fixed_num_ch,
                   sigma_fixed_clean_sample_count_by_contact=sig_fixed_den_ch,
                   sigma_fsp_power_numerator_by_contact=sig_fsp_num_ch,
                   sigma_fsp_clean_sample_count_by_contact=sig_fsp_den_ch,
                   swa_power_numerator_by_contact=swa_num_ch,
                   swa_clean_sample_count_by_contact=swa_den_ch,
                   power_samples_per_second=int(sf),
                   power_historical_minimum_clean_fraction_per_second=0.5,
                   power_support_semantics=(
                       "numerator=sum of clean squared-envelope samples in each nominal "
                       "1-s bin; denominator=count of those clean samples; historical "
                       "*_by_contact arrays require denominator >= 0.5*samples_per_second"),
                   sigma_fixed_parietal=sig_fixed_parietal,
                   sigma_fsp_parietal=sig_fsp_parietal,
                   swa_parietal=swa_parietal,
                   hr_1=hr_1, hr_4=hr_4, rr_1=rr_1, rr_4=rr_4, fs_rr=FS_RR,
                   ep_dr=ep_dr, ep_swa=ep_swa, ep_clean=ep_clean, ep_emg=ep_emg, ep_eog=ep_eog,
                   ep_dr_by_contact=ep_dr_ch,
                   ep_swa_by_contact=ep_swa_ch,
                   ep_clean_fraction_by_contact=ep_clean_ch,
                   ep_measured_fraction_by_contact=ep_measured_ch,
                   ep_valid_welch_window_mask_by_contact=ep_valid_window_ch,
                   ep_window_swa_power_by_contact=ep_window_swa_ch,
                   ep_window_total_power_by_contact=ep_window_total_ch,
                   ep_longest_clean_run_s_by_contact=ep_longest_clean_run_ch,
                   ep_valid_window_span_s_by_contact=ep_valid_window_span_ch,
                   event_3b_so_band_hz=np.asarray(SO_BAND_NAJI),
                   event_3b_negative_half_duration_s=np.asarray(
                       SO_NEGATIVE_HALF_DURATION_S),
                   event_3b_positive_half_max_s=SO_POSITIVE_HALF_MAX_S,
                   event_3b_amplitude_semantics=(
                       "Naji cites fixed Dang-Vu scalp-voltage gates, which do not "
                       "transfer to iEEG; cache retains duration-qualified down/up "
                       "amplitudes before an offline within-contact/stage percentile "
                       "sensitivity rule"),
                   ep_emg_by_channel=ep_emg_ch, ep_eog_by_channel=ep_eog_ch,
                   ep_emg_valid_window_mask_by_channel=ep_emg_valid_window_ch,
                   ep_eog_valid_window_mask_by_channel=ep_eog_valid_window_ch,
                   ep_emg_window_mean_square_by_channel=ep_emg_window_power_ch,
                   ep_eog_window_variance_by_channel=ep_eog_window_power_ch,
                   epoch_s=EPOCH, beats=beats, stage_lab=np.asarray(stage_lab, dtype="<U4"),
                   stage_lab_proxy=np.asarray(proxy_lab, dtype="<U4"),
                   stage_lab_proxy_sensitivity=np.asarray(
                       stage_lab_proxy_sensitivity, dtype="<U4"),
                   stage_lab_annotation=np.asarray(annotation_lab, dtype="<U5"),
                   stage_source=np.asarray(stage_sources, dtype="<U32"),
                   has_eog=bool(eog_indices), has_emg=bool(emg_indices),
                   stage_method=(
                       "author RESPect sleep/NREM/REM/SWS/transition/awake annotations are "
                       "primary; author-unknown sleep remains unclassified; no AASM N2/N3 "
                       "scoring; proxy-within-unknown labels stored for sensitivity only"),
                   normalized_event_annotations_json=json.dumps(
                       annotations, sort_keys=True),
                   annotation_signal_exclusion_policy=(
                       "author artifact intervals per named/all contact; seizure/stimulation "
                       "intervals globally; exact intervals set missing before filtering; "
                       "prepare_continuous_signal supplies the recorded +/-5 s missing-data pad"),
                   annotation_signal_pad_s=ANNOTATION_SIGNAL_PAD_S,
                   stage_proxy_diagnostics_json=json.dumps(stage_counts, sort_keys=True))
    for key, values in annotation_masks.items():
        payload[f"annotation_{key}_mask_1s"] = values
    for c in ctx:
        values = np.asarray(sorted(so_candidates[c]), float).reshape(-1, 4)
        payload[f"so_candidate_t_{c}"] = values[:, 0]
        payload[f"so_candidate_down_{c}"] = values[:, 1]
        payload[f"so_candidate_up_{c}"] = values[:, 2]
        payload[f"so_candidate_p2p_{c}"] = values[:, 3]
    atomic_savez(fp, **payload)
    frac = float(np.isfinite(sig_fixed).mean())
    print(f"[{subject}] cached {time.time()-t0:.0f}s | sigma_cov {frac:.0%} | {len(beats)} beats | "
          f"annotation-constrained W/R/NREM/N2/N3 = "
          f"{final_stage_counts['W']}/{final_stage_counts['R']}/"
          f"{final_stage_counts['NREM']}/{final_stage_counts['N2']}/"
          f"{final_stage_counts['N3']} of {n_ep} ep", flush=True)
    if delete_raw:
        try:
            os.remove(paths["_ieeg.eeg"])
            print(f"[{subject}] removed raw .eeg after checksum was stored", flush=True)
        except OSError:
            pass
    return "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default=",".join(SUBJECTS))
    ap.add_argument("--delete-raw", action="store_true",
                    help="delete the multi-GB .eeg only after its SHA-256 is stored; default keeps raw")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    requested = [s.strip() for s in a.subjects.split(",") if s.strip()]
    config = dict(cache_schema_version=CACHE_SCHEMA_VERSION,
                  cache_code_sha256=cache_code_sha256(ROOT),
                  source_snapshot="ds003848/1.0.3",
                  source_identity_manifest_sha256=file_sha256(
                      SNAPSHOT_IDENTITY_PATH))
    if not a.force and validated_complete_run_exists(
            OUT, pipeline="stage_ds003848", requested=requested, config=config,
            suffix=".npz", require_current_source_tree=False):
        print(f"validated existing complete cache run ({len(requested)} subjects)", flush=True)
        return
    run_id = start_run_manifest(
        OUT, pipeline="stage_ds003848", requested=requested, config=config)
    completed, skipped, failed = [], [], []
    for s in requested:
        if s not in SUBJECTS:
            error = "unknown subject"
            print(f"[{s}] {error}", flush=True)
            failed.append(dict(subject=s, error=error))
            continue
        try:
            status = run(s, delete_raw=a.delete_raw, force=a.force)
            if status == "skip":
                with np.load(os.path.join(OUT, f"{s}.npz"), allow_pickle=False) as record:
                    reason = npz_scalar_text(record, "reason", "unspecified")
                skipped.append(dict(subject=s, reason=reason))
            else:
                completed.append(s)
        except Exception as e:
            import traceback
            print(f"[{s}] ERROR {type(e).__name__}: {e}", flush=True)
            traceback.print_exc()
            failed.append(dict(subject=s, error=f"{type(e).__name__}: {e}"))
    write_run_manifest(
        OUT, pipeline="stage_ds003848", requested=requested, completed=completed,
        skipped=skipped, failed=failed,
        config=config, run_id=run_id,
        result_files_sha256={
            subject: file_sha256(os.path.join(OUT, f"{subject}.npz"))
            for subject in completed + [value["subject"] for value in skipped]
        })
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
