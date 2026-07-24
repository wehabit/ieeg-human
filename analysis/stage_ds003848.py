"""OpenNeuro ds003848 (Utrecht RESPect long-term iEEG) -> staged derived-series cache.

The independent replication cohort for 3A/3B. Unlike HUP, every subject carries EMG + EOG, allowing
rule-based REM/wake exclusion instead of the GMM-on-slow-wave proxy. This is not expert AASM scoring.

Six patients, ~1 h continuous `task-[Ss]leep` runs @ 2048 Hz, 50 Hz line. Verified against the raw
channels.tsv (2026-07): all six have iEEG (3 ECoG grid + 3 SEEG depth), ECG, EMG, EOG, and (bad)
respiration belts. Channel ROLE is taken from channels.tsv by row index, which BIDS guarantees
matches the data-column order, so naming variants (ECG+/ecg1+, emg+/EMG2, orb+/Orb+) don't matter.

Pipeline per subject: download the BrainVision triple if absent -> stream in chunks with MNE ->
derive the SAME series the HUP cache stores (so lecci_faithful_3A / event_3B_cached consume it
unchanged) PLUS EMG/EOG staging features and a scored `stage_lab` -> save data/derived/ds003848/.
The raw .eeg (~3.9 GB) is optionally deleted after caching so disk stays bounded.

STAGING (rule-based, honest about its limits). True AASM scoring needs scalp EEG, which this dataset
lacks. What EMG + EOG buy is REM/Wake exclusion, and that is what this does:
  * per 30 s epoch: submental EMG RMS (notch + 10-100 Hz), EOG movement variance (0.3-6 Hz),
    slow-wave power (0.5-4 Hz, aggregated after per-contact power extraction), all robust-z-scored
    within subject
  * Wake  = high muscle tone (EMG z > 1.0)
  * REM   = muscle atonia (EMG z < -0.3) + low SWA + phasic eye movement (EOG z > 0.5)
  * NREM  = everything else with adequate delta; split N2/N3 by a 2-component GMM on log SWA
The N2/N3 split is still SWA-driven (as in Lecci's S2/SWS axis), but now on epochs from which Wake
and REM have been removed -- which the HUP proxy could not do.

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
                             SWA_BAND_L, SO_BAND_NAJI, FS_RR, sanitize_beats,
                             _aggregate_full_night_power, MIN_SIGNAL_COVERAGE, FILTER_EDGE_S,
                             MIN_CONTACT_COVERAGE, MIN_CONTACTS,
                             MIN_CONTACT_FRACTION_PER_BIN, prepare_continuous_signal,
                             aggregate_staging_features, interpolate_tachograms,
                             staging_epoch_features)
from results_3A_tutorial_style import ied_clean_mask
from pipeline_version import (CACHE_SCHEMA_VERSION, atomic_savez, cache_code_sha256, git_is_dirty,
                              git_revision, npz_scalar_text, runtime_versions,
                              source_tree_sha256, utc_now, start_run_manifest,
                              validated_complete_run_exists, write_run_manifest)

mne.set_log_level("ERROR")
MIN_NREM_DELTA_RATIO = 0.20

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


# ---------------------------------------------------------------- download
def _dl(sub, ses, base, ext):
    url = f"{BASE}/{sub}:{ses}:ieeg:{base}{ext}"
    dst = os.path.join(RAW, f"{base}{ext}")
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        return dst
    os.makedirs(RAW, exist_ok=True)
    tmp = dst + ".part"
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as f:
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            f.write(b)
    os.replace(tmp, dst)
    return dst


def ensure_files(sub, ses, base):
    paths = {}
    for ext in ("_channels.tsv", "_ieeg.json", "_ieeg.vhdr", "_ieeg.vmrk", "_ieeg.eeg"):
        paths[ext] = _dl(sub, ses, base, ext)
    return paths


def file_sha256(path, block_size=1 << 20):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


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
    roles = channel_roles(paths["_channels.tsv"])
    if not roles["ieeg"] or not roles["ecg"] or not roles["emg"] or not roles["eog"]:
        reason = (f"required modalities missing: iEEG={bool(roles['ieeg'])}, "
                  f"ECG={bool(roles['ecg'])}, EMG={bool(roles['emg'])}, "
                  f"EOG={bool(roles['eog'])}")
        atomic_savez(fp, subject=subject, status="skip",
                     cache_schema_version=CACHE_SCHEMA_VERSION,
                     cache_code_sha256=current_cache_digest, reason=reason)
        print(f"[{subject}] SKIP {reason}", flush=True)
        return "skip"

    raw = mne.io.read_raw_brainvision(paths["_ieeg.vhdr"], preload=False, verbose="ERROR")
    sf = float(raw.info["sfreq"])
    n_samp = raw.n_times
    if len(raw.ch_names) != roles["n_rows"]:
        raise RuntimeError(
            f"channel count mismatch: raw {len(raw.ch_names)} vs TSV {roles['n_rows']}")
    mismatched = [
        (i, raw_name, tsv_name)
        for i, (raw_name, tsv_name) in enumerate(zip(raw.ch_names, roles["names"]))
        if raw_name.strip().casefold() != tsv_name.strip().casefold()
    ]
    if mismatched:
        raise RuntimeError(
            f"BrainVision/TSV channel-order mismatch; first mismatch={mismatched[0]}")
    ie = roles["ieeg"]; ecg_i = roles["ecg"][0]
    emg_i = roles["emg"][0] if roles["emg"] else None
    eog_i = roles["eog"][0] if roles["eog"] else None
    names = roles["names"]
    ctx = [names[i] for i in ie]

    total_s = int(n_samp / sf)
    n_ep = int(total_s // EPOCH)

    # Individual FSP is disabled until it is estimated from all clean NREM and manually QC'd, as
    # in Lecci. The fixed 10-15 Hz analysis is primary.
    fsp, fsp_real = 13.0, False
    print(f"[{subject}] {sf:.0f} Hz, {len(ie)} iEEG, ECG@{ecg_i}, EMG@{emg_i}, EOG@{eog_i}, "
          f"FSP {fsp:.2f} Hz, {total_s} s", flush=True)

    sig_fixed_ch = np.full((len(ie), total_s), np.nan)
    sig_fsp_ch = np.full((len(ie), total_s), np.nan)
    swa_ch = np.full((len(ie), total_s), np.nan)
    ep_dr_ch = np.full((len(ie), n_ep), np.nan)
    ep_swa_ch = np.full((len(ie), n_ep), np.nan)
    ep_clean_ch = np.full((len(ie), n_ep), np.nan)
    ep_emg = np.full(n_ep, np.nan); ep_eog = np.full(n_ep, np.nan)
    beats = []
    so_candidates = {c: [] for c in ctx}
    failed_chunks = []
    ecg_failures = []

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

    picks_all = list(ie) + [ecg_i] + ([emg_i] if emg_i is not None else []) + \
                ([eog_i] if eog_i is not None else [])
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
        col = len(ie) + 1
        x_emg = d[col] if emg_i is not None else None
        x_eog = d[col + (1 if emg_i is not None else 0)] if eog_i is not None else None

        # Keep contacts separate through filtering, artifact masking, and power extraction.
        # Raw-voltage averaging can cancel equal power with opposite polarity/phase.
        prepared = [prepare_continuous_signal(x_ie[ci], sf) for ci in range(len(ie))]
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
            ecg_filled, ecg_measured = prepare_continuous_signal(x_ecg, sf)
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
        for sos_b, dest in ((sos_fixed, sig_fixed_ch), (sos_fsp, sig_fsp_ch),
                            (sos_swa, swa_ch)):
            for ci, (xc, clean_c) in enumerate(zip(x_channels, clean_channels)):
                env2 = np.abs(signal.hilbert(signal.sosfiltfilt(sos_b, xc))) ** 2
                start = core_a
                stop = start + n_sec * k
                e2 = env2[start:stop].reshape(n_sec, k)
                cm = clean_c[start:stop].reshape(n_sec, k).astype(float)
                num, den = (e2 * cm).sum(1), cm.sum(1)
                vals_c = np.where(
                    den >= 0.5 * k, num / np.maximum(den, 1e-12), np.nan)
                sl = slice(off, min(off + n_sec, total_s))
                dest[ci, sl] = vals_c[:sl.stop - sl.start]

        # EMG / EOG epoch features
        if x_emg is not None:
            emg_filled, emg_measured = prepare_continuous_signal(x_emg, sf)
            emg_f = np.abs(signal.sosfiltfilt(sos_emg, notch50(emg_filled)))
        else:
            emg_f, emg_measured = None, None
        if x_eog is not None:
            eog_filled, eog_measured = prepare_continuous_signal(x_eog, sf)
            eog_f = signal.sosfiltfilt(sos_eog, eog_filled)
        else:
            eog_f, eog_measured = None, None

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
                dr, swa_value = staging_epoch_features(
                    seg, clean_c[epoch_a:epoch_b], measured_c[epoch_a:epoch_b], sf)
                ep_dr_ch[ci, gi] = dr
                ep_swa_ch[ci, gi] = swa_value
                ep_clean_ch[ci, gi] = float(clean_c[epoch_a:epoch_b].mean())
            if emg_f is not None and emg_measured[epoch_a:epoch_b].all():
                ep_emg[gi] = float(np.sqrt(np.mean(emg_f[epoch_a:epoch_b] ** 2)))
            if eog_f is not None and eog_measured[epoch_a:epoch_b].all():
                ep_eog[gi] = float(np.var(eog_f[epoch_a:epoch_b]))
        t += dur

    sig_fixed, sigma_contact_qc = _aggregate_full_night_power(
        sig_fixed_ch, min_contact_coverage=MIN_CONTACT_COVERAGE,
        min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN,
        return_details=True)
    eligible_contacts = sigma_contact_qc["selected_mask"]
    ep_dr, ep_swa, ep_clean, staging_contact_qc = aggregate_staging_features(
        ep_dr_ch, ep_swa_ch, ep_clean_ch, eligible_contacts)
    sig_fsp = _aggregate_full_night_power(
        sig_fsp_ch, eligible_channels=eligible_contacts,
        min_contact_coverage=MIN_CONTACT_COVERAGE, min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN)
    swa_1 = _aggregate_full_night_power(
        swa_ch, eligible_channels=eligible_contacts,
        min_contact_coverage=MIN_CONTACT_COVERAGE, min_contacts=MIN_CONTACTS,
        min_contact_fraction_per_bin=MIN_CONTACT_FRACTION_PER_BIN)

    # heart rate grids
    beats = sanitize_beats(beats)
    rr_1, rr_4, hr_1, hr_4 = interpolate_tachograms(beats, total_s)

    sigma_coverage = float(np.isfinite(sig_fixed).mean())
    hr_coverage = float(np.isfinite(hr_4).mean())
    if failed_chunks:
        raise RuntimeError(f"{len(failed_chunks)} acquisition chunks failed")
    if ecg_failures:
        raise RuntimeError(
            f"ECG detection failed in {len(ecg_failures)} chunks; publication cache fails closed")
    if sigma_coverage < MIN_SIGNAL_COVERAGE or hr_coverage < MIN_SIGNAL_COVERAGE:
        raise RuntimeError(
            f"coverage QC failed: sigma={sigma_coverage:.1%}, HR={hr_coverage:.1%}; "
            f"minimum is {MIN_SIGNAL_COVERAGE:.0%}")

    stage_lab, stage_counts = score_stages(ep_swa, ep_emg, ep_eog, ep_dr, ep_clean)

    os.makedirs(OUT, exist_ok=True)
    source_hashes = {os.path.basename(path): file_sha256(path) for path in paths.values()}
    payload = dict(status="ok", cache_schema_version=CACHE_SCHEMA_VERSION,
                   cache_code_sha256=current_cache_digest,
                   generated_at_utc=utc_now(), code_revision=git_revision(ROOT),
                   code_dirty=git_is_dirty(ROOT), source_tree_sha256=source_tree_sha256(ROOT),
                   runtime_versions_json=json.dumps(runtime_versions(), sort_keys=True),
                   source_dataset="OpenNeuro ds003848 snapshot 1.0.3",
                   source_files_sha256_json=json.dumps(source_hashes, sort_keys=True),
                   failed_chunks_json=json.dumps(failed_chunks, sort_keys=True),
                   ecg_failures_json=json.dumps(ecg_failures, sort_keys=True),
                   ecg_processing_method=(
                       "NeuroKit2 ecg_clean/ecg_peaks method=neurokit with artifact correction; "
                       "not Naji Pan-Tompkins 0.5-100 Hz; requires blinded R-peak validation"),
                   ecg_visual_validation=False,
                   sigma_coverage=sigma_coverage, hr_coverage=hr_coverage,
                   sigma_per_contact_coverage=sigma_contact_qc["per_contact_coverage"],
                   sigma_selected_contact_mask=sigma_contact_qc["selected_mask"],
                   sigma_contact_count=sigma_contact_qc["contact_count"],
                   sigma_n_selected_contacts=sigma_contact_qc["n_selected"],
                   sigma_required_contact_count=sigma_contact_qc["required_contact_count"],
                   staging_selected_contact_mask=staging_contact_qc["selected_contact_mask"],
                   staging_contact_count=staging_contact_qc["contact_count"],
                   staging_n_selected_contacts=staging_contact_qc["n_selected_contacts"],
                   staging_required_contact_count=staging_contact_qc["required_contact_count"],
                   staging_swa_normalization=staging_contact_qc["swa_normalization"],
                   subject=subject, sf=sf, night_s=0.0, hours=total_s / 3600.0,
                   cortical_chans=np.array(ctx),
                   anatomy_selection_method=(
                       "BIDS good iEEG channel type only; no homologous-region/gray-matter QC"),
                   ekg=names[ecg_i], fsp=fsp, fsp_is_real_peak=fsp_real,
                   sigma_fixed=sig_fixed, sigma_fsp=sig_fsp, swa=swa_1,
                   hr_1=hr_1, hr_4=hr_4, rr_1=rr_1, rr_4=rr_4, fs_rr=FS_RR,
                   ep_dr=ep_dr, ep_swa=ep_swa, ep_clean=ep_clean, ep_emg=ep_emg, ep_eog=ep_eog,
                   epoch_s=EPOCH, beats=beats, stage_lab=np.array(stage_lab),
                   has_eog=bool(eog_i is not None), has_emg=bool(emg_i is not None),
                   stage_method="unvalidated EMG/EOG/iEEG rule-based proxy; not AASM scored",
                   stage_proxy_diagnostics_json=json.dumps(stage_counts, sort_keys=True))
    for c in ctx:
        values = np.asarray(sorted(so_candidates[c]), float).reshape(-1, 4)
        payload[f"so_candidate_t_{c}"] = values[:, 0]
        payload[f"so_candidate_down_{c}"] = values[:, 1]
        payload[f"so_candidate_up_{c}"] = values[:, 2]
        payload[f"so_candidate_p2p_{c}"] = values[:, 3]
    atomic_savez(fp, **payload)
    frac = float(np.isfinite(sig_fixed).mean())
    print(f"[{subject}] cached {time.time()-t0:.0f}s | sigma_cov {frac:.0%} | {len(beats)} beats | "
          f"stages W/R/N2/N3 = {stage_counts['n_wake']}/{stage_counts['n_rem']}/"
          f"{int((stage_lab=='N2').sum())}/{int((stage_lab=='N3').sum())} of {n_ep} ep", flush=True)
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
                  source_snapshot="ds003848/1.0.3")
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
