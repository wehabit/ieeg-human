"""Reusable Lecci-motivated 3A estimators for sigma-infraslow coupling.

WHAT THE PREVIOUS 3A ACTUALLY TESTED, AND WHY IT IS NOT LECCI'S TEST
-------------------------------------------------------------------
The cohort 3A asked: "is magnitude-squared coherence between sigma power and heart rate significant
in the single bin at 0.0195 Hz?" Lecci never ran that test. Their human analysis is two steps:

  Step 1 (Fig 1C/G) -- does sigma power oscillate infraslow AT ALL, and at what rate?
      Morlet wavelet power time course -> a SECOND wavelet transform of that time course at 0.001 Hz
      resolution over 0.001-0.12 Hz -> one spectrum per NREM bout >=120 s -> averaged across ALL
      bouts WEIGHTED BY DURATION -> normalised to its own mean -> Gaussian fit to locate THAT
      SUBJECT'S OWN PEAK. The reported 0.019 +/- 0.001 Hz (n=27) is the mean of those fitted peaks.
  Step 2 (Fig 6) -- does heart rate track it?
      Sigma and HR resampled to a common grid, 120 s intervals z-transformed, CROSS-CORRELATION with
      heart rate as source wave, averaged within subject then across subjects. The result is a lag,
      not a coherence value.

Testing one hard-coded bin is materially different from that method. This script implements an
explicit approximation of both steps and records the remaining deviations below. It must not be
described as an exact replication of the FieldTrip pipeline.

DEVIATIONS (unavoidable, stated rather than hidden)
  * Staging: HUP uses an unvalidated GMM delta/SWA proxy. RESPect uses author coarse
    NREM/REM/SWS annotations as primary, but still lacks expert AASM N2/N3 scoring. Lecci's
    S2-versus-SWS contrast is therefore not directly replicated; pooled author NREM is primary
    where available.
  * Signal: sigma power comes from a Hilbert envelope of a Butterworth band, binned to 1 s, rather
    than a 4-cycle Morlet sampled at 0.1 s then smoothed with a 4 s moving average. The production
    path requires 4 s smoothing, but the 1 Hz Hilbert series still needs a direct benchmark against
    a reference implementation before equivalence can be claimed.
  * Region: HUP uses unvalidated lateral-contact candidates. RESPect uses non-pathological
    parietal/postcentral/precuneus iEEG contacts. Neither is Lecci's scalp C3 source; both are
    motivated adaptations.
  * Population: epilepsy patients on anti-seizure medication.

The standalone writer is withdrawn. Current endpoint execution and artifact
materialization go through ``run_qc_grid.py``.
"""
import json
import os
import numpy as np
from scipy import signal, optimize

from staging_helpers import stage_epochs, EPOCH
from spectral_gapped import (coherence_gapped, analytic_msc_threshold, fill_short_gaps,
                             contiguous_runs)
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    cache_code_sha256,
    file_sha256,
    finite_float_or_none,
    npz_scalar_text,
    runtime_versions,
    validate_completed_cache_failures,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "data", "derived", "lc_infraslow")

FS = 1.0                      # cached derived-series rate
MIN_BOUT_S = 120.0            # Lecci: bouts >= 120 s (>= 4 epochs)
F_LO, F_HI, F_STEP = 0.001, 0.120, 0.001      # Lecci: 0.001-0.12 Hz at 0.001 Hz resolution
N_CYCLES = 4.0                # Lecci baseline; seven cycles was a sensitivity analysis
PEAK_SEARCH = (0.008, 0.060)  # where an LC-type infraslow peak is allowed to be
F_LECCI = 0.02
NPERSEG = 256
ALPHA = 0.05
XCORR_WIN_S = 120.0           # Lecci: 120 s intervals, z-transformed
XCORR_MAX_LAG_S = 60.0
CORE_STUDY_WINDOW_S = 210.0 * 60.0  # Lecci core study: first 210 min from sleep onset
LECCI_XCORR_LAG_WINDOW_S = (0.0, 15.0)  # source Fig. 6: positive peak near +5 s
MIN_SURROGATE_PEAK_LOCATIONS = 20
MIN_SURROGATE_PEAK_FRACTION = 0.05
MIN_POWER_CONTACTS = 3
MIN_POWER_COVERAGE = 0.80


# ------------------------------------------------------------------ cache access
_CACHE_DIR = [CACHE]          # mutable so callers can point at an alternate cache (e.g. ds003848)


class CacheSubjectSkipped(RuntimeError):
    """Structured, expected cache exclusion rather than an analysis failure."""


def set_cache(path):
    _CACHE_DIR[0] = path


_CACHE_LINEAGE_FIELDS = (
    "cache_directory_relative",
    "cache_manifest_run_id",
    "cache_manifest_sha256",
    "cache_file_sha256",
    "cache_hours",
)
_CACHE_PIPELINES = {"cache_lc_series", "stage_ds003848"}


def _validated_cache_manifest(subject):
    """Validate the terminal cache manifest and the exact NPZ bytes for one subject."""
    cache_dir = os.path.abspath(_CACHE_DIR[0])
    manifest_path = os.path.join(cache_dir, "RUN_MANIFEST.json")
    cache_path = os.path.join(cache_dir, f"{subject}.npz")
    if not os.path.exists(manifest_path) or not os.path.exists(cache_path):
        raise RuntimeError(f"cache lineage input is missing for {subject}")
    try:
        with open(manifest_path) as handle:
            manifest = json.load(handle)
    except Exception as exc:
        raise RuntimeError(f"{manifest_path} cannot be read") from exc

    requested_values = manifest.get("requested", [])
    requested = set(requested_values)
    completed = set(manifest.get("completed", []))
    skipped_entries = manifest.get("skipped", [])
    if any(not isinstance(value, dict) or not value.get("subject") or not value.get("reason")
           for value in skipped_entries):
        raise RuntimeError(f"{manifest_path} has malformed skip entries")
    skipped = {value["subject"] for value in skipped_entries}
    if (
        manifest.get("run_state") != "complete"
        or manifest.get("pipeline") not in _CACHE_PIPELINES
        or manifest.get("analysis_version") != ANALYSIS_VERSION
        or manifest.get("cache_schema_version") != CACHE_SCHEMA_VERSION
        or manifest.get("runtime_versions") != runtime_versions()
        or manifest.get("failed")
        or len(requested_values) != len(requested)
        or completed | skipped != requested
        or completed & skipped
        or subject not in requested
    ):
        raise RuntimeError(
            f"{manifest_path} is incomplete, internally inconsistent, or not a current "
            "production cache run")
    hashes = manifest.get("result_files_sha256")
    if not isinstance(hashes, dict) or set(hashes) != requested:
        raise RuntimeError(f"{manifest_path} lacks exact per-cache file hashes")
    actual_hash = file_sha256(cache_path)
    if hashes.get(subject) != actual_hash:
        raise RuntimeError(
            f"{cache_path} bytes do not match the terminal cache manifest")
    return manifest


def cache_lineage(subject, d=None):
    """Immutable identity of the exact cache file and cache run consumed for one subject."""
    cache_dir = os.path.abspath(_CACHE_DIR[0])
    relative = os.path.relpath(cache_dir, ROOT)
    if relative == os.pardir or relative.startswith(os.pardir + os.sep):
        raise RuntimeError("production cache directory must live inside the repository root")
    cache_path = os.path.join(cache_dir, f"{subject}.npz")
    manifest_path = os.path.join(cache_dir, "RUN_MANIFEST.json")
    manifest = _validated_cache_manifest(subject)
    hours = None
    if d is not None and "hours" in getattr(d, "files", []):
        hours = float(np.asarray(d["hours"]).item())
    elif d is None:
        with np.load(cache_path, allow_pickle=False) as cached:
            if "hours" in cached.files:
                hours = float(np.asarray(cached["hours"]).item())
    return dict(
        cache_directory_relative=relative.replace(os.sep, "/"),
        cache_manifest_run_id=manifest.get("run_id"),
        cache_manifest_sha256=file_sha256(manifest_path),
        cache_file_sha256=file_sha256(cache_path),
        cache_hours=hours,
    )


def cache_lineage_entry(record):
    """Canonical manifest entry extracted from a downstream subject record."""
    return {key: record.get(key) for key in _CACHE_LINEAGE_FIELDS}


def verify_cache_lineage(record):
    """Fail if current cache artifacts differ from the immutable inputs recorded downstream."""
    expected = cache_lineage_entry(record)
    required = [key for key in _CACHE_LINEAGE_FIELDS if key != "cache_hours"]
    if any(expected.get(key) in (None, "") for key in required):
        raise RuntimeError(f"{record.get('subject')} has incomplete cache lineage")
    cache_dir = os.path.realpath(os.path.join(ROOT, expected["cache_directory_relative"]))
    if os.path.commonpath([os.path.realpath(ROOT), cache_dir]) != os.path.realpath(ROOT):
        raise RuntimeError("cache lineage path escapes the repository root")
    manifest_path = os.path.join(cache_dir, "RUN_MANIFEST.json")
    cache_path = os.path.join(cache_dir, f"{record['subject']}.npz")
    if not os.path.exists(manifest_path) or not os.path.exists(cache_path):
        raise RuntimeError(f"current cache artifact is missing for {record['subject']}")
    manifest = json.load(open(manifest_path))
    if manifest.get("run_id") != expected["cache_manifest_run_id"]:
        raise RuntimeError(f"cache manifest run changed for {record['subject']}")
    if file_sha256(manifest_path) != expected["cache_manifest_sha256"]:
        raise RuntimeError(f"cache manifest content changed for {record['subject']}")
    if file_sha256(cache_path) != expected["cache_file_sha256"]:
        raise RuntimeError(f"cache file content changed for {record['subject']}")
    with np.load(cache_path, allow_pickle=False) as current:
        cached_subject = npz_scalar_text(current, "subject")
        if cached_subject != record["subject"]:
            raise RuntimeError(
                f"cache subject {cached_subject or '<missing>'} does not match "
                f"requested {record['subject']}")
        hours = (
            float(np.asarray(current["hours"]).item())
            if "hours" in current.files else None)
    if ((hours is None) != (expected["cache_hours"] is None)
            or (hours is not None
                and not np.isclose(hours, float(expected["cache_hours"])))):
        raise RuntimeError(f"cache duration changed for {record['subject']}")
    return expected


def stages_for(d):
    """Return labels/NREM from author-constrained RESPect states or the HUP stage proxy.

    RESPect's coarse NREM/REM/SWS labels are author annotations, not expert AASM N2/N3 scoring.
    HUP uses an unvalidated GMM slow-wave proxy, whose split labels remain N2-like/N3-like.
    """
    if "stage_lab" in getattr(d, "files", []):
        raw = np.asarray(d["stage_lab"]).astype(str)
        lab = np.full(len(raw), "", dtype=object)
        lab[raw == "N2"] = "N2"
        lab[raw == "N3"] = "N3"
        lab[raw == "NREM"] = "NREM"
        # REM is not NREM, but it is a classified sleep epoch and therefore
        # must survive normalization so ``core_study_nrem_mask`` can anchor
        # the first-210-minute window at the first available sleep label.
        lab[raw == "R"] = "R"
        nrem = (lab == "N2") | (lab == "N3") | (lab == "NREM")
        return lab, nrem, None
    ep = dict(dr=d["ep_dr"], swa=d["ep_swa"], clean=d["ep_clean"])
    return stage_epochs(ep)


def load(subject):
    fp = os.path.join(_CACHE_DIR[0], f"{subject}.npz")
    if not os.path.exists(fp):
        return None
    manifest_path = os.path.join(_CACHE_DIR[0], "RUN_MANIFEST.json")
    manifest = _validated_cache_manifest(subject)
    current_cache_digest = cache_code_sha256(ROOT)
    if manifest.get("config", {}).get("cache_code_sha256") != current_cache_digest:
        raise RuntimeError(
            f"{manifest_path} was produced by different cache-building source; rebuild with --force")
    skipped_entries = manifest.get("skipped", [])
    skipped = {
        value.get("subject"): value.get("reason")
        for value in skipped_entries if isinstance(value, dict)
    }
    if subject in skipped:
        with np.load(fp, allow_pickle=False) as skipped_cache:
            cached_subject = npz_scalar_text(skipped_cache, "subject")
            status = npz_scalar_text(skipped_cache, "status")
            schema = npz_scalar_text(skipped_cache, "cache_schema_version")
            digest = npz_scalar_text(skipped_cache, "cache_code_sha256")
            reason = npz_scalar_text(
                skipped_cache, "reason", skipped[subject] or "unspecified")
        if (cached_subject != subject or status != "skip"
                or schema != CACHE_SCHEMA_VERSION
                or digest != current_cache_digest or reason != skipped[subject]):
            raise RuntimeError(f"{fp} does not match its manifest skip entry")
        raise CacheSubjectSkipped(reason)
    if subject not in set(manifest.get("completed", [])):
        raise RuntimeError(
            f"{subject} is neither completed nor explicitly skipped in {manifest_path}")
    d = np.load(fp, allow_pickle=False)
    cached_subject = npz_scalar_text(d, "subject")
    if cached_subject != subject:
        d.close()
        raise RuntimeError(
            f"{fp} embeds subject {cached_subject or '<missing>'}, not requested {subject}")
    if npz_scalar_text(d, "status") != "ok":
        d.close()
        return None
    schema = npz_scalar_text(d, "cache_schema_version")
    if schema != CACHE_SCHEMA_VERSION:
        d.close()
        raise RuntimeError(
            f"{fp} has cache schema {schema or '<missing>'}; expected {CACHE_SCHEMA_VERSION}. "
            "Rebuild the cache with --force.")
    cached_digest = npz_scalar_text(d, "cache_code_sha256")
    if cached_digest != current_cache_digest:
        d.close()
        raise RuntimeError(
            f"{fp} was produced by different cache-building source; rebuild with --force")
    try:
        validate_completed_cache_failures(
            d, source=fp, require_ecg=True)
    except Exception:
        d.close()
        raise
    if "sigma_fixed" not in d.files or "hr_1" not in d.files or "rr_4" not in d.files:
        d.close()
        return None
    return d


def nrem_bouts(nrem, min_s=MIN_BOUT_S):
    """Contiguous NREM runs, returned as (start_second, stop_second) and >= min_s long."""
    out = []
    for i0, i1 in contiguous_runs(nrem):
        s0, s1 = int(i0 * EPOCH), int(i1 * EPOCH)
        if s1 - s0 >= min_s:
            out.append((s0, s1))
    return out


def core_study_nrem_mask(labels, epoch_s=EPOCH, window_s=CORE_STUDY_WINDOW_S):
    """Restrict labels to Lecci's core-study interval without selecting on signal values.

    Lecci defined sleep onset as the first S1 epoch followed by S2.  These datasets do not
    contain expert S1/S2 scoring, so the reproducible adaptation starts at the first available
    classified sleep epoch (R/NREM/N2/N3).  The approximation and exact epoch bounds are returned
    for every participant rather than silently analysing the full seven-hour HUP window.
    """
    labels = np.asarray(labels).astype(str)
    sleep = np.isin(labels, ("R", "NREM", "N2", "N3"))
    nrem = np.isin(labels, ("NREM", "N2", "N3"))
    found = np.flatnonzero(sleep)
    if not len(found):
        return np.zeros(len(labels), bool), dict(
            start_epoch=None,
            stop_epoch=None,
            requested_window_seconds=float(window_s),
            available_window_seconds=0.0,
            onset_proxy="no classified sleep epoch",
        )
    start = int(found[0])
    n_window_epochs = int(np.ceil(float(window_s) / float(epoch_s)))
    stop = min(len(labels), start + n_window_epochs)
    within = np.zeros(len(labels), bool)
    within[start:stop] = True
    return nrem & within, dict(
        start_epoch=start,
        stop_epoch=stop,
        start_second=float(start * epoch_s),
        stop_second=float(stop * epoch_s),
        requested_window_seconds=float(window_s),
        available_window_seconds=float((stop - start) * epoch_s),
        onset_proxy=(
            "first available classified R/NREM/N2/N3 epoch; expert S1 followed by S2 "
            "is unavailable"),
    )


# ------------------------------------------------------------------ Lecci Step 1
def morlet_kernel(f, fs, n_cycles):
    sigma_t = n_cycles / (2 * np.pi * f)
    t = np.arange(-3.5 * sigma_t, 3.5 * sigma_t + 1.0 / fs, 1.0 / fs)
    w = np.exp(2j * np.pi * f * t) * np.exp(-t ** 2 / (2 * sigma_t ** 2))
    return w / np.sqrt(np.sum(np.abs(w) ** 2))


def morlet_spectrum(x, fs, freqs, n_cycles=N_CYCLES):
    """Mean Morlet power across the time steps of one bout.

    Every accepted bout contributes at every frequency, matching the paper's duration-weighted
    construction. A support-energy correction makes the zero-padded CWT edge explicit; the old
    implementation silently dropped 120-168 s bouts at 0.02 Hz and changed the bout subset at
    every frequency.
    """
    x = np.asarray(x, float)
    x = x - np.nanmean(x)
    if not np.isfinite(x).all():
        return np.full(len(freqs), np.nan)
    out = np.full(len(freqs), np.nan)
    for i, f in enumerate(freqs):
        k = morlet_kernel(f, fs, n_cycles)
        w = signal.fftconvolve(x, k, mode="same")
        support = signal.fftconvolve(np.ones(len(x)), np.abs(k) ** 2, mode="same")
        power = np.abs(w) ** 2 / np.maximum(support, 1e-12)
        out[i] = float(np.mean(power))
    return out


def subject_spectrum(sig, nrem, fs=FS, smooth_4s=True):
    """Duration-weighted mean infraslow spectrum over all NREM bouts >= 120 s (Lecci Fig 1G)."""
    freqs = np.arange(F_LO, F_HI + 1e-12, F_STEP)
    x, _, _ = fill_short_gaps(sig, fs, max_gap_s=5.0)
    if smooth_4s:
        k = int(round(4 * fs))
        if k > 1:
            valid = np.isfinite(x)
            xf = np.where(valid, x, 0.0)
            num = np.convolve(xf, np.ones(k) / k, mode="same")
            den = np.convolve(valid.astype(float), np.ones(k) / k, mode="same")
            x = np.where(den > 0.5, num / np.maximum(den, 1e-12), np.nan)
    # Split each NREM bout at any gap that survived interpolation and treat every artifact-free
    # sub-run >= MIN_BOUT_S as its own bout. Rejecting a whole bout for one long gap would select
    # only the subset of bouts that happen to be completely gap-free.
    spectra, weights = [], []
    for s0, s1 in nrem_bouts(nrem):
        finite = np.isfinite(x[s0:s1])
        for a, b in contiguous_runs(finite):
            seg = x[s0 + a:s0 + b]
            if len(seg) < MIN_BOUT_S * fs:
                continue
            amplitude_scale = max(float(np.max(np.abs(seg))), 1.0)
            if (not np.isfinite(np.std(seg))
                    or np.std(seg) <= 100 * np.finfo(float).eps * amplitude_scale):
                continue
            sp = morlet_spectrum(seg, fs, freqs)
            if np.isfinite(sp).sum() < 5:
                continue
            spectra.append(sp); weights.append(len(seg))
    if not spectra:
        return freqs, None, 0, 0.0
    S = np.array(spectra); W = np.array(weights, float)
    with np.errstate(invalid="ignore"):
        num = np.nansum(S * W[:, None], axis=0)
        den = np.nansum(np.where(np.isfinite(S), 1.0, 0.0) * W[:, None], axis=0)
    mean_spec = np.where(den > 0, num / np.maximum(den, 1e-300), np.nan)
    ok = np.isfinite(mean_spec)
    if ok.sum() < 10:
        return freqs, None, len(spectra), float(W.sum())
    normalizer = float(np.nanmean(mean_spec[ok]))
    if not np.isfinite(normalizer) or normalizer <= np.finfo(float).tiny:
        return freqs, None, len(spectra), float(W.sum())
    mean_spec = mean_spec / normalizer       # Lecci: normalise to its own mean
    if np.isfinite(mean_spec).sum() < 10:
        return freqs, None, len(spectra), float(W.sum())
    return freqs, mean_spec, len(spectra), float(W.sum())


def _scale_free(n, slope, rng):
    """Time series whose power spectrum follows f**slope (Lecci's scale-free control profile)."""
    x = rng.randn(n)
    f = np.fft.rfftfreq(n)
    S = np.fft.rfft(x)
    amp = np.zeros_like(f)
    amp[1:] = f[1:] ** (slope / 2.0)
    y = np.fft.irfft(S * amp, n)
    return y / (y.std() + 1e-12)


def peak_location_null_is_adequate(result):
    """Whether accepted surrogate locations support a nondegenerate clustering null."""
    if not result:
        return False
    peaks = result.get("surrogate_peaks", [])
    n_valid = int(result.get("n_surrogates", 0))
    fraction = len(peaks) / max(n_valid, 1)
    return (len(peaks) >= MIN_SURROGATE_PEAK_LOCATIONS
            and fraction >= MIN_SURROGATE_PEAK_FRACTION)


def peak_significance(sig, nrem, slope, observed_prominence, n_sur=200, seed=0, smooth_4s=True):
    """Lecci fig S3 G-J, per subject: is the observed peak more prominent than one produced by a
    SCALE-FREE profile with the same 1/f slope and the SAME bout-length distribution?

    HONEST LIMITATION: this per-subject control can be weak when only a few hours are available and
    must not be read as proof that a rhythm is absent. ``surrogate_peaks`` is returned so the
    cohort summary can run a subject-matched clustering test without importing hard-coded
    simulation benchmarks. This remains exploratory until the 1 Hz approximation is benchmarked
    against the paper's reference construction.

    The surrogate inherits the real series' NaN pattern, so bout count, bout lengths and gap
    structure are identical by construction; only the spectral content is replaced.
    """
    rng = np.random.RandomState(seed)
    n = len(sig)
    finite = np.isfinite(sig)
    proms, peaks = [], []
    for _ in range(n_sur):
        surr = np.where(finite, _scale_free(n, slope, rng), np.nan)
        fr, sp, _, _ = subject_spectrum(surr, nrem, smooth_4s=smooth_4s)
        if sp is None:
            continue
        pk = fit_peak(fr, sp)
        p = pk.get("prominence_over_background", np.nan)
        if np.isfinite(p):
            proms.append(float(p))
            if pk.get("peak_hz"):
                peaks.append(float(pk["peak_hz"]))
    if len(proms) < 20:
        return None
    proms = np.array(proms)
    p_value = float((1 + np.sum(proms >= observed_prominence)) / (1 + len(proms)))
    result = dict(n_surrogates=len(proms), n_surrogates_requested=int(n_sur),
                observed_prominence=float(observed_prominence),
                null_median=float(np.median(proms)), null_p95=float(np.percentile(proms, 95)),
                p_value=p_value,
                # Full surrogate peak-frequency list for the subject-matched cohort clustering test.
                surrogate_peaks=[float(p) for p in peaks],
                surrogate_peak_hz_median=(float(np.median(peaks)) if peaks else None),
                surrogate_peak_in_infraslow_frac=(
                    float(np.mean((np.array(peaks) >= 0.015) & (np.array(peaks) <= 0.025)))
                    if peaks else None))
    result["n_surrogate_peak_locations"] = int(len(peaks))
    result["surrogate_peak_acceptance_fraction"] = float(len(peaks) / len(proms))
    result["surrogate_peak_locations_adequate"] = peak_location_null_is_adequate(result)
    return result


def _gauss3(f, *p):
    return sum(p[i] * np.exp(-((f - p[i + 1]) / p[i + 2]) ** 2) for i in range(0, 9, 3))


def fit_peak(freqs, spec, search=PEAK_SEARCH):
    """Three-Gaussian peak fit with an objective ``no peak`` outcome.

    Lecci visually confirmed physiological power and reported a three-term Gaussian peak, width,
    and mean normalized power in peak +/- 0.5 spectral SD. An unattended batch analysis also needs
    a genuine-local-excess and fit-quality criterion; otherwise scale-free noise is forced to have
    a Gaussian "peak" inside the allowed range.
    """
    ok = np.isfinite(spec) & (freqs >= F_LO) & (freqs <= F_HI)
    f, y = freqs[ok], spec[ok]
    if len(f) < 30 or np.any(y <= 0):
        return dict(peak_hz=None, accepted=False, method="insufficient")

    # A local excess above a log-log aperiodic background is an automated counterpart to visual
    # confirmation. It gates the Gaussian result but is not substituted for Lecci's endpoint.
    co = np.polyfit(np.log(f), np.log(y), 1)
    resid = y - np.exp(np.polyval(co, np.log(f)))
    band = (f >= search[0]) & (f <= search[1])
    noise_values = resid[~band]
    noise = (1.4826 * np.median(np.abs(noise_values - np.median(noise_values))) + 1e-12
             if len(noise_values) > 5 else np.std(resid) + 1e-12)
    local, props = signal.find_peaks(resid[band], prominence=2.0 * noise)
    if not len(local):
        return dict(peak_hz=None, accepted=False, method="no_local_excess",
                    prominence_over_background=float(np.max(resid[band]) / noise),
                    slope_1_over_f=float(co[0]))
    best_local = int(local[np.argmax(props["prominences"])])
    detrended_peak = float(f[band][best_local])
    residual_prominence = float(props["prominences"][np.argmax(props["prominences"])] / noise)
    if resid[band][best_local] < 2.0 * noise:
        return dict(peak_hz=None, accepted=False, method="weak_local_excess",
                    detrended_peak_hz=detrended_peak,
                    prominence_over_background=residual_prominence,
                    slope_1_over_f=float(co[0]))

    try:
        a0 = float(np.max(y[band])) if band.any() else float(np.max(y))
        p0 = [a0, detrended_peak or 0.02, 0.01, a0 / 2, 0.05, 0.02, a0 / 2, 0.01, 0.05]
        lo = [0, F_LO, 1e-4] * 3
        hi = [10 * a0 + 1e-9, F_HI, 1.0] * 3
        popt, _ = optimize.curve_fit(_gauss3, f, y, p0=p0, bounds=(lo, hi), maxfev=20000)
        fitted = _gauss3(f, *popt)
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        fit_r2 = 1.0 - float(np.sum((y - fitted) ** 2)) / max(ss_tot, 1e-12)
        if not np.isfinite(fit_r2):
            raise ValueError("non-finite Gaussian fit quality")
        components = [
            (float(popt[i]), float(popt[i + 1]), float(popt[i + 2]))
            for i in range(0, 9, 3)
            if search[0] <= popt[i + 1] <= search[1]
        ]
        if not components:
            raise ValueError("no Gaussian component in the target search range")
        amplitude, peak, width = min(
            components, key=lambda value: abs(value[1] - detrended_peak))
        spectral_sd = width / np.sqrt(2.0)
        accepted = (fit_r2 > 0.25 and abs(peak - detrended_peak) <= 0.005
                    and 0.001 <= spectral_sd <= 0.03)
        if not accepted:
            return dict(peak_hz=None, accepted=False, method="gauss3_rejected",
                        detrended_peak_hz=detrended_peak, fit_r2=float(fit_r2),
                        spectral_sd_hz=float(spectral_sd),
                        prominence_over_background=residual_prominence,
                        slope_1_over_f=float(co[0]))
        peak_window = np.abs(f - peak) <= 0.5 * spectral_sd
        return dict(
            peak_hz=float(peak), accepted=True, detrended_peak_hz=detrended_peak,
            method="gauss3_validated", fit_r2=float(fit_r2),
            gaussian_amplitude=amplitude, spectral_sd_hz=float(spectral_sd),
            peak_window_mean=float(np.mean(y[peak_window])),
            prominence_over_background=residual_prominence,
            slope_1_over_f=float(co[0]))
    except Exception as exc:
        return dict(peak_hz=None, accepted=False, method="gauss3_failed",
                    detrended_peak_hz=detrended_peak,
                    prominence_over_background=residual_prominence,
                    slope_1_over_f=float(co[0]), fit_error=type(exc).__name__)


# ------------------------------------------------------------------ Lecci Step 2
def cross_correlation(sig, hr, nrem, fs=FS, win_s=XCORR_WIN_S,
                      max_lag_s=XCORR_MAX_LAG_S, smooth_4s=True,
                      minimum_windows=5):
    """Lecci Fig 6: z-transform each 120 s interval, cross-correlate with HR as source wave,
    average the cross-correlograms within subject. Positive lag => sigma FOLLOWS heart rate."""
    s, _, _ = fill_short_gaps(sig, fs, 5.0)
    h, _, _ = fill_short_gaps(hr, fs, 5.0)
    if smooth_4s:
        k = max(1, int(round(4 * fs)))
        valid = np.isfinite(s)
        num = np.convolve(np.where(valid, s, 0.0), np.ones(k), mode="same")
        den = np.convolve(valid.astype(float), np.ones(k), mode="same")
        s = np.where(den >= k / 2, num / np.maximum(den, 1), np.nan)
    w = int(win_s * fs); ml = int(max_lag_s * fs)
    acc, n = [], 0
    for s0, s1 in nrem_bouts(nrem, min_s=win_s):
        for st in range(s0, s1 - w + 1, w):
            a, b = h[st:st + w], s[st:st + w]
            if not (np.isfinite(a).all() and np.isfinite(b).all()):
                continue
            if a.std() < 1e-9 or b.std() < 1e-9:
                continue
            a = (a - a.mean()) / a.std(); b = (b - b.mean()) / b.std()
            c = signal.correlate(b, a, mode="full") / len(a)
            lags = signal.correlation_lags(len(b), len(a), mode="full")
            keep = np.abs(lags) <= ml
            acc.append(c[keep]); n += 1
    if n < int(minimum_windows):
        return None
    lags_s = signal.correlation_lags(w, w, mode="full")[np.abs(
        signal.correlation_lags(w, w, mode="full")) <= ml] / fs
    m = np.mean(acc, axis=0)
    i = int(np.argmax(np.abs(m)))
    follows = (
        (lags_s >= LECCI_XCORR_LAG_WINDOW_S[0])
        & (lags_s <= LECCI_XCORR_LAG_WINDOW_S[1]))
    i_positive = np.where(follows)[0][int(np.argmax(m[follows]))]
    i_negative = np.where(follows)[0][int(np.argmin(m[follows]))]
    return dict(lag_s=lags_s.tolist(), xcorr=m.tolist(), n_intervals=n,
                peak_r=float(m[i]), peak_lag_s=float(lags_s[i]),
                # Lecci's source-defined human direction is a positive peak at positive lag:
                # heart rate is the source wave and sigma follows.  Preserve the opposite-sign
                # extremum instead of allowing max-|r| to disguise a reversal as replication.
                lecci_direction_peak_r=float(m[i_positive]),
                lecci_direction_peak_lag_s=float(lags_s[i_positive]),
                lecci_direction_lag_window_s=list(LECCI_XCORR_LAG_WINDOW_S),
                opposite_direction_peak_r=float(m[i_negative]),
                opposite_direction_peak_lag_s=float(lags_s[i_negative]))


# ------------------------------------------------------------------ per subject
def _skip_result(subject, lineage, reason, **details):
    return dict(
        subject=subject, status="skip", analysis_version=ANALYSIS_VERSION,
        cache_schema_version=CACHE_SCHEMA_VERSION, **lineage,
        reason=reason, **details)


def _select_input_series(subject, d, band, lineage):
    if band == "fsp" and not bool(np.asarray(d["fsp_is_real_peak"]).item()):
        return None, _skip_result(
            subject, lineage,
            "no independently reliable individual fast-spindle peak")
    roi = "legacy lateral neocortical set"
    if "sigma_fixed_parietal" in d.files:
        n_roi = int(np.asarray(d["parietal_n_selected_contacts"]).item())
        roi_coverage = float(np.asarray(d["parietal_sigma_coverage"]).item())
        if n_roi < MIN_POWER_CONTACTS or roi_coverage < MIN_POWER_COVERAGE:
            return None, _skip_result(
                subject, lineage,
                (
                    f"Lecci-motivated parietal iEEG adaptation unavailable: "
                    f"{n_roi} stable contacts, "
                    f"{roi_coverage:.1%} aggregate coverage"))
        sig = (
            d["sigma_fsp_parietal"] if band == "fsp"
            else d["sigma_fixed_parietal"])
        swa = d["swa_parietal"]
        roi = "non-pathological cortical Destrieux parietal/postcentral/precuneus contacts"
    else:
        n_roi = int(np.asarray(d["sigma_n_selected_contacts"]).item())
        roi_coverage = float(np.asarray(d["sigma_coverage"]).item())
        if n_roi < MIN_POWER_CONTACTS or roi_coverage < MIN_POWER_COVERAGE:
            return None, _skip_result(
                subject, lineage,
                (
                    f"Lecci-motivated global iEEG adaptation unavailable: "
                    f"{n_roi} stable contacts, "
                    f"{roi_coverage:.1%} aggregate coverage"))
        sig = d["sigma_fsp"] if band == "fsp" else d["sigma_fixed"]
        swa = d["swa"]
    return dict(sig=sig, swa=swa, hr=d["hr_1"], roi=roi), None


def _negative_control(freqs, spec, spec_swa, peak):
    if (spec is None or spec_swa is None or not peak.get("peak_hz")
            or not np.isfinite(peak.get("spectral_sd_hz", np.nan))):
        return None
    half_width = 0.5 * float(peak["spectral_sd_hz"])
    window = np.abs(freqs - float(peak["peak_hz"])) <= half_width
    if not window.any():
        return None
    # Both bands must be evaluated in the same sigma-defined window. Requiring an accepted
    # SWA peak preferentially drops the strongest negative controls, while comparing each
    # band's own fitted window tests two different frequencies.
    return dict(
        window_source="accepted sigma peak +/- 0.5 sigma spectral SD",
        centre_hz=float(peak["peak_hz"]), half_width_hz=half_width,
        sigma_window_mean=float(np.mean(spec[window])),
        swa_same_window_mean=float(np.mean(spec_swa[window])))


def _spectral_endpoints(sig, swa, nrem, smooth_4s, n_sur):
    freqs, spec, n_bouts, tot_s = subject_spectrum(
        sig, nrem, smooth_4s=smooth_4s)
    peak = (
        fit_peak(freqs, spec)
        if spec is not None else dict(peak_hz=None, method="no spectrum"))
    _, spec_swa, _, _ = subject_spectrum(
        swa, nrem, smooth_4s=smooth_4s)  # Lecci negative control
    peak_swa = fit_peak(freqs, spec_swa) if spec_swa is not None else dict(peak_hz=None)
    negative_control = _negative_control(freqs, spec, spec_swa, peak)

    # Lecci's scale-free control: is the peak more prominent than matched 1/f noise produces?
    peak_null = None
    if spec is not None and np.isfinite(peak.get("prominence_over_background", np.nan)):
        peak_null = peak_significance(
            sig, nrem, peak.get("slope_1_over_f", 0.0),
            peak["prominence_over_background"],
            n_sur=n_sur, smooth_4s=smooth_4s)
        if peak_null:
            peak["peak_p_value"] = peak_null["p_value"]
    return dict(
        freqs=freqs, spec=spec, spec_swa=spec_swa,
        n_bouts=n_bouts, bout_seconds=tot_s,
        peak=peak, peak_swa=peak_swa,
        negative_control=negative_control, peak_null=peak_null)


def _coherence_endpoint(sig, hr, nrem, peak):
    # Coherence over all NREM, gap-aware, read at the subject's own peak.
    mask = np.zeros(len(sig), bool)
    for s0, s1 in nrem_bouts(nrem):
        mask[s0:s1] = True
    sg = np.where(mask, sig, np.nan)
    hg = np.where(mask, hr, np.nan)
    co = coherence_gapped(hg, sg, fs=FS, nperseg=NPERSEG, highpass=0.005)
    if co is None:
        return None
    crit = analytic_msc_threshold(co["K"], ALPHA)
    f, c = co["f"], co["cxy"]
    own = peak["peak_hz"] if peak.get("peak_hz") else None
    lecci_value = float(c[np.argmin(np.abs(f - F_LECCI))])
    own_value = float(c[np.argmin(np.abs(f - own))]) if own else None
    if not np.isfinite(lecci_value):
        return None
    return dict(
        K=co["K"], crit=float(crit), n_valid=co["n_valid"],
        filled_frac=co["filled_frac"],
        at_lecci=lecci_value,
        sig_at_lecci=bool(lecci_value > crit),
        at_own_peak=(
            own_value if own_value is not None and np.isfinite(own_value) else None),
        sig_at_own_peak=(
            bool(own_value > crit)
            if own_value is not None and np.isfinite(own_value) else None),
        f=f.tolist(),
        cxy=[float(value) if np.isfinite(value) else None for value in c])


def _cardiac_endpoints(sig, hr, nrem, peak, smooth_4s):
    coherence = _coherence_endpoint(sig, hr, nrem, peak)
    xcorr = cross_correlation(sig, hr, nrem, smooth_4s=smooth_4s)
    return coherence, xcorr


def _assemble_result(subject, d, band, smooth_4s, n_sur, lineage, inputs,
                     full_record_nrem, nrem, core_window, separation,
                     spectral, coherence, xcorr):
    unavailable = []
    if spectral["spec"] is None:
        unavailable.append("sigma spectrum")
    if xcorr is None:
        unavailable.append("sigma-HR cross-correlation")
    status = "ok" if not unavailable else "partial"
    reason = None if not unavailable else "unavailable primary endpoint(s): " + ", ".join(unavailable)
    peak = spectral["peak"]
    return dict(
        subject=subject, status=status, reason=reason,
        analysis_version=ANALYSIS_VERSION,
        cache_schema_version=CACHE_SCHEMA_VERSION,
        **lineage,
        anatomy_selection_method=npz_scalar_text(
            d, "anatomy_selection_method", "missing/unvalidated"),
        endpoint_roi=inputs["roi"],
        band=band, smooth_4s=smooth_4s,
        n_surrogates_requested=int(n_sur),
        n_nrem_epochs_full_record=int(full_record_nrem.sum()),
        n_nrem_epochs=int(nrem.sum()),
        core_study_window=core_window,
        n_bouts=spectral["n_bouts"], bout_seconds=spectral["bout_seconds"],
        gmm_separation=finite_float_or_none(separation),
        spectrum=dict(
            f=spectral["freqs"].tolist(),
            sigma=(None if spectral["spec"] is None else np.where(
                np.isfinite(spectral["spec"]), spectral["spec"], None).tolist()),
            swa=(None if spectral["spec_swa"] is None else np.where(
                np.isfinite(spectral["spec_swa"]), spectral["spec_swa"], None).tolist())),
        peak=peak, peak_swa=spectral["peak_swa"],
        negative_control=spectral["negative_control"],
        peak_null=spectral["peak_null"], coherence=coherence, xcorr=xcorr,
        endpoint_availability=dict(
            spectrum=spectral["spec"] is not None,
            cross_correlation=xcorr is not None,
            coherence=coherence is not None,  # compatibility alias for fixed 0.02-Hz endpoint
            coherence_fixed_0p02=coherence is not None,
            coherence_own_peak=bool(
                coherence is not None
                and coherence.get("at_own_peak") is not None
                and peak.get("peak_hz") is not None)))


def analyse(subject, band="fixed", smooth_4s=True, n_sur=200):
    try:
        d = load(subject)
    except CacheSubjectSkipped as exc:
        return _skip_result(
            subject, cache_lineage(subject), f"cache exclusion: {exc}")
    if d is None:
        return None
    lineage = cache_lineage(subject, d)
    inputs, skipped = _select_input_series(subject, d, band, lineage)
    if skipped is not None:
        return skipped

    lab, full_record_nrem, sep = stages_for(d)
    nrem, core_window = core_study_nrem_mask(lab)
    if nrem.sum() < 40:
        return _skip_result(
            subject, lineage, "insufficient NREM in the Lecci core-study window",
            n_nrem_epochs_full_record=int(full_record_nrem.sum()),
            n_nrem_epochs=int(nrem.sum()),
            core_study_window=core_window)

    spectral = _spectral_endpoints(
        inputs["sig"], inputs["swa"], nrem, smooth_4s, n_sur)
    coherence, xcorr = _cardiac_endpoints(
        inputs["sig"], inputs["hr"], nrem, spectral["peak"], smooth_4s)
    return _assemble_result(
        subject, d, band, smooth_4s, n_sur, lineage, inputs,
        full_record_nrem, nrem, core_window, sep, spectral, coherence, xcorr)


def main():
    raise SystemExit(
        "WITHDRAWN standalone 3A writer: run "
        "analysis/run_qc_grid.py to execute the profile-materialized "
        "3A endpoint")


if __name__ == "__main__":
    main()
