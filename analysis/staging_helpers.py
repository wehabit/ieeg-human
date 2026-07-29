"""Shared staging and spectral helpers used by the active LC-proxy analyses."""

import numpy as np
from scipy import signal

from infraslow_rr_sigma_coherence import notch


EPOCH = 30.0
FS_P = 1.0
CHUNK_S = 600.0
SWA_BAND = (0.5, 4.0)
SO_BAND = (0.5, 1.25)
INFRA = (0.01, 0.04)
F_TARGET = 0.02
NREM_DR = 0.90
MIN_3A_MIN = 20.0       # min contiguous minutes for a coherence estimate; short blocks give few
                        # Welch segments (high per-window threshold), but N2/N3 are length-matched
                        # within subject so the ~1/K bias is identical on both sides of the paired test
ALPHA = 0.05
POOLED_3A_CAP_MIN = 55.0  # common cap so K (and thus the ~1/K coherence floor) matches across subjects


def _fsp_candidate(x_mtl, sf):
    """Objective candidate for a fast-spindle peak in one data split."""
    ps = []
    for j in range(x_mtl.shape[1]):
        x = notch(np.nan_to_num(x_mtl[:, j].astype(float)), sf)
        f, p = signal.welch(x, sf, nperseg=int(min(8 * sf, len(x))))
        ps.append(p)
    p = np.nanmean(ps, axis=0)
    positive = (f >= 2) & (f <= 30) & (p > 0)
    # Fit the aperiodic background outside the full spindle range, so a broad real spindle bump
    # cannot pull its own baseline upward.
    background = positive & ~((f >= 8) & (f <= 18))
    if background.sum() < 20:
        return None
    co = np.polyfit(np.log(f[background]), np.log(p[background]), 1)
    resid = np.log(np.maximum(p, np.finfo(float).tiny)) - np.polyval(co, np.log(
        np.maximum(f, f[f > 0].min())))
    noise_values = resid[positive & ~((f >= 8) & (f <= 18))]
    noise = 1.4826 * np.median(np.abs(noise_values - np.median(noise_values))) + 1e-12
    band_idx = np.where((f >= 12.0) & (f <= 15.0))[0]
    peaks, props = signal.find_peaks(resid[band_idx], prominence=2.5 * noise)
    if not len(peaks):
        return None
    local = int(peaks[np.argmax(props["prominences"])])
    idx = int(band_idx[local])
    if resid[idx] < 3.0 * noise:
        return None
    return float(f[idx])


def fsp_from(x_mtl, sf):
    """Reliability-gated individual fast-spindle peak.

    Lecci visually validated FSPs from all artifact-free NREM. Automation needs an explicit
    rejection state: a local maximum in 1/f noise is not automatically a spindle. We require a
    prominent 12-15 Hz bump in the full sample and agreement within 0.75 Hz across temporal halves.
    Failure returns the fixed 13-Hz fallback with ``found_real_peak=False``; production 3A uses the
    fixed 10-15 Hz band unless a validated FSP sensitivity analysis is explicitly requested.
    """
    x_mtl = np.asarray(x_mtl)
    if x_mtl.ndim != 2 or len(x_mtl) < int(120 * sf):
        return 13.0, False
    full = _fsp_candidate(x_mtl, sf)
    mid = len(x_mtl) // 2
    first = _fsp_candidate(x_mtl[:mid], sf)
    second = _fsp_candidate(x_mtl[mid:], sf)
    reliable = (full is not None and first is not None and second is not None
                and abs(first - second) <= 0.75
                and abs(full - first) <= 0.75 and abs(full - second) <= 0.75)
    return (float(full), True) if reliable else (13.0, False)


def band_sos(band, sf, order=4):
    return signal.butter(order, list(band), btype="band", fs=sf, output="sos")


def reliable_two_state_split(values):
    """Stable high-tail Gaussian-mixture partition; not evidence for two physiological states."""
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


def nrem_mask_adaptive(dr, clean, eligibility=None):
    """NREM epochs, thresholded WITHIN subject.

    An absolute delta-ratio cutoff does not transfer: it was tuned on one subject and discarded
    most epochs in others (101 / 48 / 354 NREM epochs where ~800 were expected), collapsing the
    staging. Fit 2 classes to this subject's own delta-ratio distribution instead."""
    valid = np.isfinite(dr) & (np.asarray(clean, float) >= 0.5)
    if eligibility is not None:
        valid &= np.asarray(eligibility, bool)
    if valid.sum() < 60:
        return np.zeros(len(dr), bool)
    high, _ = reliable_two_state_split(np.asarray(dr, float)[valid])
    if high is None:
        return np.zeros(len(dr), bool)
    m = np.zeros(len(dr), bool)
    m[np.where(valid)[0]] = high
    # CONSOLIDATE. Per-epoch classification flickers in and out of NREM every few epochs, which
    # leaves no contiguous block long enough for a 0.02 Hz coherence estimate (628 NREM epochs but
    # no 20-min run). Real sleep scoring smooths over time; a 5-epoch (2.5 min) majority filter
    # fills single-epoch dropouts and removes isolated epochs.
    k = 5
    smoothed = np.convolve(m.astype(float), np.ones(k) / k, mode="same") >= 0.5
    return smoothed & valid


def stage_epochs(ep):
    """NREM epochs -> {N2-like, N3-like} by 2-component GMM on log slow-wave power."""
    dr, swa, clean = ep["dr"], ep["swa"], ep["clean"]
    stage_eligible = np.isfinite(swa) & (swa > 0)
    # Fit the delta partition on exactly the epochs eligible for downstream stage labels.  Allowing
    # SWA-missing epochs into the model can move its components and relabel otherwise unchanged data.
    nrem = nrem_mask_adaptive(dr, clean, eligibility=stage_eligible)
    lab = np.full(len(dr), "", dtype=object)
    if nrem.sum() < 40:
        return lab, nrem, None
    lab[nrem] = "NREM"
    high, diagnostic = reliable_two_state_split(np.log(swa[nrem]))
    sep = diagnostic.get("separation")
    if high is not None:
        lab[np.where(nrem)[0]] = np.where(high, "N3", "N2")
    return lab, nrem, sep


def dominant_block(sel, nrem, win=11, frac_thr=0.60):
    """Longest STAGE-DOMINANT block, not stage-pure.

    N2 and N3 alternate on a timescale of minutes, so no 40+ min block is ever purely one stage --
    demanding purity returns nothing. Coherence needs continuity, so instead smooth the stage label
    over ~5 min and keep periods where the stage is the local majority of NREM. Returns
    (start_epoch, end_epoch, minutes, purity) with purity = actual fraction of the block in-stage."""
    s = sel.astype(float); nz = nrem.astype(float)
    k = np.ones(win)
    num = np.convolve(s, k, mode="same"); den = np.convolve(nz, k, mode="same")
    frac = np.divide(num, np.maximum(den, 1e-9))
    dom = (frac >= frac_thr) & nrem
    best = (0, 0)
    i = 0
    while i < len(dom):
        if not dom[i]:
            i += 1; continue
        j = i
        while j < len(dom) and dom[j]:
            j += 1
        if j - i > best[1] - best[0]:
            best = (i, j)
        i = j
    i0, i1 = best
    if i1 <= i0:
        return None
    seg_nrem = nrem[i0:i1]
    purity = float(sel[i0:i1][seg_nrem].mean()) if seg_nrem.any() else 0.0
    return i0, i1, (i1 - i0) * EPOCH / 60.0, purity


def msc_block(hr, sig, nperseg_cap=256):
    m = np.isfinite(hr) & np.isfinite(sig)
    if m.sum() < 600:
        return None
    hr, sig = hr[m], sig[m]
    sos = signal.butter(3, 0.005, btype="high", fs=FS_P, output="sos")
    a, b = signal.sosfiltfilt(sos, hr - hr.mean()), signal.sosfiltfilt(sos, sig - sig.mean())
    nper = int(min(nperseg_cap, (len(a) // 4) // 2 * 2))
    if nper < 128:
        return None
    f, cxy = signal.coherence(a, b, fs=FS_P, nperseg=nper, noverlap=nper // 2)
    K = int((len(a) - nper // 2) // (nper // 2))
    band = (f >= INFRA[0]) & (f <= INFRA[1])
    return dict(at=float(cxy[np.argmin(np.abs(f - F_TARGET))]), bmax=float(cxy[band].max()),
                K=K, crit=1 - ALPHA ** (1 / max(K - 1, 1)), n_sec=int(len(a)), nperseg=nper)


__all__ = [
    "EPOCH",
    "FS_P",
    "CHUNK_S",
    "SWA_BAND",
    "SO_BAND",
    "INFRA",
    "F_TARGET",
    "NREM_DR",
    "MIN_3A_MIN",
    "ALPHA",
    "POOLED_3A_CAP_MIN",
    "fsp_from",
    "band_sos",
    "reliable_two_state_split",
    "nrem_mask_adaptive",
    "stage_epochs",
    "dominant_block",
    "msc_block",
]
