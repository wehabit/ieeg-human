"""Gap-aware spectral estimation for the LC-infraslow tests.

WHY THIS EXISTS
---------------
The original `msc_block()` handled missing samples like this:

    m = np.isfinite(hr) & np.isfinite(sig)
    hr, sig = hr[m], sig[m]

i.e. it DELETED bad seconds and glued the survivors together, then treated the result as an
unbroken 1 Hz recording. Measured drop rates in the HUP cohort: median 8.7%, up to 30.8%.

The damage depends entirely on HOW the samples go missing, which is why this module distinguishes
two regimes instead of applying one rule:

  * SCATTERED dropout (isolated bad seconds -- what per-second IED masking actually produces).
    Deleting every ~11th sample uniformly COMPRESSES the time axis, so a true f0 is reported at
    f0/(1-p). At the cohort's worst drop rate (30.8%) a 0.02 Hz rhythm lands at 0.029 Hz, more than
    two frequency bins from where the test looks. Simulation: detection of a real weak coupling
    falls from 100% to 4.8%.
  * BURST dropout (contiguous chunks lost). Splicing then removes whole blocks while preserving
    timing WITHIN each surviving block, so the peak frequency stays put and only coherence is
    attenuated. Simulation: the old code is essentially unbiased here.

So the fix cannot be "never interpolate" (that discards the scattered case entirely -- with 9%
scattered dropout NO 256 s window is gap-free, and the estimator returns nothing) nor "always
interpolate" (that would paper over genuine multi-minute discontinuities).

THE FIX
-------
1. Gaps up to `max_gap_s` are linearly INTERPOLATED. At 1 Hz, filling isolated 1-5 s holes is
   negligible for a 50 s oscillation and it preserves the time base exactly -- which is the whole
   problem with deletion.
2. Longer gaps are treated as true discontinuities: Welch segments are drawn only from within
   contiguous runs and never straddle one.
3. Any segment whose interpolated fraction exceeds `max_fill_frac` is rejected, so no estimate
   rests mostly on invented samples.
4. Auto/cross-spectra are accumulated over the surviving segment pool before forming coherence:
       Cxy(f) = |sum_k Pxy_k|^2 / (sum_k Pxx_k * sum_k Pyy_k)
   which is exactly Welch's estimator with a restricted segment pool. With no gaps it reproduces
   `scipy.signal.coherence` to ~1e-15.
5. The 0.005 Hz drift high-pass is applied PER RUN, never across a discontinuity.

Interpolating both series could in principle induce shared smoothness where the two signals lose
samples at the same instants. test_spectral_gapped.py therefore checks the false-positive rate
under deliberately CORRELATED gap patterns; it stays at nominal alpha.
"""
import numpy as np
from scipy import signal


def contiguous_runs(valid):
    """[(start, stop), ...] for each maximal run of True in a boolean array."""
    v = np.asarray(valid).astype(np.int8)
    if v.size == 0:
        return []
    d = np.diff(np.concatenate(([0], v, [0])))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def fill_short_gaps(x, fs, max_gap_s):
    """Linearly interpolate gaps of at most `max_gap_s`. Returns (filled, still_missing, was_filled).

    Leading/trailing NaN are never extrapolated -- they stay missing.
    """
    x = np.asarray(x, float).copy()
    finite = np.isfinite(x)
    filled = np.zeros(len(x), bool)
    if not finite.any():
        return x, ~finite, filled
    max_gap = int(round(max_gap_s * fs))
    idx = np.arange(len(x))
    first, last = np.flatnonzero(finite)[[0, -1]]
    for s, e in contiguous_runs(~finite):
        if s <= first or e > last:          # edge gap: nothing to interpolate between
            continue
        if e - s <= max_gap:
            x[s:e] = np.interp(idx[s:e], idx[finite], x[finite])
            filled[s:e] = True
    return x, ~np.isfinite(x), filled


def welch_segments(valid, nperseg, noverlap):
    """Start indices of every Welch segment lying entirely inside one gap-free run."""
    step = nperseg - noverlap
    starts = []
    for s, e in contiguous_runs(valid):
        if e - s < nperseg:
            continue
        starts.extend(range(s, e - nperseg + 1, step))
    return starts


def _prepare(x, y, fs, max_gap_s, highpass):
    """Fill short gaps in both series, intersect validity, high-pass per run."""
    xf, x_missing, x_filled = fill_short_gaps(x, fs, max_gap_s)
    yf, y_missing, y_filled = fill_short_gaps(y, fs, max_gap_s)
    valid = ~x_missing & ~y_missing
    filled = (x_filled | y_filled) & valid
    xf = np.where(valid, xf, np.nan)
    yf = np.where(valid, yf, np.nan)
    if highpass is not None:
        sos = signal.butter(3, highpass, btype="high", fs=fs, output="sos")
        padlen = 3 * (2 * len(sos) + 1)
        for s, e in contiguous_runs(valid):
            if e - s <= padlen + 1:
                valid[s:e] = False
                continue
            xf[s:e] = signal.sosfiltfilt(sos, xf[s:e] - xf[s:e].mean())
            yf[s:e] = signal.sosfiltfilt(sos, yf[s:e] - yf[s:e].mean())
    return xf, yf, valid, filled


def coherence_gapped(x, y, fs, nperseg, noverlap=None, window="hann", highpass=None,
                     detrend="constant", max_gap_s=5.0, max_fill_frac=0.25):
    """Magnitude-squared coherence that preserves the time base across missing samples.

    Returns dict(f, cxy, K, n_valid, n_runs, filled_frac) or None if fewer than 3 segments survive.
    K counts the Welch segments actually averaged -- feed it to analytic_msc_threshold.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.shape != y.shape:
        raise ValueError("x and y must be the same length")
    if noverlap is None:
        noverlap = nperseg // 2

    xf, yf, valid, filled = _prepare(x, y, fs, max_gap_s, highpass)
    starts = [s for s in welch_segments(valid, nperseg, noverlap)
              if filled[s:s + nperseg].mean() <= max_fill_frac]
    if len(starts) < 3:
        return None

    win = signal.get_window(window, nperseg)
    scale = 1.0 / (fs * (win ** 2).sum())
    Pxx = Pyy = 0.0
    Pxy = 0.0 + 0.0j
    for st in starts:
        xs, ys = xf[st:st + nperseg], yf[st:st + nperseg]
        if detrend == "constant":
            xs, ys = xs - xs.mean(), ys - ys.mean()
        elif detrend == "linear":
            xs, ys = signal.detrend(xs), signal.detrend(ys)
        X = np.fft.rfft(xs * win)
        Y = np.fft.rfft(ys * win)
        Pxx = Pxx + (X * np.conj(X)).real * scale
        Pyy = Pyy + (Y * np.conj(Y)).real * scale
        Pxy = Pxy + (X * np.conj(Y)) * scale

    freqs = np.fft.rfftfreq(nperseg, 1.0 / fs)
    non_dc = freqs > 0
    x_support = float(np.sum(Pxx[non_dc]))
    y_support = float(np.sum(Pyy[non_dc]))
    if (not np.isfinite(x_support) or not np.isfinite(y_support)
            or x_support <= np.finfo(float).tiny
            or y_support <= np.finfo(float).tiny):
        return None
    denom = Pxx * Pyy
    cxy = np.full(denom.shape, np.nan, float)
    supported = np.isfinite(denom) & (denom > 0)
    cxy[supported] = np.abs(Pxy[supported]) ** 2 / denom[supported]
    return dict(f=freqs, cxy=cxy,
                K=len(starts), n_valid=int(valid.sum()), n_runs=len(contiguous_runs(valid)),
                filled_frac=float(filled.sum() / max(valid.sum(), 1)))


def psd_gapped(x, fs, nperseg, noverlap=None, window="hann", detrend="constant",
               max_gap_s=5.0, max_fill_frac=0.25, highpass=None):
    """Welch PSD using the same gap policy as coherence_gapped."""
    xf, _, valid, filled = _prepare(x, x, fs, max_gap_s, highpass)
    if noverlap is None:
        noverlap = nperseg // 2
    starts = [s for s in welch_segments(valid, nperseg, noverlap)
              if filled[s:s + nperseg].mean() <= max_fill_frac]
    if len(starts) < 3:
        return None
    win = signal.get_window(window, nperseg)
    scale = 1.0 / (fs * (win ** 2).sum())
    Pxx = 0.0
    for st in starts:
        xs = xf[st:st + nperseg]
        if detrend == "constant":
            xs = xs - xs.mean()
        elif detrend == "linear":
            xs = signal.detrend(xs)
        X = np.fft.rfft(xs * win)
        Pxx = Pxx + (X * np.conj(X)).real * scale
    return dict(f=np.fft.rfftfreq(nperseg, 1.0 / fs), pxx=Pxx / len(starts), K=len(starts))


def analytic_msc_threshold(K, alpha=0.05):
    """Halliday 1995 confidence limit for MSC with K independent segments."""
    return 1.0 - alpha ** (1.0 / max(K - 1, 1))
