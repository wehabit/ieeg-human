"""Validation for spectral_gapped.py.

The decisive test is [3]: SCATTERED dropout, which is what per-second IED masking produces and
which is the regime where delete-and-splice uniformly compresses the time axis. Test [4] checks
that the interpolation introduced by the fix cannot manufacture significance, including under
deliberately CORRELATED gap patterns (both series losing samples at the same instants).

    .venv/bin/python analysis/test_spectral_gapped.py
"""
import numpy as np
from scipy import signal

from spectral_gapped import coherence_gapped, analytic_msc_threshold

rng = np.random.RandomState(0)
FS, N, F0 = 1.0, 3300, 0.02
NPER = 256


def old_msc(hr, sig, nperseg_cap=NPER):
    """The ORIGINAL delete-and-splice implementation, verbatim."""
    m = np.isfinite(hr) & np.isfinite(sig)
    if m.sum() < 600:
        return None
    hr, sig = hr[m], sig[m]
    sos = signal.butter(3, 0.005, btype="high", fs=FS, output="sos")
    a = signal.sosfiltfilt(sos, hr - hr.mean())
    b = signal.sosfiltfilt(sos, sig - sig.mean())
    nper = int(min(nperseg_cap, (len(a) // 4) // 2 * 2))
    if nper < 128:
        return None
    f, cxy = signal.coherence(a, b, fs=FS, nperseg=nper, noverlap=nper // 2)
    return f, cxy, int((len(a) - nper // 2) // (nper // 2))


def brown(n):
    x = rng.randn(n)
    fr = np.fft.rfftfreq(n)
    s = np.fft.rfft(x)
    s[1:] /= fr[1:]
    return np.fft.irfft(s, n)


def coupled(n, amp, f0=F0):
    t = np.arange(n) / FS
    ph = rng.uniform(0, 2 * np.pi)
    a, b = brown(n), brown(n)
    return (a / a.std() + amp * np.sin(2 * np.pi * f0 * t + ph),
            b / b.std() + amp * np.sin(2 * np.pi * f0 * t + ph + 0.4))


def scattered(n, p):
    return rng.rand(n) > p


fails = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not cond:
        fails.append(name)


print("\n[1] With no gaps, coherence_gapped must reproduce scipy.signal.coherence exactly")
a, b = coupled(N, 0.3)
f_s, c_s = signal.coherence(a, b, fs=FS, nperseg=NPER, noverlap=NPER // 2)
r = coherence_gapped(a, b, fs=FS, nperseg=NPER, highpass=None)
check("frequency grid identical", np.allclose(f_s, r["f"]))
check("coherence identical", np.allclose(c_s, r["cxy"], atol=1e-10),
      f"max |diff| = {np.abs(c_s - r['cxy']).max():.2e}")
check("K matches Welch segment count", r["K"] == (N - NPER) // (NPER // 2) + 1, f"K={r['K']}")

print("\n[2] Scattered dropout must remain ESTIMABLE (strict segmenting would return None)")
x, y = coupled(N, 0.3)
m = scattered(N, 0.09)
g = coherence_gapped(np.where(m, x, np.nan), np.where(m, y, np.nan), fs=FS, nperseg=NPER,
                     highpass=0.005)
check("9% scattered dropout still yields an estimate", g is not None,
      f"K={None if g is None else g['K']}, filled={None if g is None else round(g['filled_frac'],3)}")

print("\n[3] THE REGRESSION: scattered dropout shifts the old estimator's peak off 0.02 Hz.")
TRIALS = 300
for p, amp in ((0.31, 0.5), (0.31, 0.3), (0.19, 0.3)):
    op, npk, oh, nh, nn, on = [], [], 0, 0, 0, 0
    for _ in range(TRIALS):
        x, y = coupled(N, amp)
        mm = scattered(N, p)
        xg, yg = np.where(mm, x, np.nan), np.where(mm, y, np.nan)
        o = old_msc(xg, yg)
        if o is not None:
            fo, co, Ko = o
            on += 1
            sel = (fo >= 0.005) & (fo <= 0.06)
            op.append(fo[sel][np.argmax(co[sel])])
            oh += co[np.argmin(np.abs(fo - F0))] > analytic_msc_threshold(Ko)
        gg = coherence_gapped(xg, yg, fs=FS, nperseg=NPER, highpass=0.005)
        if gg is not None:
            nn += 1
            sel = (gg["f"] >= 0.005) & (gg["f"] <= 0.06)
            npk.append(gg["f"][sel][np.argmax(gg["cxy"][sel])])
            nh += gg["cxy"][np.argmin(np.abs(gg["f"] - F0))] > analytic_msc_threshold(gg["K"])
    op, npk = np.array(op), np.array(npk)
    print(f"\n  drop {p:.0%}, coupling {amp}:  predicted old peak = f0/(1-p) = {F0/(1-p):.4f} Hz")
    print(f"    old  median peak {np.median(op):.4f} Hz  (bias {np.median(op)-F0:+.4f})"
          f"   detection@0.02 = {oh/max(on,1):5.1%}")
    print(f"    new  median peak {np.median(npk):.4f} Hz  (bias {np.median(npk)-F0:+.4f})"
          f"   detection@0.02 = {nh/max(nn,1):5.1%}")
    check(f"[p={p:.0%},amp={amp}] new estimator unbiased at 0.02 Hz",
          abs(np.median(npk) - F0) < 0.002)
    check(f"[p={p:.0%},amp={amp}] new recovers detection power lost by splicing",
          nh / max(nn, 1) >= oh / max(on, 1),
          f"{nh/max(nn,1):.1%} vs {oh/max(on,1):.1%}")

print("\n[4] False positives must stay at nominal alpha -- including CORRELATED gaps")
for label, correlated in (("independent gaps", False), ("identical gaps in both series", True)):
    hits = tot = 0
    for _ in range(800):
        x, y = brown(N), brown(N)
        mx = scattered(N, 0.15)
        my = mx if correlated else scattered(N, 0.15)
        gg = coherence_gapped(np.where(mx, x, np.nan), np.where(my, y, np.nan),
                              fs=FS, nperseg=NPER, highpass=0.005)
        if gg is None:
            continue
        tot += 1
        hits += gg["cxy"][np.argmin(np.abs(gg["f"] - F0))] > analytic_msc_threshold(gg["K"])
    fpr = hits / max(tot, 1)
    check(f"FPR near 5% ({label})", 0.02 <= fpr <= 0.10, f"FPR = {fpr:.3f} over {tot} trials")

print("\n[5] Coherence is undefined when either input has no spectral power")
undefined = coherence_gapped(
    np.zeros(1000), rng.randn(1000), fs=FS, nperseg=NPER, highpass=0.005)
check("constant input is rejected instead of reported as zero coherence",
      undefined is None)

print("\n" + ("ALL CHECKS PASSED" if not fails else f"{len(fails)} FAILURE(S): {fails}"))
raise SystemExit(1 if fails else 0)
