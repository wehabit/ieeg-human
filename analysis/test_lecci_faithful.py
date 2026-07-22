"""Validation for lecci_faithful_3A.py -- can it recover a planted rhythm, and does it stay
quiet when there is none?

    .venv/bin/python analysis/test_lecci_faithful.py
"""
import numpy as np

from lecci_faithful_3A import subject_spectrum, fit_peak, cross_correlation, morlet_spectrum
from cohort_stages_3ABD import EPOCH

rng = np.random.RandomState(0)
FS = 1.0
TOTAL_S = 7200
N_EP = int(TOTAL_S // EPOCH)

fails = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not cond:
        fails.append(name)


def pink(n):
    x = rng.randn(n)
    f = np.fft.rfftfreq(n)
    s = np.fft.rfft(x)
    s[1:] /= np.sqrt(f[1:])
    return np.fft.irfft(s, n)


def brown(n):
    x = rng.randn(n)
    f = np.fft.rfftfreq(n)
    s = np.fft.rfft(x)
    s[1:] /= f[1:]
    return np.fft.irfft(s, n)


# NREM = one long consolidated stretch plus a couple of shorter bouts
nrem = np.zeros(N_EP, bool)
nrem[10:110] = True      # 50 min
nrem[130:180] = True     # 25 min
nrem[200:215] = True     # 7.5 min

print("\n[1] Planted rhythms must be recovered at the right frequency")
for f_true in (0.014, 0.020, 0.026, 0.033):
    t = np.arange(TOTAL_S) / FS
    x = brown(TOTAL_S)
    x = x / x.std() + 1.2 * np.sin(2 * np.pi * f_true * t + rng.uniform(0, 6.28))
    freqs, spec, nb, tot = subject_spectrum(x, nrem)
    pk = fit_peak(freqs, spec)
    err = abs(pk["peak_hz"] - f_true) if pk["peak_hz"] else float("nan")
    print(f"    true {f_true:.3f} Hz -> fitted {pk['peak_hz']:.4f} Hz "
          f"({pk['method']}, prominence {pk['prominence_over_background']:.1f}, {nb} bouts)")
    check(f"recovers {f_true} Hz within 0.003 Hz", np.isfinite(err) and err < 0.003,
          f"error {err:.4f} Hz")

print("\n[2] Scale-free noise with NO planted rhythm must not yield a prominent infraslow peak")
proms = []
for _ in range(30):
    x = brown(TOTAL_S)
    freqs, spec, _, _ = subject_spectrum(x / x.std(), nrem)
    pk = fit_peak(freqs, spec)
    if np.isfinite(pk.get("prominence_over_background", np.nan)):
        proms.append(pk["prominence_over_background"])
proms = np.array(proms)
print(f"    prominence on pure 1/f^2: median {np.median(proms):.2f}, 95th pct {np.percentile(proms,95):.2f}")

x = brown(TOTAL_S); t = np.arange(TOTAL_S) / FS
x = x / x.std() + 1.2 * np.sin(2 * np.pi * 0.02 * t)
freqs, spec, _, _ = subject_spectrum(x, nrem)
prom_signal = fit_peak(freqs, spec)["prominence_over_background"]
print(f"    prominence WITH a planted 0.02 Hz rhythm: {prom_signal:.2f}")
check("planted rhythm is more prominent than 95% of noise-only cases",
      prom_signal > np.percentile(proms, 95),
      f"{prom_signal:.2f} vs {np.percentile(proms,95):.2f}")

print("\n[3] Lecci's SWA negative control: a 1/f band must not mimic the sigma peak")
check("noise-only prominence stays modest", np.median(proms) < prom_signal / 2,
      f"median {np.median(proms):.2f} vs signal {prom_signal:.2f}")

print("\n[4] Cross-correlation must recover a known lead/lag")
for true_lag in (-8, 0, 12):
    t = np.arange(TOTAL_S) / FS
    common = np.sin(2 * np.pi * 0.02 * t)
    hr = common + 0.5 * pink(TOTAL_S) / pink(TOTAL_S).std()
    sig = np.roll(common, true_lag) + 0.5 * pink(TOTAL_S) / pink(TOTAL_S).std()
    xc = cross_correlation(sig, hr, nrem)
    print(f"    true lag {true_lag:+3d} s -> recovered {xc['peak_lag_s']:+.0f} s "
          f"(r={xc['peak_r']:+.3f}, {xc['n_intervals']} intervals)")
    check(f"recovers lag {true_lag:+d}s within 2 s", abs(xc["peak_lag_s"] - true_lag) <= 2.0,
          f"got {xc['peak_lag_s']:+.1f}")

print("\n[5] Bouts shorter than the wavelet must return NaN, not an edge artifact")
short = brown(150) / 10
sp = morlet_spectrum(short, FS, np.array([0.005, 0.02, 0.08]))
check("0.005 Hz unmeasurable in a 150 s bout", not np.isfinite(sp[0]))
check("0.08 Hz measurable in a 150 s bout", np.isfinite(sp[2]))

print("\n" + ("ALL CHECKS PASSED" if not fails else f"{len(fails)} FAILURE(S): {fails}"))
raise SystemExit(1 if fails else 0)
