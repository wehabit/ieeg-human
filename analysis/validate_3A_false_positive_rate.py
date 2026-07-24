"""Quarantined calibration of the superseded pooled single-bin 3A test.

The cohort found 7/23 subjects exceeding their own alpha=0.05 analytic coherence threshold, where
1.2 would be expected by chance (binomial p=0.0001). That inference is only valid if the analytic
threshold  crit = 1 - alpha^(1/(K-1))  actually delivers a 5% false-positive rate for THIS pipeline.
It assumes K independent Welch segments; with 50% overlapping Hann windows the segments are
correlated, so the true rate could be higher and the "finding" could be pure bias.

Test: push pairs of INDEPENDENT signals (no coupling by construction) through the identical
msc_block() code path and count how often they clear the threshold. Physiological signals are
strongly autocorrelated, so several spectral shapes are tried -- autocorrelation is exactly what
inflates coherence false positives.
"""
raise SystemExit(
    "LEGACY 3A CALIBRATION QUARANTINED: run test_coherence_calibration.py for the current "
    "gap-aware estimator.")

import numpy as np
from scipy import signal

from cohort_stages_3ABD import msc_block, FS_P, POOLED_3A_CAP_MIN

rng = np.random.RandomState(0)
N_SIM = 2000
n_sec = int(POOLED_3A_CAP_MIN * 60)          # same length as the capped pooled-3A block


def white(n):
    return rng.randn(n)


def pink(n):
    x = rng.randn(n)
    f = np.fft.rfftfreq(n)
    s = np.fft.rfft(x)
    s[1:] /= np.sqrt(f[1:])
    return np.fft.irfft(s, n)


def brown(n):
    """1/f^2 -- strongly autocorrelated, closest to slow physiological drift."""
    x = rng.randn(n)
    f = np.fft.rfftfreq(n)
    s = np.fft.rfft(x)
    s[1:] /= f[1:]
    return np.fft.irfft(s, n)


def ar1(n, rho=0.99):
    x = np.zeros(n); e = rng.randn(n)
    for i in range(1, n):
        x[i] = rho * x[i - 1] + e[i]
    return x


print(f"pooled-3A block = {POOLED_3A_CAP_MIN:.0f} min at {FS_P:.0f} Hz = {n_sec} samples")
print(f"{N_SIM} independent pairs per spectral shape, through the real msc_block()\n")
print(f"{'spectrum':10s} {'K':>4s} {'crit':>7s} {'FPR@0.02Hz':>11s} {'FPR band-max':>13s}  verdict")

for name, gen in (("white", white), ("pink 1/f", pink), ("brown 1/f^2", brown), ("AR1 r=.99", ar1)):
    hits_at = hits_bmax = 0
    K = crit = None
    for _ in range(N_SIM):
        r = msc_block(gen(n_sec), gen(n_sec))
        if r is None:
            continue
        K, crit = r["K"], r["crit"]
        hits_at += r["at"] > r["crit"]
        hits_bmax += r["bmax"] > r["crit"]
    fpr = hits_at / N_SIM
    fprb = hits_bmax / N_SIM
    verdict = "OK (~5%)" if 0.02 <= fpr <= 0.09 else ("INFLATED" if fpr > 0.09 else "conservative")
    print(f"{name:10s} {K:4d} {crit:7.3f} {fpr:11.3f} {fprb:13.3f}  {verdict}")

print("\nInterpretation: the cohort observed 7/23 = 0.304 at 0.02 Hz.")
print("If FPR ~0.05 the excess is real (binomial p=1e-4). If FPR ~0.30 it is entirely expected.")
