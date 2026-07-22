"""ADVERSARIAL calibration test for coherence_gapped, under the conditions the REAL analysis
creates -- which test_spectral_gapped.py never exercised.

test_spectral_gapped.py validated the estimator on ONE long block with scattered dropout, giving
K ~ 24 and essentially a single contiguous run. The corrected cohort analysis is nothing like that:
it pools Welch segments across ~34 separate NREM bouts and reaches K ~ 107. Three things change,
and each could inflate the false-positive rate -- which would turn the fix into a source of
spurious significance:

  A. HIGH K. crit = 1 - alpha^(1/(K-1)) assumes K INDEPENDENT segments. With 50%-overlapping Hann
     windows they are correlated. That was calibrated empirically at K ~ 24 and simply assumed to
     hold at K ~ 107.
  B. MANY RUNS. The 0.005 Hz drift high-pass is applied PER RUN. A 3rd-order Butterworth high-pass
     at 0.005 Hz has a time constant of ~32 s and settles over ~100-200 s, so every run edge injects
     a filter transient. With one run that is negligible; with 34 runs there are 68 edges. If those
     transients are similar in shape across the two signals they will manufacture coherence at
     exactly the low frequencies this study cares about.
  C. REAL BOUT STRUCTURE. Bout lengths are heavily skewed, so K is dominated by a few long bouts.

Ground truth throughout: the two signals are INDEPENDENT by construction, so every rejection is a
false positive. Anything above ~8% at 0.02 Hz means the corrected pipeline cannot be trusted.

    .venv/bin/python analysis/test_coherence_calibration.py
"""
import os
import numpy as np

from spectral_gapped import coherence_gapped, analytic_msc_threshold
from lecci_faithful_3A import load, nrem_bouts, NPERSEG, FS
from cohort_stages_3ABD import stage_epochs

rng = np.random.RandomState(0)
N_SIM = 400
F0 = 0.02

fails = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not cond:
        fails.append(name)


def brown(n):
    x = rng.randn(n)
    f = np.fft.rfftfreq(n)
    s = np.fft.rfft(x)
    s[1:] /= f[1:]
    return np.fft.irfft(s, n)


def pink(n):
    x = rng.randn(n)
    f = np.fft.rfftfreq(n)
    s = np.fft.rfft(x)
    s[1:] /= np.sqrt(f[1:])
    return np.fft.irfft(s, n)


def fpr_for(mask, gen, highpass, n_sim=N_SIM, label=""):
    """False-positive rate at 0.02 Hz for independent signals restricted to `mask`."""
    n = len(mask)
    hits = tot = 0
    Ks = []
    for _ in range(n_sim):
        x = np.where(mask, gen(n), np.nan)
        y = np.where(mask, gen(n), np.nan)
        r = coherence_gapped(x, y, fs=FS, nperseg=NPERSEG, highpass=highpass)
        if r is None:
            continue
        tot += 1
        Ks.append(r["K"])
        hits += r["cxy"][np.argmin(np.abs(r["f"] - F0))] > analytic_msc_threshold(r["K"])
    return hits / max(tot, 1), (int(np.median(Ks)) if Ks else 0), tot


# ---------------------------------------------------------------- A. high K, single run
print("\n[A] HIGH K on one long contiguous run (isolates K from run structure)")
for total_s, in ((3300,), (12000,), (24000,)):
    mask = np.ones(total_s, bool)
    fpr, K, tot = fpr_for(mask, brown, 0.005, label=f"{total_s}s")
    print(f"    {total_s:6d} s single run -> K={K:4d}  FPR@0.02Hz = {fpr:.3f}  ({tot} trials)")
    check(f"FPR calibrated at K={K}", 0.02 <= fpr <= 0.09, f"{fpr:.3f}")

# ---------------------------------------------------------------- B. many runs, edge transients
print("\n[B] MANY SEPARATE RUNS -- does the per-run high-pass inject correlated edge transients?")
for n_bouts, bout_s in ((1, 20400), (8, 2550), (34, 600), (68, 300)):
    total_s = 24000
    mask = np.zeros(total_s, bool)
    gap = (total_s - n_bouts * bout_s) // max(n_bouts, 1)
    p = 0
    for _ in range(n_bouts):
        mask[p:p + bout_s] = True
        p += bout_s + gap
    for hp, tag in ((0.005, "high-pass ON "), (None, "high-pass OFF")):
        fpr, K, tot = fpr_for(mask, brown, hp, n_sim=250)
        print(f"    {n_bouts:3d} runs x {bout_s:5d}s  {tag} -> K={K:4d}  FPR = {fpr:.3f}  ({tot} trials)")
        if hp is not None:
            check(f"FPR calibrated with {n_bouts} runs, high-pass on", 0.02 <= fpr <= 0.09,
                  f"{fpr:.3f}")

# ---------------------------------------------------------------- C. real bout structure
print("\n[C] REAL bout structure from the cached cohort (the actual analysis geometry)")
CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     "data", "derived", "lc_infraslow")
subjects = sorted(f[:-4] for f in os.listdir(CACHE) if f.endswith(".npz"))[:6] \
    if os.path.isdir(CACHE) else []
tested = 0
for s in subjects:
    d = load(s)
    if d is None:
        continue
    ep = dict(dr=d["ep_dr"], swa=d["ep_swa"], clean=d["ep_clean"])
    lab, nrem, _ = stage_epochs(ep)
    if nrem.sum() < 40:
        continue
    n = len(d["sigma_fsp"])
    mask = np.zeros(n, bool)
    bouts = nrem_bouts(nrem)
    for a, b in bouts:
        mask[a:b] = True
    # also carry that subject's real missing-sample pattern
    real_gaps = np.isfinite(d["sigma_fsp"]) & np.isfinite(d["hr_1"])
    mask &= real_gaps
    if mask.sum() < 3000:
        continue
    for gen, gname in ((brown, "1/f^2"), (pink, "1/f")):
        fpr, K, tot = fpr_for(mask, gen, 0.005, n_sim=200)
        print(f"    {s:18s} {len(bouts):3d} bouts, {int(mask.sum()):6d}s valid, {gname:6s} "
              f"-> K={K:4d}  FPR = {fpr:.3f}")
        check(f"{s} {gname} FPR calibrated", 0.02 <= fpr <= 0.09, f"{fpr:.3f}")
    tested += 1
if not tested:
    print("    (no cache yet -- section skipped)")

print("\n" + ("ALL CHECKS PASSED" if not fails else f"{len(fails)} FAILURE(S): {fails}"))
raise SystemExit(1 if fails else 0)
