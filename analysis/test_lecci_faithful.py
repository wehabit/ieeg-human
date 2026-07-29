"""Validation for lecci_faithful_3A.py -- can it recover a planted rhythm, and does it stay
quiet when there is none?

    .venv/bin/python analysis/test_lecci_faithful.py
"""
import numpy as np

import lecci_faithful_3A as lecci_module
from lecci_faithful_3A import (subject_spectrum, fit_peak, cross_correlation, morlet_spectrum,
                              peak_location_null_is_adequate)
from staging_helpers import EPOCH

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
accepted_noise = 0
for _ in range(30):
    x = brown(TOTAL_S)
    freqs, spec, _, _ = subject_spectrum(x / x.std(), nrem)
    pk = fit_peak(freqs, spec)
    accepted_noise += int(bool(pk.get("accepted")))
    if np.isfinite(pk.get("prominence_over_background", np.nan)):
        proms.append(pk["prominence_over_background"])
proms = np.array(proms)
print(f"    prominence on pure 1/f^2: median {np.median(proms):.2f}, 95th pct {np.percentile(proms,95):.2f}")
print(f"    accepted Gaussian peaks on pure 1/f^2: {accepted_noise}/30")
check("objective peak gate rarely accepts pure scale-free noise", accepted_noise <= 2,
      f"accepted {accepted_noise}/30")

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

common = np.sin(2 * np.pi * 0.02 * np.arange(TOTAL_S))
inverted = cross_correlation(-np.roll(common, 5), common, nrem)
check("opposite-sign delayed coupling cannot count as Lecci-direction support",
      inverted["lecci_direction_peak_r"] <= 0
      and inverted["opposite_direction_peak_r"] < -0.9)

print("\n[5] Every accepted bout must contribute at every frequency")
short = brown(150) / 10
sp = morlet_spectrum(short, FS, np.array([0.005, 0.02, 0.08]))
check("support-corrected CWT returns one finite value per frequency", np.isfinite(sp).all())

print("\n[6] A peak-location clustering null must not resample one accepted null peak forever")
check("one accepted surrogate peak is inadequate for a clustering null",
      not peak_location_null_is_adequate(dict(
          n_surrogates=200, surrogate_peaks=[0.02])))
check("a prespecified nondegenerate accepted-location set is adequate",
      peak_location_null_is_adequate(dict(
          n_surrogates=200, surrogate_peaks=np.linspace(0.01, 0.05, 20).tolist())))

print("\n[7] A zero-variance signal must not count as an available spectrum")
_, constant_spectrum, _, _ = subject_spectrum(
    np.ones(600), np.ones(int(600 / EPOCH), bool))
check("constant power is rejected before zero-spectrum normalization",
      constant_spectrum is None)


print("\n[8] Standalone analysis must apply one serialized 210-minute mask to every endpoint")


class SyntheticCache:
    def __init__(self, values):
        self.values = values
        self.files = list(values)

    def __getitem__(self, key):
        return self.values[key]


leading_rem_labels = np.asarray(["R"] * 20 + ["NREM"] * 500)
normalized_labels, normalized_nrem, _ = lecci_module.stages_for(
    SyntheticCache({"stage_lab": leading_rem_labels}))
leading_rem_nrem, leading_rem_window = (
    lecci_module.core_study_nrem_mask(normalized_labels))
check(
    "REM remains available to anchor the first-210-minute window without entering NREM",
    np.array_equal(normalized_labels, leading_rem_labels)
    and np.array_equal(normalized_nrem, leading_rem_labels == "NREM")
    and leading_rem_window["start_epoch"] == 0
    and leading_rem_window["stop_epoch"] == 420
    and leading_rem_nrem.sum() == 400
    and not leading_rem_nrem[:20].any()
    and not leading_rem_nrem[420:].any(),
)


long_labels = np.asarray(["W"] * 10 + ["NREM"] * 490)
long_seconds = int(len(long_labels) * EPOCH)
fake_cache = SyntheticCache({
    "sigma_n_selected_contacts": np.asarray(3),
    "sigma_coverage": np.asarray(1.0),
    "sigma_fixed": np.ones(long_seconds),
    "swa": np.ones(long_seconds),
    "hr_1": np.ones(long_seconds),
    "anatomy_selection_method": np.asarray("synthetic"),
})
originals = {
    name: getattr(lecci_module, name)
    for name in (
        "load", "cache_lineage", "stages_for", "subject_spectrum",
        "coherence_gapped", "cross_correlation",
    )
}
seen_masks = []
seen_coherence_support = []


def fake_subject_spectrum(signal_values, stage_mask, **kwargs):
    seen_masks.append(np.asarray(stage_mask, bool).copy())
    return np.arange(0.001, 0.121, 0.001), None, 0, 0.0


try:
    lecci_module.load = lambda subject: fake_cache
    lecci_module.cache_lineage = lambda subject, cache=None: {}
    lecci_module.stages_for = lambda cache: (
        long_labels, long_labels == "NREM", None)
    lecci_module.subject_spectrum = fake_subject_spectrum

    def fake_coherence(first, second, **kwargs):
        seen_coherence_support.append(
            int((np.isfinite(first) & np.isfinite(second)).sum()))
        return None

    lecci_module.coherence_gapped = fake_coherence

    def fake_cross_correlation(sig, hr, stage_mask, **kwargs):
        seen_masks.append(np.asarray(stage_mask, bool).copy())
        return None

    lecci_module.cross_correlation = fake_cross_correlation
    standalone = lecci_module.analyse("synthetic", n_sur=2)
finally:
    for name, value in originals.items():
        setattr(lecci_module, name, value)

check(
    "standalone 3A applies the same first-210-minute support to spectrum, control, and coupling",
    len(seen_masks) == 3
    and all(mask.sum() == 420 for mask in seen_masks)
    and all(not mask[430:].any() for mask in seen_masks)
    and seen_coherence_support == [420 * EPOCH]
    and standalone["n_nrem_epochs_full_record"] == 490
    and standalone["n_nrem_epochs"] == 420
    and standalone["core_study_window"]["start_epoch"] == 10
    and standalone["core_study_window"]["stop_epoch"] == 430,
)

print("\n" + ("ALL CHECKS PASSED" if not fails else f"{len(fails)} FAILURE(S): {fails}"))
raise SystemExit(1 if fails else 0)
