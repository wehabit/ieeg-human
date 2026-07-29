"""Reproducible demonstrations for the 2026-07 LC-proxy audit.

These are counterexamples to legacy estimators, not tests of the biological data.
"""
import json
import os

import numpy as np
from scipy import interpolate, signal, stats

from pipeline_version import ANALYSIS_VERSION, CACHE_SCHEMA_VERSION
from cache_lc_series import aggregate_staging_features, staging_epoch_features
from cohort_3A_cortical import delta_ratio
from staging_helpers import reliable_two_state_split, stage_epochs, SWA_BAND
from event_3d_estimators import (
    pair_one_spindle_per_so,
    participant_rotation_test,
    pooled_endpoint_passes_qc,
    spindle_rms,
    spindle_events_from_rms,
)
from lecci_faithful_3A import subject_spectrum
from spectral_gapped import coherence_gapped


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rng = np.random.default_rng(20260723)


def heading(text):
    print(f"\n{text}\n{'-' * len(text)}")


heading("E1. Averaging channel maxima inflates a participant-level 3B peak")
n_sim, n_ch, n_time = 20000, 12, 21
noise = rng.normal(size=(n_sim, n_ch, n_time))
legacy = noise.max(axis=2).mean(axis=1)       # maximise each channel, then average
correct = noise.mean(axis=1).max(axis=1)      # average curves, then maximise once
print(f"null mean, legacy estimator   : {legacy.mean():.3f}")
print(f"null mean, participant curve  : {correct.mean():.3f}")
print(f"legacy/correct inflation ratio: {legacy.mean() / correct.mean():.2f}x")
assert legacy.mean() > 2 * correct.mean()


heading("E2. Unconstrained cubic interpolation can overshoot observed HR values")
x = np.arange(4.0)
y = np.array([60.0, 180.0, 60.0, 180.0])
grid = np.linspace(0, 3, 1001)
cubic = interpolate.CubicSpline(x, y)(grid)
pchip = interpolate.PchipInterpolator(x, y)(grid)
print(f"observed range : [{y.min():.1f}, {y.max():.1f}] bpm")
print(f"CubicSpline    : [{cubic.min():.1f}, {cubic.max():.1f}] bpm")
print(f"PCHIP          : [{pchip.min():.1f}, {pchip.max():.1f}] bpm")
assert cubic.min() < y.min() and cubic.max() > y.max()
assert pchip.min() >= y.min() and pchip.max() <= y.max()


heading("E3. Chunk-wise percentiles force weak chunks to contribute false 'top-quartile' events")
weak = np.ones(100)
strong = np.full(100, 10.0)
legacy_keep = np.r_[weak >= np.percentile(weak, 75),
                    strong >= np.percentile(strong, 75)]
whole = np.r_[weak, strong]
correct_keep = whole >= np.percentile(whole, 75)
print(f"weak events retained, chunk-wise : {legacy_keep[:100].sum()}/100")
print(f"weak events retained, whole-night: {correct_keep[:100].sum()}/100")
assert legacy_keep[:100].sum() == 100 and correct_keep[:100].sum() == 0


heading("E4. Raw resultant length R is event-count biased even under uniform phases")
def uniform_r(n, sims=50000):
    phase = rng.uniform(-np.pi, np.pi, size=(sims, n))
    return np.abs(np.exp(1j * phase).mean(axis=1))

r20, r200 = uniform_r(20), uniform_r(200)
print(f"median null R with n=20 events : {np.median(r20):.3f}")
print(f"median null R with n=200 events: {np.median(r200):.3f}")
print(f"ratio                         : {np.median(r20) / np.median(r200):.2f}x")
assert np.median(r20) > 2.5 * np.median(r200)


heading("E5. Uncorrected channel tests produce subject-level false discoveries")
for channels in (6, 54, 93):
    chance_any = 1 - 0.95 ** channels
    print(f"{channels:2d} independent null channels: P(at least one p<.05) = {chance_any:.1%}")


heading(f"E6. Current QC outputs carry {ANALYSIS_VERSION}; legacy results are quarantined")
current_outputs = [
    os.path.join(ROOT, "outputs", "qc_calibration", "staging_window_calibration.json"),
]
for grid_id in (
        "coverage_oat_v1", "auxiliary_window_support_v1",
        "event_count_oat_v1", "staging_window_support_v1"):
    for cohort in ("hup", "respect"):
        current_outputs.append(os.path.join(
            ROOT, "outputs", "qc_grid_public", grid_id,
            f"{cohort}_qc_grid_summary.json"))

for path in current_outputs:
    with open(path, encoding="utf-8") as handle:
        record = json.load(handle)
    relative = os.path.relpath(path, ROOT)
    print(f"{relative:72s}: {record.get('analysis_version')}")
    assert record.get("analysis_version") == ANALYSIS_VERSION
    assert record.get("cache_schema_version") == CACHE_SCHEMA_VERSION

for directory in (
        "lecci_faithful_3A", "event_3B_cached", "event_3D_by_stage",
        "ds003848_3A", "ds003848_3B"):
    marker = os.path.join(ROOT, "outputs", directory, "LEGACY_DO_NOT_USE.md")
    with open(marker, encoding="utf-8") as handle:
        warning = handle.read().lower()
    print(f"outputs/{directory:24s}: legacy marker present")
    assert "legacy" in warning and (
        "do not use" in warning or "do not cite" in warning
    )


heading("E7. Independently randomized channel nulls are anti-conservative")
n_sim, n_events, n_channels, n_lags = 100000, 40, 12, 21
# Null observed channels share one event train, hence one noise curve. The legacy surrogate gives
# each channel an independent event train and then averages, shrinking its variance by n_channels.
shared_stat = rng.normal(scale=1 / np.sqrt(n_events), size=(n_sim, n_lags)).max(axis=1)
independent_stat = rng.normal(
    scale=1 / np.sqrt(n_events * n_channels), size=(n_sim, n_lags)).max(axis=1)
legacy_cutoff = np.percentile(independent_stat, 95)
false_positive = np.mean(shared_stat > legacy_cutoff)
print(f"nominal alpha                         : 5.0%")
print(f"legacy independent-channel null FPR  : {false_positive:.1%}")
assert false_positive > 0.50


heading("E8. Repeating correlated spindle peaks within each SO inflates Rayleigh significance")
def legacy_rayleigh_p(phases):
    n = phases.shape[-1]
    r = np.abs(np.exp(1j * phases).mean(axis=-1))
    z = n * r ** 2
    return np.exp(-z) * (1 + (2 * z - z ** 2) / (4 * n))


n_sim, n_so = 30000, 30
base_phase = rng.uniform(-np.pi, np.pi, size=(n_sim, n_so))
for repeats in (1, 2, 3, 4):
    repeated = np.repeat(base_phase, repeats, axis=1)
    p = legacy_rayleigh_p(repeated)
    fpr = np.mean(p < 0.05)
    print(f"{repeats} peak(s) per SO: false-positive rate {fpr:.1%}")
    if repeats == 1:
        assert 0.04 <= fpr <= 0.06
    else:
        assert fpr > 0.15


heading("E9. Raw-voltage averaging can cancel identical contact power")
t = np.arange(0, 30, 0.001)
a = np.sin(2 * np.pi * 13 * t)
b = -a
mean_contact_power = np.mean([np.mean(a ** 2), np.mean(b ** 2)])
power_after_voltage_mean = np.mean(((a + b) / 2) ** 2)
print(f"mean of per-contact powers : {mean_contact_power:.3f}")
print(f"power after voltage average: {power_after_voltage_mean:.3f}")
assert mean_contact_power > 0.49 and power_after_voltage_mean == 0


heading("E10. Per-chunk normalization erases genuine slower amplitude changes")
raw_power = np.r_[np.ones(600), np.full(600, 4.0)]
legacy = np.r_[raw_power[:600] / np.median(raw_power[:600]),
               raw_power[600:] / np.median(raw_power[600:])]
correct = raw_power / np.median(raw_power)
print(f"legacy chunk medians: {np.median(legacy[:600]):.1f}, {np.median(legacy[600:]):.1f}")
print(f"full-night medians : {np.median(correct[:600]):.1f}, {np.median(correct[600:]):.1f}")
assert np.allclose(legacy, 1) and not np.allclose(correct[:600], correct[600:])


heading("E11. A percentile/duration rule cannot by itself establish spindle identity")
sf = 200.0
stationary_noise = np.random.RandomState(1701).normal(size=int(3600 * sf))
sos = signal.butter(3, [12.0, 16.0], btype="band", fs=sf, output="sos")
band_noise = signal.sosfiltfilt(sos, stationary_noise)
rms = spindle_rms(band_noise, sf)
threshold = np.percentile(rms, 75)
noise_events = spindle_events_from_rms(rms, sf, threshold)
print(f"events called in 60 min stationary Gaussian noise: {len(noise_events)}")
print(f"apparent event rate                         : {len(noise_events) / 60:.2f}/min")
assert len(noise_events) > 100


heading("E12. Raw post-trough maxima are biased by event count")
local_rng = np.random.default_rng(1)
n_participant, n_lag = 30, 21
few = local_rng.normal(
    scale=1 / np.sqrt(30), size=(n_participant, n_lag)).max(axis=1)
many = local_rng.normal(
    scale=1 / np.sqrt(300), size=(n_participant, n_lag)).max(axis=1)
raw_p = stats.wilcoxon(few, many).pvalue
few_null = local_rng.normal(
    scale=1 / np.sqrt(30), size=(100000, n_lag)).max(axis=1).mean()
many_null = local_rng.normal(
    scale=1 / np.sqrt(300), size=(100000, n_lag)).max(axis=1).mean()
centered_p = stats.wilcoxon(few - few_null, many - many_null).pvalue
print(f"zero-effect mean maximum, 30 events : {few.mean():.3f}")
print(f"zero-effect mean maximum, 300 events: {many.mean():.3f}")
print(f"paired raw-max p                     : {raw_p:.3g}")
print(f"paired own-null-centered p           : {centered_p:.3f}")
assert raw_p < 0.001 and centered_p > 0.05


heading("E13. The iid R correction fails under repeated/serially dependent phases")
local_rng = np.random.default_rng(1)
blocked, independent = [], []
for _ in range(25):
    phases = (
        np.repeat(local_rng.uniform(-np.pi, np.pi, 30), 10),
        local_rng.uniform(-np.pi, np.pi, 300),
    )
    for values, destination in zip(phases, (blocked, independent)):
        n = len(values)
        R = abs(np.exp(1j * values).mean())
        destination.append(np.sqrt(max((n * R ** 2 - 1) / (n - 1), 0)))
stage_p = stats.wilcoxon(blocked, independent).pvalue
print(f"median corrected R, 30 phase blocks repeated 10x: {np.median(blocked):.3f}")
print(f"median corrected R, 300 iid phases             : {np.median(independent):.3f}")
print(f"paired stage p under equal uniform marginals    : {stage_p:.3g}")
assert stage_p < 0.001


heading("E14. Testing sigma against SWA after selecting sigma's maximum is circular")
local_rng = np.random.default_rng(44)
p_values = []
for _ in range(1000):
    sigma_null = local_rng.normal(size=(20, 50))
    swa_null = local_rng.normal(size=(20, 50))
    selected_bin = np.argmax(sigma_null, axis=1)
    row = np.arange(20)
    p_values.append(stats.wilcoxon(
        sigma_null[row, selected_bin],
        swa_null[row, selected_bin],
        alternative="greater").pvalue)
p_values = np.asarray(p_values)
print(f"nominal alpha                                      : 5.0%")
print(f"false-positive rate after sigma-bin selection      : {np.mean(p_values < .05):.1%}")
print(f"median one-sided p under identical null populations: {np.median(p_values):.3g}")
assert np.mean(p_values < 0.05) > 0.95

heading("E15. Undefined constant-signal spectra/coherence must not count as null measurements")
_, constant_spec, _, _ = subject_spectrum(
    np.ones(600), np.ones(20, bool))
constant_coherence = coherence_gapped(
    np.zeros(1000), np.random.RandomState(15).normal(size=1000),
    fs=1.0, nperseg=256, highpass=0.005)
print(f"constant power spectrum endpoint available : {constant_spec is not None}")
print(f"zero-autospectrum coherence available      : {constant_coherence is not None}")
assert constant_spec is None and constant_coherence is None


heading("E16. The former minimal 3D participant had inadequate inferential support")
former = dict(
    n_channels_tested=1,
    n_spindle_events=20,
    per_channel=[dict(
        n=20, valid_pooled_nrem_seconds=120,
        valid_pooled_nrem_fraction=1.0)],
)
supported = dict(
    n_channels_tested=3,
    n_spindle_events=200,
    per_channel=[
        dict(n=n, valid_pooled_nrem_seconds=1200,
             valid_pooled_nrem_fraction=0.80)
        for n in (67, 67, 66)
    ],
)
print(f"former 1-contact/20-event/120-s endpoint passes: "
      f"{pooled_endpoint_passes_qc(former)}")
print(f"prespecified supported boundary passes          : "
      f"{pooled_endpoint_passes_qc(supported)}")
assert not pooled_endpoint_passes_qc(former)
assert pooled_endpoint_passes_qc(supported)


heading("E17. A finite SO-centered pairing window creates phase locking under independence")
fs = 20
n = 1800 * fs
t = np.arange(n) / fs
phase = 2 * np.pi * 0.8 * t
so = np.arange(2 * fs, n - 2 * fs, 5 * fs)
participant_vectors, participant_events, participant_qc = [], [], []
for subject_i in range(25):
    contact_vectors, contact_records = [], []
    for contact_i in range(3):
        local_rng = np.random.RandomState(1000 + 3 * subject_i + contact_i)
        spindles = np.where(local_rng.rand(n) < 1 / fs)[0]
        selected, _ = pair_one_spindle_per_so(
            so, spindles, local_rng.rand(len(spindles)), 2 * fs)
        contact_vectors.append(np.mean(np.exp(1j * phase[selected])))
        contact_records.append(dict(
            n=len(selected),
            valid_pooled_nrem_seconds=1800.0,
            valid_pooled_nrem_fraction=1.0,
        ))
    participant_vectors.append(np.mean(contact_vectors))
    participant_events.append(sum(value["n"] for value in contact_records))
    participant_qc.append(pooled_endpoint_passes_qc(dict(
        per_channel=contact_records,
        n_channels_tested=3,
        n_spindle_events=participant_events[-1],
    )))
selection_group = participant_rotation_test(
    participant_vectors, n_sur=19999, rng=np.random.RandomState(0))
analytic_window_bias = abs(np.sin(2 * np.pi * 0.8 * 2) / (2 * np.pi * 0.8 * 2))
print(f"participants passing every pooled QC gate: {sum(participant_qc)}/25")
print(f"median independent paired events/subject: {np.median(participant_events):.0f}")
print(f"analytic finite-window phase bias        : {analytic_window_bias:.4f}")
print(f"simulated descriptive group R            : {selection_group['group_R']:.4f}")
print(f"invalid participant-rotation p           : {selection_group['p']:.5f}")
assert all(participant_qc)
assert len(participant_vectors) == len(participant_qc) == 25
assert selection_group["p"] < 0.001


heading("E18. A brief artifact can dominate an unmasked staging spectrum")
sf = 200.0
t = np.arange(int(30 * sf)) / sf
baseline = (
    0.1 * np.sin(2 * np.pi * 1 * t)
    + 0.1 * np.sin(2 * np.pi * 10 * t)
)
artifact = baseline.copy()
contaminated = (t >= 10) & (t < 12)
artifact[contaminated] += 100 * np.sin(2 * np.pi * 1 * t[contaminated])

def swa_welch(x):
    frequency, power = signal.welch(x, sf, nperseg=int(4 * sf))
    keep = (frequency >= SWA_BAND[0]) & (frequency < SWA_BAND[1])
    return float(np.trapezoid(power[keep], frequency[keep]))

baseline_swa = swa_welch(baseline)
unmasked_swa = swa_welch(artifact)
clean = np.ones(len(t), bool)
clean[contaminated] = False
safe_dr, safe_swa, safe_details = staging_epoch_features(
    artifact, clean, np.ones(len(t), bool), sf, return_details=True)
print(f"nominal clean fraction                    : {clean.mean():.1%}")
print(f"unmasked/baseline SWA inflation           : {unmasked_swa / baseline_swa:,.0f}x")
print(f"old unmasked delta ratio                  : {delta_ratio(artifact, sf):.6f}")
print(f"gap-aware complete windows retained       : {safe_details['n_valid_windows']}/14")
print(f"gap-aware/baseline SWA ratio              : {safe_swa / baseline_swa:.6f}")
assert unmasked_swa > 1000 * baseline_swa
assert safe_details["n_valid_windows"] == 12
assert np.isclose(safe_swa, baseline_swa, rtol=1e-6)
assert np.isfinite(safe_dr)


heading("E19. Changing contact gain/availability can manufacture stage-proxy modes")
legacy_dr = np.r_[np.full(150, 0.10), np.full(300, 0.60)]
legacy_swa = np.r_[np.full(150, 0.20), np.full(150, 0.50), np.full(150, 5.0)]
legacy_labels, legacy_nrem, legacy_sep = stage_epochs(dict(
    dr=legacy_dr, swa=legacy_swa, clean=np.ones(450)))

dr_contacts = np.full((6, 300), 0.60)
swa_contacts = np.r_[np.full((3, 300), 0.50), np.full((3, 300), 5.0)]
dr_contacts[:3, 150:] = np.nan
swa_contacts[:3, 150:] = np.nan
dr_contacts[3:, :150] = np.nan
swa_contacts[3:, :150] = np.nan
fixed_dr, fixed_swa, _, fixed_qc = aggregate_staging_features(
    dr_contacts, swa_contacts, np.ones((6, 300)), np.ones(6, bool))
print(f"legacy stationary NREM labels             : "
      f"N2={np.sum(legacy_labels == 'N2')}, N3={np.sum(legacy_labels == 'N3')}")
print(f"legacy mixture separation                 : {legacy_sep:.2f}")
print(f"corrected fixed-set requirement/support   : "
      f"{fixed_qc['required_contact_count']}/6 required, 3/6 present")
print(f"corrected manufactured epochs retained    : {np.isfinite(fixed_swa).sum()}/300")
assert legacy_nrem.sum() == 300
assert np.sum(legacy_labels == "N2") == np.sum(legacy_labels == "N3") == 150
assert np.isnan(fixed_dr).all() and np.isnan(fixed_swa).all()


heading("E20. Two-Gaussian fit gates do not prove two latent sleep states")
skewed_one_mode = np.random.RandomState(123).lognormal(0.0, 0.8, 1000)
high_tail, diagnostic = reliable_two_state_split(skewed_one_mode)
print(f"known source distribution                 : one-mode lognormal")
print(f"two-vs-one Gaussian BIC gain              : {diagnostic['bic_gain_2_vs_1']:.1f}")
print(f"standardized component separation         : {diagnostic['separation']:.2f}")
print(f"accepted high-tail fraction               : {diagnostic['high_component_fraction']:.1%}")
assert high_tail is not None
assert diagnostic["bic_gain_2_vs_1"] > 10

print("\nALL AUDIT COUNTEREXAMPLES CONFIRMED")
