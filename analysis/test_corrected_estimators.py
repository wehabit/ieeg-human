"""Synthetic regression tests for participant-level 3B and event-locked 3D helpers."""
import io

import numpy as np
from scipy import signal

from event_3B_mednick import subject_so_triggered, rr_baseline_hr, FS_RR
from event_3B_cached import coverage_eligible_contacts
from event_3D_by_stage import (detect_so_events, detect_spindle_events, bh_fdr, rayleigh,
                               pair_one_spindle_per_so, so_event_candidates, select_so_events,
                               spindle_rms, spindle_events_from_rms,
                               participant_rotation_test, same_stage_pair_mask,
                               pooled_endpoint_passes_qc,
                               PRODUCTION_3D_INFERENCE_ENABLED)
from cache_lc_series import (detect_so_candidates, threshold_so_candidates, sanitize_beats,
                             _write_multichannel_power, _aggregate_full_night_power,
                             aggregate_staging_features, interpolate_tachograms,
                             staging_epoch_features)
from cohort_stages_3ABD import (
    band_sos, fsp_from, nrem_mask_adaptive,
    reliable_two_state_split as mixture_high_tail_split, stage_epochs,
)
from results_3A_tutorial_style import ied_clean_mask
from stage_ds003848 import (
    score_stages, channel_roles_from_rows, selected_channel_name_mismatches,
)
from spectral_gapped import fill_short_gaps


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


rng = np.random.RandomState(7)

# 3B: two channels have true peaks at different lags. The participant estimator must report the
# peak of the averaged curves, not the average of independently maximised channel values.
total_s = 1800
grid = np.arange(0, total_s, 1 / FS_RR)
hr = 70.0 + 0.15 * rng.randn(len(grid))
troughs_a = np.arange(30, total_s - 30, 15.0)
troughs_b = troughs_a + 4.0
for t in troughs_a:
    hr += 0.35 * np.exp(-0.5 * ((grid - (t + 1.0)) / 0.5) ** 2)
for t in troughs_b:
    hr += 0.35 * np.exp(-0.5 * ((grid - (t + 3.0)) / 0.5) ** 2)
pool = np.arange(len(hr))
res = subject_so_triggered(
    hr, [troughs_a, troughs_b], float(hr.mean()), pool, n_sur=199,
    rng=np.random.RandomState(8))
check("3B returns exactly one participant-level estimate", res is not None and res["n_channels"] == 2)
check("3B former stage-shift p is retained only as a bounded diagnostic",
      res["p_upper"] is None
      and 0 < res["stage_shift_p_upper_diagnostic"] <= 1
      and res["inference_status"].startswith("disabled"))
check("3B peak belongs to the averaged curve", 0 <= res["peak_lag_s"] < 5)

# A contact can have >=30 stored candidates yet lose them after complete finite in-stage RR windows
# are enforced. The production two-contact requirement must apply after that final filter.
one_valid_channel = np.arange(30, total_s - 30, 30.0)
boundary_only_channel = np.full(40, 1.0)
post_window_qc = subject_so_triggered(
    hr, [one_valid_channel, boundary_only_channel], float(hr.mean()), pool,
    n_sur=19, rng=np.random.RandomState(81), minimum_channels=2,
    channel_ids=["eligible", "boundary-only"])
check("3B minimum-contact gate is applied after complete-window eligibility",
      post_window_qc is None)


class _SyntheticCacheContacts:
    files = ["sigma_selected_contact_mask"]

    def __getitem__(self, key):
        if key != "sigma_selected_contact_mask":
            raise KeyError(key)
        return np.array([True, True, True, False])


eligible_contacts, excluded_contacts = coverage_eligible_contacts(
    _SyntheticCacheContacts(), ["A", "B", "C", "low-coverage"])
check("3B excludes contacts rejected by the cache's stable coverage gate",
      eligible_contacts == ["A", "B", "C"]
      and excluded_contacts == ["low-coverage"])

# A perfectly shared event train must have the same null as one channel. Independently randomizing
# every channel would make the multichannel null much narrower and is known to inflate FPR to 78%.
single = subject_so_triggered(
    hr, [troughs_a], float(hr.mean()), pool, n_sur=199, rng=np.random.RandomState(18))
shared = subject_so_triggered(
    hr, [troughs_a] * 6, float(hr.mean()), pool, n_sur=199, rng=np.random.RandomState(18))
check("3B shared-shift null preserves perfectly correlated channels",
      np.isclose(single["null_mean_pct"], shared["null_mean_pct"])
      and np.isclose(single["null_sd_pct"], shared["null_sd_pct"])
      and single["stage_shift_p_upper_diagnostic"]
      == shared["stage_shift_p_upper_diagnostic"])

# Paper-aligned RR tachogram: HR bursts are RR minima, and electrode peak times are averaged.
rr = np.full_like(grid, 60.0 / 70.0)
for t0 in troughs_a:
    rr -= 0.025 * np.exp(-0.5 * ((grid - (t0 + 2.0)) / 0.35) ** 2)
rr_result = subject_so_triggered(
    rr, [troughs_a, troughs_a], 70.0, pool, n_sur=99,
    rng=np.random.RandomState(19), domain="rr")
check("3B RR-domain estimator finds RR-minimum/HR-maximum timing",
      rr_result["tachogram_domain"] == "rr"
      and abs(rr_result["peak_lag_s"] - 2.0) <= 0.25)

# A smooth HR bump with SOs merely clustered near it has no event-locked transient. Whole-stage
# shifts nevertheless call the high local level significant; this diagnostic must not be exposed
# as inference, and the local event curve should reveal that it is flat.
trend_n = int(3600 * FS_RR)
trend_t = np.arange(trend_n) / FS_RR
trend_hr = 65 + 10 * np.exp(-0.5 * ((trend_t - 3000) / 400) ** 2)
trend_rr = 60 / trend_hr
trend_troughs = np.arange(2950, 3050, 2.5)
trend_result = subject_so_triggered(
    trend_rr, [trend_troughs], rr_baseline_hr(trend_rr), np.arange(trend_n),
    n_sur=499, rng=np.random.RandomState(0), domain="rr")
check("3B disables a nonstationarity-sensitive stage-shift false positive",
      trend_result["stage_shift_p_upper_diagnostic"] < 0.01
      and trend_result["p_upper"] is None
      and abs(trend_result["event_locked_local_change_pct"]) < 0.05)

# The RR-domain curve must use 60/mean(RR) as its baseline.  The arithmetic mean of instantaneous
# HR is a different estimand (Jensen's inequality) and creates a nonzero "effect" with no event.
heterogeneous_rr = np.tile([0.5, 1.5], 100)
check("3B RR baseline uses the same average-then-invert denominator as its event curve",
      np.isclose(rr_baseline_hr(heterogeneous_rr), 60.0)
      and not np.isclose(rr_baseline_hr(heterogeneous_rr),
                         np.mean(60.0 / heterogeneous_rr)))

# A >5-s beat-knot discontinuity must remain longer than the generic <=5-s filling allowance.
# Otherwise a 5.1-s acquisition gap can be reopened after 4-Hz discretisation.
gap_beats = np.r_[np.arange(0.0, 31.0), np.arange(35.1, 70.1)]
_, rr_gap, _, _ = interpolate_tachograms(gap_beats, 70)
rr_refilled, _, _ = fill_short_gaps(rr_gap, FS_RR, max_gap_s=5.0)
gap_grid = np.arange(0, 70, 1 / FS_RR)
long_gap = (gap_grid >= 30.0) & (gap_grid <= 36.1)
check("3B >5-s ECG gap cannot be reopened by downstream short-gap filling",
      np.isnan(rr_gap[long_gap]).all() and np.isnan(rr_refilled[long_gap]).all())

# 3D: construct complete 1-Hz SO cycles and 13-Hz spindle bursts of valid duration.
sf = 200.0
t = np.arange(0, 120, 1 / sf)
so = np.sin(2 * np.pi * t)
sp = np.zeros_like(t)
for centre in np.arange(5, 115, 5):
    win = np.abs(t - centre) <= 0.4
    sp[win] += 3.0 * np.sin(2 * np.pi * 13 * t[win]) * np.hanning(win.sum())
so_events = detect_so_events(so, sf)
sp_events = detect_spindle_events(sp + 0.02 * rng.randn(len(sp)), sf)
check("3D detects complete SO events", len(so_events) > 20)
check("3D spindle detector enforces event duration and finds bursts", len(sp_events) >= 15)

selected_sp, selected_so = pair_one_spindle_per_so(
    [100, 300], [90, 105, 110, 295, 305], [1, 5, 3, 2, 4], max_distance_samples=30)
check("3D retains one maximum-amplitude spindle per independent SO",
      selected_sp.tolist() == [105, 305] and selected_so.tolist() == [100, 300])
check("3D stage analyses reject SO-spindle pairs that cross a stage boundary",
      not same_stage_pair_mask([1], [0], {1})[0]
      and same_stage_pair_mask([1], [1], {1})[0])

check("Rayleigh probability is bounded for perfect locking",
      0 < rayleigh(np.zeros(20))[1] <= 1)
uniform_p = np.array([
    rayleigh(rng.uniform(-np.pi, np.pi, 30))[1] for _ in range(3000)
])
check("Rayleigh test is calibrated on independent uniform phases",
      0.035 <= np.mean(uniform_p < 0.05) <= 0.065)

# One event per SO does not make nearby SO phases independent. Repeating 30 independent null
# phase blocks five times makes the ordinary Rayleigh p spuriously tiny. Cohort inference must use
# one vector per participant rather than treating those repeated phases as extra observations.
block_phase = np.random.RandomState(0).uniform(-np.pi, np.pi, 30)
repeated_phase = np.repeat(block_phase, 5)
check("3D ordinary event-level Rayleigh can be anti-conservative under serial dependence",
      rayleigh(repeated_phase)[1] < 0.05)

# The mathematical rotation helper is calibrated only when participant directions are truly
# uniform. Production use is disabled because the finite pairing window violates that assumption.
participant_null = 0.2 * np.exp(
    1j * np.random.RandomState(11).uniform(-np.pi, np.pi, 25))
null_group = participant_rotation_test(
    participant_null, n_sur=4999, rng=np.random.RandomState(12))
participant_locked = 0.2 * np.exp(1j * np.full(25, 0.4))
locked_group = participant_rotation_test(
    participant_locked, n_sur=4999, rng=np.random.RandomState(13))
check("rotation helper does not reject explicitly scattered participant directions",
      null_group["p"] > 0.05)
check("rotation helper detects explicitly common directions",
      locked_group["group_R"] > 0.19 and locked_group["p"] < 0.01)
group_p = []
for seed in range(300):
    local = np.random.RandomState(1000 + seed)
    null_vectors = local.uniform(0.02, 0.25, 25) * np.exp(
        1j * local.uniform(-np.pi, np.pi, 25))
    group_p.append(participant_rotation_test(
        null_vectors, n_sur=499, rng=np.random.RandomState(5000 + seed))["p"])
check("rotation helper is calibrated for explicitly uniform-direction cohorts",
      0.025 <= np.mean(np.asarray(group_p) < 0.05) <= 0.075)
check("3D production inference remains disabled pending a pairing-aware null",
      not PRODUCTION_3D_INFERENCE_ENABLED)

check("3D pooled descriptive endpoint rejects the former one-contact/20-event record",
      not pooled_endpoint_passes_qc(dict(
          n_channels_tested=1, n_channels_passing_coverage=1,
          n_spindle_events=20,
          per_channel=[dict(
              n=20, valid_pooled_nrem_seconds=120,
              valid_pooled_nrem_fraction=1.0)])))
check("3D pooled descriptive endpoint accepts only the prespecified supported record",
      pooled_endpoint_passes_qc(dict(
          n_channels_tested=3, n_channels_passing_coverage=3,
          n_spindle_events=200,
          per_channel=[
              dict(n=n, valid_pooled_nrem_seconds=1200,
                   valid_pooled_nrem_fraction=0.80)
              for n in (67, 67, 66)
          ])))

decision = bh_fdr([0.001, 0.01, 0.04, 0.5])
check("BH-FDR rejects only supported channel tests", decision.tolist() == [True, True, False, False])

# 3B channel-night threshold: a weak chunk and a strong chunk must be ranked together, not each
# forced to contribute its own top quartile.
weak = [(float(i), 1.0, 1.0, 2.0) for i in range(100)]
strong = [(float(i + 100), 10.0, 10.0, 20.0) for i in range(100)]
selected = threshold_so_candidates(weak + strong)
check("3B percentile is applied once across the channel-night",
      len(selected) == 100 and np.all(selected >= 100))
check("R-peak refractory filtering compares with the last retained peak",
      np.allclose(sanitize_beats([0.0, 0.2, 0.4, 0.8]), [0.0, 0.4, 0.8]))

# Per-contact power must be invariant to flipping one contact's polarity. Raw-voltage averaging
# would cancel the second construction.
sf_power = 200.0
t_power = np.arange(0, 180, 1 / sf_power)
base = np.sin(2 * np.pi * 13 * t_power)
sos = band_sos((12, 14), sf_power)
same = np.column_stack([base, base])
opposite = np.column_stack([base, -base])
dest_same = np.full((2, 180), np.nan)
dest_opposite = np.full((2, 180), np.nan)
_write_multichannel_power(same, sos, sf_power, 180, 180, 0, dest_same)
_write_multichannel_power(opposite, sos, sf_power, 180, 180, 0, dest_opposite)
check("per-contact power is invariant to contact polarity",
      np.allclose(_aggregate_full_night_power(dest_same),
                  _aggregate_full_night_power(dest_opposite), atol=1e-8, equal_nan=True))

# Missing raw samples may be filled only as a numerical scaffold for filtering. They must remain
# missing in derived power, rather than becoming a plausible trough that inflates coverage.
with_gap = base.copy()
with_gap[int(60 * sf_power):int(70 * sf_power)] = np.nan
dest_gap = np.full((1, 180), np.nan)
_write_multichannel_power(
    with_gap[:, None], sos, sf_power, 180, 180, 0, dest_gap)
check("raw NaN gaps remain ineligible in derived sigma power",
      np.isnan(dest_gap[0, 60:70]).all())

# A brief high-amplitude artifact can dominate a 30-s Welch spectrum even if most samples remain
# clean. Contact-level staging therefore fails closed rather than using the unmasked epoch.
stage_sf = 200.0
stage_t = np.arange(int(30 * stage_sf)) / stage_sf
clean_stage_signal = np.sin(2 * np.pi * stage_t) + 0.05 * rng.randn(len(stage_t))
clean_stage_mask = np.ones(len(stage_t), bool)
measured_stage_mask = np.ones(len(stage_t), bool)
clean_dr, clean_swa = staging_epoch_features(
    clean_stage_signal, clean_stage_mask, measured_stage_mask, stage_sf)
artifact_stage_signal = clean_stage_signal.copy()
artifact_region = (stage_t >= 10) & (stage_t < 12)
artifact_stage_signal[artifact_region] += 100 * np.sin(
    2 * np.pi * stage_t[artifact_region])
artifact_clean_mask = clean_stage_mask.copy()
artifact_clean_mask[artifact_region] = False
artifact_dr, artifact_swa = staging_epoch_features(
    artifact_stage_signal, artifact_clean_mask, measured_stage_mask, stage_sf)
check("artifact-tainted contact epochs cannot contribute delta/SWA staging power",
      np.isfinite(clean_dr) and np.isfinite(clean_swa)
      and np.isnan(artifact_dr) and np.isnan(artifact_swa))

# Contact gain/availability must not manufacture a two-mode epoch series. In this counterexample,
# three low-gain contacts are available in one half and three high-gain contacts in the other.
# The fixed-set 80% support rule makes both halves unavailable rather than inventing N2/N3.
changing_dr = np.full((6, 300), 0.6)
changing_swa = np.r_[
    np.full((3, 300), 0.5),
    np.full((3, 300), 5.0),
]
changing_clean = np.ones((6, 300))
changing_dr[:3, 150:] = np.nan
changing_swa[:3, 150:] = np.nan
changing_dr[3:, :150] = np.nan
changing_swa[3:, :150] = np.nan
fixed_dr, fixed_swa, _, fixed_qc = aggregate_staging_features(
    changing_dr, changing_swa, changing_clean, np.ones(6, bool))
check("changing contact gain/availability cannot manufacture stage-proxy modes",
      fixed_qc["required_contact_count"] == 5
      and np.isnan(fixed_dr).all() and np.isnan(fixed_swa).all())

# Values and missingness on a contact excluded by full-night QC must have no effect.
base_dr = np.tile(np.r_[np.full(50, 0.1), np.full(50, 0.6)], (4, 1))
base_swa = np.tile(np.r_[np.ones(50), np.full(50, 4.0)], (4, 1))
base_clean = np.ones_like(base_dr)
selected_mask = np.array([True, True, True, False])
aggregate_a = aggregate_staging_features(
    base_dr, base_swa, base_clean, selected_mask)[:2]
base_dr[-1] = 100
base_swa[-1] = np.nan
aggregate_b = aggregate_staging_features(
    base_dr, base_swa, base_clean, selected_mask)[:2]
check("full-night-excluded contacts cannot alter staging features",
      all(np.allclose(a, b, equal_nan=True) for a, b in zip(aggregate_a, aggregate_b)))

disjoint = np.full((6, 600), np.nan)
for channel in range(6):
    disjoint[channel, channel * 100:(channel + 1) * 100] = 1.0
qualified, contact_qc = _aggregate_full_night_power(
    disjoint, min_contact_coverage=0.80, min_contacts=3,
    min_contact_fraction_per_bin=0.80, return_details=True)
check("disjoint low-coverage contacts cannot masquerade as 100% aggregate coverage",
      np.isnan(qualified).all() and contact_qc["n_selected"] == 0)

# Filter/Hilbert context must make a streamed result effectively invariant to an internal 600-s
# boundary. Without overlap, the synthetic boundary discrepancy is about 7%.
duration_s = 1200
t_long = np.arange(0, duration_s, 1 / sf_power)
modulated = (1.0 + 0.4 * np.sin(2 * np.pi * 0.013 * t_long)) * np.sin(
    2 * np.pi * 13 * t_long)
x_long = np.column_stack([modulated, 0.8 * modulated])
dest_full = np.full((2, duration_s), np.nan)
_write_multichannel_power(
    x_long, sos, sf_power, duration_s, duration_s, 0, dest_full)
dest_stream = np.full((2, duration_s), np.nan)
edge = int(30 * sf_power)
_write_multichannel_power(
    x_long[:int(630 * sf_power)], sos, sf_power, 600, duration_s, 0, dest_stream)
_write_multichannel_power(
    x_long[int(570 * sf_power):], sos, sf_power, 600, duration_s, 600, dest_stream,
    core_start_sample=edge)
full_power = _aggregate_full_night_power(dest_full)
stream_power = _aggregate_full_night_power(dest_stream)
check("30-s filter context removes the internal chunk-boundary artifact",
      np.nanmax(np.abs(full_power - stream_power)) < 1e-3)

# Artifact-only coupling counterexample: padded IED masks must reject filter ringing.
x = 0.001 * rng.randn(int(240 * sf))
x[(np.arange(10, 230, 10) * sf).astype(int)] += 1000
clean = ied_clean_mask(x, sf, pad_s=2.5)
so_art = signal.sosfiltfilt(band_sos((0.16, 1.25), sf, 3), x)
sp_art = signal.sosfiltfilt(band_sos((12, 16), sf), x)
candidates = so_event_candidates(so_art, sf)
if len(candidates):
    candidates = candidates[clean[candidates[:, 0].astype(int)]]
so_art_events = select_so_events(candidates)
rms_art = spindle_rms(sp_art, sf)
if clean.any():
    threshold = np.percentile(rms_art[clean], 75)
    sp_art_events, sp_art_amp = spindle_events_from_rms(
        rms_art, sf, threshold, valid=clean, return_amplitudes=True)
else:
    sp_art_events, sp_art_amp = np.array([], int), np.array([])
paired_art, _ = pair_one_spindle_per_so(
    so_art_events, sp_art_events, sp_art_amp, int(2.5 * sf))
check("artifact/IED padding rejects impulse-only false coupling", len(paired_art) == 0)

# EMG/EOG can exclude obvious wake/REM but cannot make every remaining epoch NREM. A flat,
# near-zero delta ratio must yield no NREM labels; a separable high-delta component may yield
# conservative NREM candidates.
n_stage = 100
flat_lab, flat_counts = score_stages(
    np.ones(n_stage), np.ones(n_stage), np.ones(n_stage),
    np.zeros(n_stage), np.ones(n_stage))
check("RESPect stage proxy does not default low-delta epochs to NREM",
      flat_counts["n_nrem"] == 0 and not np.isin(flat_lab, ["N2", "N3"]).any())
unimodal_rng = np.random.RandomState(31)
unimodal_lab, unimodal_counts = score_stages(
    np.exp(unimodal_rng.normal(0, 0.05, n_stage)),
    np.exp(unimodal_rng.normal(0, 0.05, n_stage)),
    np.exp(unimodal_rng.normal(0, 0.05, n_stage)),
    unimodal_rng.normal(0.30, 0.015, n_stage), np.ones(n_stage))
check("RESPect narrow near-Gaussian delta example is left unclassified",
      unimodal_counts["n_nrem"] == 0
      and not np.isin(unimodal_lab, ["N2", "N3"]).any())
separated_lab, separated_counts = score_stages(
    np.r_[np.ones(50), np.full(50, 4.0)], np.ones(n_stage), np.ones(n_stage),
    np.r_[np.full(50, 0.10), np.full(50, 0.60)], np.ones(n_stage))
check("RESPect stage proxy uses a separable high-delta component",
      separated_counts["n_nrem"] >= 40
      and np.isin(separated_lab, ["NREM", "N2", "N3"]).any())
stage_buffer = io.BytesIO()
np.savez_compressed(stage_buffer, stage_lab=separated_lab)
stage_buffer.seek(0)
with np.load(stage_buffer, allow_pickle=False) as staged:
    reloaded_stage_labels = staged["stage_lab"]
check("RESPect stage labels remain readable with pickle disabled",
      reloaded_stage_labels.dtype.kind == "U"
      and np.array_equal(reloaded_stage_labels, separated_lab))

roles = channel_roles_from_rows([
    {"name": "ecg_bad", "type": "ECG", "status": "bad"},
    {"name": "ecg_good", "type": "ECG", "status": "good"},
    {"name": "seeg", "type": "SEEG", "status": "good"},
    {"name": "emg", "type": "EMG", "status": "good"},
    {"name": "eog", "type": "EOG", "status": "good"},
])
check("RESPect channel parser excludes a bad ECG before selecting the usable ECG",
      roles["ecg"] == [1])
check("RESPect order QC ignores an MNE-renamed unused duplicate placeholder",
      selected_channel_name_mismatches(
          ["ecg", "seeg", ".....-1"], ["ecg", "seeg", "....."], [0, 1]) == [])
check("RESPect order QC still rejects a mismatch in a selected modality",
      selected_channel_name_mismatches(
          ["ecg", "wrong", ".....-1"], ["ecg", "seeg", "....."], [0, 1])
      == [(1, "wrong", "seeg")])

# Excluded epochs cannot set the robust centres/MADs used to classify retained epochs.
clean_rng = np.random.RandomState(77)
n_clean, n_dirty = 100, 20
clean_flag = np.r_[np.ones(n_clean), np.zeros(n_dirty)]
clean_swa = np.r_[np.ones(50), np.full(50, 4.0)]
clean_emg = np.exp(clean_rng.normal(0, 0.4, n_clean))
clean_eog = np.exp(clean_rng.normal(0, 0.4, n_clean))
dr_all = np.r_[np.full(50, 0.10), np.full(50, 0.60), np.full(n_dirty, 0.10)]
labels_missing_dirty, _ = score_stages(
    np.r_[clean_swa, np.full(n_dirty, np.nan)],
    np.r_[clean_emg, np.full(n_dirty, np.nan)],
    np.r_[clean_eog, np.full(n_dirty, np.nan)],
    dr_all, clean_flag)
labels_extreme_dirty, _ = score_stages(
    np.r_[clean_swa, np.full(n_dirty, 1e6)],
    np.r_[clean_emg, np.full(n_dirty, 1e-6)],
    np.r_[clean_eog, np.full(n_dirty, 1e6)],
    dr_all, clean_flag)
check("dirty RESPect epochs cannot alter robust thresholds or labels of clean epochs",
      np.array_equal(
          labels_missing_dirty[:n_clean], labels_extreme_dirty[:n_clean]))

# Epochs missing a required modality can never receive a RESPect label and must not move the
# delta-model reference distribution used for eligible epochs.
reference_swa = np.r_[np.ones(50), np.full(50, 4.0)]
reference_dr = np.r_[np.full(50, 0.05), np.full(50, 0.50)]
reference_emg = np.ones(100)
reference_eog = np.ones(100)
reference_clean = np.ones(100)
reference_labels, _ = score_stages(
    reference_swa, reference_emg, reference_eog, reference_dr, reference_clean)
extended_labels, _ = score_stages(
    np.r_[reference_swa, np.ones(100)],
    np.r_[reference_emg, np.full(100, np.nan)],
    np.r_[reference_eog, np.full(100, np.nan)],
    np.r_[reference_dr, np.full(100, 5.0)],
    np.ones(200))
check("RESPect modality-ineligible epochs cannot relabel eligible epochs",
      np.array_equal(reference_labels, extended_labels[:100]))

# HUP staging must not let temporal majority smoothing re-admit an explicitly dirty epoch, and a
# unimodal SWA distribution must remain pooled NREM rather than being forced into N2/N3.
hup_dr = np.r_[np.full(60, 0.15), np.full(60, 0.55)]
hup_clean = np.ones(120)
hup_clean[80] = 0.0
hup_nrem = nrem_mask_adaptive(hup_dr, hup_clean)
check("HUP NREM smoothing never re-admits a dirty epoch", not hup_nrem[80])
hup_swa_rng = np.random.RandomState(41)
hup_lab, hup_nrem_all, _ = stage_epochs(dict(
    dr=hup_dr, swa=np.exp(hup_swa_rng.normal(0, 0.05, 120)),
    clean=np.ones(120)))
check("HUP staging does not force unimodal SWA into N2/N3",
      hup_nrem_all.sum() >= 40 and (hup_lab == "NREM").sum() >= 40
      and not np.isin(hup_lab, ["N2", "N3"]).any())

hup_reference = dict(
    dr=np.r_[np.full(50, 0.05), np.full(50, 0.50)],
    swa=np.r_[np.ones(50), np.full(50, 4.0)],
    clean=np.ones(100),
)
hup_reference_lab, hup_reference_nrem, _ = stage_epochs(hup_reference)
hup_extended_lab, hup_extended_nrem, _ = stage_epochs(dict(
    dr=np.r_[hup_reference["dr"], np.full(100, 5.0)],
    swa=np.r_[hup_reference["swa"], np.full(100, np.nan)],
    clean=np.ones(200),
))
check("HUP SWA-ineligible epochs cannot relabel eligible epochs",
      np.array_equal(hup_reference_lab, hup_extended_lab[:100])
      and np.array_equal(hup_reference_nrem, hup_extended_nrem[:100]))

# Density-fit gates can prefer two Gaussians for a one-mode skew distribution. This deliberately
# demonstrates why these classes remain high-tail proxies rather than validated sleep states.
skewed_unimodal = np.random.RandomState(123).lognormal(0.0, 0.8, 1000)
skew_partition, skew_diagnostic = mixture_high_tail_split(skewed_unimodal)
check("GMM fit gates do not establish physiological states in skewed unimodal data",
      skew_partition is not None
      and skew_diagnostic["bic_gain_2_vs_1"] > 10
      and 0.10 <= skew_diagnostic["high_component_fraction"] <= 0.90)

# FSP automation must be allowed to say "no peak" on aperiodic data.
def pink_noise(n, local_rng):
    white = local_rng.randn(n)
    freq = np.fft.rfftfreq(n)
    spectrum = np.fft.rfft(white)
    spectrum[1:] /= np.sqrt(freq[1:])
    return np.fft.irfft(spectrum, n)

noise = np.column_stack([pink_noise(int(300 * sf), np.random.RandomState(seed))
                         for seed in range(4)])
_, is_real_noise = fsp_from(noise, sf)
check("FSP detector rejects aperiodic 1/f data", not is_real_noise)
planted = noise + 4.0 * np.sin(2 * np.pi * 13.25 * np.arange(len(noise))[:, None] / sf)
fsp_hz, is_real_planted = fsp_from(planted, sf)
check("FSP detector recovers a stable planted fast-spindle peak",
      is_real_planted and abs(fsp_hz - 13.25) <= 0.5)

print("ALL CHECKS PASSED")
