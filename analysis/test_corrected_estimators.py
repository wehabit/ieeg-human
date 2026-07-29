"""Synthetic regression tests for participant-level 3B and event-locked 3D helpers."""
import hashlib
import io
import os
import tempfile
import threading
import time

import numpy as np
from scipy import signal

import cohort_3A_cortical as cortical_module
import hup_portal as portal_module
from event_3B_mednick import (
    subject_so_triggered, so_triggered, rr_baseline_hr, FS_RR)
from event_3B_cached import (
    coverage_eligible_contacts,
    stable_stage_epoch_indices,
    stage_diagnostic_rng,
)
from event_3d_estimators import (
    EVENT_FS,
    PRODUCTION_3D_INFERENCE_ENABLED,
    bh_fdr,
    epoch_set_sample_mask,
    intervals_wholly_inside,
    pair_one_spindle_per_so,
    participant_rotation_test,
    pooled_endpoint_passes_qc,
    rayleigh,
    select_so_events,
    so_event_candidates,
    spindle_events_from_rms,
    spindle_rms,
    stage_event_pairs,
)
from cache_lc_series import (detect_so_candidates, threshold_so_candidates, sanitize_beats,
                             _write_multichannel_power, _aggregate_full_night_power,
                             aggregate_staging_features, interpolate_tachograms,
                             staging_epoch_features, _binned_power_values,
                             power_from_binned_support,
                             empty_channel_activity_extrema,
                             update_channel_activity_extrema,
                             finalize_channel_activity_qc,
                             finalize_ecg_cache_qc,
                             require_complete_acquisition)
from staging_helpers import (
    band_sos, fsp_from, nrem_mask_adaptive,
    reliable_two_state_split as mixture_high_tail_split, stage_epochs,
)
from cohort_3A_cortical import HUP_SOURCE_PINS, find_night, verify_hup_source_identity
from infraslow_rr_sigma_coherence import configure_http_session
from signal_qc import dilate_boolean_mask, ied_clean_mask
from run_qc_grid import analyse_3b
from qc_profiles import load_qc_profile
from stage_ds003848 import (
    score_stages, channel_roles_from_rows, selected_channel_name_mismatches,
    electrode_eligibility_from_rows, event_annotations_from_rows,
    event_exclusion_mask, annotation_epoch_context, annotation_second_masks,
    constrain_proxy_to_annotations, aggregate_auxiliary_epoch_features,
    auxiliary_epoch_window_features, SNAPSHOT_FILES, verify_pinned_snapshot_file,
)
from spectral_gapped import fill_short_gaps


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


rng = np.random.RandomState(7)


stage_n3_first = stage_diagnostic_rng(
    "HUPTEST_phaseII", "N3").rand(8)
stage_diagnostic_rng(
    "HUPTEST_phaseII", "N2").rand(1000)
stage_n3_after_n2 = stage_diagnostic_rng(
    "HUPTEST_phaseII", "N3").rand(8)
stage_n2 = stage_diagnostic_rng(
    "HUPTEST_phaseII", "N2").rand(8)
check(
    "3B diagnostic RNG is reproducible and isolated by subject/stage",
    np.array_equal(stage_n3_first, stage_n3_after_n2)
    and not np.array_equal(stage_n3_first, stage_n2),
)


low_coverage_warnings = finalize_ecg_cache_qc([], hr_coverage=0.25)
check(
    "completed ECG detection permits low coverage and records a warning",
    len(low_coverage_warnings) == 1
    and "HR coverage 25.0%" in low_coverage_warnings[0],
)
try:
    finalize_ecg_cache_qc(
        [{"start_s": 0.0, "duration_s": 60.0, "error": "SyntheticError: detector failed"}],
        hr_coverage=0.95,
    )
except RuntimeError as exc:
    ecg_failure_rejected = (
        "ECG detector chunks failed" in str(exc)
        and "refusing status='ok' cache" in str(exc)
    )
else:
    ecg_failure_rejected = False
check("ECG detector exceptions are fatal even when nominal coverage is high",
      ecg_failure_rejected)

try:
    require_complete_acquisition([
        {
            "start_s": 0.0,
            "duration_s": 600.0,
            "error": "ValueError: synthetic portal failure",
        },
    ])
except RuntimeError as exc:
    acquisition_failure_visible = (
        "1 acquisition chunks failed" in str(exc)
        and "start_s=0.0" in str(exc)
        and "duration_s=600.0" in str(exc)
        and "ValueError: synthetic portal failure" in str(exc)
    )
else:
    acquisition_failure_visible = False
check(
    "fatal acquisition errors retain chunk geometry and the underlying exception",
    acquisition_failure_visible,
)


class _FakeHttp:
    def __init__(self):
        self.calls = []
        self.mounts = {}

    def request(self, method, url, **kwargs):
        self.calls.append((method, url, kwargs))
        return "response"

    def mount(self, prefix, adapter):
        self.mounts[prefix] = adapter


fake_http = _FakeHttp()
configure_http_session(fake_http)
check("iEEG requests receive a bounded connect/read timeout by default",
      fake_http.request("GET", "https://example.invalid") == "response"
      and fake_http.calls[-1][2]["timeout"] == (10.0, 90.0))
fake_http.request("GET", "https://example.invalid", timeout=(1.0, 2.0))
check("an explicit iEEG request timeout is preserved",
      fake_http.calls[-1][2]["timeout"] == (1.0, 2.0))


class _FakePortalNode:
    def __init__(self, label, data_check):
        self.values = {"channelLabel": label, "dataCheck": data_check}

    def findtext(self, key):
        return self.values[key]


class _FakePortalDetail:
    def __init__(self, revision_id):
        self.portal_id = revision_id
        self.start_time = 1
        self.end_time = 2
        self.duration = 1.0
        self.number_of_samples = 100
        self.sample_rate = 100.0


class _FakePinnedPortalDataset:
    name = "HUP116_phaseII"

    def __init__(self, snapshot_id):
        expected_revision, expected_check = HUP_SOURCE_PINS[self.name]["channels"]["EKG1"]
        self.snap_id = snapshot_id
        self.ts_array = [_FakePortalNode("EKG1", expected_check)]
        self._detail = _FakePortalDetail(expected_revision)

    def get_channel_labels(self):
        return ["EKG1"]

    def get_time_series_details(self, label):
        if label != "EKG1":
            raise KeyError(label)
        return self._detail


pinned_hup = HUP_SOURCE_PINS["HUP116_phaseII"]
verified_hup = verify_hup_source_identity(
    _FakePinnedPortalDataset(pinned_hup["snapshot_id"]), [], "EKG1")
check("HUP portal input matching the pinned snapshot/revision identity is accepted",
      verified_hup["snapshot_id"] == pinned_hup["snapshot_id"])
try:
    verify_hup_source_identity(
        _FakePinnedPortalDataset("retargeted-snapshot"), [], "EKG1")
    retargeted_hup_rejected = False
except RuntimeError:
    retargeted_hup_rejected = True
check("a HUP dataset name retargeted to another snapshot is rejected",
      retargeted_hup_rejected)

scalp_mixed_labels = [
    "LF1", "LF12", "RF1", "RF12",
    "F8", "F03", "C3", "C04", "CZ", "A01", "M2", "EKG1",
]
selected_intracranial = cortical_module.cortical_channels(scalp_mixed_labels)
check("conventional and zero-padded 10-20 labels cannot become cortical contacts",
      selected_intracranial == ["LF12", "RF12"]
      and all(cortical_module.is_standard_scalp_eeg_label(label)
              for label in ("F8", "F03", "C3", "C04", "CZ", "A01", "M2")))


class _NoCandidateCache:
    files = ()


masked_3b = analyse_3b(
    _NoCandidateCache(),
    {
        "subject": "synthetic",
        "cohort": "HUP",
        "rr_4": np.asarray([], float),
        "stage_lab": np.asarray([], dtype="<U5"),
        "contacts": np.asarray(["flat", "active"]),
        "frontal_contact_mask": np.asarray([False, True]),
        "hr_coverage": 0.0,
        "hr_meets_profile": False,
        "staging_qc": {},
    },
    load_qc_profile("overlap11_endpoint_local"),
)
check("HUP 3B cannot re-admit a contact removed by raw activity QC",
      masked_3b["eligible_contact_ids"] == ["active"])


class _SyntheticNightDataset:
    """Thread-safe deterministic probe source used to verify parallel night selection."""

    def __init__(self):
        self._lock = threading.Lock()
        self.active = 0
        self.max_active = 0

    def get_data(self, start_us, duration_us, channels):
        del duration_us, channels
        with self._lock:
            self.active += 1
            self.max_active = max(self.max_active, self.active)
        try:
            time.sleep(0.005)
            start_h = start_us / 3.6e9
            local_rng = np.random.RandomState(int(round(start_h * 2)) + 100)
            t_probe = np.arange(600) / 100.0
            delta_amp = 4.0 if 2.0 <= start_h < 5.0 else 0.5
            x = (
                delta_amp * np.sin(2 * np.pi * 1.0 * t_probe)
                + np.sin(2 * np.pi * 10.0 * t_probe)
                + 0.01 * local_rng.randn(len(t_probe))
            )
            return x[:, None]
        finally:
            with self._lock:
                self.active -= 1


night_ds = _SyntheticNightDataset()
night_sequential = find_night(
    night_ds, {"CTX": 0}, "CTX", 100.0, total_h=8.0,
    required_h=3.0, probe_workers=1)
night_parallel, night_diagnostics = find_night(
    night_ds, {"CTX": 0}, "CTX", 100.0, total_h=8.0,
    required_h=3.0, probe_workers=4, return_diagnostics=True)
check("parallel sparse-probe night selection is identical to sequential selection",
      night_parallel == night_sequential)
check("night probes execute concurrently instead of serially",
      night_ds.max_active > 1)
check("night-selection provenance records complete successful probe coverage",
      night_diagnostics["n_probes"] == 16
      and night_diagnostics["n_finite_scores"] == 16
      and night_diagnostics["n_probe_failures"] == 0)

original_probe_get = portal_module.get


def _too_many_failed_probes(ds, idx, start_s, dur_s):
    if int(round(start_s / (30 * 60))) < 4:
        raise ConnectionError("synthetic portal failure")
    return ds.get_data(int(start_s * 1e6), int(dur_s * 1e6), idx)


try:
    portal_module.get = _too_many_failed_probes
    try:
        find_night(
            _SyntheticNightDataset(), {"CTX": 0}, "CTX", 100.0, total_h=8.0,
            required_h=3.0, probe_workers=1)
        excessive_probe_failures_rejected = False
    except RuntimeError:
        excessive_probe_failures_rejected = True
finally:
    portal_module.get = original_probe_get
check("night selection fails closed when more than 20% of probes fail",
      excessive_probe_failures_rejected)

pinned_bytes = b"pinned snapshot bytes"
with tempfile.NamedTemporaryFile(delete=False) as pinned_test_file:
    pinned_test_file.write(pinned_bytes)
    pinned_test_path = pinned_test_file.name
try:
    pinned_test_name = "_synthetic_snapshot_identity_test"
    SNAPSHOT_FILES[pinned_test_name] = {
        "size": len(pinned_bytes),
        "sha256": hashlib.sha256(pinned_bytes).hexdigest(),
        "s3_version_id": "synthetic-version",
    }
    check("a local OpenNeuro input matching the pinned size and SHA-256 is accepted",
          verify_pinned_snapshot_file(pinned_test_name, pinned_test_path)
          == SNAPSHOT_FILES[pinned_test_name])
    SNAPSHOT_FILES[pinned_test_name]["sha256"] = "0" * 64
    try:
        verify_pinned_snapshot_file(pinned_test_name, pinned_test_path)
        corrupted_snapshot_rejected = False
    except RuntimeError:
        corrupted_snapshot_rejected = True
    check("a nonempty but wrong local OpenNeuro input is rejected",
          corrupted_snapshot_rejected)
finally:
    SNAPSHOT_FILES.pop("_synthetic_snapshot_identity_test", None)
    os.remove(pinned_test_path)

for mask_length, half_width in ((1, 0), (8, 1), (31, 4), (100, 17)):
    raw_mask = rng.rand(mask_length) < 0.2
    expected_dilation = (
        np.convolve(
            raw_mask.astype(float), np.ones(2 * half_width + 1), mode="same") > 0
        if half_width else raw_mask
    )
    check(
        f"linear-time artifact dilation matches centred convolution ({mask_length}, {half_width})",
        np.array_equal(dilate_boolean_mask(raw_mask, half_width), expected_dilation))

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

right_edge_series = np.linspace(1.0, 1.2, 61)
right_edge_result = subject_so_triggered(
    right_edge_series,
    [[(len(right_edge_series) - int(5 * FS_RR)) / FS_RR]],
    60.0 / right_edge_series.mean(),
    np.arange(len(right_edge_series)),
    n_sur=2,
    rng=np.random.RandomState(70),
    domain="rr",
    minimum_events_per_channel=1,
    minimum_surrogate_pool_samples=1,
)
check(
    "3B accepts a half-open event window ending exactly at the final sample",
    right_edge_result is not None
    and right_edge_result["n_so_total"] == 1,
)


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

compatibility_result = so_triggered(
    hr, troughs_a, float(hr.mean()), pool, n_sur=99,
    rng=np.random.RandomState(82))
check(
    "withdrawn single-channel 3B helper routes to descriptive estimator with z/p disabled",
    compatibility_result is not None
    and compatibility_result["z"] is None
    and compatibility_result["p_upper"] is None
    and compatibility_result["inference_status"].startswith(
        "descriptive only"),
)

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
so_events = select_so_events(so_event_candidates(so, sf))
spindle_test_rms = spindle_rms(sp + 0.02 * rng.randn(len(sp)), sf)
sp_events = spindle_events_from_rms(
    spindle_test_rms, sf,
    np.percentile(spindle_test_rms[np.isfinite(spindle_test_rms)], 75))
check("3D detects complete SO events", len(so_events) > 20)
check("3D spindle detector enforces event duration and finds bursts", len(sp_events) >= 15)

selected_sp, selected_so = pair_one_spindle_per_so(
    [100, 300], [90, 105, 110, 295, 305], [1, 5, 3, 2, 4], max_distance_samples=30)
check("3D retains one maximum-amplitude spindle per independent SO",
      selected_sp.tolist() == [105, 305] and selected_so.tolist() == [100, 300])
samples_per_epoch = int(EVENT_FS * 30)
stage_one_mask = epoch_set_sample_mask(
    {1}, 3 * samples_per_epoch, EVENT_FS)
check(
    "3D stage membership requires complete event extents, not only centres",
    intervals_wholly_inside(
        stage_one_mask,
        [samples_per_epoch + 10, samples_per_epoch - 2],
        [samples_per_epoch + 20, samples_per_epoch + 10],
    ).tolist() == [True, False],
)

competition_record = {
    "total_samples": 3 * samples_per_epoch,
    "eligible_so_indices": np.asarray([samples_per_epoch + 50]),
    "eligible_so_start_samples": np.asarray([samples_per_epoch + 20]),
    "eligible_so_stop_samples_exclusive": np.asarray([
        samples_per_epoch + 80]),
    "eligible_spindle_indices": np.asarray([
        samples_per_epoch + 40, samples_per_epoch + 60]),
    "eligible_spindle_amplitudes": np.asarray([10.0, 5.0]),
    "eligible_spindle_start_samples": np.asarray([
        samples_per_epoch - 10, samples_per_epoch + 50]),
    "eligible_spindle_stop_samples_exclusive": np.asarray([
        samples_per_epoch + 50, samples_per_epoch + 70]),
    "eligible_spindle_phases": np.asarray([0.1, 0.9]),
}
stage_specific_pairs = stage_event_pairs(competition_record, {1})
check(
    "3D repeats pairing after extent-based stage restriction",
    stage_specific_pairs["indices"].tolist()
    == [samples_per_epoch + 60]
    and np.allclose(stage_specific_pairs["phases"], [0.9]),
)

bounded_rms = np.ones(100)
bounded_rms[20:32] = 3.0
bounded_peak, bounded_amplitude, bounded_start, bounded_stop = (
    spindle_events_from_rms(
        bounded_rms, 20.0, 2.0,
        return_amplitudes=True, return_bounds=True,
        require_complete_valid_extent=True))
check(
    "3D spindle detector returns exact half-open onset/offset bounds",
    bounded_peak.tolist() == [20]
    and bounded_amplitude.tolist() == [3.0]
    and bounded_start.tolist() == [20]
    and bounded_stop.tolist() == [32],
)

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
duration_probe = np.asarray(
    [1.0] * 5 + [-1.0] * 12 + [1.0] * 2 + [-1.0] * 5 + [1.0] * 5)
check(
    "3B uses Naji/Dang-Vu asymmetric half-wave durations rather than the old 0.3-1.0 s rule",
    len(detect_so_candidates(duration_probe, sf=10, edge_s=0)) == 1,
)
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

# The neutral cache must retain enough information to vary the clean-sample rule without filtering
# raw data again, while the 0.5 reconstruction remains exactly the historical array.
support_env2 = np.arange(1.0, 41.0)
support_clean = np.zeros(40, bool)
support_clean[:4] = True
support_clean[10:15] = True
support_clean[20:27] = True
support_values, support_num, support_den = _binned_power_values(
    support_env2, support_clean, 10, 4, return_support=True)
legacy_support_values = np.where(
    support_den >= 0.5 * 10,
    support_num / np.maximum(support_den, 1e-12),
    np.nan)
check("per-second power support exactly reconstructs the historical 50% rule",
      np.array_equal(np.isnan(support_values), np.isnan(legacy_support_values))
      and np.array_equal(
          support_values[np.isfinite(support_values)],
          legacy_support_values[np.isfinite(legacy_support_values)])
      and np.array_equal(
          support_values[np.isfinite(support_values)],
          power_from_binned_support(
              support_num, support_den, 10, 0.5)[np.isfinite(support_values)])
      and np.array_equal(support_den, [4, 5, 7, 0])
      and np.isnan(support_values[[0, 3]]).all()
      and np.isfinite(support_values[[1, 2]]).all())
support_70 = power_from_binned_support(support_num, support_den, 10, 0.7)
check("per-second clean support threshold is reversible offline",
      np.isnan(support_70[[0, 1, 3]]).all()
      and np.isclose(support_70[2], support_num[2] / support_den[2]))

# Missing raw samples may be filled only as a numerical scaffold for filtering. They must remain
# missing in derived power, rather than becoming a plausible trough that inflates coverage.
with_gap = base.copy()
with_gap[int(60 * sf_power):int(70 * sf_power)] = np.nan
dest_gap = np.full((1, 180), np.nan)
dest_gap_num = np.full((1, 180), np.nan)
dest_gap_den = np.zeros((1, 180), dtype=np.uint32)
_write_multichannel_power(
    with_gap[:, None], sos, sf_power, 180, 180, 0, dest_gap,
    numerator_dest=dest_gap_num, clean_sample_count_dest=dest_gap_den)
check("raw NaN gaps remain ineligible in derived sigma power",
      np.isnan(dest_gap[0, 60:70]).all()
      and np.allclose(
          dest_gap,
          power_from_binned_support(dest_gap_num, dest_gap_den, int(sf_power), 0.5),
          equal_nan=True))

# A brief high-amplitude artifact can dominate a 30-s Welch spectrum even if most samples remain
# clean. Gap-aware staging must reject every overlapping 4-s periodogram while retaining clean,
# temporally separate windows from the same epoch.
stage_sf = 200.0
stage_t = np.arange(int(30 * stage_sf)) / stage_sf
clean_stage_signal = np.sin(2 * np.pi * stage_t) + 0.05 * rng.randn(len(stage_t))
clean_stage_mask = np.ones(len(stage_t), bool)
measured_stage_mask = np.ones(len(stage_t), bool)
clean_dr, clean_swa, clean_stage_details = staging_epoch_features(
    clean_stage_signal, clean_stage_mask, measured_stage_mask, stage_sf,
    return_details=True)
welch_f, welch_p = signal.welch(
    clean_stage_signal, stage_sf, nperseg=int(4 * stage_sf))
welch_swa_band = (welch_f >= 0.5) & (welch_f < 4.0)
welch_total_band = (welch_f >= 0.5) & (welch_f < 25.0)
welch_swa = float(np.trapezoid(
    welch_p[welch_swa_band], welch_f[welch_swa_band]))
welch_dr = welch_swa / float(np.trapezoid(
    welch_p[welch_total_band], welch_f[welch_total_band]))
check("fourteen clean staging windows reproduce scipy Welch",
      clean_stage_details["n_valid_windows"] == 14
      and np.isclose(clean_swa, welch_swa, rtol=1e-12, atol=1e-12)
      and np.isclose(clean_dr, welch_dr, rtol=1e-12, atol=1e-12))
artifact_stage_signal = clean_stage_signal.copy()
artifact_region = (stage_t >= 10) & (stage_t < 12)
artifact_stage_signal[artifact_region] += 100 * np.sin(
    2 * np.pi * stage_t[artifact_region])
artifact_clean_mask = clean_stage_mask.copy()
artifact_clean_mask[artifact_region] = False
artifact_dr, artifact_swa, artifact_stage_details = staging_epoch_features(
    artifact_stage_signal, artifact_clean_mask, measured_stage_mask, stage_sf,
    return_details=True)
check("masked artifact samples never enter a retained staging periodogram",
      artifact_stage_details["n_valid_windows"] == 12
      and not artifact_stage_details["valid_window_mask"][4]
      and not artifact_stage_details["valid_window_mask"][5])
check("clean subwindows recover staging power instead of discarding the whole epoch",
      np.isfinite(clean_dr) and np.isfinite(clean_swa)
      and np.isfinite(artifact_dr) and np.isfinite(artifact_swa)
      and abs(np.log(artifact_swa / clean_swa)) < 0.02
      and abs(artifact_dr - clean_dr) < 0.01)
dirty_dr, dirty_swa, dirty_details = staging_epoch_features(
    artifact_stage_signal, np.zeros_like(artifact_clean_mask),
    measured_stage_mask, stage_sf, return_details=True)
check("an epoch without a complete clean Welch window remains unavailable",
      dirty_details["n_valid_windows"] == 0
      and np.isnan(dirty_dr) and np.isnan(dirty_swa))

# Applying the historical reference support after neutral feature extraction reproduces the former
# all-clean eligibility without destroying the recoverable per-window data in the cache.
support_counts = np.array([[14, 14], [14, 12], [14, 14]])
reference_dr, reference_swa, _, reference_qc = aggregate_staging_features(
    np.ones((3, 2)), np.ones((3, 2)), np.ones((3, 2)), np.ones(3, bool),
    valid_window_count_by_contact=support_counts, min_valid_windows=14,
    min_contact_feature_coverage=0.5, min_contacts=3,
    min_contact_fraction_per_epoch=1.0)
check("historical all-clean support is an offline profile rather than destructive extraction",
      reference_qc["minimum_valid_welch_windows"] == 14
      and np.isfinite(reference_dr[0]) and np.isnan(reference_dr[1])
      and np.isfinite(reference_swa[0]) and np.isnan(reference_swa[1]))

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
      fixed_qc["n_selected_contacts"] == 0
      and np.isnan(fixed_dr).all() and np.isnan(fixed_swa).all())

# A contact with only one valid 30-s feature must not enter the fixed staging denominator merely
# because its independent 1-s sigma coverage happened to pass.
sparse_dr = np.ones((6, 100))
sparse_swa = np.ones((6, 100))
sparse_dr[3:, 1:] = np.nan
sparse_swa[3:, 1:] = np.nan
qualified_dr, qualified_swa, _, qualified_qc = aggregate_staging_features(
    sparse_dr, sparse_swa, np.ones_like(sparse_dr), np.ones(6, bool))
check("staging-sparse contacts cannot invalidate a stable full-night contact set",
      qualified_qc["n_selected_contacts"] == 3
      and qualified_qc["required_contact_count"] == 3
      and np.isfinite(qualified_dr).all() and np.isfinite(qualified_swa).all()
      and np.allclose(qualified_qc["per_contact_feature_coverage"][:3], 1.0)
      and np.allclose(qualified_qc["per_contact_feature_coverage"][3:], 0.01))

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

# Finite derived power is not evidence that the source voltage varied: detrending/filter roundoff
# can turn a constant raw channel into tiny positive band power.  The raw activity mask must be
# inherited by both sigma aggregation and staging-contact selection.
activity_state = empty_channel_activity_extrema(3)
update_channel_activity_extrema(
    activity_state,
    np.asarray([
        np.zeros(100),
        1e6 + np.linspace(0.0, 1e-10, 100),
        np.sin(np.linspace(0.0, 4 * np.pi, 100)),
    ]))
activity_qc = finalize_channel_activity_qc(activity_state)
finite_power = np.ones((3, 100))
nonflat_power, nonflat_power_qc = _aggregate_full_night_power(
    finite_power, eligible_channels=activity_qc["nonflat_mask"],
    min_contacts=1, return_details=True)
nonflat_dr, nonflat_swa, _, nonflat_staging_qc = aggregate_staging_features(
    finite_power, finite_power, finite_power, activity_qc["nonflat_mask"],
    min_contacts=1, min_contact_feature_coverage=0.0)
check("numerically flat raw channels cannot pass sigma or staging contact selection",
      activity_qc["nonflat_mask"].tolist() == [False, False, True]
      and nonflat_power_qc["selected_mask"].tolist() == [False, False, True]
      and nonflat_staging_qc["selected_contact_mask"].tolist()
      == [False, False, True]
      and np.isfinite(nonflat_power).all()
      and np.isfinite(nonflat_dr).all()
      and np.isfinite(nonflat_swa).all())

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
    complete_clean = np.asarray([
        clean[int(start):int(stop)].all()
        for start, stop in candidates[:, 2:4]
    ])
    candidates = candidates[complete_clean]
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

electrode_qc = electrode_eligibility_from_rows([
    {"name": ".....", "group": "grid"},
    {"name": ".....", "group": "grid"},
    {"name": "A1", "group": "grid", "soz": "no", "resected": "no", "edge": "no",
     "Destrieux_label_text": "G_front_middle"},
    {"name": "A2", "group": "depth", "soz": "yes", "graymatter": "yes",
     "Destrieux_label_text": "G_parietal_sup"},
], ["A1", "A2"])
check("unused duplicate electrode placeholders do not invalidate the exact requested join",
      np.array_equal(electrode_qc["selected_mask"], [True, False]))
check("RESPect ROI masks are derived only after pathology exclusions",
      np.array_equal(electrode_qc["frontal_mask"], [True, False])
      and not electrode_qc["parietal_mask"].any()
      and electrode_qc["exclusion_reasons"]["A2"] == ["soz"])

synthetic_events = event_annotations_from_rows([
    {"onset": "0", "duration": "60", "offset": "60", "sample_start": "0",
     "sample_end": "600", "trial_type": "sleep", "sub_type": "NREM",
     "electrodes_involved_onset": "all", "electrodes_involved_offset": "all"},
    {"onset": "60", "duration": "60", "offset": "120", "sample_start": "600",
     "sample_end": "1200", "trial_type": "sleep", "sub_type": "REM",
     "electrodes_involved_onset": "all", "electrodes_involved_offset": "all"},
    {"onset": "120", "duration": "30", "offset": "150", "sample_start": "1200",
     "sample_end": "1500", "trial_type": "sleep-wake transition", "sub_type": "unknown",
     "electrodes_involved_onset": "all", "electrodes_involved_offset": "all"},
    {"onset": "150", "duration": "30", "offset": "180", "sample_start": "1500",
     "sample_end": "1800", "trial_type": "sleep", "sub_type": "unknown",
     "electrodes_involved_onset": "all", "electrodes_involved_offset": "all"},
    {"onset": "10", "duration": "2", "offset": "12", "sample_start": "100",
     "sample_end": "120", "trial_type": "artefact", "sub_type": "n/a",
     "electrodes_involved_onset": "A1", "electrodes_involved_offset": "A1"},
], total_s=180, sf=10)
annotation_labels, annotation_sources = annotation_epoch_context(
    synthetic_events, 6, epoch_s=30)
primary_labels, _ = constrain_proxy_to_annotations(
    np.array(["N2", "N3", "N2", "N3", "N2", "N3"]),
    annotation_labels, annotation_sources)
sensitivity_labels, _ = constrain_proxy_to_annotations(
    np.array(["N2", "N3", "N2", "N3", "N2", "N3"]),
    annotation_labels, annotation_sources, allow_proxy_unknown_sleep=True)
check("author disturbances break otherwise author-defined stable sleep epochs",
      np.array_equal(primary_labels[:5], ["", "NREM", "R", "R", ""]))
check("author-unknown sleep stays unclassified in the primary stage array",
      primary_labels[5] == "" and sensitivity_labels[5] == "N3")

selection_annotations = [
    dict(onset_s=0.25, stop_s=30.25, category="sws_selection", contacts=["all"]),
]
selection_labels, selection_sources = annotation_epoch_context(
    selection_annotations, 2, epoch_s=30)
selection_masks = annotation_second_masks(selection_annotations, total_s=60)
check("a curated SWS selection is N3 and cannot simultaneously be reported as awake",
      selection_labels[0] == "N3"
      and selection_sources[0] == "author_sws_selection"
      and not selection_masks["awake"][:31].any())

check("contact-specific artifact annotations mask only that requested contact",
      event_exclusion_mask(
          synthetic_events, 0, 1800, 10, channel="A1", pad_s=0).sum() == 20
      and not event_exclusion_mask(
          synthetic_events, 0, 1800, 10, channel="A2", pad_s=0).any())

auxiliary, auxiliary_qc = aggregate_auxiliary_epoch_features(np.array([
    [1.0, 2.0, 4.0, 8.0],
    [10.0, 20.0, 40.0, 80.0],
    [1.0, np.nan, np.nan, np.nan],
]), min_channel_coverage=0.75)
check("all stable EMG/EOG channels contribute after within-channel scale normalization",
      auxiliary_qc["n_selected_channels"] == 2
      and np.allclose(auxiliary / auxiliary[0], [1.0, 2.0, 4.0, 8.0]))

# EMG/EOG now use the same complete 4-s, 2-s-stride support geometry as iEEG staging.  Fully
# measured epochs must retain the exact former full-epoch definition, while a gap rejects only
# overlapping windows and remains available to a profile that permits partial support.
aux_sf = 100.0
aux_t = np.arange(int(30 * aux_sf)) / aux_sf
aux_emg_signal = (1.0 + 0.2 * np.sin(2 * np.pi * 0.07 * aux_t)) * np.sin(
    2 * np.pi * 20 * aux_t)
aux_eog_signal = np.sin(2 * np.pi * 1.2 * aux_t) + 0.1 * np.sin(
    2 * np.pi * 0.4 * aux_t)
aux_measured = np.ones(len(aux_t), bool)
clean_emg_value, clean_emg_details = auxiliary_epoch_window_features(
    aux_emg_signal, aux_measured, aux_sf, "emg_rms", return_details=True)
clean_eog_value, clean_eog_details = auxiliary_epoch_window_features(
    aux_eog_signal, aux_measured, aux_sf, "eog_variance", return_details=True)
check("all-clean auxiliary windows exactly preserve former full-epoch features",
      clean_emg_details["n_valid_windows"] == 14
      and clean_eog_details["n_valid_windows"] == 14
      and np.isclose(clean_emg_value, np.sqrt(np.mean(aux_emg_signal ** 2)),
                     rtol=0, atol=1e-15)
      and np.isclose(clean_eog_value, np.var(aux_eog_signal), rtol=0, atol=1e-15))
aux_gap = (aux_t >= 10) & (aux_t < 12)
partial_measured = aux_measured.copy()
partial_measured[aux_gap] = False
contaminated_emg = aux_emg_signal.copy()
contaminated_eog = aux_eog_signal.copy()
contaminated_emg[aux_gap] = 1e6
contaminated_eog[aux_gap] = 1e6
partial_emg_value, partial_emg_details = auxiliary_epoch_window_features(
    contaminated_emg, partial_measured, aux_sf, "emg_rms", return_details=True)
partial_eog_value, partial_eog_details = auxiliary_epoch_window_features(
    contaminated_eog, partial_measured, aux_sf, "eog_variance", return_details=True)
check("auxiliary gaps reject only overlapping complete windows",
      partial_emg_details["n_valid_windows"] == 12
      and partial_eog_details["n_valid_windows"] == 12
      and not partial_emg_details["valid_window_mask"][4:6].any()
      and not partial_eog_details["valid_window_mask"][4:6].any()
      and np.isfinite(partial_emg_value) and np.isfinite(partial_eog_value)
      and np.isnan(partial_emg_details["window_power"][4:6]).all()
      and np.isnan(partial_eog_details["window_power"][4:6]).all())
aux_values = np.array([
    [clean_emg_value, partial_emg_value],
    [2 * clean_emg_value, 2 * partial_emg_value],
])
aux_window_counts = np.array([[14, 12], [14, 12]])
aux_reference, _ = aggregate_auxiliary_epoch_features(
    aux_values, min_channel_coverage=0.0, min_channel_fraction_per_epoch=1.0,
    valid_window_count_by_channel=aux_window_counts, min_valid_windows=14)
aux_partial, aux_partial_qc = aggregate_auxiliary_epoch_features(
    aux_values, min_channel_coverage=0.0, min_channel_fraction_per_epoch=1.0,
    valid_window_count_by_channel=aux_window_counts, min_valid_windows=12)
check("auxiliary complete-window and channel gates are selectable offline",
      np.isfinite(aux_reference[0]) and np.isnan(aux_reference[1])
      and np.isfinite(aux_partial).all()
      and aux_partial_qc["minimum_valid_windows"] == 12
      and aux_partial_qc["minimum_channel_coverage"] == 0.0)

isolated_stage_labels = np.array(["N2", "W"] * 20)
stable_stage_labels = np.array(["W", "N2", "N2", "N2", "N2", "N2", "N2", "W"])
check("isolated 30-s stage labels cannot satisfy Naji's stable 3-min-bin rule",
      len(stable_stage_epoch_indices(isolated_stage_labels, "N2")) == 0)
check("an uninterrupted 3-min stage run satisfies Naji's stable-bin rule",
      np.array_equal(
          stable_stage_epoch_indices(stable_stage_labels, "N2"),
          np.arange(1, 7)))
interrupted_stage_labels = np.array(["N2"] * 6 + [""] + ["N2"] * 5)
check("an author-disturbed epoch breaks rather than bridges a Naji stable-stage run",
      np.array_equal(
          stable_stage_epoch_indices(interrupted_stage_labels, "N2"),
          np.arange(6)))

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
