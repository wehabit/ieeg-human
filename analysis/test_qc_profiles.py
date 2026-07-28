"""Regression checks for immutable QC profiles and deterministic sensitivity grids."""
from __future__ import annotations

import copy

import numpy as np

from lecci_faithful_3A import core_study_nrem_mask
from materialize_qc_cache import _materialize_auxiliary, materialize
from pipeline_version import CACHE_SCHEMA_VERSION
from qc_profiles import (
    expand_qc_grid,
    load_qc_profile,
    load_profile_set,
    qc_profile_sha256,
    validate_qc_profile,
    validated_staging_calibration,
)
from run_qc_grid import analyse_3b, materialization_diagnostics


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


profile_set = load_profile_set()
calibration_provenance = validated_staging_calibration()
audit80 = load_qc_profile("audit80")
overlap11 = load_qc_profile("overlap11_endpoint_local")
check(
    "staging calibration pins the exact purpose-specific source files",
    set(calibration_provenance["calibration_source_files_sha256"]) == {
        "analysis/calibrate_staging_windows.py",
        "analysis/pipeline_version.py",
    },
)
check(
    "historical audit profile records the exact former 80% power gates",
    audit80["power"]["minimum_contact_coverage"] == 0.8
    and audit80["power"]["minimum_contact_fraction_per_bin"] == 0.8
    and audit80["power"]["minimum_aggregate_coverage"] == 0.8,
)
check(
    "historical all-clean epoch rule is represented as all 14 fixed Welch windows",
    audit80["staging"]["minimum_valid_welch_windows"] == 14,
)
check(
    "historical auxiliary staging gates are represented explicitly",
    audit80["staging"]["auxiliary_aggregation"]
    == "fixed_channel_geometric_median"
    and audit80["staging"]["minimum_auxiliary_channel_coverage"] == 0.8
    and audit80["staging"][
        "minimum_auxiliary_channel_fraction_per_epoch"] == 0.8
    and audit80["staging"]["minimum_valid_auxiliary_windows"] == 14,
)
check(
    "conservative endpoint-local profile records eleven-window staging support",
    overlap11["staging"]["minimum_valid_welch_windows"] == 11
    and overlap11["staging"]["minimum_valid_auxiliary_windows"] == 14
    and overlap11["power"]["aggregation"] == "overlap_connected_median_polish"
    and overlap11["staging"]["aggregation"] == "overlap_connected_median_polish",
)

first = expand_qc_grid("coverage_oat_v1")
second = expand_qc_grid("coverage_oat_v1")
check(
    "one-factor coverage grid expansion is deterministic",
    first == second
    and [qc_profile_sha256(value) for value in first]
    == [qc_profile_sha256(value) for value in second],
)
coverage_grid = profile_set["grids"]["coverage_oat_v1"]
expected_n = 1 + sum(
    len(set(axis.get("levels", coverage_grid["levels"])) - {
        audit80[axis["path"].split(".")[0]][axis["path"].split(".")[1]]
    })
    for axis in coverage_grid["axes"]
)
check(
    "coverage grid changes one axis at a time and skips the duplicated 80% baseline",
    len(first) == expected_n
    and sum(value.get("changed_axis") is None for value in first) == 1,
)
for value in first[1:]:
    changed_path = value["changed_path"]
    changed_top = changed_path.split(".")[0]
    unchanged_sections = set(audit80) - {
        "profile_id",
        "role",
        "selection_blinded_to_endpoint_values",
        changed_top,
    }
    check(
        f"{value['profile_id']} leaves every unrelated profile section unchanged",
        all(value[key] == audit80[key] for key in unchanged_sections),
    )

boundary = copy.deepcopy(audit80)
boundary["profile_id"] = "boundary"
boundary["hr"]["minimum_coverage"] = 0.7499
check(
    "fractional threshold values are preserved exactly rather than rounded",
    validate_qc_profile(boundary)["hr"]["minimum_coverage"] == 0.7499,
)

invalid = copy.deepcopy(audit80)
invalid["profile_id"] = "invalid"
invalid["power"]["minimum_contact_coverage"] = 1.01
try:
    validate_qc_profile(invalid)
except ValueError:
    invalid_rejected = True
else:
    invalid_rejected = False
check("out-of-range coverage thresholds are rejected", invalid_rejected)

tampered = copy.deepcopy(audit80)
before = qc_profile_sha256(tampered)
tampered["endpoint_3b"]["minimum_so_per_contact"] = 29
check(
    "changing one gate changes the canonical profile hash",
    qc_profile_sha256(tampered) != before,
)

event_first = expand_qc_grid("event_count_oat_v1")
event_second = expand_qc_grid("event_count_oat_v1")
event_grid = profile_set["grids"]["event_count_oat_v1"]
expected_event_n = 1 + sum(
    len(set(axis["levels"]) - {
        overlap11[axis["path"].split(".")[0]][axis["path"].split(".")[1]]
    })
    for axis in event_grid["axes"]
)
check(
    "axis-specific event-count grid expansion is deterministic and complete",
    event_first == event_second
    and len(event_first) == expected_event_n == 25
    and len({value["profile_id"] for value in event_first}) == len(event_first),
)
check(
    "event-count grid retains its calibrated endpoint-local base exactly",
    event_first[0] == overlap11,
)
for value in event_first[1:]:
    changed_path = value["changed_path"]
    changed_top = changed_path.split(".")[0]
    unchanged_sections = set(overlap11) - {
        "profile_id",
        "role",
        "selection_blinded_to_endpoint_values",
        changed_top,
    }
    check(
        f"{value['profile_id']} changes only its declared event-grid section",
        all(value[key] == overlap11[key] for key in unchanged_sections)
        and value["grid_id"] == "event_count_oat_v1",
    )


class ArrayCache:
    """Small NPZ-like object for endpoint-selection regressions."""

    def __init__(self, values):
        self.values = values
        self.files = list(values)

    def __getitem__(self, key):
        return self.values[key]


# The overlap profile's staging endpoint must start from every anatomy-qualified
# contact, even when a disjoint power endpoint selects a different contact.
n_epochs = 60
endpoint_cache = ArrayCache({
    "cache_schema_version": np.asarray(CACHE_SCHEMA_VERSION),
    "subject": np.asarray("synthetic"),
    "cortical_chans": np.asarray(["power_only", "stage_a", "stage_b"]),
    "cortical_signal_nonflat_mask": np.ones(3, bool),
    "sigma_fixed_by_contact": np.r_[
        np.ones((1, 180)), np.full((2, 180), np.nan)],
    "swa_by_contact": np.r_[
        np.ones((1, 180)), np.full((2, 180), np.nan)],
    "sigma_fixed_power_numerator_by_contact": np.r_[
        np.full((1, 180), 10.0), np.full((2, 180), np.nan)],
    "sigma_fixed_clean_sample_count_by_contact": np.r_[
        np.full((1, 180), 10), np.zeros((2, 180), int)],
    "swa_power_numerator_by_contact": np.r_[
        np.full((1, 180), 10.0), np.full((2, 180), np.nan)],
    "swa_clean_sample_count_by_contact": np.r_[
        np.full((1, 180), 10), np.zeros((2, 180), int)],
    "power_samples_per_second": np.asarray(10),
    "power_historical_minimum_clean_fraction_per_second": np.asarray(0.5),
    "ep_dr_by_contact": np.r_[
        np.full((1, n_epochs), np.nan),
        np.full((1, n_epochs), 0.35),
        np.full((1, n_epochs), 0.45),
    ],
    "ep_swa_by_contact": np.r_[
        np.full((1, n_epochs), np.nan),
        np.ones((1, n_epochs)),
        np.full((1, n_epochs), 2.0),
    ],
    "ep_clean_fraction_by_contact": np.ones((3, n_epochs)),
    "ep_valid_welch_window_mask_by_contact": np.r_[
        np.zeros((1, n_epochs, 14), bool),
        np.ones((2, n_epochs, 14), bool),
    ],
    "hr_1": np.ones(180),
    "rr_4": np.ones(720),
})
endpoint_materialized = materialize(endpoint_cache, overlap11)
check(
    "overlap staging selection is endpoint-local rather than inherited from power",
    endpoint_materialized["global_power_qc"]["sigma"][
        "selected_contact_mask"].tolist() == [True, False, False]
    and endpoint_materialized["staging_qc"][
        "selected_contact_mask"].tolist() == [False, True, True]
    and np.isfinite(endpoint_materialized["ep_dr"]).all(),
)

legacy_v8_values = copy.deepcopy(endpoint_cache.values)
legacy_v8_values["cache_schema_version"] = np.asarray(
    "2026-07-neutral-per-contact-gap-aware-source-pin-v8")
legacy_v8_values.pop("cortical_signal_nonflat_mask")
legacy_v8_materialized = materialize(ArrayCache(legacy_v8_values), overlap11)
check(
    "legacy v8 caches remain readable with explicit missing-flat-QC provenance",
    legacy_v8_materialized["cortical_signal_nonflat_mask"].all()
    and "not applied" in legacy_v8_materialized[
        "cortical_signal_activity_qc_provenance"],
)

# Per-second support must be rematerialized from numerator/denominator rather
# than inheriting the historical 50%-clean compatibility array.
power_threshold_cache = copy.deepcopy(endpoint_cache.values)
power_threshold_cache["sigma_fixed_power_numerator_by_contact"][0, 0] = 4.0
power_threshold_cache["sigma_fixed_clean_sample_count_by_contact"][0, 0] = 4
power_threshold_cache["swa_power_numerator_by_contact"][0, 0] = 4.0
power_threshold_cache["swa_clean_sample_count_by_contact"][0, 0] = 4
power_threshold_cache["sigma_fixed_by_contact"][0, 0] = np.nan
power_threshold_cache["swa_by_contact"][0, 0] = np.nan
low_support_profile = copy.deepcopy(overlap11)
low_support_profile["profile_id"] = "power_support_25"
low_support_profile["power"]["minimum_clean_fraction_per_second"] = 0.25
low_support_result = materialize(
    ArrayCache(power_threshold_cache), low_support_profile)
check(
    "power clean-sample threshold is applied offline from reversible support",
    np.isfinite(low_support_result["sigma_global"][0]),
)

aux_window_mask = np.zeros((2, 3, 14), bool)
aux_window_mask[:, 0, :] = True
aux_window_mask[:, 1, :12] = True
aux_window_mask[:, 2, :10] = True
aux_cache = ArrayCache({
    "ep_emg_by_channel": np.asarray([
        [1.0, 2.0, 3.0],
        [10.0, 20.0, 30.0],
    ]),
    "ep_emg_valid_window_mask_by_channel": aux_window_mask,
})
fixed_aux_profile = copy.deepcopy(audit80)
fixed_aux_profile["staging"]["minimum_auxiliary_channel_coverage"] = 0.0
fixed_aux_profile[
    "staging"]["minimum_auxiliary_channel_fraction_per_epoch"] = 1.0
fixed_aux, fixed_aux_qc = _materialize_auxiliary(
    aux_cache, "emg", fixed_aux_profile, 3)
partial_aux_profile = copy.deepcopy(overlap11)
partial_aux_profile[
    "staging"]["minimum_valid_auxiliary_windows"] = 11
overlap_aux, overlap_aux_qc = _materialize_auxiliary(
    aux_cache, "emg", partial_aux_profile, 3)
check(
    "RESP auxiliary features are rematerialized at the profile's complete-window threshold",
    np.isfinite(fixed_aux[0])
    and np.isnan(fixed_aux[1:]).all()
    and np.isfinite(overlap_aux[:2]).all()
    and np.isnan(overlap_aux[2])
    and fixed_aux_qc["minimum_valid_auxiliary_windows"] == 14
    and overlap_aux_qc["minimum_valid_auxiliary_windows"] == 11,
)
diagnostic = materialization_diagnostics({
    "stage_lab": np.asarray(["NREM", ""]),
    "stage_lab_proxy": np.asarray(["NREM", "N2"]),
    "ep_emg": np.asarray([1.0, np.nan]),
    "ep_eog": np.asarray([1.0, 2.0]),
    "auxiliary_qc": {"status": "ok"},
    "qc_profile": partial_aux_profile,
})
check(
    "profile outputs expose compact auxiliary and final-label diagnostics",
    diagnostic["n_finite_emg_epochs"] == 1
    and diagnostic["n_finite_eog_epochs"] == 2
    and diagnostic["n_final_labels_different_from_proxy"] == 1
    and diagnostic["minimum_valid_auxiliary_windows"] == 11,
)

long_labels = np.asarray(["W"] * 10 + ["NREM"] * 490)
core_nrem, core_window = core_study_nrem_mask(long_labels)
check(
    "3A is restricted to Lecci's first 210 minutes from the available sleep-onset proxy",
    core_nrem.sum() == 420
    and core_window["start_epoch"] == 10
    and core_window["stop_epoch"] == 430
    and not core_nrem[430:].any(),
)

# Pooled exploratory NREM is the physiological union, including transitions
# between exact labels.  It must not mean only epochs literally named NREM.
pooled_labels = np.asarray(
    ["N2"] * 3 + ["N3"] * 3 + ["NREM"] * 3 + [""] * 2)
pooled_rr = np.ones(len(pooled_labels) * 30 * 4)
empty_event_cache = ArrayCache({})
pooled_materialized = {
    "rr_4": pooled_rr,
    "stage_lab": pooled_labels,
    "contacts": np.asarray(["A"]),
    "frontal_contact_mask": np.asarray([True]),
    "cohort": "RESPect",
    "hr_coverage": 1.0,
    "hr_meets_profile": True,
    "subject": "synthetic",
}
pooled_result = analyse_3b(
    empty_event_cache, pooled_materialized, overlap11)["stages"]["NREM"]
check(
    "exploratory pooled NREM includes contiguous N2/N3/NREM transitions",
    pooled_result["raw_epochs"] == 9
    and pooled_result["stable_epochs"] == 9,
)

rr_gate_profile = copy.deepcopy(overlap11)
rr_gate_profile["profile_id"] = "rr_gate_regression"
rr_gate_profile["endpoint_3b"]["minimum_finite_rr_samples"] = (
    len(pooled_rr) + 1)
rr_gate_result = analyse_3b(
    empty_event_cache, pooled_materialized, rr_gate_profile)["stages"]["NREM"]
check(
    "3B enforces its configured whole-record finite-RR support gate",
    not rr_gate_result["available_under_profile"]
    and any(
        "finite RR samples" in reason
        for reason in rr_gate_result["support_reasons"]),
)
