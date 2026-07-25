"""Focused regressions for neutral, offline 3D event reconstruction."""
from __future__ import annotations

import copy

import numpy as np

from cache_lc_series import (
    EVENT_3D_IED_PAD_S,
    EVENT_3D_SAMPLING_HZ,
    EVENT_3D_SO_BAND,
    EVENT_3D_SO_DURATION_S,
    EVENT_3D_SPINDLE_BAND,
    event_3d_so_candidates,
    event_3d_spindle_rms,
)
from event_3d_cache_support import (
    NEUTRAL_3D_REQUIRED_FIELDS,
    analyse_3d_cache_support,
)
from event_3D_by_stage import (
    EVENT_FS,
    IED_PAD_S,
    SO_BAND,
    SO_DUR,
    SPINDLE_BAND,
    so_event_candidates,
    spindle_rms,
)
from qc_profiles import load_qc_profile


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


def contains_non_null_p_value(value):
    if isinstance(value, dict):
        return any(
            (
                ("p_value" in str(key) or str(key) == "p")
                and item is not None
            )
            or contains_non_null_p_value(item)
            for key, item in value.items()
        )
    if isinstance(value, (list, tuple)):
        return any(contains_non_null_p_value(item) for item in value)
    return False


profile = load_qc_profile("overlap11_endpoint_local")
materialized_hup = {
    "cohort": "HUP",
    "contacts": np.asarray(["A1", "A2", "A3"]),
    "stage_lab": np.asarray(["W", "NREM", "N2", "N3", "NREM"]),
}
current_power_cache = {
    "sigma_fixed_by_contact": np.ones((3, 150)),
    # This deliberately resembles the legacy Naji/3B candidate payload.
    "so_candidate_t_A1": np.asarray([30.0, 60.0]),
}
result = analyse_3d_cache_support(
    current_power_cache, materialized_hup, profile)

check(
    "legacy neutral power cache cannot reconstruct the 3D phase endpoint",
    result["status"] == "not_reconstructable"
    and result["support_passes_profile"] is None
    and result["descriptive_effect"] is None,
)
check(
    "missing neutral event fields are enumerated rather than treated as QC failures",
    result["neutral_event_cache_schema"]["missing_fields"]
    == list(NEUTRAL_3D_REQUIRED_FIELDS)
    and not any(result["threshold_evaluability"].values()),
)
check(
    "coarse sigma and 3B SO candidates are explicitly rejected as 3D substitutes",
    any("one-second 10-15-Hz power is not a substitute" in reason
        for reason in result["support_reasons"])
    and any("different 3B Naji" in reason
            for reason in result["support_reasons"]),
)
check(
    "unreconstructable result contains no inferential p value",
    result["inference_enabled"] is False
    and not contains_non_null_p_value(result),
)

changed = copy.deepcopy(profile)
changed["endpoint_3d"]["minimum_paired_events"] = 50
changed_result = analyse_3d_cache_support(
    current_power_cache, materialized_hup, changed)
check(
    "locked threshold is recorded exactly even when its effect is not evaluable",
    changed_result["profile_thresholds"]["minimum_paired_events"] == 50
    and changed_result["support_passes_profile"] is None,
)

materialized_respect = {
    "cohort": "RESPect",
    "contacts": np.asarray(["F1", "F2"]),
    "stage_lab": np.asarray(["W", "N2", "N3"]),
}
respect_result = analyse_3d_cache_support(
    current_power_cache, materialized_respect, profile)
check(
    "RESPect is labeled outside the specified direct-stream 3D cohort",
    respect_result["status"] == "not_in_scope"
    and respect_result["cohort_in_scope"] is False,
)


# The cache producer retains fixed-duration candidates before the amplitude
# percentile.  A 1-Hz cosine has one 1-s nonnegative-to-negative cycle per
# second, all of which are scientifically eligible before thresholding.
helper_fs = 100.0
helper_t = np.arange(int(6 * helper_fs)) / helper_fs
helper_so = np.cos(2 * np.pi * helper_t)
helper_candidates = event_3d_so_candidates(helper_so, helper_fs)
check(
    "neutral cache constants match the direct-stream 3D estimator",
    EVENT_3D_SAMPLING_HZ == EVENT_FS
    and EVENT_3D_SO_BAND == SO_BAND
    and EVENT_3D_SPINDLE_BAND == SPINDLE_BAND
    and EVENT_3D_SO_DURATION_S == SO_DUR
    and EVENT_3D_IED_PAD_S == IED_PAD_S,
)
check(
    "cache SO helper retains duration-qualified candidates before percentile thresholding",
    len(helper_candidates) >= 4
    and np.allclose(
        np.diff(helper_candidates[:, 0]) / helper_fs, 1.0, atol=0.02),
)
check(
    "cache SO candidate producer matches the direct-stream 3D detector",
    np.array_equal(
        helper_candidates,
        so_event_candidates(helper_so, helper_fs),
    ),
)
constant_rms = event_3d_spindle_rms(np.ones(200), helper_fs)
check(
    "cache spindle helper uses a 200-ms RMS window",
    np.allclose(constant_rms[20:-20], 1.0),
)
helper_spindle = np.sin(2 * np.pi * 13 * helper_t)
check(
    "cache spindle RMS producer matches the direct-stream 3D detector",
    np.array_equal(
        event_3d_spindle_rms(helper_spindle, helper_fs),
        spindle_rms(helper_spindle, helper_fs),
    ),
)


# Fully synthetic 3-contact night.  Candidate amplitudes are unique, so the
# channel-night top quartile is deterministic.  Every selected SO has one
# 0.6-s duration-qualified RMS event at a constant SO phase.
n_contacts = 3
duration_s = 180
n_samples = int(duration_s * EVENT_3D_SAMPLING_HZ)
rms = np.ones((n_contacts, n_samples), float)
phase = np.full((n_contacts, n_samples), np.pi / 3, float)
valid = np.ones((n_contacts, n_samples), np.uint8)
candidate_seconds = np.arange(20, 161, 4)
candidate_samples = np.round(
    candidate_seconds * EVENT_3D_SAMPLING_HZ).astype(int)
candidate_amplitudes = np.arange(1, len(candidate_samples) + 1, dtype=float)
for contact in range(n_contacts):
    for sample in candidate_samples:
        # The peak is offset from the SO but remains within the fixed +/-2-s
        # pairing window. Twelve samples at 20 Hz gives a 0.6-s event.
        rms[contact, sample + 4:sample + 16] = 3.0

synthetic_cache = {
    "event_3d_rms_12_16_by_contact": rms,
    "event_3d_so_phase_0p16_1p25_by_contact": phase,
    "event_3d_valid_sample_mask_by_contact": valid,
    "event_3d_so_candidate_contact_index": np.repeat(
        np.arange(n_contacts), len(candidate_samples)),
    "event_3d_so_candidate_sample": np.tile(
        candidate_samples, n_contacts),
    "event_3d_so_candidate_amplitude": np.tile(
        candidate_amplitudes, n_contacts),
    "event_3d_sampling_hz": EVENT_3D_SAMPLING_HZ,
    "ep_clean_fraction_by_contact": np.ones((n_contacts, 6)),
    "ep_measured_fraction_by_contact": np.ones((n_contacts, 6)),
}
synthetic_materialized = {
    "cohort": "HUP",
    "contacts": np.asarray(["A1", "A2", "A3"]),
    "stage_lab": np.asarray(["N2", "N2", "N2", "N3", "N3", "N3"]),
}
synthetic_profile = copy.deepcopy(profile)
synthetic_profile["endpoint_3d"].update({
    "minimum_events_per_contact": 5,
    "minimum_paired_events": 15,
    "minimum_contact_acquisition_fraction": 0.5,
    "minimum_contact_event_valid_fraction": 0.5,
    "minimum_valid_nrem_fraction_per_contact": 0.8,
    "minimum_valid_nrem_seconds_per_contact": 120,
    "minimum_contacts": 2,
})
synthetic_result = analyse_3d_cache_support(
    synthetic_cache, synthetic_materialized, synthetic_profile)
pooled = synthetic_result["stage_descriptive_support"]["pooled_NREM"]
check(
    "neutral fields make every locked 3D support threshold evaluable",
    synthetic_result["neutral_event_cache_schema"]["schema_validated"]
    and all(synthetic_result["threshold_evaluability"].values()),
)
check(
    "offline detector applies event/contact/valid-NREM support and reconstructs an effect",
    synthetic_result["support_passes_profile"] is True
    and synthetic_result["descriptive_effect_available"] is True
    and pooled["n_contacts_meeting_all_contact_thresholds"] == 3
    and pooled["n_paired_events_across_qualified_contacts"] >= 15,
)
check(
    "equal-contact descriptive vector recovers the injected phase without inference",
    np.isclose(
        synthetic_result["descriptive_effect"][
            "participant_preferred_phase_deg"],
        60.0,
    )
    and np.isclose(
        synthetic_result["descriptive_effect"]["participant_R"], 1.0)
    and synthetic_result["inference_enabled"] is False
    and not contains_non_null_p_value(synthetic_result),
)

strict_profile = copy.deepcopy(synthetic_profile)
strict_profile["endpoint_3d"]["minimum_paired_events"] = 1000
strict_result = analyse_3d_cache_support(
    synthetic_cache, synthetic_materialized, strict_profile)
check(
    "a stricter support profile withholds the effect rather than calling it null",
    strict_result["status"] == "support_below_profile"
    and strict_result["support_passes_profile"] is False
    and strict_result["descriptive_effect"] is None
    and strict_result["threshold_evaluability"]["minimum_paired_events"],
)

malformed_cache = dict(synthetic_cache)
malformed_cache["event_3d_so_phase_0p16_1p25_by_contact"] = phase[:, :-1]
malformed_result = analyse_3d_cache_support(
    malformed_cache, synthetic_materialized, synthetic_profile)
check(
    "malformed neutral arrays fail closed as unreconstructable",
    malformed_result["status"] == "not_reconstructable"
    and malformed_result["support_passes_profile"] is None
    and malformed_result["neutral_event_cache_schema"]["validation_errors"],
)

fraction_cache = dict(synthetic_cache)
measured_fraction = np.ones((n_contacts, 6), float)
measured_fraction[0] = np.asarray([1.0, 0.5, 0.0, 0.0, 0.75, 0.0])
fraction_cache["ep_measured_fraction_by_contact"] = measured_fraction
strict_acquisition = copy.deepcopy(synthetic_profile)
strict_acquisition["endpoint_3d"][
    "minimum_contact_acquisition_fraction"] = 0.4
lenient_acquisition = copy.deepcopy(synthetic_profile)
lenient_acquisition["endpoint_3d"][
    "minimum_contact_acquisition_fraction"] = 0.3
strict_acquisition_result = analyse_3d_cache_support(
    fraction_cache, synthetic_materialized, strict_acquisition)
lenient_acquisition_result = analyse_3d_cache_support(
    fraction_cache, synthetic_materialized, lenient_acquisition)
check(
    "3D acquisition fraction is sample-weighted rather than finite-placeholder coverage",
    np.isclose(
        strict_acquisition_result["per_contact_acquisition_fraction"][0],
        0.375,
    )
    and not strict_acquisition_result[
        "per_contact_passes_acquisition_event_valid_support"][0]
    and lenient_acquisition_result[
        "per_contact_passes_acquisition_event_valid_support"][0],
)

unmeasured_cache = dict(synthetic_cache)
unmeasured_cache["ep_measured_fraction_by_contact"] = np.zeros(
    (n_contacts, 6), float)
unmeasured_result = analyse_3d_cache_support(
    unmeasured_cache, synthetic_materialized, synthetic_profile)
check(
    "all-unmeasured contacts have zero acquisition support rather than 100 percent",
    unmeasured_result["per_contact_acquisition_fraction"]
    == [0.0, 0.0, 0.0]
    and not any(
        unmeasured_result[
            "per_contact_passes_acquisition_event_valid_support"]),
)

out_of_range_cache = dict(synthetic_cache)
out_of_range_cache["ep_measured_fraction_by_contact"] = np.full(
    (n_contacts, 6), 1.1)
out_of_range_result = analyse_3d_cache_support(
    out_of_range_cache, synthetic_materialized, synthetic_profile)
check(
    "malformed measured acquisition support fails closed",
    out_of_range_result["status"] == "not_reconstructable"
    and out_of_range_result["support_passes_profile"] is None
    and not out_of_range_result["threshold_evaluability"][
        "minimum_contact_acquisition_fraction"],
)
