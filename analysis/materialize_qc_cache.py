"""Apply a validated QC profile to one neutral v8 derived cache offline."""
from __future__ import annotations

import argparse
import json
import os
import warnings

import numpy as np

from cache_lc_series import (
    _aggregate_full_night_power,
    aggregate_staging_features,
    power_from_binned_support,
)
from staging_helpers import stage_epochs
from overlap_aggregate import (
    overlap_connected_aggregate,
    overlap_connected_staging,
)
from pipeline_version import (
    CACHE_SCHEMA_VERSION,
    atomic_savez,
    npz_scalar_text,
)
from qc_profiles import (
    load_qc_profile,
    manifest_qc_config,
    validate_qc_profile,
)
from stage_ds003848 import (
    aggregate_auxiliary_epoch_features,
    constrain_proxy_to_annotations,
    score_stages,
)


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEGACY_CACHE_SCHEMA_V8 = "2026-07-neutral-per-contact-gap-aware-source-pin-v8"


def _required(cache, key):
    if key not in cache.files:
        raise RuntimeError(
            f"neutral cache lacks {key!r}; rebuild with schema {CACHE_SCHEMA_VERSION}")
    return np.asarray(cache[key])


def _fixed_power(values, eligible, profile):
    power = profile["power"]
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="All-NaN slice encountered",
            category=RuntimeWarning,
        )
        aggregate, details = _aggregate_full_night_power(
            values,
            eligible_channels=eligible,
            min_contact_coverage=power["minimum_contact_coverage"],
            min_contacts=power["minimum_contacts"],
            min_contact_fraction_per_bin=power[
                "minimum_contact_fraction_per_bin"],
            return_details=True,
        )
    # The historical fixed-contact mean is a direct calculation, not an
    # iterative fit.  Mark convergence as not applicable/satisfied so the
    # fail-closed gate reserved for median polish does not withhold it merely
    # because the fixed estimator has no ``fit`` object.
    details.update(
        fit_required=False,
        support_passes_fit_convergence=True,
        fit_convergence_status="not_applicable_closed_form",
    )
    return aggregate, details


def _power_contact_values(cache, band, profile):
    """Reconstruct one contact-by-second band at the profile's support threshold."""
    numerator = _required(
        cache, f"{band}_power_numerator_by_contact").astype(float)
    denominator = _required(
        cache, f"{band}_clean_sample_count_by_contact").astype(float)
    samples_per_second = int(
        _required(cache, "power_samples_per_second").reshape(()).item())
    minimum = float(
        profile["power"]["minimum_clean_fraction_per_second"])
    values = power_from_binned_support(
        numerator,
        denominator,
        samples_per_second,
        minimum_clean_fraction_per_second=minimum,
    )

    # v8 keeps the historical 50%-clean arrays only as a compatibility copy.  Fail
    # closed if a freshly rebuilt cache cannot reproduce that copy exactly; otherwise
    # a supposedly neutral threshold sweep would not be changing just one QC choice.
    legacy_key = f"{band}_by_contact"
    historical_minimum = (
        float(np.asarray(cache[
            "power_historical_minimum_clean_fraction_per_second"]).item())
        if "power_historical_minimum_clean_fraction_per_second" in cache.files
        else 0.5
    )
    compatibility_checked = False
    if minimum == historical_minimum and legacy_key in cache.files:
        legacy = np.asarray(cache[legacy_key], float)
        compatibility_checked = True
        if not np.array_equal(values, legacy, equal_nan=True):
            raise RuntimeError(
                f"{band} support reconstruction disagrees with its historical "
                f"{historical_minimum:g} compatibility array")
    return values, dict(
        minimum_clean_fraction_per_second=minimum,
        samples_per_second=samples_per_second,
        compatibility_copy_checked=compatibility_checked,
        support_source=(
            f"{band}_power_numerator_by_contact and "
            f"{band}_clean_sample_count_by_contact"),
    )


def _power_matrices(cache, profile):
    sigma_values, sigma_support = _power_contact_values(
        cache, "sigma_fixed", profile)
    swa_values, swa_support = _power_contact_values(cache, "swa", profile)
    if sigma_values.shape != swa_values.shape:
        raise RuntimeError("neutral sigma/SWA contact matrices do not align")
    return sigma_values, swa_values, sigma_support, swa_support


def _power_pair(cache, eligible, profile, prepared=None):
    if prepared is None:
        prepared = _power_matrices(cache, profile)
    sigma_values, swa_values, sigma_support, swa_support = prepared
    if profile["power"]["aggregation"] == "fixed_contact_mean":
        # Select and aggregate on joint sigma/SWA observations.  Calling the fixed aggregator a
        # second time on raw SWA could silently drop a sigma-selected contact or use a different
        # contact subset in a time bin, defeating the negative-control comparison.
        joint = (
            np.isfinite(sigma_values) & (sigma_values > 0)
            & np.isfinite(swa_values) & (swa_values > 0)
        )
        sigma, sigma_qc = _fixed_power(
            np.where(joint, sigma_values, np.nan), eligible, profile)
        swa, swa_qc = _fixed_power(
            np.where(joint, swa_values, np.nan), eligible, profile)
        if not np.array_equal(
                sigma_qc["selected_mask"], swa_qc["selected_mask"]):
            raise RuntimeError("joint fixed sigma/SWA contact selection diverged")
        sigma_qc["per_second_support"] = sigma_support
        swa_qc["per_second_support"] = swa_support
        return sigma, swa, dict(sigma=sigma_qc, swa=swa_qc)
    joint = (
        np.isfinite(sigma_values) & (sigma_values > 0)
        & np.isfinite(swa_values) & (swa_values > 0)
    )
    sigma, sigma_qc = overlap_connected_aggregate(
        sigma_values, eligible_contacts=eligible,
        observation_mask=joint, transform="log",
        minimum_contacts=profile["power"]["minimum_contacts"])
    # Force SWA onto the exact sigma observation component.
    selected_contacts = sigma_qc["selected_contact_mask"]
    selected_times = sigma_qc["selected_time_mask"]
    component_mask = joint & selected_contacts[:, None] & selected_times[None, :]
    swa, swa_qc = overlap_connected_aggregate(
        swa_values, eligible_contacts=selected_contacts,
        observation_mask=component_mask, transform="log",
        minimum_contacts=profile["power"]["minimum_contacts"])
    sigma_qc["per_second_support"] = sigma_support
    swa_qc["per_second_support"] = swa_support
    return sigma, swa, dict(sigma=sigma_qc, swa=swa_qc)


def _materialize_staging(cache, candidate_contacts, profile):
    dr = _required(cache, "ep_dr_by_contact").astype(float)
    swa = _required(cache, "ep_swa_by_contact").astype(float)
    clean = _required(cache, "ep_clean_fraction_by_contact").astype(float)
    window_mask = _required(
        cache, "ep_valid_welch_window_mask_by_contact").astype(bool)
    window_count = window_mask.sum(axis=2)
    staging = profile["staging"]
    if staging["aggregation"] == "fixed_contact_median":
        ep_dr, ep_swa, ep_clean, details = aggregate_staging_features(
            dr,
            swa,
            clean,
            candidate_contacts,
            min_contacts=staging["minimum_contacts"],
            min_contact_feature_coverage=staging[
                "minimum_contact_feature_coverage"],
            min_contact_fraction_per_epoch=staging[
                "minimum_contact_fraction_per_epoch"],
            valid_window_count_by_contact=window_count,
            min_valid_windows=staging["minimum_valid_welch_windows"],
        )
        details.update(
            fit_required=False,
            support_passes_fit_convergence=True,
            fit_convergence_status="not_applicable_closed_form",
        )
    else:
        ep_dr, ep_swa, ep_clean, details = overlap_connected_staging(
            dr,
            swa,
            clean,
            candidate_contacts,
            valid_window_count_by_contact=window_count,
            minimum_valid_windows=staging["minimum_valid_welch_windows"],
            minimum_contacts=staging["minimum_contacts"],
        )
    details["valid_window_count_by_contact"] = window_count
    return ep_dr, ep_swa, ep_clean, details


def _materialize_auxiliary(cache, modality, profile, n_epochs):
    """Aggregate one RESP EMG/EOG modality from reversible complete-window support."""
    values = _required(cache, f"ep_{modality}_by_channel").astype(float)
    window_mask = _required(
        cache, f"ep_{modality}_valid_window_mask_by_channel").astype(bool)
    if values.ndim != 2 or values.shape[1] != int(n_epochs):
        raise RuntimeError(
            f"{modality.upper()} channel features do not align with staging epochs")
    if window_mask.ndim != 3 or window_mask.shape[:2] != values.shape:
        raise RuntimeError(
            f"{modality.upper()} complete-window support does not align with features")
    window_count = window_mask.sum(axis=2)
    staging = profile["staging"]
    minimum_windows = int(staging["minimum_valid_auxiliary_windows"])
    minimum_coverage = float(
        staging["minimum_auxiliary_channel_coverage"])
    minimum_fraction = float(
        staging["minimum_auxiliary_channel_fraction_per_epoch"])

    if staging["auxiliary_aggregation"] == "fixed_channel_geometric_median":
        aggregate, details = aggregate_auxiliary_epoch_features(
            values,
            min_channel_coverage=minimum_coverage,
            min_channel_fraction_per_epoch=minimum_fraction,
            valid_window_count_by_channel=window_count,
            min_valid_windows=minimum_windows,
        )
    else:
        observation = (
            np.isfinite(values)
            & (values > 0)
            & (window_count >= minimum_windows)
        )
        per_channel_coverage = observation.mean(axis=1)
        eligible = per_channel_coverage >= minimum_coverage
        aggregate, details = overlap_connected_aggregate(
            values,
            eligible_contacts=eligible,
            observation_mask=observation,
            transform="log",
            minimum_contacts=1,
        )
        selected_n = int(details["selected_contact_mask"].sum())
        required = (
            max(1, int(np.ceil(minimum_fraction * selected_n)))
            if selected_n
            else 0
        )
        if required:
            aggregate = np.where(
                np.asarray(details["contact_count"]) >= required,
                aggregate,
                np.nan,
            )
        details.update(
            per_channel_coverage=per_channel_coverage,
            minimum_channel_coverage=minimum_coverage,
            minimum_channel_fraction_per_epoch=minimum_fraction,
            required_channel_count=required,
            n_selected_channels=selected_n,
        )
    details.update(
        modality=modality.upper(),
        minimum_valid_auxiliary_windows=minimum_windows,
        valid_window_count_by_channel=window_count,
    )
    return aggregate, details


def materialize(cache, profile):
    """Return profile-specific arrays and structured support without modifying the cache."""
    profile = validate_qc_profile(profile)
    schema = npz_scalar_text(cache, "cache_schema_version")
    if schema not in {CACHE_SCHEMA_VERSION, LEGACY_CACHE_SCHEMA_V8}:
        raise RuntimeError(
            f"cache schema {schema!r} does not match {CACHE_SCHEMA_VERSION!r}")
    contacts = [str(value) for value in _required(cache, "cortical_chans")]
    n_contacts = len(contacts)
    if "cortical_signal_nonflat_mask" in cache.files:
        nonflat_contacts = np.asarray(
            cache["cortical_signal_nonflat_mask"], bool).ravel()
        activity_qc_provenance = (
            "raw-signal numerical flat-line mask stored by the cache")
    elif schema == LEGACY_CACHE_SCHEMA_V8:
        # v8 predates raw-voltage flat-line metadata.  Preserve its prior behaviour explicitly so
        # pinned results remain inspectable, but never claim that legacy all-True eligibility is a
        # completed activity check.
        nonflat_contacts = np.ones(n_contacts, bool)
        activity_qc_provenance = (
            "legacy v8 compatibility: activity mask unavailable; all contacts provisionally "
            "eligible and numerical flat-line exclusion not applied")
    else:
        raise RuntimeError(
            "current neutral cache lacks the required raw-signal activity mask")
    if len(nonflat_contacts) != n_contacts:
        raise RuntimeError("raw-signal activity mask does not align with neutral contact matrices")
    if "parietal_contact_mask" in cache.files:
        parietal = np.asarray(cache["parietal_contact_mask"], bool) & nonflat_contacts
        frontal = np.asarray(cache["frontal_contact_mask"], bool) & nonflat_contacts
        cohort = "RESPect"
    else:
        parietal = nonflat_contacts.copy()
        frontal = nonflat_contacts.copy()
        cohort = "HUP"
    if len(parietal) != n_contacts or len(frontal) != n_contacts:
        raise RuntimeError("ROI masks do not align with neutral contact matrices")

    prepared_power = _power_matrices(cache, profile)
    sigma_global, swa_global, global_qc = _power_pair(
        cache, nonflat_contacts, profile, prepared=prepared_power)
    if np.array_equal(parietal, nonflat_contacts):
        sigma_parietal = sigma_global
        swa_parietal = swa_global
        parietal_qc = global_qc
    else:
        sigma_parietal, swa_parietal, parietal_qc = _power_pair(
            cache, parietal, profile, prepared=prepared_power)
    # Preserve the historical fixed-set reference exactly, but do not make a staging feature
    # depend on whether that contact happened to pass an unrelated sigma-power coverage gate.
    # The overlap-connected profile starts from every anatomy-qualified cortical contact and
    # resolves support on the staging endpoint's own observation graph.
    candidate_contacts = (
        np.asarray(global_qc["sigma"]["selected_mask"], bool)
        if profile["power"]["aggregation"] == "fixed_contact_mean"
        else nonflat_contacts.copy()
    )
    ep_dr, ep_swa, ep_clean, staging_qc = _materialize_staging(
        cache, candidate_contacts, profile)

    if cohort == "HUP":
        stage_lab, nrem, separation = stage_epochs(
            dict(dr=ep_dr, swa=ep_swa, clean=ep_clean))
        stage_proxy = stage_lab.copy()
        stage_proxy_sensitivity = stage_lab.copy()
        stage_diagnostics = {"gmm_separation": separation}
        ep_emg = np.full(len(ep_dr), np.nan)
        ep_eog = np.full(len(ep_dr), np.nan)
        auxiliary_qc = {
            "status": "not_applicable",
            "reason": "HUP staging proxy does not use RESP EMG/EOG channels",
        }
    else:
        ep_emg, emg_qc = _materialize_auxiliary(
            cache, "emg", profile, len(ep_dr))
        ep_eog, eog_qc = _materialize_auxiliary(
            cache, "eog", profile, len(ep_dr))
        auxiliary_qc = {"status": "ok", "emg": emg_qc, "eog": eog_qc}
        proxy, stage_diagnostics = score_stages(
            ep_swa, ep_emg, ep_eog, ep_dr, ep_clean)
        annotation = np.asarray(cache["stage_lab_annotation"]).astype(str)
        sources = np.asarray(cache["stage_source"]).astype(str)
        stage_lab, _ = constrain_proxy_to_annotations(
            proxy, annotation, sources, allow_proxy_unknown_sleep=False)
        stage_proxy_sensitivity, _ = constrain_proxy_to_annotations(
            proxy, annotation, sources, allow_proxy_unknown_sleep=True)
        stage_proxy = proxy
        nrem = np.isin(stage_lab, ("NREM", "N2", "N3"))
        separation = None

    hr = np.asarray(cache["hr_1"], float)
    rr = np.asarray(cache["rr_4"], float)
    hr_coverage = float(np.isfinite(rr).mean())
    power_gate = profile["power"]["minimum_aggregate_coverage"]
    parietal_coverage = float(np.isfinite(sigma_parietal).mean())
    return dict(
        subject=npz_scalar_text(cache, "subject"),
        cohort=cohort,
        contacts=np.asarray(contacts, dtype="<U96"),
        cortical_signal_nonflat_mask=nonflat_contacts,
        cortical_signal_activity_qc_provenance=activity_qc_provenance,
        parietal_contact_mask=parietal,
        frontal_contact_mask=frontal,
        sigma_global=sigma_global,
        swa_global=swa_global,
        sigma_parietal=sigma_parietal,
        swa_parietal=swa_parietal,
        hr_1=hr,
        rr_4=rr,
        ep_dr=ep_dr,
        ep_swa=ep_swa,
        ep_clean=ep_clean,
        ep_emg=ep_emg,
        ep_eog=ep_eog,
        stage_lab=np.asarray(stage_lab, dtype="<U5"),
        stage_lab_proxy=np.asarray(stage_proxy, dtype="<U5"),
        stage_lab_proxy_sensitivity=np.asarray(
            stage_proxy_sensitivity, dtype="<U5"),
        nrem=np.asarray(nrem, bool),
        gmm_separation=separation,
        global_power_qc=global_qc,
        parietal_power_qc=parietal_qc,
        staging_qc=staging_qc,
        auxiliary_qc=auxiliary_qc,
        stage_diagnostics=stage_diagnostics,
        hr_coverage=hr_coverage,
        hr_meets_profile=bool(
            hr_coverage >= profile["hr"]["minimum_coverage"]),
        parietal_power_coverage=parietal_coverage,
        parietal_power_meets_profile=bool(parietal_coverage >= power_gate),
        qc_profile=profile,
        qc_manifest=manifest_qc_config(profile),
    )


def _json_safe(value):
    """Return strict-JSON data: no NumPy objects and no bare NaN/Infinity tokens."""
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cache")
    parser.add_argument("--profile", default="audit80")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    profile = load_qc_profile(args.profile)
    with np.load(args.cache, allow_pickle=False) as cache:
        result = materialize(cache, profile)
    # Materialized caches are local analysis intermediates, not publication outputs.
    payload = {
        key: value
        for key, value in result.items()
        if key not in {
            "global_power_qc",
            "parietal_power_qc",
            "staging_qc",
            "auxiliary_qc",
            "stage_diagnostics",
            "qc_profile",
            "qc_manifest",
        }
    }
    if payload.get("gmm_separation") is None:
        payload["gmm_separation"] = np.asarray(np.nan)
    payload["support_json"] = json.dumps(
        _json_safe({
            "global_power_qc": result["global_power_qc"],
            "parietal_power_qc": result["parietal_power_qc"],
            "staging_qc": result["staging_qc"],
            "auxiliary_qc": result["auxiliary_qc"],
            "stage_diagnostics": result["stage_diagnostics"],
            **result["qc_manifest"],
        }),
        sort_keys=True,
        allow_nan=False,
    )
    atomic_savez(os.path.abspath(args.output), **payload)


if __name__ == "__main__":
    main()
