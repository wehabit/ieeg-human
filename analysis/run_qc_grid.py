"""Run locked, endpoint-local 3A/3B QC sensitivity profiles from neutral caches.

The grid is diagnostic: no profile is promoted because it produces an LC-like result.  Scientific
method constants (Lecci's 120-s bouts and Naji's 180-s stable stages) stay fixed while one QC axis
changes at a time.  The HUP neutral cache also permits offline reconstruction of the current hybrid
3D event detector and pairing given the common profile-specific stage labels. Endpoint-specific
3D acquisition/event-valid gates filter effect contacts after staging; they do not redefine sleep
stages. Phase vectors remain descriptive because a pairing-aware null and validated HUP sleep
stages are unavailable.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os

import numpy as np

import lecci_faithful_3A as lecci
from event_3B_cached import stable_stage_epoch_indices, stage_so_times
from event_3B_mednick import subject_so_triggered, rr_baseline_hr, FS_RR
from event_3d_cache_support import analyse_3d_cache_support
from materialize_qc_cache import materialize
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    atomic_json_dump,
    cache_code_sha256,
    file_sha256,
    npz_scalar_text,
    runtime_versions,
    source_tree_sha256,
    utc_now,
)
from qc_profiles import (
    expand_qc_grid,
    manifest_qc_config,
    validated_staging_calibration,
)
from spectral_gapped import coherence_gapped, analytic_msc_threshold


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_OUTPUT = os.path.join(ROOT, "outputs", "qc_grid")
GRID_SURROGATES_3B = 199


def _seed(*parts):
    digest = hashlib.sha256("|".join(map(str, parts)).encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "little")


def _cache_key(*values):
    """Deterministic in-process memoization key for profile subsections."""
    return json.dumps(
        values, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _selected_count(details):
    if "n_selected" in details:
        return int(details["n_selected"])
    if "selected_contact_mask" in details:
        return int(np.asarray(details["selected_contact_mask"], bool).sum())
    return 0


def materialization_diagnostics(materialized):
    """Compact profile-specific staging/auxiliary support diagnostics.

    Endpoint results can remain unchanged when an auxiliary threshold changes because RESPect
    author annotations override the affected proxy labels, and HUP has no EMG/EOG inputs.  Keeping
    these intermediate counts in every profile makes that invariance auditable instead of making a
    correctly wired sensitivity axis look inactive.
    """
    labels = np.asarray(materialized["stage_lab"]).astype(str)
    proxy = np.asarray(materialized["stage_lab_proxy"]).astype(str)

    def counts(values):
        return {
            label: int(np.sum(values == label))
            for label in ("W", "R", "NREM", "N2", "N3", "")
            if np.any(values == label)
        }

    auxiliary_qc = materialized.get("auxiliary_qc", {})
    return dict(
        auxiliary_status=str(auxiliary_qc.get("status", "unknown")),
        auxiliary_reason=auxiliary_qc.get("reason"),
        n_finite_emg_epochs=int(
            np.isfinite(np.asarray(materialized["ep_emg"], float)).sum()),
        n_finite_eog_epochs=int(
            np.isfinite(np.asarray(materialized["ep_eog"], float)).sum()),
        final_stage_label_counts=counts(labels),
        proxy_stage_label_counts=counts(proxy),
        n_final_labels_different_from_proxy=int(np.sum(labels != proxy)),
        minimum_valid_auxiliary_windows=int(
            materialized["qc_profile"]["staging"][
                "minimum_valid_auxiliary_windows"]),
        minimum_valid_welch_windows=int(
            materialized["qc_profile"]["staging"][
                "minimum_valid_welch_windows"]),
    )


def analyse_3a(materialized, profile):
    sig = np.asarray(materialized["sigma_parietal"], float)
    swa = np.asarray(materialized["swa_parietal"], float)
    hr = np.asarray(materialized["hr_1"], float)
    lab = np.asarray(materialized["stage_lab"]).astype(str)
    nrem, analysis_window = lecci.core_study_nrem_mask(lab)
    power_qc = materialized["parietal_power_qc"]["sigma"]
    power_fit_converged = bool(
        power_qc.get("support_passes_fit_convergence", True))
    n_contacts = _selected_count(power_qc)
    aggregate_coverage = float(np.isfinite(sig).mean())
    min_nrem = int(profile["endpoint_3a"]["minimum_nrem_epochs"])
    n_nrem = int(nrem.sum())
    support_pass = bool(
        power_fit_converged
        and n_contacts >= int(profile["power"]["minimum_contacts"])
        and aggregate_coverage >= float(
            profile["power"]["minimum_aggregate_coverage"])
        and n_nrem >= min_nrem
    )

    freqs = None
    spectrum = None
    spectrum_swa = None
    n_bouts = 0
    bout_seconds = 0.0
    if support_pass:
        freqs, spectrum, n_bouts, bout_seconds = lecci.subject_spectrum(
            sig, nrem)
        _, spectrum_swa, _, _ = lecci.subject_spectrum(swa, nrem)
    peak = (
        lecci.fit_peak(freqs, spectrum)
        if spectrum is not None
        else {
            "peak_hz": None,
            "accepted": False,
            "method": (
                "no spectrum"
                if support_pass
                else "withheld: support below locked profile"
            ),
        }
    )
    spectrum_available = spectrum is not None
    negative_control = None
    if (
        spectrum is not None
        and spectrum_swa is not None
        and peak.get("peak_hz") is not None
        and np.isfinite(peak.get("spectral_sd_hz", np.nan))
    ):
        half_width = 0.5 * float(peak["spectral_sd_hz"])
        window = np.abs(freqs - float(peak["peak_hz"])) <= half_width
        if window.any():
            negative_control = dict(
                centre_hz=float(peak["peak_hz"]),
                half_width_hz=half_width,
                sigma_window_mean=float(np.mean(spectrum[window])),
                swa_same_window_mean=float(np.mean(spectrum_swa[window])),
            )

    cardiac_pass = bool(materialized["hr_meets_profile"])
    coherence = None
    xcorr = None
    if support_pass and cardiac_pass:
        second_mask = np.zeros(len(sig), bool)
        for start, stop in lecci.nrem_bouts(nrem):
            second_mask[start:stop] = True
        co = coherence_gapped(
            np.where(second_mask, hr, np.nan),
            np.where(second_mask, sig, np.nan),
            fs=lecci.FS,
            nperseg=lecci.NPERSEG,
            highpass=0.005,
        )
        if co is not None:
            index = int(np.argmin(np.abs(co["f"] - lecci.F_LECCI)))
            value = float(co["cxy"][index])
            if np.isfinite(value):
                threshold = float(analytic_msc_threshold(co["K"], lecci.ALPHA))
                coherence = dict(
                    at_0p02_hz=value,
                    analytic_threshold=threshold,
                    above_analytic_threshold=bool(value > threshold),
                    K=int(co["K"]),
                    n_valid=int(co["n_valid"]),
                )
        xcorr = lecci.cross_correlation(
            sig,
            hr,
            nrem,
            minimum_windows=profile["endpoint_3a"][
                "minimum_cross_correlation_windows"],
        )
        if xcorr is not None:
            xcorr = {
                key: value
                for key, value in xcorr.items()
                if key not in ("lag_s", "xcorr")
            }
    reasons = []
    if not power_fit_converged:
        reasons.append("overlap-connected power fit did not converge")
    if n_contacts < int(profile["power"]["minimum_contacts"]):
        reasons.append(
            f"{n_contacts} selected contacts < "
            f"{profile['power']['minimum_contacts']}")
    if aggregate_coverage < profile["power"]["minimum_aggregate_coverage"]:
        reasons.append(
            f"aggregate coverage {aggregate_coverage:.3f} < "
            f"{profile['power']['minimum_aggregate_coverage']:.3f}")
    if n_nrem < min_nrem:
        reasons.append(f"{n_nrem} NREM epochs < {min_nrem}")
    if not cardiac_pass:
        reasons.append(
            f"RR coverage {materialized['hr_coverage']:.3f} < "
            f"{profile['hr']['minimum_coverage']:.3f}; cardiac endpoints only")
    return dict(
        n_selected_contacts=n_contacts,
        selected_contact_ids=[
            str(contact)
            for contact, keep in zip(
                materialized["contacts"],
                np.asarray(
                    power_qc.get(
                        "selected_mask",
                        power_qc.get("selected_contact_mask", [])),
                    bool,
                ),
            )
            if keep
        ],
        aggregate_coverage=aggregate_coverage,
        hr_coverage=float(materialized["hr_coverage"]),
        n_nrem_epochs=n_nrem,
        analysis_window=analysis_window,
        n_bouts=int(n_bouts),
        bout_seconds=float(bout_seconds),
        support_passes_profile=support_pass,
        power_fit_converged=power_fit_converged,
        spectrum_computed=bool(spectrum_available),
        spectrum_available_under_profile=bool(support_pass and spectrum_available),
        peak=_jsonable(peak),
        negative_control=negative_control,
        coherence=coherence,
        cross_correlation=xcorr,
        endpoint_availability=dict(
            spectrum=bool(support_pass and spectrum_available),
            fixed_0p02_coherence=coherence is not None,
            cross_correlation=xcorr is not None,
        ),
        support_reasons=reasons,
    )


def _stage_pool(rr, stable_epochs):
    stage = np.zeros(len(rr), bool)
    for epoch in stable_epochs:
        start = int(epoch * lecci.EPOCH * FS_RR)
        stop = int((epoch + 1) * lecci.EPOCH * FS_RR)
        stage[start:stop] = True
    return np.where(stage & np.isfinite(rr))[0]


def analyse_3b(cache, materialized, profile):
    rr = np.asarray(materialized["rr_4"], float)
    lab = np.asarray(materialized["stage_lab"]).astype(str)
    contacts = [str(value) for value in materialized["contacts"]]
    roi_mask = np.asarray(
        materialized["frontal_contact_mask"]
        if materialized["cohort"] == "RESPect"
        else np.ones(len(contacts), bool),
        bool,
    )
    eligible_contacts = [
        contact for contact, keep in zip(contacts, roi_mask) if keep
    ]
    result = dict(
        endpoint_contact_selection=(
            "anatomy/ROI plus stage-local SO and complete RR-window support; "
            "does not inherit the sigma full-night coverage mask"),
        eligible_contact_ids=eligible_contacts,
        hr_coverage=float(materialized["hr_coverage"]),
        hr_meets_profile=bool(materialized["hr_meets_profile"]),
        n_finite_rr_samples=int(np.isfinite(rr).sum()),
        staging_fit_converged=bool(
            materialized.get("staging_qc", {}).get(
                "support_passes_fit_convergence", True)),
        stages={},
    )
    stages = ("N2", "N3", "NREM")
    for stage_name in stages:
        if stage_name == "NREM":
            # Exploratory pooling means the union of NREM/N2/N3 and preserves continuity across
            # an N2-to-N3 transition.  Matching only the literal "NREM" label silently excluded
            # author-scored N2/N3 epochs from an endpoint named pooled NREM.
            pooled_labels = np.where(
                np.isin(lab, ("NREM", "N2", "N3")), "NREM", "")
            stable = stable_stage_epoch_indices(pooled_labels, "NREM")
            raw_epochs = int(np.isin(lab, ("NREM", "N2", "N3")).sum())
        else:
            stable = stable_stage_epoch_indices(lab, stage_name)
            raw_epochs = int((lab == stage_name).sum())
        stage_result = dict(
            raw_epochs=raw_epochs,
            stable_epochs=int(len(stable)),
            stable_seconds=int(len(stable) * lecci.EPOCH),
            exploratory=(stage_name == "NREM"),
        )
        reasons = []
        if not result["staging_fit_converged"]:
            reasons.append("overlap-connected staging fit did not converge")
        if not materialized["hr_meets_profile"]:
            reasons.append("cardiac coverage below profile")
        if np.isfinite(rr).sum() < profile["endpoint_3b"]["minimum_finite_rr_samples"]:
            reasons.append(
                f"{int(np.isfinite(rr).sum())} finite RR samples < "
                f"{profile['endpoint_3b']['minimum_finite_rr_samples']}")
        if len(stable) < 6:
            reasons.append("no uninterrupted 180-s stage run")
        pool = _stage_pool(rr, stable)
        if len(pool) < profile["endpoint_3b"]["minimum_finite_stage_samples"]:
            reasons.append(
                f"{len(pool)} finite stage samples < "
                f"{profile['endpoint_3b']['minimum_finite_stage_samples']}")
        keep_epochs = set(stable.tolist())
        event_counts = {}
        troughs = []
        event_contacts = []
        amplitude_percentile = float(
            profile["endpoint_3b"]["so_amplitude_percentile"])
        for contact in eligible_contacts:
            times = stage_so_times(
                cache, contact, keep_epochs,
                percentile=amplitude_percentile)
            event_counts[contact] = int(len(times))
            if len(times) >= profile["endpoint_3b"]["minimum_so_per_contact"]:
                troughs.append(times)
                event_contacts.append(contact)
        stage_result["candidate_so_counts_by_contact"] = event_counts
        stage_result["so_amplitude_percentile"] = amplitude_percentile
        stage_result["so_amplitude_rule_provenance"] = (
            "intracranial adaptation/sensitivity choice; Naji cites fixed "
            "Dang-Vu scalp-amplitude criteria and does not specify this percentile")
        stage_result["contacts_meeting_event_count"] = event_contacts
        if len(event_contacts) < profile["endpoint_3b"]["minimum_contacts"]:
            reasons.append(
                f"{len(event_contacts)} event-supported contacts < "
                f"{profile['endpoint_3b']['minimum_contacts']}")
        estimate = None
        if not reasons:
            baseline = rr_baseline_hr(rr[pool])
            estimate = subject_so_triggered(
                rr,
                troughs,
                baseline,
                pool,
                n_sur=GRID_SURROGATES_3B,
                rng=np.random.RandomState(
                    _seed(materialized["subject"], stage_name, "3B-grid")),
                domain="rr",
                minimum_channels=profile["endpoint_3b"]["minimum_contacts"],
                channel_ids=event_contacts,
                minimum_events_per_channel=profile["endpoint_3b"][
                    "minimum_so_per_contact"],
                minimum_surrogate_pool_samples=profile["endpoint_3b"][
                    "minimum_finite_stage_samples"],
            )
            if estimate is None:
                reasons.append(
                    "complete in-stage RR windows leave insufficient event/contact support")
            else:
                estimate = {
                    key: value
                    for key, value in estimate.items()
                    if key not in ("curve", "rr_curve", "lag_s")
                }
                estimate["stage_mean_hr"] = float(baseline)
        stage_result["estimate"] = _jsonable(estimate)
        stage_result["available_under_profile"] = estimate is not None
        stage_result["support_reasons"] = reasons
        result["stages"][stage_name] = stage_result
    return result


def _cache_manifest(cache_dir):
    path = os.path.join(cache_dir, "RUN_MANIFEST.json")
    with open(path) as handle:
        manifest = json.load(handle)
    if manifest.get("run_state") != "complete":
        raise RuntimeError(f"{path} is not a completed cache run")
    expected_pipeline = (
        "stage_ds003848"
        if os.path.basename(os.path.abspath(cache_dir)) == "ds003848"
        else "cache_lc_series"
    )
    if manifest.get("pipeline") != expected_pipeline:
        raise RuntimeError(
            f"{path} pipeline {manifest.get('pipeline')!r} is not "
            f"{expected_pipeline!r}")
    if manifest.get("analysis_version") != ANALYSIS_VERSION:
        raise RuntimeError(
            f"{path} analysis version differs from {ANALYSIS_VERSION}")
    if manifest.get("cache_schema_version") != CACHE_SCHEMA_VERSION:
        raise RuntimeError(
            f"{path} cache schema differs from {CACHE_SCHEMA_VERSION}")
    if manifest.get("config", {}).get("cache_code_sha256") != cache_code_sha256(ROOT):
        raise RuntimeError(f"{path} cache builder digest is stale")
    if manifest.get("failed"):
        raise RuntimeError(f"{path} contains failed participants")
    if manifest.get("runtime_versions") != runtime_versions():
        raise RuntimeError(
            f"{path} runtime versions differ from the current analysis runtime")
    requested = list(manifest.get("requested", []))
    completed = list(manifest.get("completed", []))
    skipped_entries = list(manifest.get("skipped", []))
    if any(
        not isinstance(value, dict)
        or not value.get("subject")
        or not value.get("reason")
        for value in skipped_entries
    ):
        raise RuntimeError(f"{path} has malformed skipped-subject records")
    skipped = [value["subject"] for value in skipped_entries]
    for label, values in (
            ("requested", requested), ("completed", completed), ("skipped", skipped)):
        if any(not isinstance(value, str) or not value for value in values):
            raise RuntimeError(f"{path} has invalid {label} subject identifiers")
        if len(values) != len(set(values)):
            raise RuntimeError(f"{path} has duplicate {label} subject identifiers")
    if set(completed) & set(skipped):
        raise RuntimeError(f"{path} lists a subject as both completed and skipped")
    if set(requested) != set(completed) | set(skipped):
        raise RuntimeError(
            f"{path} requested subjects are not exactly completed or explicitly skipped")

    expected_hashes = manifest.get("result_files_sha256")
    if (
        not isinstance(expected_hashes, dict)
        or set(expected_hashes) != set(requested)
    ):
        raise RuntimeError(
            f"{path} lacks exact hashes for all terminal subject caches")
    current_digest = cache_code_sha256(ROOT)
    verified_hashes = {}
    for subject in completed + skipped:
        cache_path = os.path.join(cache_dir, f"{subject}.npz")
        if not os.path.exists(cache_path):
            raise RuntimeError(f"{path} references missing cache {cache_path}")
        actual_hash = file_sha256(cache_path)
        if expected_hashes.get(subject) != actual_hash:
            raise RuntimeError(f"{cache_path} does not match its manifest SHA-256")
        with np.load(cache_path, allow_pickle=False) as cache:
            expected_status = "ok" if subject in completed else "skip"
            if npz_scalar_text(cache, "subject") != subject:
                raise RuntimeError(f"{cache_path} embeds a different subject")
            if npz_scalar_text(cache, "status") != expected_status:
                raise RuntimeError(
                    f"{cache_path} status does not match the manifest")
            if npz_scalar_text(cache, "cache_schema_version") != CACHE_SCHEMA_VERSION:
                raise RuntimeError(f"{cache_path} has a stale cache schema")
            if npz_scalar_text(cache, "cache_code_sha256") != current_digest:
                raise RuntimeError(f"{cache_path} has a stale cache builder digest")
        verified_hashes[subject] = actual_hash
    return manifest, path, verified_hashes


def run_grid(cache_dir, grid_id):
    manifest, manifest_path, cache_hashes = _cache_manifest(cache_dir)
    calibration = validated_staging_calibration()
    if (
        calibration["analysis_version"] != ANALYSIS_VERSION
        or calibration["cache_schema_version"] != CACHE_SCHEMA_VERSION
    ):
        raise RuntimeError(
            "pinned staging calibration analysis/cache version is stale")
    if manifest["pipeline"] == calibration["cache_pipeline"]:
        if (
            manifest["run_id"] != calibration["cache_run_id"]
            or file_sha256(manifest_path)
                != calibration["cache_manifest_sha256"]
            or cache_hashes != calibration["cache_files_sha256"]
        ):
            raise RuntimeError(
                "RESPect cache bytes no longer match the pinned staging calibration")
        calibration_application = (
            "same RESPect cache cohort used for outcome-blind calibration")
    else:
        calibration_application = calibration["hup_transport_status"]
    subjects = list(manifest.get("completed", []))
    profiles = expand_qc_grid(grid_id)
    subject_results_by_profile = [[] for _ in profiles]

    # A count/event grid changes many endpoint thresholds while leaving the
    # expensive neutral materialization identical.  Process one subject at a
    # time and memoize only by the profile subsections that can change each
    # intermediate/result.  This is computational reuse, not profile
    # selection: every locked profile still receives an explicit result.
    for subject in subjects:
        path = os.path.join(cache_dir, f"{subject}.npz")
        with np.load(path, allow_pickle=False) as cache:
            materialized_cache = {}
            result_3a_cache = {}
            result_3b_cache = {}
            result_3d_cache = {}
            for profile_index, profile in enumerate(profiles):
                materialized_key = _cache_key(
                    profile["hr"], profile["power"], profile["staging"])
                if materialized_key not in materialized_cache:
                    materialized_cache[materialized_key] = materialize(
                        cache, profile)
                materialized = materialized_cache[materialized_key]

                key_3a = _cache_key(
                    materialized_key, profile["endpoint_3a"])
                if key_3a not in result_3a_cache:
                    result_3a_cache[key_3a] = analyse_3a(
                        materialized, profile)

                stage_key = hashlib.sha256(
                    np.asarray(
                        materialized["stage_lab"], dtype="<U5"
                    ).tobytes()
                ).hexdigest()
                key_3b = _cache_key(
                    stage_key,
                    bool(materialized["hr_meets_profile"]),
                    profile["endpoint_3b"],
                )
                if key_3b not in result_3b_cache:
                    result_3b_cache[key_3b] = analyse_3b(
                        cache, materialized, profile)
                key_3d = _cache_key(stage_key, profile["endpoint_3d"])
                if key_3d not in result_3d_cache:
                    result_3d_cache[key_3d] = analyse_3d_cache_support(
                        cache, materialized, profile)

                subject_results_by_profile[profile_index].append(dict(
                    subject=subject,
                    support_diagnostics=materialization_diagnostics(
                        materialized),
                    result_3a=result_3a_cache[key_3a],
                    result_3b=result_3b_cache[key_3b],
                    result_3d=result_3d_cache[key_3d],
                ))

    profile_results = []
    for profile, subject_results in zip(
            profiles, subject_results_by_profile):
        counts = dict(
            n_subjects=len(subject_results),
            n_3a_spectrum=sum(
                value["result_3a"]["endpoint_availability"]["spectrum"]
                for value in subject_results),
            n_3a_coherence=sum(
                value["result_3a"]["endpoint_availability"][
                    "fixed_0p02_coherence"]
                for value in subject_results),
            n_3a_cross_correlation=sum(
                value["result_3a"]["endpoint_availability"][
                    "cross_correlation"]
                for value in subject_results),
            n_3b_N2=sum(
                value["result_3b"]["stages"]["N2"]["available_under_profile"]
                for value in subject_results),
            n_3b_N3=sum(
                value["result_3b"]["stages"]["N3"]["available_under_profile"]
                for value in subject_results),
            n_3b_pooled_NREM_exploratory=sum(
                value["result_3b"]["stages"]["NREM"][
                    "available_under_profile"]
                for value in subject_results),
            n_3d_descriptive_effect=sum(
                value["result_3d"]["descriptive_effect_available"]
                for value in subject_results),
            n_3d_support_evaluable=sum(
                value["result_3d"]["support_passes_profile"] is not None
                for value in subject_results),
        )
        minimum_participants = int(profile["cohort"]["minimum_participants"])
        endpoint_counts = {
            key: value for key, value in counts.items() if key.startswith("n_3")
        }
        counts["minimum_participants"] = minimum_participants
        counts["meets_minimum_participants"] = {
            key.removeprefix("n_"): bool(value >= minimum_participants)
            for key, value in endpoint_counts.items()
        }
        profile_results.append(dict(
            **manifest_qc_config(profile),
            subjects=subject_results,
            summary=counts,
        ))
    return dict(
        analysis_version=ANALYSIS_VERSION,
        cache_schema_version=CACHE_SCHEMA_VERSION,
        generated_at_utc=utc_now(),
        source_tree_sha256=source_tree_sha256(ROOT),
        grid_id=grid_id,
        grid_role="post-audit sensitivity; no result-based profile selection",
        cache_directory_relative=os.path.relpath(
            cache_dir, ROOT).replace(os.sep, "/"),
        cache_manifest_run_id=manifest["run_id"],
        cache_manifest_sha256=file_sha256(manifest_path),
        cache_pipeline=manifest["pipeline"],
        cache_runtime_versions=manifest["runtime_versions"],
        cache_files_sha256=cache_hashes,
        staging_calibration_provenance={
            key: value
            for key, value in calibration.items()
            if key != "artifact_path"
        },
        staging_calibration_application_to_this_cache=(
            calibration_application),
        runtime_versions=runtime_versions(),
        requested_subjects=list(manifest.get("requested", [])),
        analysed_subjects=subjects,
        explicitly_skipped_subjects=list(manifest.get("skipped", [])),
        profiles=profile_results,
        limitations=[
            "3B NREM is an explicitly exploratory pooled-state endpoint, not Naji's N2/SWS contrast",
            "HUP proxy stages are not PSG/expert validated",
            "3D is a Staresina/Helfrich hybrid adaptation: it first detects a duration-qualified "
            "spindle and then retains at most one spindle within +/-2 s of each SO",
            "RESPect is outside the specified 3D cohort; only HUP neutral caches contain the "
            "pre-threshold event fields",
            "all 3D inference remains disabled; no 3D p value is produced",
            "3A uses the first 210 minutes from the first available classified sleep epoch; "
            "expert S1-followed-by-S2 sleep onset is unavailable",
            "the 3B iEEG amplitude percentile is an explicit adaptation/sensitivity choice, "
            "not Naji's fixed scalp-voltage criterion",
            "profile grids are retrospective sensitivity analyses, not prospective preregistration",
        ],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument(
        "--grid",
        choices=(
            "coverage_oat_v1",
            "staging_window_support_v1",
            "auxiliary_window_support_v1",
            "event_count_oat_v1",
        ),
        default="coverage_oat_v1",
    )
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    cache_dir = os.path.abspath(args.cache_dir)
    result = run_grid(cache_dir, args.grid)
    cohort = (
        "respect"
        if os.path.basename(cache_dir) == "ds003848"
        else "hup"
    )
    output = os.path.join(
        os.path.abspath(args.output_dir),
        args.grid,
        f"{cohort}_qc_grid.json",
    )
    atomic_json_dump(result, output)
    print(
        f"wrote {output}: {len(result['profiles'])} profiles, "
        f"{len(result['requested_subjects'])} subjects",
        flush=True,
    )


if __name__ == "__main__":
    main()
