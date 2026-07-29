"""Compare simultaneous scalp EEG with the frozen iEEG LC-proxy adaptations.

The comparison holds participant, portal snapshot, seven-hour interval, cached
ECG/RR, and profile-materialized sleep labels fixed.  For 3A the primary scalp
sensor is C3/C03.  For 3B, unilateral F3 and Fz are explicitly exploratory:
Naji derived a cardiac curve for each referenced F3 and F4 electrode and then
averaged the electrode-specific HR-maximum/RR-minimum times, whereas the portal
reference is undocumented and Fz was not a Naji sensor. F4 is in the HUP138
acquisition plan for the separate scalp-only bilateral exploratory module, but
is not analyzed by this paired scalp–iEEG pipeline.

Before accepting any paired result this script:

* verifies both cache files byte-for-byte against terminal manifests;
* requires the current cache schema/source digest and stored raw-signal
  activity mask; no legacy cache is patched in memory;
* recomputes the original iEEG 3A and 3B result objects and requires an exact
  match to the locked ``overlap11_endpoint_local`` profile snapshot;
* computes spectra, peaks, and SWA controls on common finite-positive EEG
  support, then computes coherence and cross-correlation on the stricter
  subset that also has finite shared HR;
* hashes both endpoint-local support masks plus the shared RR, HR, and
  stage-label arrays used on both modality arms.

Outputs are strict JSON (no NaN/Infinity), a normalized subject-stage CSV, an
explicit role-pair CSV, and PNG/SVG figures.

Usage
-----
    .venv/bin/python analysis/paired_scalp_ieeg_comparison.py
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os

import numpy as np

from artifact_contracts import (
    PAIRED_RESULT_SCHEMA,
    SCALP_INVENTORY_SCHEMA,
)
from cache_lc_series import power_from_binned_support
from cache_paired_scalp import (
    DEFAULT_IEEG_CACHE,
    DEFAULT_OUTPUT as DEFAULT_SCALP_CACHE,
    SCALP_CHANNEL_PLAN,
    validate_pinned_ieeg_cache,
    validate_scalp_sidecar_payload,
    validate_terminal_sidecar_manifest,
)
from event_3B_cached import stable_stage_epoch_indices, stage_so_times
from event_3b_estimators import FS_RR, rr_baseline_hr, subject_so_triggered
from materialize_qc_cache import materialize
from overlap_aggregate import overlap_connected_aggregate
from paired_artifact_validation import _validate_classifications
from paired_reporting import (
    EXACT_SIGNFLIP_MAX_PAIRS,
    MONTE_CARLO_SIGNFLIP_DRAWS,
    STAGES_3B,
    TARGET_COMPATIBLE_PEAK_BAND_HZ,
    _group_summary,
    _ieeg_3b_csv_row,
    _joined_support_reasons,
    _make_3a_figure,
    _make_3b_figure,
    _metric,
    _negative_control_ratio,
    _normalized_csv_rows,
    _pair_summary,
    _paired_coherence_ratio_plot,
    _paired_plot,
    _role_pair_csv_rows,
    _write_csv,
)
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    atomic_json_dump,
    file_sha256,
    npz_scalar_text,
    runtime_versions,
    utc_now,
)
from release_provenance import (
    require_clean_release_provenance,
    require_release_source_unchanged,
    validate_recorded_release_provenance,
)
from qc_profiles import (
    load_qc_profile,
    profile_file_sha256,
    qc_profile_sha256,
)
from run_qc_grid import analyse_3a, analyse_3b
import lecci_faithful_3A as lecci


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_QC_GRID = os.path.join(
    ROOT, "outputs", "qc_grid_public", "locked",
    "overlap11_endpoint_local__hup.json")
DEFAULT_RESULTS = os.path.join(ROOT, "outputs", "paired_scalp_ieeg")
DEFAULT_SCALP_INVENTORY = os.path.join(
    ROOT, "outputs", "paired_scalp_inventory",
    "hup_scalp_channel_inventory.json")
PROFILE_ID = "overlap11_endpoint_local"
PIPELINE = "paired_scalp_ieeg_comparison"
RESULT_SCHEMA = PAIRED_RESULT_SCHEMA
N_SURROGATES_3B = 199
POLARITY_SENSITIVITY_STATUS = (
    "not_computable_from_outcome_neutral_sidecar"
)
NAJI_EXACT_BLOCKING_REASONS = (
    (
        "F4 is not analyzed and electrode-specific F3/F4 cardiac curves are "
        "not combined as in Naji"
    ),
    (
        "online reference is undocumented, so F3/F4 labels do not establish "
        "Naji's F3/A2 and F4/A1 derivations"
    ),
    (
        "sleep stages are held-fixed iEEG-derived proxy labels, not "
        "independently visually scored scalp PSG bins"
    ),
)


def _json_safe(value):
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


def _canonical_json(value):
    return json.dumps(
        _json_safe(value),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _array_sha256(value):
    array = np.ascontiguousarray(np.asarray(value))
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(b"\0")
    digest.update(str(array.shape).encode("ascii"))
    digest.update(b"\0")
    digest.update(array.tobytes())
    return digest.hexdigest()


def _validate_scalp_inventory(
        path,
        selected_subjects=None,
        *,
        qc_grid_path=DEFAULT_QC_GRID,
        ieeg_cache_dir=DEFAULT_IEEG_CACHE,
        scalp_cache_dir=DEFAULT_SCALP_CACHE):
    """Validate the terminal all-cohort inventory and selected intersection."""
    path = os.path.abspath(path)
    manifest_path = os.path.join(os.path.dirname(path), "RUN_MANIFEST.json")
    if not os.path.isfile(path) or not os.path.isfile(manifest_path):
        raise FileNotFoundError(
            "paired comparison requires the terminal all-HUP scalp inventory")
    with open(path) as handle:
        payload = json.load(handle)
    with open(manifest_path) as handle:
        manifest = json.load(handle)
    for value, label in ((payload, "inventory"), (manifest, "manifest")):
        if (
            value.get("schema_version") != SCALP_INVENTORY_SCHEMA
            or value.get("pipeline") != "audit_hup_scalp_inventory"
            or value.get("run_state") != "complete"
            or value.get("analysis_version") != ANALYSIS_VERSION
            or value.get("cache_schema_version") != CACHE_SCHEMA_VERSION
        ):
            raise RuntimeError(
                f"scalp inventory {label} is not a terminal current-schema run")
    if manifest.get("inventory_file") != os.path.basename(path):
        raise RuntimeError("scalp inventory manifest names a different artifact")
    for key in (
        "requested_subjects",
        "n_requested",
        "classification_counts",
        "classification_subjects",
        "query_error_subjects",
        "code_revision",
        "code_dirty",
        "source_tree_sha256",
        "audit_script_sha256",
    ):
        if manifest.get(key) != payload.get(key):
            raise RuntimeError(
                f"scalp inventory payload/manifest disagree at {key}")
    actual_sha256 = file_sha256(path)
    if manifest.get("inventory_file_sha256") != actual_sha256:
        raise RuntimeError("scalp inventory bytes differ from its terminal manifest")
    audit_script = os.path.join(ROOT, "analysis", "audit_hup_scalp_inventory.py")
    if payload.get("audit_script_sha256") != file_sha256(audit_script):
        raise RuntimeError("scalp inventory was produced by a different audit script")
    inventory_provenance = validate_recorded_release_provenance(
        payload, ROOT, label="scalp inventory")
    manifest_provenance = validate_recorded_release_provenance(
        manifest, ROOT, label="scalp inventory manifest")
    if inventory_provenance != manifest_provenance:
        raise RuntimeError(
            "scalp inventory payload/manifest provenance differs")
    expected_cache_manifest = os.path.join(
        os.path.abspath(ieeg_cache_dir), "RUN_MANIFEST.json")
    expected_sidecar_manifest = os.path.join(
        os.path.abspath(scalp_cache_dir), "RUN_MANIFEST.json")
    lineage_checks = {
        "frozen_cache_manifest_sha256": file_sha256(
            expected_cache_manifest),
        "qc_grid_sha256": file_sha256(os.path.abspath(qc_grid_path)),
        "sidecar_manifest_sha256": file_sha256(
            expected_sidecar_manifest),
        "source_pin_sha256": file_sha256(os.path.join(
            ROOT, "analysis", "hup_ieeg_source_pin.json")),
        "current_locked_qc_profile_sha256": qc_profile_sha256(
            load_qc_profile(PROFILE_ID)),
    }
    for key, expected in lineage_checks.items():
        if payload.get(key) != expected:
            raise RuntimeError(
                f"scalp inventory lineage differs at {key}")
    _validate_classifications(payload)
    records = payload["subjects"]
    eligible = payload.get("classification_subjects", {}).get(
        "paired_3a_eligible", [])
    if (
        selected_subjects is not None
        and list(selected_subjects) != list(eligible)
    ):
        raise RuntimeError(
            "paired subjects do not exactly equal the all-cohort audited "
            "3A-eligible intersection in its frozen order: "
            f"selected={list(selected_subjects)!r}, "
            f"eligible={list(eligible)!r}")
    unique_f3_f4_labels = []
    for record in records:
        roles = record.get("roles", {})
        if (
            len(roles.get("f3", {}).get("matches", [])) == 1
            and len(roles.get("f4", {}).get("matches", [])) == 1
        ):
            unique_f3_f4_labels.append(record["subject"])
    return {
        "path_relative": os.path.relpath(path, ROOT).replace(os.sep, "/"),
        "sha256": actual_sha256,
        "manifest_path_relative": os.path.relpath(
            manifest_path, ROOT).replace(os.sep, "/"),
        "manifest_sha256": file_sha256(manifest_path),
        "n_frozen_hup_participants": int(payload["n_requested"]),
        "classification_counts": payload["classification_counts"],
        "classification_subjects": payload["classification_subjects"],
        "paired_3a_eligible_subjects": eligible,
        "participants_with_unique_f3_and_f4_labels": unique_f3_f4_labels,
        "reference_warning": payload["reference_warning"],
        "selection_was_blind_to_scalp_endpoint_values": payload[
            "selection_was_blind_to_scalp_endpoint_values"],
    }


def _seed(*parts):
    digest = hashlib.sha256(
        "|".join(map(str, parts)).encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "little")


def _validate_scalp_cache(subject, scalp_dir, pinned_ieeg):
    manifest_path = os.path.join(scalp_dir, "RUN_MANIFEST.json")
    cache_path = os.path.join(scalp_dir, f"{subject}.npz")
    if not os.path.isfile(manifest_path) or not os.path.isfile(cache_path):
        raise FileNotFoundError(f"missing paired-scalp cache for {subject}")
    with open(manifest_path) as handle:
        manifest = json.load(handle)
    cache_dir = os.path.dirname(pinned_ieeg["manifest_path"])
    manifest_pin = validate_terminal_sidecar_manifest(
        manifest,
        source=manifest_path,
        output_dir=scalp_dir,
        cache_dir=cache_dir,
    )
    current_dependency = manifest_pin["cache_dependency_sha256"]
    completed = set(manifest_pin["completed"]) | set(
        manifest_pin["reused"])
    if subject not in completed:
        raise RuntimeError(f"{subject} is not complete in the scalp manifest")
    expected = manifest_pin["result_files_sha256"].get(subject)
    actual = file_sha256(cache_path)
    if not expected or actual != expected:
        raise RuntimeError(
            f"{subject} scalp cache bytes differ from its manifest")
    with np.load(cache_path, allow_pickle=False) as scalp:
        validated_payload = validate_scalp_sidecar_payload(
            scalp,
            path=cache_path,
            subject=subject,
            pinned=pinned_ieeg,
            dependency_digest=current_dependency,
        )
    return {
        "path": cache_path,
        "sha256": actual,
        "manifest_path": manifest_path,
        "manifest_sha256": file_sha256(manifest_path),
        "cache_dependency_sha256": current_dependency,
        **validated_payload,
    }


def _load_frozen_profile_record(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"frozen QC grid is missing: {path}")
    with open(path) as handle:
        payload = json.load(handle)
    matching = [
        record for record in payload.get("profiles", [])
        if record.get("qc_profile_id") == PROFILE_ID
        and record.get("qc_profile_sha256")
        == qc_profile_sha256(load_qc_profile(PROFILE_ID))
    ]
    if len(matching) != 1:
        raise RuntimeError(
            f"{path} contains {len(matching)} exact frozen {PROFILE_ID} profiles")
    profile_record = matching[0]
    subjects = {}
    for value in profile_record.get("subjects", []):
        subject = value.get("subject")
        if not subject or subject in subjects:
            raise RuntimeError("frozen QC-grid subject records are malformed")
        subjects[subject] = value
    return payload, profile_record, subjects


def _assert_frozen_ieeg_match(subject, result_3a, result_3b, frozen):
    if subject not in frozen:
        raise RuntimeError(f"{subject} is absent from the frozen QC grid")
    expected = frozen[subject]
    for endpoint, actual in (
        ("result_3a", result_3a),
        ("result_3b", result_3b),
    ):
        if endpoint not in expected:
            raise RuntimeError(
                f"frozen record for {subject} lacks {endpoint}")
        actual_json = _canonical_json(actual)
        expected_json = _canonical_json(expected[endpoint])
        if actual_json != expected_json:
            # Give a bounded diagnostic rather than writing mismatched outputs.
            raise RuntimeError(
                f"recomputed {subject} {endpoint} differs from the frozen "
                "overlap11 QC-grid result; comparison aborted")
    return {
        "status": "exact_match",
        "checked_fields": ["result_3a", "result_3b"],
    }


def _contact_power(cache, band):
    numerator = np.asarray(
        cache[f"{band}_power_numerator_by_contact"], float)
    denominator = np.asarray(
        cache[f"{band}_clean_sample_count_by_contact"], float)
    samples_per_second = int(np.asarray(
        cache["power_samples_per_second"]).item())
    return power_from_binned_support(
        numerator,
        denominator,
        samples_per_second,
        minimum_clean_fraction_per_second=0.5,
    )


def _scalp_power(cache, band):
    numerator = np.asarray(
        cache[f"{band}_power_numerator_by_channel"], float)
    denominator = np.asarray(
        cache[f"{band}_clean_sample_count_by_channel"], float)
    samples_per_second = int(np.asarray(
        cache["power_samples_per_second"]).item())
    return power_from_binned_support(
        numerator,
        denominator,
        samples_per_second,
        minimum_clean_fraction_per_second=0.5,
    )


def _selected_mask(details, n_contacts):
    value = details.get(
        "selected_mask", details.get("selected_contact_mask"))
    if value is None:
        return np.zeros(n_contacts, bool)
    result = np.asarray(value, bool).ravel()
    if len(result) != n_contacts:
        raise RuntimeError("power selection mask does not align with contacts")
    return result


def _single_scalp_power_materialization(scalp, role):
    channels = [str(value) for value in scalp["scalp_chans"]]
    roles = json.loads(npz_scalar_text(scalp, "channel_roles_json"))
    if role not in roles:
        return None
    channel = str(roles[role])
    if channel not in channels:
        raise RuntimeError(f"scalp role {role} does not align with channels")
    index = channels.index(channel)
    nonflat = np.asarray(
        scalp["scalp_signal_nonflat_mask"], bool).ravel()
    if len(nonflat) != len(channels):
        raise RuntimeError("scalp flat-line mask does not align with channels")
    sigma_matrix = _scalp_power(scalp, "sigma_fixed")[[index]]
    swa_matrix = _scalp_power(scalp, "swa")[[index]]
    joint = (
        np.isfinite(sigma_matrix) & (sigma_matrix > 0)
        & np.isfinite(swa_matrix) & (swa_matrix > 0)
    )
    eligible = np.asarray([nonflat[index]], bool)
    sigma, sigma_qc = overlap_connected_aggregate(
        sigma_matrix,
        eligible_contacts=eligible,
        observation_mask=joint,
        transform="log",
        minimum_contacts=1,
    )
    selected_contacts = np.asarray(
        sigma_qc["selected_contact_mask"], bool)
    selected_times = np.asarray(sigma_qc["selected_time_mask"], bool)
    component = (
        joint
        & selected_contacts[:, None]
        & selected_times[None, :]
    )
    swa, swa_qc = overlap_connected_aggregate(
        swa_matrix,
        eligible_contacts=selected_contacts,
        observation_mask=component,
        transform="log",
        minimum_contacts=1,
    )
    raw_sigma = sigma_matrix[0]
    raw_swa = swa_matrix[0]
    return {
        "channel": channel,
        "channel_index": index,
        "nonflat": bool(nonflat[index]),
        "sigma": sigma,
        "swa": swa,
        "sigma_qc": sigma_qc,
        "swa_qc": swa_qc,
        "raw_sigma_median": (
            float(np.nanmedian(raw_sigma))
            if np.isfinite(raw_sigma).any() else None),
        "raw_swa_median": (
            float(np.nanmedian(raw_swa))
            if np.isfinite(raw_swa).any() else None),
    }


def _scalp_3a_materialized(base, scalp_power):
    result = dict(base)
    selected = np.asarray(
        scalp_power["sigma_qc"]["selected_contact_mask"], bool)
    sigma_qc = dict(scalp_power["sigma_qc"])
    sigma_qc["selected_mask"] = selected
    sigma_qc["n_selected"] = int(selected.sum())
    swa_qc = dict(scalp_power["swa_qc"])
    swa_qc["selected_mask"] = np.asarray(
        swa_qc["selected_contact_mask"], bool)
    swa_qc["n_selected"] = int(
        np.asarray(swa_qc["selected_contact_mask"], bool).sum())
    result.update(
        contacts=np.asarray([scalp_power["channel"]], dtype="<U16"),
        cortical_signal_nonflat_mask=np.asarray(
            [scalp_power["nonflat"]], bool),
        cortical_signal_activity_qc_provenance=(
            "raw full-interval numerical flat-line mask in scalp sidecar"),
        parietal_contact_mask=np.ones(1, bool),
        frontal_contact_mask=np.ones(1, bool),
        sigma_global=scalp_power["sigma"],
        swa_global=scalp_power["swa"],
        sigma_parietal=scalp_power["sigma"],
        swa_parietal=scalp_power["swa"],
        global_power_qc={"sigma": sigma_qc, "swa": swa_qc},
        parietal_power_qc={"sigma": sigma_qc, "swa": swa_qc},
        parietal_power_coverage=float(
            np.isfinite(scalp_power["sigma"]).mean()),
    )
    return result


def _endpoint_local_3a_materializations(ieeg, scalp, profile):
    """Build matched EEG-only and EEG-plus-HR inputs for paired 3A.

    Missing HR is irrelevant to the sigma spectrum, its fitted peak, and the
    same-window SWA control. Applying the cardiac mask to those EEG-only
    endpoints can split or remove otherwise valid 120-s EEG bouts. Coherence
    and cross-correlation do require HR, so they use the strict subset.
    """
    keys = ("sigma_parietal", "swa_parietal", "hr_1")
    arrays = {
        ("ieeg", key): np.asarray(ieeg[key], float)
        for key in keys
    }
    arrays.update({
        ("scalp", key): np.asarray(scalp[key], float)
        for key in keys
    })
    shapes = {value.shape for value in arrays.values()}
    if len(shapes) != 1 or len(next(iter(shapes))) != 1:
        raise RuntimeError("paired 3A arrays do not share one one-dimensional time base")
    if not np.array_equal(
            arrays[("ieeg", "hr_1")],
            arrays[("scalp", "hr_1")],
            equal_nan=True):
        raise RuntimeError("paired 3A arms do not share the exact HR array")
    eeg_common = np.ones(
        next(iter(shapes)), dtype=bool)
    for arm in ("ieeg", "scalp"):
        for key in ("sigma_parietal", "swa_parietal"):
            values = arrays[(arm, key)]
            eeg_common &= np.isfinite(values) & (values > 0)
    cardiac_common = (
        eeg_common & np.isfinite(arrays[("ieeg", "hr_1")])
    )
    minimum_hr_coverage = float(profile["hr"]["minimum_coverage"])
    if (
        not np.isfinite(minimum_hr_coverage)
        or not 0 <= minimum_hr_coverage <= 1
    ):
        raise RuntimeError(
            "paired 3A profile has an invalid HR coverage threshold")

    outputs = {}
    for arm, source in (("ieeg", ieeg), ("scalp", scalp)):
        endpoint_materializations = {}
        for endpoint, support in (
                ("eeg", eeg_common),
                ("cardiac", cardiac_common)):
            result = dict(source)
            result["sigma_parietal"] = np.where(
                support, arrays[(arm, "sigma_parietal")], np.nan)
            result["swa_parietal"] = np.where(
                support, arrays[(arm, "swa_parietal")], np.nan)
            result["hr_1"] = (
                np.asarray(arrays[(arm, "hr_1")], float).copy()
                if endpoint == "eeg"
                else np.where(
                    support, arrays[(arm, "hr_1")], np.nan)
            )
            if endpoint == "eeg":
                # analyse_3a currently combines EEG and cardiac endpoints.
                # Disable its unused cardiac pass here; the assembled result
                # receives those endpoints only from the cardiac-common pass.
                result["hr_coverage"] = 0.0
                result["hr_meets_profile"] = False
            else:
                shared_hr_coverage = float(
                    np.isfinite(result["hr_1"]).mean())
                result["hr_coverage"] = shared_hr_coverage
                result["hr_meets_profile"] = bool(
                    shared_hr_coverage >= minimum_hr_coverage)
            # These aliases are not read by analyse_3a today, but keeping them
            # coherent prevents a later diagnostic from silently switching
            # back to an independently supported series.
            result["sigma_global"] = result["sigma_parietal"]
            result["swa_global"] = result["swa_parietal"]
            result["parietal_power_coverage"] = float(support.mean())
            result["paired_3a_support_arm"] = arm
            result["paired_3a_support_endpoint"] = endpoint
            endpoint_materializations[endpoint] = result
        outputs[arm] = endpoint_materializations

    labels = np.asarray(ieeg["stage_lab"]).astype(str)
    nrem, analysis_window = lecci.core_study_nrem_mask(labels)
    second_nrem = np.repeat(nrem, int(lecci.EPOCH))
    if len(second_nrem) != len(eeg_common):
        raise RuntimeError(
            "paired 3A second-level arrays do not align exactly with "
            "30-s stage labels")
    support = {
        "rule": {
            "eeg_common": (
                "finite positive iEEG sigma+SWA AND finite positive scalp "
                "sigma+SWA; used for spectra, peaks, and negative controls"),
            "cardiac_common": (
                "EEG-common support AND finite shared HR; used for fixed "
                "0.02-Hz coherence and 120-s cross-correlation"),
        },
        "n_total_seconds": int(len(eeg_common)),
        "eeg_common": {
            "n_seconds": int(eeg_common.sum()),
            "fraction": float(eeg_common.mean()),
            "n_seconds_in_core_nrem": int(
                (eeg_common & second_nrem).sum()),
            "support_mask_sha256": _array_sha256(eeg_common),
        },
        "cardiac_common": {
            "n_seconds": int(cardiac_common.sum()),
            "fraction": float(cardiac_common.mean()),
            "minimum_hr_coverage": minimum_hr_coverage,
            "passes_hr_coverage_profile": bool(
                cardiac_common.mean() >= minimum_hr_coverage),
            "n_seconds_in_core_nrem": int(
                (cardiac_common & second_nrem).sum()),
            "support_mask_sha256": _array_sha256(cardiac_common),
        },
        "cardiac_common_is_subset_of_eeg_common": bool(
            np.all(~cardiac_common | eeg_common)),
        "analysis_window": analysis_window,
    }
    return outputs, support


def _xcorr_retained_window_starts(materialized):
    """Return the exact 120-s starts admitted by Lecci cross-correlation."""
    sig = np.asarray(materialized["sigma_parietal"], float)
    hr = np.asarray(materialized["hr_1"], float)
    labels = np.asarray(materialized["stage_lab"]).astype(str)
    nrem, _ = lecci.core_study_nrem_mask(labels)
    smoothed, _, _ = lecci.fill_short_gaps(sig, lecci.FS, 5.0)
    heart, _, _ = lecci.fill_short_gaps(hr, lecci.FS, 5.0)
    kernel = max(1, int(round(4 * lecci.FS)))
    valid = np.isfinite(smoothed)
    numerator = np.convolve(
        np.where(valid, smoothed, 0.0), np.ones(kernel), mode="same")
    denominator = np.convolve(
        valid.astype(float), np.ones(kernel), mode="same")
    smoothed = np.where(
        denominator >= kernel / 2,
        numerator / np.maximum(denominator, 1),
        np.nan,
    )
    window = int(lecci.XCORR_WIN_S * lecci.FS)
    starts = []
    for start, stop in lecci.nrem_bouts(
            nrem, min_s=lecci.XCORR_WIN_S):
        for candidate in range(start, stop - window + 1, window):
            heart_window = heart[candidate:candidate + window]
            signal_window = smoothed[candidate:candidate + window]
            if not (
                np.isfinite(heart_window).all()
                and np.isfinite(signal_window).all()
            ):
                continue
            if heart_window.std() < 1e-9 or signal_window.std() < 1e-9:
                continue
            starts.append(int(candidate))
    return starts


def _assemble_endpoint_local_3a_result(eeg_result, cardiac_result):
    """Keep EEG-only outputs independent of the stricter cardiac support."""
    result = copy.deepcopy(eeg_result)
    eeg_support_reasons = [
        reason for reason in eeg_result["support_reasons"]
        if not str(reason).endswith("; cardiac endpoints only")
    ]
    result["support_reasons"] = eeg_support_reasons
    result["hr_coverage"] = cardiac_result["hr_coverage"]
    result["coherence"] = copy.deepcopy(cardiac_result["coherence"])
    result["cross_correlation"] = copy.deepcopy(
        cardiac_result["cross_correlation"])
    result["endpoint_availability"] = {
        "spectrum": bool(
            eeg_result["endpoint_availability"]["spectrum"]),
        "fixed_0p02_coherence": bool(
            cardiac_result["endpoint_availability"][
                "fixed_0p02_coherence"]),
        "cross_correlation": bool(
            cardiac_result["endpoint_availability"][
                "cross_correlation"]),
    }
    result["endpoint_support_diagnostics"] = {
        "spectrum_peak_negative_control": {
            "aggregate_coverage": eeg_result["aggregate_coverage"],
            "support_passes_profile": eeg_result[
                "support_passes_profile"],
            "n_bouts": eeg_result["n_bouts"],
            "bout_seconds": eeg_result["bout_seconds"],
            "support_reasons": copy.deepcopy(eeg_support_reasons),
        },
        "coherence_cross_correlation": {
            "aggregate_coverage": cardiac_result["aggregate_coverage"],
            "hr_coverage": cardiac_result["hr_coverage"],
            "support_passes_profile": cardiac_result[
                "support_passes_profile"],
            "n_bouts": cardiac_result["n_bouts"],
            "bout_seconds": cardiac_result["bout_seconds"],
            "support_reasons": copy.deepcopy(
                cardiac_result["support_reasons"]),
        },
    }
    return result


def _assert_matched_3a_geometry(
        subject, ieeg_eeg_result, scalp_eeg_result,
        ieeg_cardiac_result, scalp_cardiac_result,
        ieeg_cardiac_materialized, scalp_cardiac_materialized):
    """Fail if either endpoint-local paired support geometry diverges."""
    for key in ("n_bouts", "bout_seconds"):
        if ieeg_eeg_result.get(key) != scalp_eeg_result.get(key):
            raise RuntimeError(
                f"{subject} paired 3A EEG-common support produced different "
                f"spectrum {key}")
    left_spectrum = ieeg_eeg_result[
        "endpoint_availability"]["spectrum"]
    right_spectrum = scalp_eeg_result[
        "endpoint_availability"]["spectrum"]
    if bool(left_spectrum) != bool(right_spectrum):
        raise RuntimeError(
            f"{subject} paired 3A spectrum availability differs on "
            "EEG-common support")

    left_co = ieeg_cardiac_result.get("coherence")
    right_co = scalp_cardiac_result.get("coherence")
    if (left_co is None) != (right_co is None):
        raise RuntimeError(
            f"{subject} paired 3A coherence availability differs on "
            "cardiac-common support")
    if left_co is not None:
        for key in ("K", "n_valid", "analytic_threshold"):
            if not np.isclose(
                    float(left_co[key]), float(right_co[key]), rtol=0, atol=1e-12):
                raise RuntimeError(
                    f"{subject} paired 3A coherence geometry differs at {key}")
    left_xc = ieeg_cardiac_result.get("cross_correlation")
    right_xc = scalp_cardiac_result.get("cross_correlation")
    if (left_xc is None) != (right_xc is None):
        raise RuntimeError(
            f"{subject} paired 3A xcorr availability differs on "
            "cardiac-common support")
    if (
        left_xc is not None
        and int(left_xc["n_intervals"]) != int(right_xc["n_intervals"])
    ):
        raise RuntimeError(
            f"{subject} paired 3A xcorr interval counts differ on shared support")
    left_starts = _xcorr_retained_window_starts(
        ieeg_cardiac_materialized)
    right_starts = _xcorr_retained_window_starts(
        scalp_cardiac_materialized)
    if left_starts != right_starts:
        raise RuntimeError(
            f"{subject} paired 3A xcorr retained different 120-s windows")
    starts_hash = _array_sha256(np.asarray(left_starts, dtype=np.int64))
    return {
        "status": "exact_endpoint_local_geometry",
        "spectrum_peak_negative_control": {
            "status": "exact_shared_eeg_geometry",
            "checked_fields": [
                "endpoint_availability.spectrum",
                "n_bouts",
                "bout_seconds",
            ],
        },
        "coherence_cross_correlation": {
            "status": "exact_shared_cardiac_geometry",
            "cross_correlation_retained_window_starts_sha256": starts_hash,
            "cross_correlation_retained_window_count": len(left_starts),
            "checked_fields": [
                "coherence.K",
                "coherence.n_valid",
                "coherence.analytic_threshold",
                "cross_correlation.n_intervals",
                "cross_correlation.retained_window_start_indices",
            ],
        },
    }


def _stage_epochs(labels, stage):
    labels = np.asarray(labels).astype(str)
    if stage == "NREM":
        pooled = np.where(
            np.isin(labels, ("NREM", "N2", "N3")), "NREM", "")
        return (
            stable_stage_epoch_indices(pooled, "NREM"),
            int(np.isin(labels, ("NREM", "N2", "N3")).sum()),
        )
    return (
        stable_stage_epoch_indices(labels, stage),
        int((labels == stage).sum()),
    )


def _stage_pool(rr, epochs):
    mask = np.zeros(len(rr), bool)
    for epoch in epochs:
        start = int(epoch * lecci.EPOCH * FS_RR)
        stop = int((epoch + 1) * lecci.EPOCH * FS_RR)
        mask[start:stop] = True
    return np.where(mask & np.isfinite(rr))[0]


def _analyse_scalp_3b(scalp, base, role, profile):
    roles = json.loads(npz_scalar_text(scalp, "channel_roles_json"))
    if role not in roles:
        return {
            "role": role,
            "channel": None,
            "available": False,
            "reason": "sensor absent",
            "stages": {},
        }
    channel = str(roles[role])
    channels = [str(value) for value in scalp["scalp_chans"]]
    channel_index = channels.index(channel)
    nonflat = bool(np.asarray(
        scalp["scalp_signal_nonflat_mask"], bool)[channel_index])
    rr = np.asarray(base["rr_4"], float)
    labels = np.asarray(base["stage_lab"]).astype(str)
    endpoint = profile["endpoint_3b"]
    result = {
        "role": role,
        "channel": channel,
        "available": True,
        "sensor_nonflat": nonflat,
        "exploratory": True,
        "profile_id": profile["profile_id"],
        "minimum_contacts": 1,
        "original_profile_minimum_contacts": 2,
        "method_status": (
            "unilateral F3 sensitivity; Naji derived a cardiac curve per "
            "referenced F3/F4 electrode and averaged the electrode-specific "
            "HR-maximum/RR-minimum times"
            if role == "f3"
            else "Fz exploratory sensitivity; Fz was not a Naji sensor"),
        "reference_status": (
            "portal channel reference is undocumented; not known to be "
            "F3/A2, F4/A1, or linked mastoids"),
        "polarity_sensitivity": {
            "status": POLARITY_SENSITIVITY_STATUS,
            "reason": (
                "the outcome-neutral sidecar stores candidates detected in the "
                "recorded polarity, not the filtered waveform or the complementary "
                "positive-to-negative half-wave candidate set; opposite-polarity "
                "detection requires a prespecified re-stream"),
        },
        "stages": {},
    }
    for stage in STAGES_3B:
        stable, raw_epochs = _stage_epochs(labels, stage)
        pool = _stage_pool(rr, stable)
        reasons = []
        if not nonflat:
            reasons.append("scalp sensor is numerically flat")
        if not bool(base.get("hr_meets_profile", False)):
            reasons.append("cardiac coverage below profile")
        if np.isfinite(rr).sum() < endpoint["minimum_finite_rr_samples"]:
            reasons.append(
                f"{int(np.isfinite(rr).sum())} finite RR samples < "
                f"{endpoint['minimum_finite_rr_samples']}")
        if len(stable) < 6:
            reasons.append("no uninterrupted 180-s stage run")
        if len(pool) < endpoint["minimum_finite_stage_samples"]:
            reasons.append(
                f"{len(pool)} finite stage samples < "
                f"{endpoint['minimum_finite_stage_samples']}")
        keep = set(stable.tolist())
        times = stage_so_times(
            scalp,
            channel,
            keep,
            percentile=float(endpoint["so_amplitude_percentile"]),
        )
        if len(times) < endpoint["minimum_so_per_contact"]:
            reasons.append(
                f"{len(times)} eligible SOs < "
                f"{endpoint['minimum_so_per_contact']}")
        estimate = None
        if not reasons:
            baseline = rr_baseline_hr(rr[pool])
            estimate = subject_so_triggered(
                rr,
                [times],
                baseline,
                pool,
                n_sur=N_SURROGATES_3B,
                rng=np.random.RandomState(
                    _seed(base["subject"], stage, role, "paired-scalp-3B")),
                domain="rr",
                minimum_channels=1,
                channel_ids=[channel],
                minimum_events_per_channel=endpoint[
                    "minimum_so_per_contact"],
                minimum_surrogate_pool_samples=endpoint[
                    "minimum_finite_stage_samples"],
            )
            if estimate is None:
                reasons.append(
                    "complete in-stage RR windows leave insufficient SO support")
            else:
                estimate["stage_mean_hr"] = float(baseline)
        result["stages"][stage] = {
            "raw_epochs": raw_epochs,
            "stable_epochs": int(len(stable)),
            "stable_seconds": int(len(stable) * lecci.EPOCH),
            "pooled_nrem_exploratory": stage == "NREM",
            "candidate_so_count": int(len(times)),
            "so_amplitude_percentile": float(
                endpoint["so_amplitude_percentile"]),
            "so_amplitude_rule_provenance": (
                "within-sensor/stage percentile adaptation; Naji cites fixed "
                "Dang-Vu scalp-voltage gates and does not specify this percentile"),
            "estimate": estimate,
            "available_under_exploratory_profile": estimate is not None,
            "support_reasons": reasons,
        }
    return result


def _available_3b_stages(result):
    """Return stages with a valid current 3B estimate in canonical order."""
    stages = result.get("stages", {}) if isinstance(result, dict) else {}
    return [
        stage for stage in STAGES_3B
        if bool(stages.get(stage, {}).get("available_under_profile"))
    ]


def _available_exploratory_scalp_stages(record, role):
    """Return valid stages for one explicitly exploratory scalp role."""
    if record is None:
        return []
    result = record.get("scalp", {}).get("result_3b", {}).get(role, {})
    stages = result.get("stages", {}) if isinstance(result, dict) else {}
    return [
        stage for stage in STAGES_3B
        if bool(
            stages.get(stage, {}).get(
                "available_under_exploratory_profile"
            )
        )
    ]


def _naji_method_scope():
    """Describe the implemented 3B scope without implying exact transfer."""
    return {
        "exact_comparison_available": False,
        "implemented_scalp_estimators": [
            "unilateral F3 exploratory sensitivity",
            "Fz exploratory sensitivity (not a Naji sensor)",
        ],
        "f4_status": (
            "planned for a separate bilateral analysis; not analyzed here"
        ),
        "blocking_reasons": list(NAJI_EXACT_BLOCKING_REASONS),
    }


def _naji_label_inventory_status(subject, frozen_subjects, records):
    """Separate channel-label inventory from the analyses actually performed.

    A valid staged iEEG endpoint cannot make the scalp comparison exact. The
    current scalp estimator analyzes unilateral F3 and non-Naji Fz only. F4 may
    be present in a validated terminal sidecar for the separate bilateral
    sensitivity, but it is not analyzed here and the portal references are
    undocumented.
    """
    record_by_subject = {
        record["subject"]: record
        for record in records
    }
    record = record_by_subject.get(subject)
    frozen_3b = frozen_subjects.get(subject, {}).get("result_3b", {})
    acquired_roles = sorted(
        record.get("lineage", {}).get("scalp_sidecar_roles", {})
        if record is not None
        else []
    )
    return {
        "subject": subject,
        "channel_label_inventory": {
            "unique_f3_label_present": True,
            "unique_f4_label_present": True,
            "labels_establish_naji_reference_montage": False,
        },
        "locked_ieeg_3b_available_stages": _available_3b_stages(
            frozen_3b
        ),
        "participant_in_current_paired_3a_intersection": record is not None,
        "terminal_sidecar_roles_acquired": acquired_roles,
        "current_exploratory_scalp_3b": {
            "unilateral_f3": {
                "analyzed": record is not None and "f3" in acquired_roles,
                "available_stages": _available_exploratory_scalp_stages(
                    record, "f3"
                ),
                "interpretation": (
                    "unilateral/reference-incomplete sensitivity only"
                ),
            },
            "fz": {
                "analyzed": record is not None and "fz" in acquired_roles,
                "available_stages": _available_exploratory_scalp_stages(
                    record, "fz"
                ),
                "interpretation": (
                    "exploratory sensitivity; Fz was not a Naji sensor"
                ),
            },
            "f4": {
                "acquired": "f4" in acquired_roles,
                "analyzed": False,
                "available_stages": [],
                "interpretation": (
                    (
                        "validated terminal sidecar role reserved for the "
                        "separate bilateral analysis; excluded from this "
                        "estimator and group claims"
                    )
                    if "f4" in acquired_roles
                    else (
                        "channel plan/label inventory does not prove terminal "
                        "acquisition; excluded from this estimator and claims"
                    )
                ),
            },
        },
        "current_exact_naji_comparison_available": False,
        "exact_naji_blocking_reasons": list(
            NAJI_EXACT_BLOCKING_REASONS
        ),
    }


def _raw_ieeg_power_summary(cache, materialized):
    sigma = _contact_power(cache, "sigma_fixed")
    swa = _contact_power(cache, "swa")
    selected = _selected_mask(
        materialized["parietal_power_qc"]["sigma"], sigma.shape[0])
    sigma_medians = np.asarray([
        np.nanmedian(row) if np.isfinite(row).any() else np.nan
        for row in sigma
    ])
    swa_medians = np.asarray([
        np.nanmedian(row) if np.isfinite(row).any() else np.nan
        for row in swa
    ])
    return {
        "selected_contact_ids": [
            str(value) for value, keep in zip(
                materialized["contacts"], selected) if keep
        ],
        "selected_contact_sigma_medians": sigma_medians[selected],
        "selected_contact_swa_medians": swa_medians[selected],
        "median_selected_contact_sigma_power": (
            float(np.nanmedian(sigma_medians[selected]))
            if np.isfinite(sigma_medians[selected]).any() else None),
        "median_selected_contact_swa_power": (
            float(np.nanmedian(swa_medians[selected]))
            if np.isfinite(swa_medians[selected]).any() else None),
        "interpretation": (
            "descriptive recorded-unit band power only; reference, impedance, "
            "gain, geometry, and spatial scale differ, so the ratio is not a "
            "universal scalp-versus-iEEG amplitude law"),
    }


def analyse_subject(subject, *, profile, exploratory_profile,
                    frozen_subjects, ieeg_cache_dir, scalp_cache_dir):
    pinned = validate_pinned_ieeg_cache(subject, ieeg_cache_dir)
    scalp_lineage = _validate_scalp_cache(
        subject, scalp_cache_dir, pinned)
    with (
        np.load(pinned["path"], allow_pickle=False) as ieeg,
        np.load(scalp_lineage["path"], allow_pickle=False) as scalp,
    ):
        materialized = materialize(ieeg, profile)
        mask = np.asarray(
            materialized["cortical_signal_nonflat_mask"], bool)
        activity_provenance = materialized[
            "cortical_signal_activity_qc_provenance"]
        frozen_result_3a_ieeg = analyse_3a(materialized, profile)
        result_3b_ieeg = analyse_3b(ieeg, materialized, profile)
        frozen_match = _assert_frozen_ieeg_match(
            subject, frozen_result_3a_ieeg, result_3b_ieeg, frozen_subjects)

        c3 = _single_scalp_power_materialization(scalp, "c3")
        if c3 is None:
            raise RuntimeError(f"{subject} has no C3/C03 primary channel")
        scalp_materialized = _scalp_3a_materialized(materialized, c3)
        # These hashes prove that the EEG-arm comparison did not change ECG or
        # staging. Modality, location, reference, and aggregation still differ.
        for key in ("hr_1", "rr_4", "stage_lab"):
            if not np.array_equal(
                    np.asarray(scalp_materialized[key]),
                    np.asarray(materialized[key]),
                    equal_nan=True):
                raise RuntimeError(
                    f"scalp 3A adapter changed shared input {key}")
        independent_result_3a_scalp = analyse_3a(
            scalp_materialized, profile)
        endpoint_materializations, shared_3a_support = (
            _endpoint_local_3a_materializations(
                materialized, scalp_materialized, profile)
        )
        ieeg_eeg_result = analyse_3a(
            endpoint_materializations["ieeg"]["eeg"], profile)
        scalp_eeg_result = analyse_3a(
            endpoint_materializations["scalp"]["eeg"], profile)
        ieeg_cardiac_result = analyse_3a(
            endpoint_materializations["ieeg"]["cardiac"], profile)
        scalp_cardiac_result = analyse_3a(
            endpoint_materializations["scalp"]["cardiac"], profile)
        result_3a_ieeg = _assemble_endpoint_local_3a_result(
            ieeg_eeg_result, ieeg_cardiac_result)
        result_3a_scalp = _assemble_endpoint_local_3a_result(
            scalp_eeg_result, scalp_cardiac_result)
        shared_3a_geometry = _assert_matched_3a_geometry(
            subject,
            ieeg_eeg_result,
            scalp_eeg_result,
            ieeg_cardiac_result,
            scalp_cardiac_result,
            endpoint_materializations["ieeg"]["cardiac"],
            endpoint_materializations["scalp"]["cardiac"],
        )
        result_3b_scalp = {
            role: _analyse_scalp_3b(
                scalp, materialized, role, exploratory_profile)
            for role in ("f3", "fz")
        }
        ieeg_raw = _raw_ieeg_power_summary(ieeg, materialized)
        shared = {
            "stage_lab_sha256": _array_sha256(
                materialized["stage_lab"]),
            "hr_1_sha256": _array_sha256(materialized["hr_1"]),
            "rr_4_sha256": _array_sha256(materialized["rr_4"]),
            "n_stage_epochs": int(len(materialized["stage_lab"])),
            "n_hr_seconds": int(len(materialized["hr_1"])),
            "n_rr_4hz_samples": int(len(materialized["rr_4"])),
            "paired_3a_endpoint_support": shared_3a_support,
            "paired_3a_geometry_verification": shared_3a_geometry,
            "statement": (
                "identical ECG/RR and sleep labels are used for both arms; "
                "3A spectra/peaks/SWA controls use identical finite-positive "
                "EEG support, while coherence/cross-correlation use its "
                "finite-shared-HR subset; ECG was not redetected and sleep "
                "was not restaged"),
        }
        raw_ratio = None
        if (
            ieeg_raw["median_selected_contact_sigma_power"] is not None
            and c3["raw_sigma_median"] is not None
            and ieeg_raw["median_selected_contact_sigma_power"] > 0
        ):
            raw_ratio = float(
                c3["raw_sigma_median"]
                / ieeg_raw["median_selected_contact_sigma_power"])
        return {
            "subject": subject,
            "lineage": {
                "ieeg_cache_sha256": pinned["sha256"],
                "ieeg_cache_manifest_sha256": pinned["manifest_sha256"],
                "ieeg_cache_manifest_run_id": pinned["manifest_run_id"],
                "scalp_cache_sha256": scalp_lineage["sha256"],
                "scalp_cache_manifest_sha256": scalp_lineage[
                    "manifest_sha256"],
                "scalp_sidecar_roles": dict(scalp_lineage["roles"]),
                "night_s": pinned["night_s"],
                "hours": pinned["hours"],
                "sf": pinned["sf"],
            },
            "frozen_ieeg_verification": frozen_match,
            "ieeg_activity_qc": {
                "mask": mask,
                "provenance": activity_provenance,
                "selected_contact_ids": [
                    str(value) for value, keep in zip(
                        ieeg["cortical_chans"], mask) if keep
                ],
            },
            "shared_inputs": shared,
            "ieeg": {
                "result_3a": result_3a_ieeg,
                "frozen_independent_support_result_3a": (
                    frozen_result_3a_ieeg),
                "result_3b": result_3b_ieeg,
                "raw_power": ieeg_raw,
                "profile_id": profile["profile_id"],
            },
            "scalp": {
                "result_3a": result_3a_scalp,
                "independent_support_result_3a": (
                    independent_result_3a_scalp),
                "result_3b": result_3b_scalp,
                "raw_power": {
                    "channel": c3["channel"],
                    "raw_sigma_median": c3["raw_sigma_median"],
                    "raw_swa_median": c3["raw_swa_median"],
                    "scalp_to_ieeg_median_sigma_power_ratio": raw_ratio,
                    "interpretation": ieeg_raw["interpretation"],
                },
                "profile_3a_id": profile["profile_id"],
                "profile_3b_id": exploratory_profile["profile_id"],
            },
        }


def _parse_subjects(value):
    if value is None or not str(value).strip():
        return None
    result = []
    for item in str(value).split(","):
        token = item.strip()
        if not token:
            continue
        subject = (
            token if token.startswith("HUP")
            else f"HUP{int(token)}_phaseII")
        if subject.startswith("HUP") and not subject.endswith("_phaseII"):
            subject += "_phaseII"
        if subject not in SCALP_CHANNEL_PLAN:
            raise ValueError(f"{subject} is outside the paired-scalp plan")
        result.append(subject)
    if len(result) != len(set(result)):
        raise ValueError("subjects contain duplicates")
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Run pinned simultaneous scalp-versus-iEEG comparisons")
    parser.add_argument("--subjects")
    parser.add_argument("--ieeg-cache", default=DEFAULT_IEEG_CACHE)
    parser.add_argument("--scalp-cache", default=DEFAULT_SCALP_CACHE)
    parser.add_argument(
        "--scalp-inventory", default=DEFAULT_SCALP_INVENTORY)
    parser.add_argument("--qc-grid", default=DEFAULT_QC_GRID)
    parser.add_argument("--output-dir", default=DEFAULT_RESULTS)
    args = parser.parse_args()
    release_provenance = require_clean_release_provenance(ROOT)
    requested_subjects = _parse_subjects(args.subjects)
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    profile = load_qc_profile(PROFILE_ID)
    exploratory_profile = copy.deepcopy(profile)
    exploratory_profile["profile_id"] = (
        "overlap11_endpoint_local__single_scalp_channel_exploratory")
    exploratory_profile["role"] = (
        "post-audit simultaneous single-scalp-sensor exploratory comparison")
    exploratory_profile["selection_blinded_to_endpoint_values"] = True
    exploratory_profile["endpoint_3b"]["minimum_contacts"] = 1

    frozen_payload, frozen_profile, frozen_subjects = (
        _load_frozen_profile_record(os.path.abspath(args.qc_grid)))
    scalp_inventory = _validate_scalp_inventory(
        args.scalp_inventory,
        requested_subjects,
        qc_grid_path=args.qc_grid,
        ieeg_cache_dir=args.ieeg_cache,
        scalp_cache_dir=args.scalp_cache,
    )
    subjects = list(scalp_inventory["paired_3a_eligible_subjects"])
    records = []
    for subject in subjects:
        print(f"[{subject}] paired comparison", flush=True)
        records.append(analyse_subject(
            subject,
            profile=profile,
            exploratory_profile=exploratory_profile,
            frozen_subjects=frozen_subjects,
            ieeg_cache_dir=os.path.abspath(args.ieeg_cache),
            scalp_cache_dir=os.path.abspath(args.scalp_cache),
        ))
    summary = _group_summary(records)
    naji_label_status = [
        _naji_label_inventory_status(subject, frozen_subjects, records)
        for subject in scalp_inventory[
            "participants_with_unique_f3_and_f4_labels"
        ]
    ]
    metadata = {
        "schema_version": RESULT_SCHEMA,
        "pipeline": PIPELINE,
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "generated_at_utc": utc_now(),
        **release_provenance,
        "runtime_versions": runtime_versions(),
        "requested_subjects": subjects,
        "completed_subjects": [value["subject"] for value in records],
        "profile_id": PROFILE_ID,
        "profile_sha256": qc_profile_sha256(profile),
        "profile_file_sha256": profile_file_sha256(),
        "exploratory_single_channel_profile": exploratory_profile,
        "exploratory_single_channel_profile_sha256": qc_profile_sha256(
            exploratory_profile),
        "frozen_qc_grid_path_relative": os.path.relpath(
            os.path.abspath(args.qc_grid), ROOT).replace(os.sep, "/"),
        "frozen_qc_grid_sha256": file_sha256(
            os.path.abspath(args.qc_grid)),
        "frozen_qc_grid_profile_sha256": frozen_profile[
            "qc_profile_sha256"],
        "frozen_qc_grid_cache_manifest_sha256": frozen_payload.get(
            "cache_manifest_sha256"),
        "scalp_inventory": scalp_inventory,
        "naji_label_inventory_status": naji_label_status,
        "naji_method_scope": _naji_method_scope(),
        "comparison_design": (
            "within-participant simultaneous EEG-arm comparison: same portal "
            "snapshot, frozen night, duration, sample rate, cached ECG/RR, "
            "and profile-materialized stage labels. Spectra, peaks, and SWA "
            "controls use exact common finite-positive EEG support; coherence "
            "and cross-correlation use its exact finite-shared-HR subset. The "
            "EEG arms differ jointly in modality, location, reference, and "
            "contact aggregation, so this does not isolate a pure modality "
            "effect"),
        "3a_sensor": (
            "single C3/C03 scalp channel, Lecci-motivated location adaptation "
            "but "
            "unknown online reference; compared with the frozen iEEG aggregate"),
        "3b_sensors": (
            "single F3 (unilateral/reference-incomplete) and Fz (not a Naji "
            "sensor), both exploratory one-channel profile clones"),
        "sleep_label_provenance": (
            "HUP sleep labels are an iEEG-derived delta/SWA GMM proxy without "
            "independent visual PSG scoring or EOG/EMG. Holding the same labels "
            "fixed controls timing but makes this comparison asymmetric and "
            "does not validate a scalp-paper staging pipeline."),
        "lc_specificity": (
            "none: these are downstream cortical/cardiac observables. Neither "
            "arm records LC neurons, norepinephrine, nor an independently "
            "validated LC signal, and neither selectively perturbs LC."),
        "csv_artifacts": {
            "paired_metrics.csv": (
                "normalized observations; unique by "
                "subject/question/stage/modality, including exactly one iEEG "
                "row per participant-stage for 3B"
            ),
            "role_pair_metrics.csv": (
                "explicit scalp-role pairs; a 3B iEEG observation is "
                "intentionally repeated under F3 and Fz when it contributes "
                "to both contrasts, so do not pool this table across roles"
            ),
        },
    }
    full = {
        **metadata,
        "subjects": records,
        "group_summary": summary,
    }
    require_release_source_unchanged(ROOT, release_provenance)
    atomic_json_dump(
        _json_safe(full),
        os.path.join(output_dir, "subject_results.json"),
    )
    atomic_json_dump(
        _json_safe({**metadata, "group_summary": summary}),
        os.path.join(output_dir, "group_summary.json"),
    )
    role_pair_rows = _role_pair_csv_rows(records)
    normalized_rows = _normalized_csv_rows(records, role_pair_rows)
    _write_csv(
        os.path.join(output_dir, "paired_metrics.csv"),
        normalized_rows,
    )
    _write_csv(
        os.path.join(output_dir, "role_pair_metrics.csv"),
        role_pair_rows,
    )
    _make_3a_figure(records, output_dir)
    _make_3b_figure(records, output_dir)
    require_release_source_unchanged(ROOT, release_provenance)
    manifest = {
        **metadata,
        "run_state": "complete",
        "result_files_sha256": {
            name: file_sha256(os.path.join(output_dir, name))
            for name in (
                "subject_results.json",
                "group_summary.json",
                "paired_metrics.csv",
                "role_pair_metrics.csv",
                "paired_3A_C3.png",
                "paired_3A_C3.svg",
                "paired_3B_F3_Fz.png",
                "paired_3B_F3_Fz.svg",
            )
        },
    }
    atomic_json_dump(
        _json_safe(manifest),
        os.path.join(output_dir, "RUN_MANIFEST.json"),
    )
    print(
        f"{len(records)} paired participants -> {output_dir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
