"""Offline, outcome-neutral reconstruction of the HUP 3D event endpoint.

The HUP neutral cache retains 20-Hz, per-contact 12--16-Hz RMS and
0.16--1.25-Hz SO phase, a sample-validity mask, and every clean
duration-qualified SO candidate before amplitude thresholding.  This is enough
to repeat stage restriction, channel-night percentile thresholds, spindle
detection, SO/spindle pairing, and every locked support profile without
streaming the raw night again.

Only descriptive equal-contact participant vectors are reported.  The finite
SO-centred pairing window can itself align phases under independent event
trains, so neither Rayleigh nor cohort p values are valid without a null that
shifts or block-resamples complete spindle trains and repeats the full pairing
and contact-aggregation procedure.
"""
from __future__ import annotations

import numpy as np

from event_3d_estimators import (
    EVENT_FS,
    channel_night_events,
    circular_summary,
    stage_event_pairs,
)


NEUTRAL_3D_REQUIRED_FIELDS = (
    "event_3d_rms_12_16_by_contact",
    "event_3d_so_phase_0p16_1p25_by_contact",
    "event_3d_valid_sample_mask_by_contact",
    "event_3d_so_candidate_contact_index",
    "event_3d_so_candidate_sample",
    "event_3d_so_candidate_amplitude",
    "event_3d_so_candidate_cycle_start_sample",
    "event_3d_so_candidate_cycle_stop_sample_exclusive",
    "event_3d_sampling_hz",
)

_THRESHOLD_KEYS = (
    "minimum_events_per_contact",
    "minimum_paired_events",
    "minimum_contact_acquisition_fraction",
    "minimum_contact_event_valid_fraction",
    "minimum_valid_nrem_fraction_per_contact",
    "minimum_valid_nrem_seconds_per_contact",
    "minimum_contacts",
)


def _cache_files(cache):
    files = getattr(cache, "files", None)
    if files is None:
        files = cache.keys()
    return set(files)


def _profile_thresholds(profile):
    endpoint = profile["endpoint_3d"]
    return {
        "minimum_events_per_contact": int(
            endpoint["minimum_events_per_contact"]),
        "minimum_paired_events": int(endpoint["minimum_paired_events"]),
        "minimum_contact_acquisition_fraction": float(
            endpoint["minimum_contact_acquisition_fraction"]),
        "minimum_contact_event_valid_fraction": float(
            endpoint["minimum_contact_event_valid_fraction"]),
        "minimum_valid_nrem_fraction_per_contact": float(
            endpoint["minimum_valid_nrem_fraction_per_contact"]),
        "minimum_valid_nrem_seconds_per_contact": int(
            endpoint["minimum_valid_nrem_seconds_per_contact"]),
        "minimum_contacts": int(endpoint["minimum_contacts"]),
    }


def _validated_neutral_payload(cache, n_contacts):
    """Load and shape-check neutral fields, returning errors instead of guessing."""
    files = _cache_files(cache)
    missing = [
        field for field in NEUTRAL_3D_REQUIRED_FIELDS if field not in files
    ]
    if missing:
        return None, missing, []

    errors = []
    # Preserve the cache's compact float32 representation across a sensitivity
    # grid.  The event detector promotes one contact at a time to float, so a
    # second full-night float64 copy would add memory without changing results.
    rms = np.asarray(
        cache["event_3d_rms_12_16_by_contact"], dtype=np.float32)
    phase = np.asarray(
        cache["event_3d_so_phase_0p16_1p25_by_contact"], dtype=np.float32)
    valid = np.asarray(
        cache["event_3d_valid_sample_mask_by_contact"], bool)
    contact_index = np.asarray(
        cache["event_3d_so_candidate_contact_index"])
    candidate_sample = np.asarray(cache["event_3d_so_candidate_sample"])
    candidate_amplitude = np.asarray(
        cache["event_3d_so_candidate_amplitude"], float)
    candidate_start = np.asarray(
        cache["event_3d_so_candidate_cycle_start_sample"])
    candidate_stop = np.asarray(
        cache["event_3d_so_candidate_cycle_stop_sample_exclusive"])
    sampling_hz_values = np.asarray(
        cache["event_3d_sampling_hz"], float).reshape(-1)

    if rms.ndim != 2 or rms.shape[0] != n_contacts:
        errors.append(
            "event RMS must be a contact-by-sample matrix aligned to anatomy contacts")
    if phase.shape != rms.shape or valid.shape != rms.shape:
        errors.append("event RMS, SO phase, and validity mask must align")
    if (
        contact_index.ndim != 1
        or candidate_sample.ndim != 1
        or candidate_amplitude.ndim != 1
        or candidate_start.ndim != 1
        or candidate_stop.ndim != 1
        or not (
            len(contact_index)
            == len(candidate_sample)
            == len(candidate_amplitude)
            == len(candidate_start)
            == len(candidate_stop)
        )
    ):
        errors.append(
            "SO candidate identity, amplitude, and cycle-bound arrays must "
            "be aligned one-dimensional vectors")
    if len(sampling_hz_values) != 1 or not np.isclose(
            sampling_hz_values[0], EVENT_FS):
        errors.append(
            f"event sampling rate must be the validated {EVENT_FS:g} Hz")

    if not errors:
        contact_index_float = np.asarray(contact_index, float)
        candidate_sample_float = np.asarray(candidate_sample, float)
        candidate_start_float = np.asarray(candidate_start, float)
        candidate_stop_float = np.asarray(candidate_stop, float)
        if (
            not np.isfinite(contact_index_float).all()
            or not np.equal(
                contact_index_float, np.floor(contact_index_float)).all()
        ):
            errors.append("SO candidate contact indices must be finite integers")
        if (
            not np.isfinite(candidate_sample_float).all()
            or not np.equal(
                candidate_sample_float, np.floor(candidate_sample_float)).all()
        ):
            errors.append("SO candidate samples must be finite integers")
        if (
            not np.isfinite(candidate_start_float).all()
            or not np.equal(
                candidate_start_float, np.floor(candidate_start_float)).all()
            or not np.isfinite(candidate_stop_float).all()
            or not np.equal(
                candidate_stop_float, np.floor(candidate_stop_float)).all()
        ):
            errors.append("SO candidate cycle bounds must be finite integers")
        if not np.isfinite(candidate_amplitude).all():
            errors.append("SO candidate amplitudes must be finite")

    if not errors:
        contact_index = contact_index.astype(int)
        candidate_sample = candidate_sample.astype(int)
        candidate_start = candidate_start.astype(int)
        candidate_stop = candidate_stop.astype(int)
        if (
            (contact_index < 0).any()
            or (contact_index >= n_contacts).any()
        ):
            errors.append("SO candidate contact index is out of bounds")
        if rms.shape[1] == 0 or (
            (candidate_sample < 0).any()
            or (candidate_sample >= rms.shape[1]).any()
        ):
            errors.append("SO candidate sample is out of bounds")
        if (
            (candidate_start < 0).any()
            or (candidate_stop > rms.shape[1]).any()
            or (candidate_stop <= candidate_start).any()
            or (candidate_sample < candidate_start).any()
            or (candidate_sample >= candidate_stop).any()
        ):
            errors.append(
                "SO candidate cycle bounds are invalid or do not contain "
                "the trough")

    if errors:
        return None, missing, errors

    # The explicit validity mask is authoritative.  Nonfinite values beneath a
    # true mask remain unavailable and are handled fail-closed by the detector.
    jointly_valid = valid & np.isfinite(rms) & np.isfinite(phase)
    payload = {
        "rms": np.where(jointly_valid, rms, np.nan),
        "phase": np.where(jointly_valid, phase, np.nan),
        "valid": jointly_valid,
        "candidate_contact_index": contact_index,
        "candidate_sample": candidate_sample,
        "candidate_amplitude": candidate_amplitude,
        "candidate_start": candidate_start,
        "candidate_stop": candidate_stop,
        "sampling_hz": float(sampling_hz_values[0]),
    }
    return payload, missing, errors


def _staging_contact_support(cache, payload, n_contacts, n_epochs):
    """Reconstruct outcome-blind staging-candidate support."""
    files = _cache_files(cache)
    required = {
        "ep_clean_fraction_by_contact",
        "ep_measured_fraction_by_contact",
    }
    missing = sorted(required - files)
    if missing:
        return None, None, [
            f"cache lacks per-contact staging support fields: {missing}"
        ]
    clean = np.asarray(cache["ep_clean_fraction_by_contact"], float)
    measured = np.asarray(cache["ep_measured_fraction_by_contact"], float)
    if (
        clean.ndim != 2
        or clean.shape != measured.shape
        or clean.shape != (n_contacts, n_epochs)
    ):
        return None, None, [
            "per-contact epoch clean/measured fractions do not align with "
            "anatomy contacts and materialized stage epochs"
        ]
    finite_measured = np.isfinite(measured)
    if (
        ((measured[finite_measured] < 0)
         | (measured[finite_measured] > 1)).any()
    ):
        return None, None, [
            "per-contact epoch measured fractions must lie within [0, 1]"
        ]
    if payload is None:
        return None, None, [
            "validated event sample support is unavailable for the staging clean-fraction gate"
        ]
    # Equal 30-s epochs make this mean exactly the observed sample fraction
    # over the nominal night.  A nonfinite epoch is unavailable (zero), not a
    # fully observed epoch merely because a placeholder clean fraction exists.
    observed = np.where(finite_measured, measured, 0.0).mean(axis=1)
    # The clean quantity uses the event detector's +/-2.5-s IED/missing-data
    # mask.  The ordinary staging cache uses a different 0.5-s mask, so its
    # clean values cannot be substituted here.
    clean_fraction = payload["valid"].mean(axis=1)
    return observed, clean_fraction, []


def _stage_summary(
        stage_name,
        labels,
        contacts,
        events,
        base_contact_support,
        thresholds,
):
    """Reduce one stage to equal-contact descriptive vectors and locked support."""
    if stage_name == "pooled_NREM":
        keep_epochs = set(np.where(
            np.isin(labels, ("NREM", "N2", "N3"))
        )[0].tolist())
    else:
        keep_epochs = set(np.where(labels == stage_name)[0].tolist())

    per_contact = []
    qualifying_vectors = []
    qualifying_event_count = 0
    for contact_index, contact in enumerate(contacts):
        record = events[contact_index]
        stage_pairs = stage_event_pairs(record, keep_epochs)
        phases = np.asarray(stage_pairs["phases"], float)
        circular = circular_summary(phases)
        n_events = int(len(phases))
        meets_events = (
            n_events >= thresholds["minimum_events_per_contact"])
        qualifies = bool(base_contact_support[contact_index] and meets_events)
        if qualifies:
            vector = complex(np.mean(np.exp(1j * phases)))
            qualifying_vectors.append(vector)
            qualifying_event_count += n_events
        per_contact.append({
            "contact": contact,
            "n_paired_events": n_events,
            "n_selected_so_candidates_pooled_nrem": int(record["n_so"]),
            "n_duration_qualified_spindles_pooled_nrem": int(
                record["n_spindle"]),
            "valid_pooled_nrem_seconds": float(
                record["valid_nrem_seconds"]),
            "valid_pooled_nrem_fraction": float(
                record["valid_nrem_fraction"]),
            "meets_base_contact_support": bool(
                base_contact_support[contact_index]),
            "meets_minimum_events_per_contact": bool(meets_events),
            "qualifies_for_stage_vector": qualifies,
            "descriptive_circular_summary": circular,
        })

    n_qualified = int(len(qualifying_vectors))
    passes = bool(
        n_qualified >= thresholds["minimum_contacts"]
        and qualifying_event_count >= thresholds["minimum_paired_events"]
    )
    participant_effect = None
    if passes:
        participant_vector = complex(np.mean(qualifying_vectors))
        participant_effect = {
            "n_equal_weight_contacts": n_qualified,
            "n_paired_events_across_qualified_contacts": int(
                qualifying_event_count),
            "participant_vector_real": float(participant_vector.real),
            "participant_vector_imag": float(participant_vector.imag),
            "participant_R": float(abs(participant_vector)),
            "participant_preferred_phase_rad": float(
                np.angle(participant_vector)),
            "participant_preferred_phase_deg": float(
                np.degrees(np.angle(participant_vector))),
            "aggregation": (
                "complex mean within contact, then equal-contact complex mean"),
        }
    return {
        "stage": stage_name,
        "n_stage_epochs": int(
            np.isin(labels, ("NREM", "N2", "N3")).sum()
            if stage_name == "pooled_NREM"
            else (labels == stage_name).sum()
        ),
        "n_contacts_meeting_all_contact_thresholds": n_qualified,
        "n_paired_events_across_qualified_contacts": int(
            qualifying_event_count),
        "support_passes_profile": passes,
        "descriptive_effect": participant_effect,
        "per_contact": per_contact,
    }


def analyse_3d_cache_support(cache, materialized, profile):
    """Apply one locked profile and reconstruct descriptive 3D vectors offline."""
    files = _cache_files(cache)
    cohort = str(materialized["cohort"])
    labels = np.asarray(materialized["stage_lab"]).astype(str)
    contacts = [str(value) for value in materialized["contacts"]]
    nrem = np.isin(labels, ("NREM", "N2", "N3"))
    thresholds = _profile_thresholds(profile)
    in_scope = cohort == "HUP"
    staging_fit_converged = bool(
        materialized.get("staging_qc", {}).get(
            "support_passes_fit_convergence", False))
    payload, missing, validation_errors = _validated_neutral_payload(
        cache, len(contacts))
    observed, clean_fraction, staging_errors = _staging_contact_support(
        cache, payload, len(contacts), len(labels))
    validation_errors.extend(staging_errors)
    reconstructable = bool(
        in_scope
        and staging_fit_converged
        and payload is not None
        and not validation_errors)

    reasons = []
    if not in_scope:
        reasons.append(
            "the specified 3D endpoint is defined for HUP; "
            "RESPect has no neutral 3D event cache")
    if not staging_fit_converged:
        reasons.append(
            "common overlap-connected staging fit did not converge")
    if missing:
        reasons.append(
            "neutral cache lacks pre-threshold 12-16-Hz RMS, 0.16-1.25-Hz "
            "SO phase/candidates, and/or their sample-validity support")
    reasons.extend(validation_errors)
    if "sigma_fixed_by_contact" in files and missing:
        reasons.append(
            "one-second 10-15-Hz power is not a substitute for 12-16-Hz "
            "duration-qualified spindle times and SO phase")
    if (
        any(name.startswith("so_candidate_t_") for name in files)
        and missing
    ):
        reasons.append(
            "stored SO candidates implement the different 3B Naji "
            "0.15-4-Hz half-wave definition and are not 3D candidates")

    common = {
        "endpoint": "3D_SO_spindle_phase",
        "cohort_in_scope": in_scope,
        "profile_thresholds": thresholds,
        "cache_prerequisite_support": {
            "n_anatomy_candidate_contacts": int(len(contacts)),
            "n_N2_like_epochs": int(np.sum(labels == "N2")),
            "n_N3_like_epochs": int(np.sum(labels == "N3")),
            "n_pooled_nrem_epochs": int(nrem.sum()),
            "n_prethreshold_duration_qualified_so_candidates": (
                None
                if payload is None
                else int(len(payload["candidate_sample"]))
            ),
            "valid_event_samples_by_contact": (
                None
                if payload is None
                else payload["valid"].sum(axis=1).astype(int).tolist()
            ),
            "stage_source": (
                "unvalidated HUP intracranial power proxy"
                if in_scope
                else "RESPect annotation-constrained proxy"),
        },
        "neutral_event_cache_schema": {
            "required_fields": list(NEUTRAL_3D_REQUIRED_FIELDS),
            "missing_fields": missing,
            "schema_present": not missing,
            "validation_errors": validation_errors,
            "schema_validated": bool(not missing and not validation_errors),
        },
        "inference_enabled": False,
        "inferential_p_value": None,
        "staging_fit_converged": staging_fit_converged,
    }

    if not reconstructable:
        return {
            **common,
            "threshold_evaluability": {
                key: False for key in _THRESHOLD_KEYS
            },
            "support_passes_profile": None,
            "descriptive_effect": None,
            "descriptive_effect_available": False,
            "stage_descriptive_support": None,
            "status": "not_in_scope" if not in_scope else "not_reconstructable",
            "support_reasons": reasons,
            "required_next_step": (
                "rebuild the HUP neutral cache with the validated 3D fields"
                if in_scope
                else "none for this cohort unless a RESPect 3D endpoint is specified"
            ),
        }

    candidate_contact = payload["candidate_contact_index"]
    candidate_sample = payload["candidate_sample"]
    candidate_amplitude = payload["candidate_amplitude"]
    candidate_start = payload["candidate_start"]
    candidate_stop = payload["candidate_stop"]
    events = []
    for contact_index in range(len(contacts)):
        contact_candidates = candidate_contact == contact_index
        candidates = np.column_stack((
            candidate_sample[contact_candidates],
            candidate_amplitude[contact_candidates],
            candidate_start[contact_candidates],
            candidate_stop[contact_candidates],
        ))
        record = channel_night_events(
            payload["rms"][contact_index],
            payload["phase"][contact_index],
            candidates,
            labels,
        )
        events.append(record)

    contact_availability_pass = (
        (observed >= thresholds[
            "minimum_contact_acquisition_fraction"])
        & (clean_fraction >= thresholds[
            "minimum_contact_event_valid_fraction"])
    )
    valid_seconds = np.asarray([
        value["valid_nrem_seconds"] for value in events
    ], float)
    valid_fraction = np.asarray([
        value["valid_nrem_fraction"] for value in events
    ], float)
    valid_nrem_pass = (
        (valid_seconds >= thresholds[
            "minimum_valid_nrem_seconds_per_contact"])
        & (valid_fraction >= thresholds[
            "minimum_valid_nrem_fraction_per_contact"])
    )
    base_contact_support = contact_availability_pass & valid_nrem_pass

    stages = {
        stage_name: _stage_summary(
            stage_name,
            labels,
            contacts,
            events,
            base_contact_support,
            thresholds,
        )
        for stage_name in ("N2", "N3", "pooled_NREM")
    }
    pooled = stages["pooled_NREM"]
    support_passes = bool(pooled["support_passes_profile"])
    if not support_passes:
        if int(contact_availability_pass.sum()) < thresholds["minimum_contacts"]:
            reasons.append(
                f"{int(contact_availability_pass.sum())} contacts pass acquisition/event-valid "
                f"support < {thresholds['minimum_contacts']}")
        if int(valid_nrem_pass.sum()) < thresholds["minimum_contacts"]:
            reasons.append(
                f"{int(valid_nrem_pass.sum())} contacts pass valid-NREM support "
                f"< {thresholds['minimum_contacts']}")
        if pooled[
                "n_contacts_meeting_all_contact_thresholds"
        ] < thresholds["minimum_contacts"]:
            reasons.append(
                f"{pooled['n_contacts_meeting_all_contact_thresholds']} "
                f"event-supported contacts < {thresholds['minimum_contacts']}")
        if pooled[
                "n_paired_events_across_qualified_contacts"
        ] < thresholds["minimum_paired_events"]:
            reasons.append(
                f"{pooled['n_paired_events_across_qualified_contacts']} paired "
                f"events < {thresholds['minimum_paired_events']}")

    return {
        **common,
        "threshold_evaluability": {
            key: True for key in _THRESHOLD_KEYS
        },
        "per_contact_acquisition_fraction": (
            observed.tolist()),
        "per_contact_event_valid_fraction": (
            clean_fraction.tolist()),
        "per_contact_passes_acquisition_event_valid_support": (
            contact_availability_pass.tolist()),
        "per_contact_passes_valid_nrem_support": valid_nrem_pass.tolist(),
        "support_passes_profile": support_passes,
        "descriptive_effect": (
            pooled["descriptive_effect"] if support_passes else None),
        "descriptive_effect_available": bool(support_passes),
        "stage_descriptive_support": stages,
        "status": (
            "descriptive_available"
            if support_passes
            else "support_below_profile"
        ),
        "support_reasons": reasons,
        "required_next_step": (
            "validate HUP proxy stages and implement a complete-train "
            "time-shift/block null before any inferential 3D claim"
        ),
    }
