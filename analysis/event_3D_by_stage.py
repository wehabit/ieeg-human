"""Withdrawn direct-stream entry point for the 3D SO-to-spindle endpoint.

The direct portal runner had its own staging aggregation and could disagree
with the locked, profile-materialized analysis.  Current analysis must build
the neutral cache with ``cache_lc_series.py`` and evaluate it through
``run_qc_grid.py``.  Pure event estimators live in
``event_3d_estimators.py``.

Scientific constants and helpers are re-exported temporarily for compatibility
with historical audit code.  No direct-stream analysis is executed here.
"""
from __future__ import annotations

from cache_lc_series import (
    HUP_ANATOMY_SELECTION_METHOD,
    MIN_CONTACT_COVERAGE,
    MIN_CONTACT_FRACTION_PER_BIN,
    MIN_CONTACTS,
)
from event_3d_estimators import (
    EVENT_FS,
    EVENT_PERCENTILE,
    IED_PAD_S,
    MIN_EVENTS,
    MIN_POOLED_CONTACTS,
    MIN_POOLED_EVENTS,
    MIN_POOLED_NREM_COVERAGE,
    MIN_POOLED_VALID_S_PER_CONTACT,
    PAIR_WINDOW_S,
    PRODUCTION_3D_INFERENCE_ENABLED,
    SO_BAND,
    SO_DUR,
    SPINDLE_BAND,
    SP_DUR,
    bh_fdr,
    channel_night_events,
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
from hup_portal import HUP_SOURCE_PIN_SCHEMA_VERSION
from pipeline_version import ANALYSIS_VERSION


def production_config(hours=7.0):
    """Return the frozen configuration of the withdrawn direct-stream run.

    This remains available only to validate already-generated historical
    artifacts.  It is not the configuration entry point for current analysis.
    """
    return {
        "hours": float(hours),
        "analysis_version": ANALYSIS_VERSION,
        "source_pin_schema_version": HUP_SOURCE_PIN_SCHEMA_VERSION,
        "so_band_hz": list(SO_BAND),
        "spindle_band_hz": list(SPINDLE_BAND),
        "so_duration_s": list(SO_DUR),
        "spindle_duration_s": list(SP_DUR),
        "event_percentile": EVENT_PERCENTILE,
        "event_sampling_hz": EVENT_FS,
        "ied_mask_padding_s": IED_PAD_S,
        "minimum_events_per_contact": MIN_EVENTS,
        "threshold_scope": "channel-night pooled NREM",
        "pairing": (
            "duration-qualified hybrid; one spindle per SO within "
            f"+/-{PAIR_WINDOW_S:g} s"),
        "staging_qc": {
            "minimum_contacts": MIN_CONTACTS,
            "minimum_contact_feature_coverage": MIN_CONTACT_COVERAGE,
            "minimum_contact_fraction_per_epoch":
                MIN_CONTACT_FRACTION_PER_BIN,
            "minimum_candidate_observed_fraction": 0.80,
            "minimum_candidate_clean_fraction": 0.80,
            "swa_normalization": (
                "divide each fixed contact by its full-night median "
                "clean-epoch SWA"),
        },
        "anatomy_selection_method": HUP_ANATOMY_SELECTION_METHOD,
        "pooled_qc": {
            "minimum_contacts": MIN_POOLED_CONTACTS,
            "minimum_paired_events": MIN_POOLED_EVENTS,
            "minimum_valid_nrem_seconds_per_contact":
                MIN_POOLED_VALID_S_PER_CONTACT,
            "minimum_valid_nrem_fraction_per_contact":
                MIN_POOLED_NREM_COVERAGE,
        },
        "production_inference_enabled":
            PRODUCTION_3D_INFERENCE_ENABLED,
        "group_inference": (
            "disabled; participant vectors descriptive pending a "
            "time-shift/block null that repeats event pairing and contact "
            "aggregation"),
    }


def main():
    raise SystemExit(
        "WITHDRAWN direct 3D runner: build neutral caches with "
        "analysis/cache_lc_series.py, then run "
        "analysis/run_qc_grid.py --cache-dir data/derived/lc_infraslow")


if __name__ == "__main__":
    main()
