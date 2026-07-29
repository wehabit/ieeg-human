"""WITHDRAWN legacy 3A entry point.

Current HUP portal/source helpers live in :mod:`hup_portal`.  They are
re-exported here only so historical scripts remain importable.  The superseded
raw-stream/FSP estimators were removed rather than left as an unreachable second
scientific pipeline.
"""
from __future__ import annotations

import numpy as np

from hup_portal import (
    COHORT,
    HUP_SOURCE_PINS,
    HUP_SOURCE_PIN_PATH,
    HUP_SOURCE_PIN_SCHEMA_VERSION,
    NIGHT_PROBE_WORKERS,
    STANDARD_SCALP_EEG_LABELS,
    NightProbeSourceMismatch,
    PortalSampleCountMismatch,
    cortical_channels,
    delta_ratio,
    expected_portal_sample_count,
    find_night,
    is_standard_scalp_eeg_label,
    pull_continuous_exact,
    validate_hup_series_geometry,
    verify_hup_source_identity,
)


def tort_mi(phase, amplitude, nbins=18):
    """Legacy compatibility helper used by the withdrawn cohort-stage script."""
    phase = np.asarray(phase, float)
    amplitude = np.asarray(amplitude, float)
    bins = np.linspace(-np.pi, np.pi, int(nbins) + 1)
    means = np.asarray([
        np.mean(amplitude[(phase >= bins[index]) & (phase < bins[index + 1])])
        if np.any((phase >= bins[index]) & (phase < bins[index + 1]))
        else np.nan
        for index in range(int(nbins))
    ])
    if np.isnan(means).any() or means.sum() <= 0:
        return np.nan
    probability = means / means.sum()
    return float(
        (
            np.log(nbins)
            + np.sum(probability * np.log(probability + 1e-12))
        )
        / np.log(nbins)
    )


def main():
    raise SystemExit(
        "WITHDRAWN: use cache_lc_series.py followed by the current offline "
        "3A/3B/3D endpoint modules")


if __name__ == "__main__":
    main()
