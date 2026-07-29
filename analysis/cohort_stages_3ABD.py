"""Withdrawn cohort analysis retained only as a staging-helper compatibility shim.

The historical executable used superseded estimators and must not be run. Active code should
import the shared helpers from :mod:`staging_helpers`.
"""

from staging_helpers import (
    ALPHA,
    CHUNK_S,
    EPOCH,
    FS_P,
    F_TARGET,
    INFRA,
    MIN_3A_MIN,
    NREM_DR,
    POOLED_3A_CAP_MIN,
    SO_BAND,
    SWA_BAND,
    band_sos,
    dominant_block,
    fsp_from,
    msc_block,
    nrem_mask_adaptive,
    reliable_two_state_split,
    stage_epochs,
)


__all__ = [
    "EPOCH",
    "FS_P",
    "CHUNK_S",
    "SWA_BAND",
    "SO_BAND",
    "INFRA",
    "F_TARGET",
    "NREM_DR",
    "MIN_3A_MIN",
    "ALPHA",
    "POOLED_3A_CAP_MIN",
    "fsp_from",
    "band_sos",
    "reliable_two_state_split",
    "nrem_mask_adaptive",
    "stage_epochs",
    "dominant_block",
    "msc_block",
]


def main():
    raise SystemExit(
        "WITHDRAWN: use the versioned cache and corrected analyses; "
        "import shared staging/spectral helpers from staging_helpers.py.")


if __name__ == "__main__":
    main()
