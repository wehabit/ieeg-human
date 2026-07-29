"""Withdrawn 3A analysis retained only as a signal-QC compatibility shim.

The historical executable used superseded estimators and must not be run. Active code should
import the shared helpers from :mod:`signal_qc`.
"""

from signal_qc import dilate_boolean_mask, ied_clean_mask, robust_z


__all__ = ["robust_z", "dilate_boolean_mask", "ied_clean_mask"]


def main():
    raise SystemExit(
        "WITHDRAWN: use cache_lc_series.py followed by lecci_faithful_3A.py; "
        "import shared signal-QC helpers from signal_qc.py.")


if __name__ == "__main__":
    main()
