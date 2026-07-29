"""Withdrawn direct-stream 3B entry point.

The reusable, publication-facing SO/RR estimators now live in
``event_3b_estimators.py``.  Names are re-exported here so historical imports
remain compatible without maintaining a second estimator implementation.

Current endpoint execution and artifact materialization go through
``run_qc_grid.py``.
"""
from event_3b_estimators import (
    FS_RR,
    HALF_WIN,
    rr_baseline_hr,
    so_triggered,
    subject_so_triggered,
)

__all__ = [
    "FS_RR",
    "HALF_WIN",
    "rr_baseline_hr",
    "so_triggered",
    "subject_so_triggered",
]


def main():
    raise SystemExit(
        "LEGACY/WITHDRAWN direct 3B entry point: run "
        "analysis/run_qc_grid.py to execute the profile-materialized "
        "3B endpoint")


if __name__ == "__main__":
    main()
