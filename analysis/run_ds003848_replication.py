"""Withdrawn ds003848 endpoint writer.

This script previously called the standalone 3A and 3B analysis functions and
created a second artifact contract. Current RESPect analysis is materialized
from the neutral ds003848 cache by ``run_qc_grid.py`` so the same locked
profiles and endpoint implementations are used for both cohorts.
"""


def main():
    raise SystemExit(
        "WITHDRAWN ds003848 3A/3B writer: first build "
        "data/derived/ds003848 with analysis/stage_ds003848.py, then run "
        "analysis/run_qc_grid.py --cache-dir data/derived/ds003848")


if __name__ == "__main__":
    main()
