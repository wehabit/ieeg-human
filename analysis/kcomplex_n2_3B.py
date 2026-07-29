"""Withdrawn K-complex follow-up.

The historical script depended on obsolete whole-night SO arrays and exposed a
nonstationarity-sensitive event-locking statistic.  It is retained only as a
clear pointer for old commands; no scientific estimator lives here.
"""


def main():
    raise SystemExit(
        "WITHDRAWN: the former K-complex follow-up is not a valid production "
        "analysis. Use cache_lc_series.py plus run_qc_grid.py for the current "
        "descriptive 3B endpoint."
    )


if __name__ == "__main__":
    main()
