"""Shared schema names for checked-in LC-proxy evidence artifacts.

Keeping these constants in a dependency-light module lets producers and
clean-checkout validators enforce the same contract without importing portal
clients, plotting code, or scientific estimators.
"""

SCALP_INVENTORY_SCHEMA = (
    "2026-07-hup-scalp-channel-inventory-v2-versioned")
PAIRED_RESULT_SCHEMA = (
    "2026-07-paired-scalp-ieeg-results-v4-versioned-evidence")

QC_PUBLIC_SUMMARY_SCHEMA = "2026-07-qc-grid-public-summary-v1"
QC_LOCKED_SNAPSHOT_SCHEMA = (
    "2026-07-qc-grid-locked-profile-snapshot-v1")
QC_PUBLIC_MANIFEST_SCHEMA = "2026-07-qc-grid-public-manifest-v1"
QC_PUBLIC_PIPELINE = "compact_qc_grid_artifacts"
