# Analysis artifact policy

This repository separates source code, private derived data, and public result
evidence.  The goal is a clean checkout that is small enough to review but still
fails closed when a claimed result file is missing, modified, or stale.

The checked-in numerical evidence remains the frozen v8 baseline until the
ordered v9 rebuild below completes. Current v9 readers must reject those bytes
rather than relabeling them.

## What belongs in Git

Commit:

- analysis code, locked configuration, source pins, and tests;
- compact QC summaries under `outputs/qc_grid_public/`;
- the portable, hash-pinned v8-to-v9 comparison summary under
  `outputs/rebuild_comparison/`;
- the one full locked-profile snapshot used by the paired scalp–iEEG
  exact-match check;
- terminal result manifests, participant/group tables, and figures that are
  cited by the documentation.

Do not commit:

- raw EEG, portal downloads, credentials, or derived NPZ caches under `data/`;
- full profile-grid JSON under `outputs/qc_grid/`;
- raw field-level rebuild-comparison reports, which embed local input paths;
- temporary logs, interrupted-run files, local plotting caches, or obsolete
  free-text summaries from withdrawn writers.

Legacy `.txt` cohort summaries are not an authoritative result format. A
current claim must resolve to a terminal manifest and its hashed structured
participant/group tables or figures; narrative explanation belongs in the
versioned documentation. Historical text output can be recovered from Git
history when needed and must not be retained merely because an old script once
generated it.

The eight former full QC grids repeated spectra and event diagnostics for every
profile and occupied about 58.5 MB (1.35 million lines).  The tracked public
replacement retains every profile definition, every cohort count, and every
participant’s endpoint-availability booleans.  A separate locked snapshot
retains the complete `overlap11_endpoint_local` HUP records needed by the
paired comparison.  Full grids remain deterministic local build products.

Removing files from the current tree does not shrink existing Git history.  A
history rewrite would affect every clone and is therefore a separate,
explicitly coordinated operation; this analysis does not perform one.

## Rebuild order for a future dataset or source change

Create the pinned environment, then rebuild in dependency order:

```bash
python -m venv .venv
.venv/bin/python -m pip install -r env/requirements.txt

# Private, ignored neutral caches.
.venv/bin/python analysis/cache_lc_series.py --force
.venv/bin/python analysis/stage_ds003848.py --force
.venv/bin/python analysis/calibrate_staging_windows.py \
  --update-profile-pin

# Full local grids (ignored).
for grid in coverage_oat_v1 staging_window_support_v1 \
            auxiliary_window_support_v1 event_count_oat_v1; do
  .venv/bin/python analysis/run_qc_grid.py \
    --cache-dir data/derived/lc_infraslow --grid "$grid"
  .venv/bin/python analysis/run_qc_grid.py \
    --cache-dir data/derived/ds003848 --grid "$grid"
done

# Reviewable public grid evidence plus the full locked-profile snapshot.
.venv/bin/python analysis/compact_qc_grid_artifacts.py

# Simultaneous scalp sidecars, all-cohort inventory, and paired results.
.venv/bin/python analysis/cache_paired_scalp.py --force
.venv/bin/python analysis/audit_hup_scalp_inventory.py
# Replace the placeholder with the inventory's exact ordered paired_3a_eligible list.
.venv/bin/python analysis/paired_scalp_ieeg_comparison.py \
  --subjects PAIRED_3A_ELIGIBLE_COMMA_LIST
```

The calibration pin update changes provenance only; it does not tune a QC
threshold. It fails if the newly measured recommendation differs from the
locked recommendation. Such a difference requires an explicit scientific
method decision before grids are run.

Do not run the withdrawn standalone `lecci_faithful_3A.py`,
`event_3B_cached.py`, `event_3D_by_stage.py`, or
`run_ds003848_replication.py` writers. They intentionally exit. The two neutral
cache builders plus `run_qc_grid.py` are the single production path.

Do not manually replace a stored digest after changing code. A
`cache_code_sha256` mismatch requires the neutral caches to be regenerated. A
downstream-only edit does not invalidate those expensive cache bytes, but its
dependent grid/result artifacts must be regenerated with the new
`source_tree_sha256`. When a cache-producing change requires portal data, the
honest interim state is **rebuild pending**, not a relabeled historical artifact.

## Validation modes

Run publication validation on a workstation that has the private derived
manifests:

```bash
.venv/bin/python analysis/test_paired_scalp.py \
  --artifact-validation publication
.venv/bin/python analysis/test_coherence_calibration.py --require-real-cache
```

Clean-checkout CI uses:

```bash
.venv/bin/python analysis/test_paired_scalp.py \
  --artifact-validation offline
```

Offline mode still requires and hashes every tracked inventory, summary,
locked-profile snapshot, result table, and figure.  It also checks their
cross-links and the current source-tree digest.  It may skip only the absent
ignored `data/derived` manifests, and reports those skips explicitly.  Deleting
a tracked evidence file is therefore a test failure, never a silent pass.

Paired validation also reconstructs both CSV exports from
`subject_results.json` and requires exact serialized equality. This prevents a
secondary flat file from being modified or duplicated and then legitimized by
updating only its hash.

`outputs/qc_grid_public/RUN_MANIFEST.json` is written last and hashes the exact
eight summaries plus the locked snapshot. Validators additionally recompute
profile hashes, require every profile prescribed by `qc_profiles_v1.json`, and
derive the locked profile's compact availability record from its full subject
records. Rehashing a modified or incomplete artifact therefore cannot bypass
semantic validation. When `compact_qc_grid_artifacts.py --output-root PATH` is
used, the summaries, locked snapshot, and terminal manifest all remain under
that path; a validation/test build cannot overwrite the repository default.

Publication mode adds exact validation of the live cache and paired-sidecar
manifests.  It must remain red while a required current-data rebuild is pending.
