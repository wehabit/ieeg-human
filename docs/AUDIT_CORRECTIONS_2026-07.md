# Method and code audit — 2026-07

## Status

All committed numerical cohort results are legacy artifacts. Current readers require
`analysis_version = 2026-07-corrected-v5` and
`cache_schema_version = 2026-07-per-contact-power-pchip-night-threshold-v5`.
No committed cohort result currently passes those gates, so biological conclusions are pending a
raw-data rerun and QC.

In this document, **FIXED** means that the corresponding implementation change is present in the
working tree; the issue register gives its executable regression where available and otherwise a
specific acceptance check. **OPEN** means that the concern still requires code, real-data
validation, expert review, or a new measurement. A code-level fix does not convert staging,
anatomy, raw detector quality, or LC construct validity into a validated result.

The complete severity/evidence/fix/acceptance-test matrix is in
[ISSUE_REGISTER_2026-07.md](ISSUE_REGISTER_2026-07.md).

## What SWA means

Slow-wave activity (SWA) is broadband EEG/iEEG power at approximately 0.5–4 Hz. It is used here as
a sleep-depth feature and Lecci negative-control band. It is not a discrete slow-oscillation event
and is not an LC measurement.

## Implemented corrections

### Signal/cache layer

- Computes band power independently per contact before contact aggregation.
- Stores per-contact full-night power and normalizes once over the complete night.
- Pulls 30 s of context around each 600 s cache chunk and writes only the core, removing internal
  filter/Hilbert boundary artifacts.
- Retains the original finite-sample mask with a 5 s exclusion margin. Interpolated/filled values
  are filter scaffolding and cannot enter power, staging, events, or ECG peaks.
- Uses human SWA 0.5–4 Hz.
- Stores RR and HR at 1/4 Hz after shape-preserving RR interpolation. Gaps longer than 5 s,
  including both endpoint grid samples, remain missing so a 5.1 s gap cannot be reopened later.
- Removes duplicate/near-duplicate R peaks sequentially.
- Stores clean, unthresholded SO candidates and applies amplitude thresholds only inside the
  relevant staged NREM data.
- Rejects IED/artifact-contaminated SO windows.
- Excludes RESPect BIDS channels marked `bad` and requires good iEEG/ECG/EMG/EOG modalities.
- Makes both HUP and RESPect staging proxies fail closed when high-delta or N2-like/N3-like
  mixture partitions are unstable; this is a high-tail algorithm, not evidence for two latent
  physiological states.
- Computes staging spectra only from fully artifact-free contact-epochs, fits every staging model
  on its exact downstream-eligible reference set, and aggregates a fixed full-night contact set
  after within-contact SWA normalization with an 80% per-epoch support requirement.
- Requires at least 80% sigma and RR/HR coverage.
- Adds cache schema, cache-producer source digest, code revision, runtime versions, timestamps,
  source selection/checksums where available, error logs, and atomic writes. Cache acceptance is
  tied to the exact cache-producing code, not a schema label or broad repository hash alone.
- Requires the embedded cache subject to match its filename/request and propagates documented cache
  exclusions as structured skip records. Every downstream result pins the exact cache manifest and
  NPZ SHA-256; summaries reject swapped, changed, or lineage-incomplete inputs.
- Terminal manifests hash every output artifact. Existing results are reusable only when the prior
  complete manifest, requested set, configuration, runtime, and every byte hash match. Strict JSON
  rejects NaN/Infinity tokens, and any ECG detector exception fails the cache.

### 3A

- Fixed 10–15 Hz sigma is primary; unvalidated one-window FSP selection is disabled.
- Uses mandatory 4 s smoothing, four-cycle second-stage Morlets, 0.001–0.12 Hz, and common bout
  support across frequencies.
- Gaussian fits now have objective local-excess/fit/width criteria and can return “no peak.”
- Constant/zero-variance traces and nonpositive normalization scales are unavailable rather than
  counted as all-NaN spectra; coherence is unavailable unless both inputs have positive spectral
  support and the target value is finite.
- Records the Gaussian spectral width and normalized peak-window endpoint separately.
- Evaluates the paired SWA negative control in the accepted sigma peak window rather than selecting
  a different SWA peak/window, but keeps the sigma-selected comparison descriptive because an
  ordinary paired p value would condition on the selection.
- Uses subject-matched scale-free surrogates and a group max-lag sign-flip test.
- Removes invalid independent-frequency-bin and independently selected-lag inference.
- Calls the method a Lecci-aligned approximation, pending validation against FieldTrip/raw data.

### 3B

- Event averages are computed in RR space, matching Naji; HR is used for presentation.
- Takes a peak time separately per channel and averages those times for Naji’s latency endpoint.
- Uses one participant-level magnitude from the average channel RR curve.
- Computes the stage baseline in RR space and converts once to HR, matching the participant curve
  and avoiding `60/mean(RR)` versus `mean(60/RR)` denominator mixing.
- Preserves multichannel event dependence with one shared circular shift in eligible stage-time.
- Uses stage-specific clean SO amplitude thresholds and 1,000 diagnostic shifts.
- Restricts events to the cache's fixed coverage-qualified contact set and requires at least two
  contacts to retain ≥30 complete finite in-stage RR windows after final window eligibility.
- Disables 3B inferential p/z reporting because whole-stage shifts do not preserve local
  nonstationary HR trends or event-density clustering; raw/local magnitudes remain descriptive.

### 3D

- Uses fixed SO 0.16–1.25 Hz and spindle 12–16 Hz bands.
- Applies channel-night thresholds to clean pooled NREM, not separately to each 600 s chunk.
- Pads IED/artifact masks by ±2.5 s.
- Assigns each spindle to its nearest SO and retains at most the maximum-amplitude spindle per SO.
- Retains bounded/calibrated per-contact Rayleigh values as diagnostics only; no contact-significance
  count is used as cohort evidence.
- Uses raw per-contact complex phase means and averages contacts equally into one descriptive
  participant vector for pooled NREM. Participant rotations are not production inference: the
  finite SO-centered pairing window itself aligns directions under independent event trains.
- Requires a pooled descriptive endpoint to have at least three included contacts and 200 paired events from
  those contacts, with at least 1,200 s and 80% valid pooled NREM per included contact; the cohort
  summary revalidates the per-contact evidence rather than trusting a status flag.
- Disables all 3D inference pending a pairing-aware time-shift/block null. The N2-like/N3-like
  comparison additionally needs matched contact sets and event-count-controlled validation.

### Reproducibility/reporting

- Hard-stops the withdrawn 3A/3B/3D generators, follow-ups, figures, combined pipeline, and legacy
  summarizer while leaving only shared helper functions importable where the corrected path needs
  them.
- Ignores all legacy JSON/text results.
- Creates an `in_progress` manifest before work. Caught subject failures are finalized and every
  requested subject is completed, skipped with a reason, or failed; runner `None` values cannot
  disappear. An unexpected process interruption intentionally remains `in_progress` and is rejected.
- Records endpoint availability and a partial-status reason for unavailable required 3A/3B
  endpoints. Available endpoints from a partial record still enter their own endpoint-specific
  denominator; each inferential endpoint has its own minimum sample-size gate. Summaries reject
  non-production pipelines/configurations, byte-hash mismatches, mixed, stale, unclassified,
  lineage-mismatched, failed, or manifest-free runs.
- Describes the HUP selector accurately as a sparse-probe high-delta candidate. It searches the
  complete record, preserves failed probes in physical time, includes the last legal window, and
  requires enough remaining recording for the requested duration.
- Pins the direct scientific dependencies and iEEG client revision; adds CI synthetic tests.
- Redacts the portal account identifier/history from setup documentation.

## Required rerun

```bash
python -m venv .venv
.venv/bin/python -m pip install -r env/requirements.txt

.venv/bin/python analysis/cache_lc_series.py --force
.venv/bin/python analysis/lecci_faithful_3A.py --force
.venv/bin/python analysis/event_3B_cached.py --force
.venv/bin/python analysis/event_3D_by_stage.py --force

.venv/bin/python analysis/stage_ds003848.py --force
.venv/bin/python analysis/run_ds003848_replication.py --force
```

Then run all synthetic checks and
`analysis/test_coherence_calibration.py --require-real-cache`.

## Human/raw-data blockers

- Expert AASM/R&K scoring or validation against expert hypnograms.
- Coordinate/tissue/region/SOZ/bad-contact localization and prespecified homologous regions.
- Blinded ECG, SO, spindle, artifact, missing-data-boundary, and staging validation.
- Benchmarking the 1 Hz Hilbert/Butterworth 3A power series against Lecci’s FieldTrip Morlet
  implementation.
- A pairing-aware null that shifts/block-resamples complete spindle trains relative to SOs and
  repeats event pairing/contact aggregation; no 3D p value is currently valid. A stage contrast
  additionally requires matched contacts and event-count control.
- Verified HUP lights-off/sleep-onset timing; sparse high-delta probes do not establish either.
- A complete current-version real-data rerun with endpoint and manifest QC.
- An independent LC/NE-sensitive measurement or intervention for any LC-specific claim.
