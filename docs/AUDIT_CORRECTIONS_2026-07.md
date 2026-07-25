# Method and code audit — 2026-07

## Status

Legacy numerical cohort results are historical artifacts. Current readers require
`analysis_version = 2026-07-qc-sensitivity-v8` and
`cache_schema_version = 2026-07-neutral-per-contact-gap-aware-source-pin-v8`.
The former v7 fixed-80 result is retained only as historical `audit80`. Under the outcome-blind
endpoint-local base, RESPect has 3/3/3 available 3A records and 0/1/2 available 3B N2/N3/pooled
records; HUP has 19/14/15 available 3A records, 6/5/15 available 3B records, and 7 descriptive 3D
records among 24 completed caches plus one explicit skip. Availability is not hypothesis support;
RESPect remains below the default cohort minimum and all 3B/3D inference remains disabled. See
[QC_SENSITIVITY_RESULTS_2026-07.md](QC_SENSITIVITY_RESULTS_2026-07.md).

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
- Loads, validates, and hashes RESPect `events.tsv`. Author NREM/REM/SWS and outside-sleep context
  are primary; author-unknown sleep remains unclassified, while proxy labels within it are stored
  only for sensitivity. Author artifact/seizure/stimulation intervals are excluded before
  filtering and break stable-stage runs.
- Joins RESPect iEEG channels exactly to `electrodes.tsv`; excludes documented
  pathological/non-cortical contacts; stores Destrieux labels; and constructs conservative
  parietal/postcentral/precuneus and frontal endpoint intersections.
- Uses all coverage-qualified good EMG/EOG channels after within-channel normalization and robust
  aggregation rather than selecting the first channel.
- Makes the HUP staging proxy fail closed when high-delta or N2-like/N3-like mixture partitions are
  unstable; this is a high-tail algorithm, not evidence for two latent physiological states.
- Computes staging spectra from complete artifact-free 4-second windows, fits every staging model
  on its exact downstream-eligible reference set, and retains joint delta/SWA observation support
  for historical fixed-set and endpoint-local overlap-connected aggregation.
- Validates the RESPect BrainVision header rate against the BIDS nominal rate within 10 ppm, stores
  both, and uses the nominal BIDS rate for exact sample/event boundaries.
- Retains reversible support and applies power, RR/HR, staging, contact, and event rules at their
  own endpoints rather than as one cache-wide participant qualification.
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
- Uses subject-matched scale-free surrogates. Reports cohort peak clustering as a tightness-only
  statistic, not evidence that peaks are near Lecci's 0.019 Hz.
- Separates the prespecified positive 0–15 s Lecci-direction cross-correlation test (sigma follows
  HR) from signed max-absolute and opposite-sign sensitivity summaries.
- Removes invalid independent-frequency-bin and independently selected-lag inference.
- Uses the conservative RESPect parietal/postcentral/precuneus intersection while calling the
  method a motivated iEEG adaptation, pending scalp-source and FieldTrip/raw validation.

### 3B

- Event averages are computed in RR space, matching Naji; HR is used for presentation.
- Takes a peak time separately per channel and averages those times for Naji’s latency endpoint.
- Uses one participant-level magnitude from the average channel RR curve.
- Computes the stage baseline in RR space and converts once to HR, matching the participant curve
  and avoiding `60/mean(RR)` versus `mean(60/RR)` denominator mixing.
- Preserves multichannel event dependence with one shared circular shift in eligible stage-time.
- Uses stage-specific clean SO amplitude thresholds and 1,000 diagnostic shifts.
- Requires an uninterrupted 180 s run of the same stage, matching Naji's stable-bin rule.
- Restricts events to the cache's fixed coverage-qualified contact set, plus the RESPect frontal
  intersection when available, and requires at least two contacts to retain ≥30 complete finite
  in-stage RR windows after final window eligibility.
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
- Uses three contacts, 200 paired events, and 1,200 valid pooled-NREM seconds per included contact
  as the base descriptive sensitivity, then varies contact/event/support rules offline; the cohort
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
- Uses atomic, reasoned `status=skip` records only for structural/cache-wide exclusions;
  endpoint-local support failure does not discard unrelated endpoints. Acquisition-chunk,
  ECG-detector, and runtime failures remain fatal.
- Records endpoint availability and a partial-status reason for unavailable required 3A/3B
  endpoints. Available endpoints from a partial record still enter their own endpoint-specific
  denominator; each inferential endpoint has its own minimum sample-size gate. Summaries reject
  non-production pipelines/configurations, byte-hash mismatches, mixed, stale, unclassified,
  lineage-mismatched, failed, or manifest-free runs.
- Validates and reports 3A and 3B independently; failing one analysis's sample-size gate cannot
  suppress a valid endpoint from the other.
- Describes the HUP selector accurately as a sparse-probe high-delta candidate. It searches the
  complete record, preserves failed probes in physical time, includes the last legal window, and
  requires enough remaining recording for the requested duration.
- Pins every HUP dataset snapshot and used-channel revision/`dataCheck` identity and fails closed
  if the portal snapshot or channel set differs.
- Bounds portal requests with 10 s connect and 90 s read timeouts, limits automatic HTTP retries to
  transient idempotent metadata GETs, and closes each subject/worker session.
- Scores HUP sparse night probes concurrently with four independent sessions while preserving
  deterministic probe-time result order. Records probe failures and fails closed when more than
  20% of probes fail.
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

- Expert AASM/R&K N2/N3 scoring or validation against expert hypnograms. RESPect author annotations
  are coarse/incomplete, and strict primary use may leave too few estimable participants.
- HUP coordinate/tissue/region/SOZ/bad-contact localization. RESPect's filtered Destrieux
  parietal/frontal iEEG intersections remain adaptations, not validated homologues of the cited
  scalp sources and not direct LC measurements.
- Blinded ECG, SO, spindle, artifact, missing-data-boundary, and staging validation.
- A clinically justified RESPect postictal sensitivity exclusion; exact seizure intervals are
  masked, but possible longer cardiac/sleep effects are not thereby removed.
- Benchmarking the 1 Hz Hilbert/Butterworth 3A power series against Lecci’s FieldTrip Morlet
  implementation.
- A pairing-aware null that shifts/block-resamples complete spindle trains relative to SOs and
  repeats event pairing/contact aggregation; no 3D p value is currently valid. A stage contrast
  additionally requires matched contacts and event-count control.
- Verified HUP lights-off/sleep-onset timing; sparse high-delta probes do not establish either.
- New or longer appropriately staged, anatomically validated data for robust inference. The v8
  grids recover estimates but show substantial profile/staging sensitivity; do not select a
  profile after seeing the biological direction.
- An independent LC/NE-sensitive measurement or intervention for any LC-specific claim.
