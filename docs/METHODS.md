# Audited methods — 3A, 3B, and 3D

## Status and construct

These are methods for a pending corrected rerun. No committed cohort output currently satisfies the
required analysis/cache versions.

3A is an LC-motivated candidate signature. 3B is a cortical–autonomic timing measure. 3D is generic
SO–spindle nesting. None is a validated human LC measurement. Direct LC/NE validation would require
an independent LC/NE-sensitive signal or intervention.

The current production contract is `analysis_version = 2026-07-corrected-v5` and
`cache_schema_version = 2026-07-per-contact-power-pchip-night-threshold-v5`.

## Data

- HUP phaseII: multi-day clinical iEEG and ECG, no scalp PSG/EOG/EMG.
- OpenNeuro ds003848 RESPect: iEEG, ECG, EMG, and EOG; approximately one-hour sleep runs.

HUP contacts are currently lateral-contact candidates selected from contact numbering. That is not
anatomic validation. Coordinates, gray-matter/region labels, bad-contact/SOZ exclusions, and
clinical review are required before publication interpretation.

The HUP interval is selected from 6 s delta-ratio probes sampled every 30 min. Candidate 3 h windows
are evaluated in their true physical positions across the complete record: failed probes stay
missing, the final legal window is included, and the selected start must leave enough recording for
the requested analysis duration. This is a **sparse-probe high-delta candidate interval**, not a
dense whole-window delta estimate, verified night, lights-off time, or sleep onset.

## Staging

- HUP uses an unvalidated delta-ratio/SWA Gaussian-mixture proxy. A split is accepted only with
  BIC and held-out improvement, stable component means, adequate separation, and nondegenerate
  component sizes. Temporal smoothing cannot re-admit explicitly dirty epochs. If an N2-like/N3-like
  split is unsupported, epochs remain pooled `NREM` rather than being forced into two stages.
- RESPect uses an unvalidated EMG/EOG/iEEG rule-based proxy. Production caching requires good
  iEEG, ECG, EMG, and EOG channels and excludes BIDS channels marked `bad`. Wake/REM rules only
  exclude candidates; NREM additionally requires a reproducible high-delta component. An unstable
  NREM or N2-like/N3-like split fails closed.

RESPect robust SWA/EMG/EOG centers and scales are estimated only from jointly finite clean epochs;
dirty epochs cannot move the thresholds applied to retained data.

Labels must be called N2-like/N3-like. Neither cohort has expert AASM/R&K labels in this pipeline.
Pooled proxy-NREM is primary; stage contrasts are exploratory.

SWA means broadband 0.5–4 Hz slow-wave activity. It is not an individual SO event.

## Versioned cache and preprocessing

Each current cache carries a schema version, timestamp, code revision, runtime versions, source
identifier/checksum where possible, selected interval, error logs, and coverage. It also carries a
SHA-256 digest of the exact cache-producing source and environment-specification files. Loaders
compare that cache-specific digest and recorded runtime versions with the current environment; a
broad Git revision or downstream source-tree digest alone is insufficient. Atomic replacement
prevents truncated caches. At least 80% sigma and cardiac coverage is required; acquisition-chunk
failures or any ECG detector exception invalidate a production cache.

Every cache also embeds its participant ID, which must match the requested filename. A downstream
subject result pins the exact cache-manifest run ID/hash and NPZ SHA-256; publication summaries
re-hash those inputs. A documented cache exclusion is a subject-identified skip artifact whose
reason/status/schema/producer digest must agree with its manifest, then propagates as a structured
downstream skip rather than a missing record.

Every terminal manifest also stores the SHA-256 of every output NPZ/JSON. Existing artifacts are
reused only when their prior complete manifest, exact requested set, configuration, runtime, and
byte hashes all match before any new manifest is written. Publication summaries re-hash result
JSON as well as cache inputs and reject non-production configurations. JSON output is strict:
non-standard NaN/Infinity tokens fail atomically.

For each contact:

1. Pull each 600 s core with 30 s of filter/Hilbert context on either side; write only the core.
2. Preserve the original finite-sample mask. Filling is permitted only as numerical filter
   scaffolding; missing samples and a 5 s dilation around them remain ineligible.
3. Detrend and line-notch the raw signal.
4. Detect/pad broadband IED and extreme-amplitude artifacts.
5. Compute band power per contact; never average raw voltages before power.
6. Bin clean, measured envelope-squared power to 1 Hz.
7. Normalize each contact once over the full night, then average contact powers.

The aggregate uses a fixed coverage-qualified contact set: each selected contact must cover at
least 80% of the interval, at least three contacts must qualify, and a reported second requires at
least 80% of that selected set. Per-contact coverage, the selected mask, and the time-resolved
contact count are stored; mutually disjoint low-coverage contacts cannot manufacture a
fully-covered aggregate.

The same measured-data mask excludes contaminated staging epochs, SO/spindle candidates, and ECG
peaks. A contact-epoch contributes delta ratio/SWA only when every sample is measured and
artifact-clean. Staging fixes one full-night contact set, divides SWA by each contact's full-night
clean-epoch median, and requires at least 80% of that exact set (and at least three contacts) in an
epoch. Staging models are fit only on epochs eligible for their downstream labels. This closes
zero-fill, brief-artifact, changing-contact, and excluded-reference leaks; blinded validation
against real data remains required.

ECG R peaks are deduplicated with a sequential refractory rule. Physiological RR intervals
0.33–1.5 s are interpolated with PCHIP inside continuous runs only. For a gap longer than 5 s, both
endpoint grid samples and everything between them are explicitly missing; this prevents a 5.1 s
beat gap from later appearing as a fillable 5.0 s internal NaN run. Both RR and derived HR are
stored at 1 and 4 Hz.

The current cache uses NeuroKit2 `ecg_clean`/`ecg_peaks` with `method="neurokit"` and artifact
correction. This is a documented deviation from Naji's 0.5–100 Hz preprocessing, Pan–Tompkins
detector, and visual confirmation. Cache metadata records the algorithm and `visual_validation =
false`; blinded manual validation and detector sensitivity analysis remain required.

## 3A: infraslow sigma and cardiac dynamics

Primary sigma is fixed 10–15 Hz; SWA 0.5–4 Hz is the negative-control band. Individual FSP analysis
is disabled until FSP can be estimated from all artifact-free NREM and manually quality-controlled.

For every proxy-NREM run at least 120 s:

1. Apply the mandatory symmetric 4 s smoothing.
2. Compute a four-cycle Morlet spectrum from 0.001–0.12 Hz at 0.001 Hz spacing.
3. Let every accepted bout contribute at every frequency using explicit support correction.
4. Average bout spectra weighted by duration and normalize the participant spectrum to its mean.
   Constant/zero-variance bouts and nonpositive or nonfinite normalization scales are unavailable,
   not measured zero spectra.
5. Fit three Gaussian terms. Accept a peak only when a local aperiodic-background excess, fit
   quality, location agreement, and spectral-width criteria pass. Otherwise report no peak.
6. Record fitted location, within-spectrum SD, and normalized power averaged over
   peak ±0.5 spectral SD.
7. Evaluate the SWA negative control in that **same sigma-defined window**. An independently
   selected SWA peak/window is not the paired control. The same-window values are descriptive:
   because the window was selected for high/accepted sigma, an ordinary sigma>SWA p value is
   circular unless peak selection is repeated inside a joint null.

Scale-free surrogates preserve each participant’s bout/gap pattern. Cohort peak clustering uses one
matched surrogate peak per participant per iteration.

For cardiac coupling, 120 s intervals are z-scored and cross-correlated with HR as the source wave.
The group correlogram is tested with a sign-flip maximum statistic across lags. Gap-aware coherence
at 0.02 Hz and at an accepted personal peak is a secondary analysis. Coherence is unavailable when
either autospectrum lacks positive non-DC support or the target-bin value is nonfinite.

This remains an approximation: the cached 1 Hz Hilbert/Butterworth series is not Lecci’s 0.1 s
FieldTrip Morlet series and must be benchmarked on identical raw input.

## 3B: SO–RR/HR timing

Clean per-contact signals are zero-phase filtered at 0.15–4 Hz. Complete negative/down and
positive/up half-waves must each last 0.3–1.0 s. After staging, up-state and peak-to-peak amplitude
thresholds are computed separately within each channel and stage at the 75th percentile.

For each channel, the 4 Hz RR tachogram is averaged in a ±5 s window around SO down-state troughs.
The post-trough RR minimum defines the HR-burst time. Channel-specific times are averaged, matching
Naji’s electrode timing endpoint. A participant magnitude is obtained from the average channel RR
curve. Its stage baseline is calculated in the same RR domain and converted once to HR, so the
percentage does not mix `60 / mean(RR)` with `mean(60 / RR)`. Observed and surrogate effects use
that identical denominator. This is a prespecified consistency choice because the paper does not
fully disambiguate the order of averaging and HR conversion for its stage baseline.

The former null applies one common circular shift in eligible stage-time to the complete
multichannel event ensemble, preserving cross-channel SO synchrony and fixing the older
stage-offset/independent-channel defects. It remains a **diagnostic only**: a whole-stage shift can
destroy local slow HR trends and SO-density clustering. The output therefore exposes no
event-locking p/z claim. It reports the raw Naji magnitude, local post-peak versus pre-event mean,
and shift-null diagnostic while a local-trend/dependence-preserving null is developed and
validated.

This is still an iEEG adaptation: Naji used visually scored stable sleep, F3/F4 scalp, absolute
Dang-Vu criteria, visually checked R peaks, and a behavioral timing endpoint that is absent here.
The production adaptation uses only the cache's fixed coverage-qualified contacts and requires at
least two contacts to retain 30 complete finite in-stage RR windows after final window eligibility.
Intracranial voltage polarity is not yet oriented to Naji's negative scalp downstate, so downstate
timing remains a source-transfer limitation rather than a validated homologous marker.

## 3D: independent-SO spindle phase

Per contact:

- SO detection: 0.16–1.25 Hz, complete 0.8–2.0 s cycles, top-quartile amplitude.
- Spindle detection: fixed 12–16 Hz, 200 ms RMS, 75th-percentile threshold, 0.5–3.0 s duration.
- Thresholds use one artifact-free pooled-NREM channel-night distribution.
- IED/artifact samples are expanded by ±2.5 s; SO troughs and spindle RMS/event samples inside
  that expanded mask are ineligible.
- Each spindle is assigned to its nearest SO; at most the maximum-amplitude spindle is retained for
  one SO epoch.

Per-contact phase nonuniformity and its finite-sample Rayleigh probability are retained only as
diagnostics under an independence assumption. They are not used to declare significant contacts or
to count participant evidence.

For the pooled-NREM descriptive endpoint, each contact contributes its raw complex phase mean and
contacts are averaged with equal weight within participant. The usual algebraic event-count
correction assumes independent event phases and is retained only as a diagnostic. Participant
rotations are also diagnostic only: selecting a spindle inside a finite SO-centered window induces
a common SO-phase direction even when spindle and SO trains are independent. Production inference
is disabled. A valid null must shift or block-resample complete spindle trains relative to SOs and
repeat event selection, pairing, contact aggregation, and the cohort statistic.

Before that pooled endpoint is admitted, at least three included contacts must each provide at least
1,200 s and 80% valid pooled NREM, and those exact included contacts must contribute at least 200
paired SO–spindle events. Events from contacts excluded by the per-contact ≥20-event rule do not
inflate that total. The cohort summary reconstructs this gate from the per-contact records.

All 3D inference is **disabled/open**. The N2-like versus N3-like contrast has additional problems:
separate stage estimates can contain different contact sets and very different event counts. Raw
vector magnitudes remain count- and dependence-sensitive, while the algebraic correction is invalid
under serial dependence. A valid contrast must use matched contacts and a dependence-preserving
within-participant stage/block null before it can be reported.

## Reproducibility rules

- Current outputs require the shared versions in `analysis/pipeline_version.py`.
- A run creates `RUN_MANIFEST.json` as `in_progress` before subject work. Caught subject-level
  failures are recorded in the terminal manifest; an unexpected process interruption intentionally
  leaves `run_state = in_progress`, which is visible and rejected by summaries.
- Required-subject failures cause a nonzero exit.
- Every requested subject is classified as completed, explicitly skipped with a reason, or failed.
  A runner returning `None` is a classified failure/skip, never an omitted participant.
- Subject outputs carry explicit endpoint-availability fields; when any prespecified endpoint is
  unavailable, the subject is marked `partial` with endpoint-specific reasons. Its other valid
  endpoints remain eligible only for their own denominators.
- Summaries reject missing manifests, failures, unexpected subjects, non-production pipelines or
  configurations, result-byte hash mismatches, mixed schema/code-lineage versions, unclassified
  endpoints, and legacy files.
- Publication summaries apply a minimum of five estimable participants separately to each
  inferential endpoint/stage by default; a correctly accounted all-unavailable run cannot exit
  successfully as an `n=0` cohort result.
- Withdrawn 3A/3B/3D command-line generators, figures, follow-ups, and combined scripts are
  hard-stopped; corrected modules may still import explicitly retained helper functions.

## Open scientific and validation limitations

The implementation defects above are fixed in the working tree. Expert/validated staging,
electrode anatomy and pathology review, blinded ECG/SO/spindle/artifact and missing-boundary QC,
the Lecci reference-implementation benchmark, verified HUP sleep timing, and a corrected real-data
rerun remain open. None of those limitations is resolved by a passing synthetic test. All 3D
inference remains disabled pending a pairing-aware null; stage contrasts have the additional
matched-contact/count-control requirement above.

See [ISSUE_REGISTER_2026-07.md](ISSUE_REGISTER_2026-07.md) for evidence and acceptance criteria.
