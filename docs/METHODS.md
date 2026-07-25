# Audited methods — 3A, 3B, and 3D

## Status and construct

These are the methods for the neutral-cache v8 QC-sensitivity analyses. The prior fixed-80
production rerun remains reproducible as a historical profile, but it is not treated as a
paper-derived or preregistered qualification rule. The authoritative numerical report is
[QC_SENSITIVITY_RESULTS_2026-07.md](QC_SENSITIVITY_RESULTS_2026-07.md).

3A is an LC-motivated candidate signature. 3B is a cortical–autonomic timing measure. 3D is generic
SO–spindle nesting. None is a validated human LC measurement. Direct LC/NE validation would require
an independent LC/NE-sensitive signal or intervention.

The current contract is `analysis_version = 2026-07-qc-sensitivity-v8` and
`cache_schema_version = 2026-07-neutral-per-contact-gap-aware-source-pin-v8`.

The complete six-subject RESPect cache and 25-request HUP cache pass exact cache-code,
input-identity, result-file, and runtime checks. HUP has 24 completed neutral caches and one
structured skip (HUP116, no cortical-contact candidates). Under the outcome-blind endpoint-local
base sensitivity, RESPect has 3/3/3 available 3A spectrum/coherence/cross-correlation records and
0/1/2 available 3B N2/N3/pooled records. HUP has 19/14/15 available 3A records, 6/5/15 available
3B records, and 7 descriptive 3D records. These are availability counts, not effect support.
RESPect remains below the default cohort minimum; 3B and 3D inference is disabled.

## Data

- HUP phaseII: multi-day clinical iEEG and ECG, no scalp PSG/EOG/EMG.
- OpenNeuro ds003848 RESPect: iEEG, ECG, EMG, and EOG; approximately one-hour sleep runs.

HUP contacts are currently lateral-contact candidates selected from contact numbering. That is not
anatomic validation. Coordinates, gray-matter/region labels, bad-contact/SOZ exclusions, and
clinical review are required before publication interpretation.

For RESPect, good iEEG channels are joined by exact name to the subject `electrodes.tsv`. Contacts
documented as SOZ, resected, edge, silicon/screw, CSF/white matter, lesion/gliosis, non-gray depth,
or without a cortical atlas label are excluded. The remaining Destrieux labels define a
parietal/postcentral/precuneus intersection for the 3A adaptation and a frontal intersection for
3B. These are conservative metadata filters and motivated iEEG adaptations; they do not validate
homology to Lecci's C3/parietal or Naji's F3/F4 scalp sensors.

The HUP interval is selected from 6 s delta-ratio probes sampled every 30 min. Candidate 3 h windows
are evaluated in their true physical positions across the complete record: failed probes stay
missing, the final legal window is included, and the selected start must leave enough recording for
the requested analysis duration. This is a **sparse-probe high-delta candidate interval**, not a
dense whole-window delta estimate, verified night, lights-off time, or sleep onset.

The HUP source manifest pins each dataset snapshot plus the revision ID and `dataCheck` identity of
every used lateral-contact candidate and cardiac channel. Cache construction verifies that identity
before selecting an interval. Portal requests use a 10 s connect timeout and 90 s read timeout;
transient idempotent metadata GETs are retried, while data pulls retain their explicit bounded retry
path. Four independent portal sessions score night probes concurrently. Ordered result collection
makes selection deterministic with respect to probe time rather than worker completion order, and
worker snapshot IDs must match the parent snapshot. The cache records the probe count, finite-score
count, failure count, and failure times; the search fails closed if more than 20% of probes fail.

## Staging

- HUP uses an unvalidated delta-ratio/SWA Gaussian-mixture proxy. A split is accepted only with
  BIC and held-out improvement, stable component means, adequate separation, and nondegenerate
  component sizes. Temporal smoothing cannot re-admit explicitly dirty epochs. If an N2-like/N3-like
  split is unsupported, epochs remain pooled `NREM` rather than being forced into two stages.
- RESPect author `events.tsv` sleep annotations are primary. Epochs covered for at least 29 of
  30 seconds by author NREM, REM, or curated SWS are labeled `NREM`, `R`, or `N3`; epochs outside
  sleep/transition are `W`;
  transition/boundary-conflict epochs remain unclassified. Author-unknown sleep also remains
  unclassified in the primary analysis. The EMG/EOG/iEEG proxy inside author-unknown sleep is
  retained only as a sensitivity analysis, not silently promoted to primary NREM.

RESPect requires good iEEG, ECG, EMG, and EOG modalities and excludes BIDS channels marked `bad`.
Features from all coverage-qualified good EMG/EOG channels are normalized within channel and
robustly aggregated. SWA/EMG/EOG centers and scales are estimated only from jointly finite clean
epochs; dirty epochs cannot move the thresholds applied to retained data.

HUP mixture labels and any RESPect proxy-sensitivity subdivision must be called N2-like/N3-like.
Neither cohort has expert AASM/R&K N2/N3 labels in this pipeline. Strict author annotations can
leave too few RESPect participants estimable; this is an honest endpoint-availability result, not a
reason to substitute proxy labels or lower the cohort gate.

SWA means broadband 0.5–4 Hz slow-wave activity. It is not an individual SO event.

## Versioned cache and preprocessing

Each current neutral cache carries a schema version, timestamp, code revision, runtime versions, source
identifier/checksum where possible, selected interval, error logs, and coverage. It also carries a
SHA-256 digest of the exact cache-producing source and environment-specification files. Loaders
compare that cache-specific digest and recorded runtime versions with the current environment; a
broad Git revision or downstream source-tree digest alone is insufficient. Atomic replacement
prevents truncated caches. Power, staging, cardiac, contact, and event support is retained in
reversible arrays and enforced at the endpoint rather than as one global participant
disqualification. Acquisition-chunk failures or any ECG detector exception still invalidate a
production cache.

RESPect caches pin and hash the raw signal plus `channels.tsv`, `events.tsv`, and `electrodes.tsv`.
Event onset/duration/offset/sample fields are cross-checked at the nominal BIDS sampling rate.
Author artifact intervals are applied to their named/all contacts; seizure and stimulation
intervals are global. Exact intervals are set missing before filtering, after which the standard
5 s missing-data dilation applies. The BrainVision header rate and BIDS nominal rate must agree
within 10 ppm; both are stored, and the nominal BIDS rate is used for exact boundary alignment.

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
6. Bin clean, measured envelope-squared power to 1 Hz while retaining the numerator and
   clean-sample denominator.
7. Normalize each contact once over the full night, then aggregate contact powers offline.

The historical `audit80` profile reconstructs its fixed coverage-qualified contact set exactly.
The outcome-blind endpoint-local profile instead finds an overlap-connected observation component,
requires the profile's minimum contacts, and uses median polish to remove stable contact gain.
Disconnected contact/time islands cannot manufacture a full series, and a nonconverged fit fails
closed. Sigma and SWA use the same joint support.

The same measured-data mask excludes contaminated staging windows, SO/spindle candidates, and ECG
peaks. Staging retains fourteen possible complete 4-second Hann/Welch windows at 2-second overlap
for every contact/epoch. The base endpoint-local sensitivity requires 11 iEEG windows, chosen by an
outcome-blind RESPect reconstruction calibration; the 1–14 grid is reported. This is not a paper
rule, stage validation, or calibrated HUP threshold. EMG/EOG retain separate complete-window
support and require 14 in the base profile. Staging models are fit only on epochs eligible for
their downstream labels. Profile outputs serialize auxiliary, proxy-label, final-label, and
disagreement counts.

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

Primary sigma is fixed 10–15 Hz; SWA 0.5–4 Hz is the negative-control band. For RESPect, both use
the coverage-qualified non-pathological Destrieux parietal/postcentral/precuneus intersection.
Individual FSP analysis is disabled until FSP can be estimated from all artifact-free NREM and
manually quality-controlled.

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
matched surrogate peak per participant per iteration. This statistic tests only whether accepted
peaks are unusually tight; it does not test whether their location is compatible with Lecci's
0.019 Hz. Any accepted-subset location interval is explicitly conditional.

For cardiac coupling, 120 s intervals are z-scored and cross-correlated with HR as the source wave.
The primary paper-direction summary is the maximum positive group correlation at 0–15 s, where
positive lag means sigma follows HR; it is tested with a one-sided sign-flip maximum statistic in
that window. Signed max-absolute correlation across all lags is an omnibus sensitivity, and the
opposite-sign extremum is reported separately so a reversal cannot count as paper-aligned support.
Gap-aware coherence at 0.02 Hz and at an accepted personal peak is a secondary analysis. Coherence
is unavailable when either autospectrum lacks positive non-DC support or the target-bin value is
nonfinite.

This remains an approximation: the cached 1 Hz Hilbert/Butterworth series is not Lecci’s 0.1 s
FieldTrip Morlet series and must be benchmarked on identical raw input.

## 3B: SO–RR/HR timing

Clean per-contact signals are zero-phase filtered at 0.15–4 Hz. Complete negative/down and
positive/up half-waves follow the cited Dang-Vu timing: negative-to-positive 0.3–1.5 s, followed by
a positive-to-negative crossing in less than 1.0 s. Only epochs inside uninterrupted runs of the
same stage lasting at least six 30 s epochs (180 s) are eligible; author disturbance intervals
break a run. Dang-Vu's absolute scalp-amplitude cutoffs cannot be transferred directly to iEEG, so
up-state and peak-to-peak amplitude percentiles are explicit adaptations tested at 0/50/75/90,
with 75 as the base sensitivity.

For each channel, the 4 Hz RR tachogram is averaged in a ±5 s window around SO down-state troughs.
The post-trough RR minimum defines the HR-burst time. Channel-specific times are averaged, matching
Naji’s electrode timing endpoint. A participant magnitude is obtained from the average channel RR
curve. Its stage baseline is calculated in the same RR domain and converted once to HR, so the
percentage does not mix `60 / mean(RR)` with `mean(60 / RR)`. Observed and diagnostic-shift effects
use that identical denominator. This is a documented consistency choice because the paper does not
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
For RESPect, the production adaptation intersects the cache's fixed coverage-qualified contacts
with the conservative Destrieux frontal ROI and requires at least two contacts to retain 30
complete finite in-stage RR windows after final window eligibility. Intracranial voltage polarity
is not yet oriented to Naji's negative scalp downstate, so downstate timing remains a
source-transfer limitation rather than a validated homologous marker.

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

The base descriptive sensitivity requires at least three included contacts, at least 1,200 valid
pooled-NREM seconds per included contact, and at least 200 paired SO–spindle events from those exact
contacts. Events from contacts excluded by the per-contact ≥20-event rule do not inflate that
total. Acquisition, event-valid, contact, duration, and event-count choices are stored and varied
offline. They are repository QC sensitivities, not thresholds supplied by Staresina or Helfrich.

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
- A cache-wide acquisition impossibility is written atomically as a subject-identified
  `status=skip` cache with measured support and reason. Endpoint-local support failure remains an
  endpoint-specific unavailable record rather than disqualifying unrelated analyses.
- Subject outputs carry explicit endpoint-availability fields; when any locked endpoint is
  unavailable, the subject is marked `partial` with endpoint-specific reasons. Its other valid
  endpoints remain eligible only for their own denominators.
- Summaries reject missing manifests, failures, unexpected subjects, non-production pipelines or
  configurations, result-byte hash mismatches, mixed schema/code-lineage versions, unclassified
  endpoints, and legacy files.
- Publication summaries apply a minimum of five estimable participants separately to each
  inferential endpoint/stage by default; a correctly accounted all-unavailable run cannot exit
  successfully as an `n=0` cohort result.
- The combined 3A/3B summary validates and reports the two analyses independently. Failure of the
  3A sample-size gate cannot prevent inspection of a valid 3B endpoint, and vice versa.
- Withdrawn 3A/3B/3D command-line generators, figures, follow-ups, and combined scripts are
  hard-stopped; corrected modules may still import explicitly retained helper functions.

## Open scientific and validation limitations

The implementation defects above are fixed in the working tree. Expert/validated N2/N3 staging,
HUP anatomy/pathology review, blinded ECG/SO/spindle/artifact and missing-boundary QC, the Lecci
reference-implementation benchmark, and verified HUP sleep timing remain open. The v8 endpoint
counts are highly profile-sensitive: endpoint-local materialization recovers many HUP estimates
and three RESPect 3A estimates, whereas the historical fixed-80 stack does not. That sensitivity
must be reported rather than bypassed or described as a biological finding. The
author-marked seizure intervals are excluded, but a clinically justified postictal sensitivity
exclusion remains necessary. The RESPect parietal/frontal iEEG ROIs remain scalp-source
adaptations, not direct LC measurements. None of those limitations is resolved by a passing
synthetic test. All 3D inference remains disabled pending a pairing-aware null; stage contrasts
have the additional matched-contact/count-control requirement above.

See [ISSUE_REGISTER_2026-07.md](ISSUE_REGISTER_2026-07.md) for evidence and acceptance criteria.
