# What the audit fixed

[← Back to the main README](../README.md)

This document records the implementation, integrity, and reporting defects corrected by the
July 2026 audit.

- 3A power is computed per contact before averaging, so phase/polarity cancellation cannot erase
  power.
- Per-contact power is normalized once over the complete selected analysis interval, not
  independently in each 10-minute chunk.
- Streamed cache filtering uses 30 s of context on both sides of every 10-minute chunk and writes
  only the core, preventing internal filter/Hilbert edge artifacts.
- Original missing-sample masks are retained and dilated by 5 s; filled values are numerical filter
  scaffolding only and cannot contribute to power bins, epochs, events, or R peaks.
- Neutral caches retain per-contact power numerators, clean-sample denominators, and observation
  masks. The historical fixed-contact 80% rules remain reproducible as `audit80`; endpoint-local
  overlap-connected aggregation and 70/75/80/90 one-axis sensitivities are evaluated offline.
- Human SWA is 0.5–4 Hz; fixed 10–15 Hz sigma is primary. The SWA negative control is measured in
  the same sigma-defined frequency window, not at an independently selected SWA peak; this
  sigma-selected comparison is descriptive unless selection is repeated inside a joint null.
- The second Morlet analysis uses 0.001–0.12 Hz, four cycles, mandatory four-second smoothing, and
  common bout support across frequency.
- Gaussian peak fitting may return “no peak” and records fit quality, spectral width, and the
  paper-defined peak-window value.
- Constant/zero-variance power and zero-autospectrum coherence are treated as unavailable endpoints,
  not as measured null results.
- RR interpolation is shape-preserving and explicitly masks both endpoints of gaps longer than
  5 s, so a 5.1 s beat gap cannot be reopened downstream as a nominal 5 s NaN run. The Naji event
  curve and its stage denominator use one consistent RR-to-HR definition.
- RESPect validates and uses its author `events.tsv` annotations: coarse NREM/REM/SWS and wake
  context is primary, author-unknown sleep remains unclassified, and artifact/seizure/stimulation
  intervals are excluded before filtering. Proxy labels inside unknown sleep are stored only as a
  sensitivity analysis.
- RESPect joins the selected channels to `electrodes.tsv`, conservatively excludes documented
  pathological/non-cortical contacts, and stores Destrieux ROI masks. The 3A adaptation uses a
  parietal/postcentral/precuneus intersection; 3B uses a frontal intersection.
- All coverage-qualified EMG and EOG channels are normalized and robustly combined rather than
  choosing the first channel. Robust feature centers/scales are fit only on jointly finite clean
  epochs, so discarded data cannot move thresholds applied to retained epochs.
- Staging spectra use complete, artifact-free 4-second Hann/Welch windows at 2-second overlap.
  An outcome-blind RESPect reconstruction calibration selected 11 of 14 iEEG windows as an
  engineering sensitivity; it is not a paper rule or validated HUP threshold. Sigma and SWA share
  one observation graph/contact set, and median-polish aggregation fails closed on nonconvergence.
- 3A reports the signed max-absolute cross-correlation only as a sensitivity and separately tests
  Lecci's source-defined direction: a positive association at 0–15 s with sigma following HR.
  Opposite-sign coupling cannot count as paper-aligned support. The current artifacts do not
  publish a cohort peak-clustering test; the tightness-only calculation survives only in a
  legacy/unshipped summarizer and is not a test that peaks lie near 0.019 Hz.
- 3B thresholds clean SO candidates within stage, averages Naji derivation/contact peak times, and
  uses only uninterrupted stage runs of at least 180 s and cache/ROI-qualified contacts, with at
  least two contacts still contributing 30 complete in-stage RR windows. One shared within-stage
  circular shift preserves cross-channel dependence, but is diagnostic only because it does not
  preserve local nonstationary HR trends/event clustering; 3B event-locking inference is disabled
  pending a valid local null.
- 3D uses fixed 0.16–1.25 Hz SO and 12–16 Hz spindle bands, whole-NREM channel-night thresholds,
  padded IED/artifact rejection, and at most one maximum-amplitude spindle per SO.
- Per-contact Rayleigh values and pooled participant vectors in 3D are descriptive only. The former
  participant-rotation test is disabled because selecting the nearest spindle inside a finite
  SO-centered window creates a shared phase direction even for independent event trains.
- The historical 3D contact/event/coverage rules are represented explicitly and varied offline.
  The base descriptive profile requires three contacts, 200 paired events, and 1,200 valid NREM
  seconds per included contact; these are repository sensitivity choices, not paper-derived
  qualification cliffs.
- Caches store a digest of the cache-producing source and environment specification, not only a
  broad Git/tree label. Cache files also embed and verify participant identity. Results pin the
  exact cache manifest and NPZ hashes; terminal manifests pin every output JSON/NPZ byte hash;
  reuse verifies the prior terminal hashes before touching a manifest. Non-standard NaN/Infinity
  JSON is rejected. Results record explicit endpoint availability and mark a subject partial with
  a reason when a required endpoint is not estimable.
- Every runner creates an `in_progress` manifest before subject work. Caught subject failures are
  finalized explicitly; an unexpected process interruption deliberately leaves the manifest
  `in_progress`, which summaries reject. Missing/`None` results are classified rather than
  disappearing; failed runs are rejected. A partial subject contributes only to endpoints that are
  actually available, with endpoint-specific denominators, reasons, and minimum sample-size gates.
- Only true cache-wide acquisition requirements create cache skips. Endpoint-local power, staging,
  RR, contact, and event support is resolved downstream and recorded per endpoint; acquisition,
  ECG-detector, and runtime failures remain fatal.
- The combined summary validates 3A and 3B independently, so an unavailable 3A endpoint cannot
  suppress a valid 3B report. Within 3A, spectra/peak fitting use shared finite EEG support while
  coherence/cross-correlation additionally require finite HR; missing cardiac samples no longer
  erase otherwise valid spectral evidence.
- The HUP interval selector now preserves failed sparse probes in physical time, scans the complete
  recording, includes the final legal window, and labels its result as a sparse-probe high-delta
  candidate—not a densely measured night or verified sleep onset.
- Every HUP portal analysis pull now requires the exact requested sample/channel geometry and
  records requested-versus-returned counts. All selected cortical and ECG channels must share one
  pinned time base; a short response is a fatal acquisition failure rather than high apparent
  coverage.
- HUP inputs are pinned to the portal dataset snapshot plus every used channel revision and
  `dataCheck` identity; a changed snapshot/channel set fails closed before acquisition. Portal
  requests have bounded 10 s connect and 90 s read timeouts. Four independent sessions score
  sparse night probes concurrently in deterministic probe order, record probe failures, and abort
  the search if more than 20% of probes fail.
- Withdrawn 3A/3B/3D command-line generators, follow-ups, figures, the combined pipeline, and its
  summarizer are hard-stopped.

Synthetic regression tests demonstrate many of the corresponding failure modes and fixes; the
[issue register](ISSUE_REGISTER_2026-07.md) states additional acceptance checks. Neither
substitutes for blinded raw-data QC. Compact RESPect/HUP v9 evidence and the paired scalp–iEEG
outputs have been regenerated and publication-validated. The tracked v8 snapshot is retained only
to make every v8-to-v9 change machine-auditable; old v8 bytes are not silently accepted as current.

See the [unresolved limitations](UNRESOLVED_LIMITATIONS_2026-07.md) for scientific and validation
questions that remain open.
