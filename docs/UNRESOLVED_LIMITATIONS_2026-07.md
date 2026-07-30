# Important unresolved limitations

[← Back to the main README](../README.md)

The implementation defects listed in [What the audit fixed](AUDIT_FIXES_2026-07.md) are corrected.
The following scientific and validation limitations remain **open** and cannot be repaired by
relabeling an output:

- HUP N2-like/N3-like labels are unvalidated slow-wave/GMM proxies, not expert AASM/R&K stages.
  RESPect author annotations are now primary but are coarse and incomplete: explicit NREM/REM
  exists for only part of the cohort, and author-unknown sleep is intentionally not promoted to
  NREM. Strict handling can leave too few estimable participants; the sample-size gate must not be
  lowered to recover a finding.
- RESPect's sidecars now support conservative pathology exclusions and motivated parietal/frontal
  ROI intersections, but those iEEG ROIs are adaptations of Lecci's C3/parietal source and Naji's
  F3/A2 and F4/A1 scalp derivations—not validated homologous sensors and not direct LC
  measurements. HUP still needs coordinates, gray-matter/region labels, bad-contact/SOZ
  exclusions, and clinical review.
- The 1 Hz Hilbert/Butterworth sigma series is an approximation to Lecci’s 0.1 s FieldTrip Morlet
  construction and must be benchmarked against a validated reference implementation.
- The chosen HUP interval is a high-delta candidate interval, not known lights-off/sleep-onset.
- Missing-data exclusion is now explicit in code, but it and the NeuroKit R peaks, SOs, spindles,
  artifacts, and stage labels still need blinded raw-data validation.
- Current ECG preprocessing/detection uses NeuroKit's default automated method without manual
  validation, not Naji's Pan–Tompkins detector plus visual R-wave confirmation. The detector and
  validation deviation is stored in cache metadata and requires blinded validation before a
  paper-fidelity claim.
- Author-marked seizure intervals are now excluded, but RESP0779 and RESP0800 still require a
  clinically justified postictal sensitivity exclusion; masking only ictal samples does not remove
  possible postictal cardiac or sleep effects.
- All current 3D inference is disabled. A valid null must shift/block-resample complete spindle
  trains relative to SOs and repeat pairing, contact aggregation, and the cohort statistic. The
  N2-like versus N3-like comparison additionally needs matched contact sets and event-count control.
- Intracranial voltage polarity/reference is not calibrated to Naji's negative scalp downstate, and
  a sign flip rotates SO phase by π. Until polarity is anatomically/physiologically oriented or a
  validated sign-invariant sensitivity analysis is prespecified, 3B timing and 3D preferred phase
  are not source-comparable.
- The v9 rerun materially differs from the original fixed-80 analysis: RESPect contributes three
  3A records but remains below the cohort minimum; HUP contributes many subject-level estimates.
  The fixed-80 reference still produces mostly/all unavailable endpoints. This profile
  sensitivity must be reported and cannot be converted into evidence for or against human LC
  tracking.
- A nonsignificant result means “not detected,” not “absent,” unless an equivalence analysis
  excludes a prespecified scientifically meaningful effect.

These limitations define work that remains necessary before any direct human LC-tracking or
paper-fidelity claim.
