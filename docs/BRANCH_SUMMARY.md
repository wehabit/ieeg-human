# Archived branch summary — all legacy results withdrawn

> **QUARANTINED DOCUMENT.** Earlier versions of this page described numerical 3A, 3B, and 3D
> results as corrected, validated, negative, weak, or positive. Those descriptions were based on
> superseded caches and estimators and have been removed. Do not recover or cite them as findings.
> No corrected real-data cohort result has been generated.

The authoritative current status is the
[issue register](ISSUE_REGISTER_2026-07.md). The current analysis contract is in
[METHODS.md](METHODS.md), and the distinction between code corrections and unresolved scientific
validation is in [AUDIT_CORRECTIONS_2026-07.md](AUDIT_CORRECTIONS_2026-07.md).

## Historical purpose of the branch

This branch was created to examine three LC-motivated sleep-physiology questions:

| Label | Question | Permitted interpretation |
|---|---|---|
| 3A | Do infraslow sigma-power dynamics covary with RR/heart-rate dynamics? | Candidate sigma–cardiac physiology; not a human LC measurement |
| 3B | Does RR/heart rate change around cortical slow-oscillation down-states? | CNS–autonomic timing; not LC-specific |
| 3D | At what slow-oscillation phase do spindle events occur? | Generic SO–spindle nesting and a physiology/QC endpoint; not an LC proxy |

The papers motivate these questions, but only the animal work directly measures or manipulates LC
or noradrenergic mechanisms. Human iEEG, sigma, SWA, SO events, spindles, RR, and heart rate do not
by themselves establish LC specificity.

## Why the former methods and results are not retained here

The audit identified defects affecting preprocessing, missing-data eligibility, cache lineage,
stage-proxy assignment, cardiac interpolation and normalization, endpoint accounting, and the
inferential unit for 3D. It also identified unresolved anatomy, staging, raw-QC, and construct
validity limitations. The specific failure mechanisms, executable counterexamples, fixes, and open
acceptance criteria are recorded in the [issue register](ISSUE_REGISTER_2026-07.md).

A banner alone would not prevent old search snippets from presenting withdrawn numbers as current,
so the obsolete result tables and method-verification claims were deliberately removed from this
file.

## Current status

- All committed cohort numbers are historical artifacts and must not be interpreted.
- HUP N2-like/N3-like values are algorithmic proxies, not expert-scored sleep stages.
- RESPect rule-based exclusions and proxy labels are not validated expert staging.
- Per-contact Rayleigh results and pooled participant vectors are descriptive only. The former
  rotation null is invalid because the finite pairing window creates phase alignment under
  independent event trains; all 3D inference is disabled pending a pairing-aware null.
- The N2-like versus N3-like 3D comparison additionally needs matched-contact, count-controlled
  participant-level validation.
- A complete version-gated raw-data rerun, blinded detector/raw-boundary QC, and the publication
  acceptance checks in the issue register are still required.

Until those conditions are met, this repository supports a corrected implementation under test,
not a biological result and not validation of heart rate as an LC proxy.
