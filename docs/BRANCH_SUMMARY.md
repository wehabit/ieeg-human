# Branch summary — corrected cohort endpoints unavailable

> **QUARANTINED DOCUMENT.** Earlier versions of this page described numerical 3A, 3B, and 3D
> results as corrected, validated, negative, weak, or positive. Those descriptions were based on
> superseded caches and estimators and have been removed. Do not recover or cite them as findings.
> The v7 RESPect rerun is now complete, but none of its prespecified 3A/3B endpoints is estimable.

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

- Legacy cohort numbers are historical artifacts and must not be interpreted.
- The complete six-subject v7 RESPect cache and downstream manifests passed lineage checks. 3A
  spectrum, cross-correlation, fixed-0.02-Hz coherence, and own-peak coherence were each estimable
  in 0/6; 3B N2 and N3 were each estimable in 0/6. This is an endpoint-unavailable result, not a
  null effect and not evidence for or against LC tracking.
- The v7 HUP cache accounted for 25/25 requested participants: 17 completed, 8 structured skips,
  and 0 failures. All terminal hashes validated. Every 3A endpoint and both 3B stages were
  estimable in 0/25; the direct-stream 3D pooled endpoint was also estimable in 0/25, with
  inference disabled.
- HUP N2-like/N3-like values are algorithmic proxies, not expert-scored sleep stages.
- RESPect author annotations are now primary but remain coarse/incomplete and do not provide
  expert AASM/R&K N2/N3 scoring; author-unknown sleep stays unclassified.
- RESPect sidecars now support conservative pathology exclusions and parietal/frontal Destrieux
  intersections, but those iEEG ROIs are adaptations of the cited scalp sources—not validated
  homologues and not LC measurements. HUP anatomy remains unverified.
- Exact author-marked seizure intervals are excluded; possible longer postictal cardiac/sleep
  effects remain an open sensitivity analysis.
- Per-contact Rayleigh results and pooled participant vectors are descriptive only. The former
  rotation null is invalid because the finite pairing window creates phase alignment under
  independent event trains; all 3D inference is disabled pending a pairing-aware null.
- The N2-like versus N3-like 3D comparison additionally needs matched-contact, count-controlled
  participant-level validation.
- Blinded detector/raw-boundary QC and the remaining publication acceptance checks in the issue
  register are still required.

The corrected reruns materially withdraw the former numerical conclusions because none survives as
an estimable current endpoint. They do not show that the biological effects are absent, and they do
not validate or refute heart rate as an LC proxy.
