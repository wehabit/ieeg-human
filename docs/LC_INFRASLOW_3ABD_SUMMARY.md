# Archived LC-infraslow 3A/3B/3D summary — legacy results withdrawn

> **QUARANTINED DOCUMENT.** Earlier versions of this page contained current-looking numerical
> results and statements that 3A was negative, 3B was weak or present, 3D was positive, staging was
> “real,” and the analyses were validated or faithfully corrected. Those statements are withdrawn.
> They came from superseded caches and estimators. This archive must not be used for the v9
> availability counts or estimates; see
> [QC_SENSITIVITY_RESULTS_2026-07.md](QC_SENSITIVITY_RESULTS_2026-07.md).

Do not cite an earlier revision of this file for a biological finding. The authoritative record of
what failed, why it is a real concern, what has been fixed in code, and what remains open is the
[LC-proxy audit issue register](ISSUE_REGISTER_2026-07.md). See [METHODS.md](METHODS.md) for the
current production contract and [AUDIT_CORRECTIONS_2026-07.md](AUDIT_CORRECTIONS_2026-07.md) for
rerun requirements.

## Scientific question retained from the legacy summary

The branch asks whether human sleep physiology contains:

- an infraslow association between sigma power and cardiac timing (3A);
- cardiac timing changes around cortical slow-oscillation down-states (3B); and
- slow-oscillation/spindle event nesting (3D).

These are distinct endpoints. They must not be collapsed into one “LC proxy.” The cited human papers
measure sleep electrophysiology and cardiac timing; the direct LC/noradrenergic evidence motivating
the hypothesis comes from animal work. Without an independent LC/NE-sensitive measurement or
intervention, a human association would remain LC-motivated rather than LC-specific.

## Withdrawn claim categories

The old results and method sections were removed, rather than left below a banner, because search
snippets and partial quotations could otherwise present them as current.

| Former claim category | Current status |
|---|---|
| Numerical 3A cohort verdicts, spectral peaks, coherence, cross-correlation, and p-values | **LEGACY VALUES WITHDRAWN.** V9 recovers individual estimates, but the former cohort claims/p-values are not current. |
| Numerical 3B effect sizes, stage contrasts, significant-subject counts, and p-values | **LEGACY VALUES WITHDRAWN.** V9 reports descriptive estimates with inference explicitly disabled. |
| Numerical 3D Rayleigh fractions, vector lengths, participant counts, and p-values | **WITHDRAWN.** Per-contact/event results are diagnostic, not cohort inference. |
| “Positive,” “negative,” “weak,” “present,” “confirmed,” or “reproduced” biological verdicts | **WITHDRAWN.** V9 availability and descriptive direction are not a validated cohort or LC result. |
| “Validated implementation,” “faithful method,” or proof based only on synthetic helper tests | **WITHDRAWN.** Tests establish bounded code behavior, not end-to-end biological validity. |
| HUP N2/N3 labels | **PROXY ONLY.** Use N2-like/N3-like and disclose that they are not expert-scored stages. |
| RESPect “real staging” or validated wake/REM exclusion | **WITHDRAWN.** Author annotations are now primary but are coarse/incomplete and do not provide expert AASM/R&K N2/N3 validation. |
| N2-like versus N3-like 3D inference | **DISABLED / OPEN.** Contact sets and event counts are unmatched and residual finite-sample bias remains. |
| Heart rate, sigma, SWA, SOs, or spindle nesting as human LC-specific measures | **NOT VALIDATED.** Construct validity remains open. |

## Current implementation status is not a result

The working tree contains code-level corrections for the defects marked **FIXED** in the issue
register, including per-contact power computation, full-night normalization, chunk overlap,
same-window SWA control, missing-mask preservation, cache-specific lineage, conservative proxy
staging, RESPect author-annotation and anatomy-sidecar integration, all-good-channel EMG/EOG
aggregation, stable 3-minute 3B stage runs, direction-preserving 3A cross-correlation,
long-RR-gap preservation, a consistent RR/HR denominator, bad-ECG exclusion, endpoint availability,
fail-closed byte-hashed manifests, coverage skips distinct from execution failures, runner return
accounting, sparse-probe night selection, and descriptive participant-level 3D vectors. The former
rotation inference is disabled.

Those corrections establish only that specified implementation defects were addressed. They do not
resolve the **OPEN** limitations: LC construct validity, HUP contact anatomy, expert N2/N3 staging,
scalp-to-iEEG source homology, blinded raw detector and acquisition-gap QC, postictal/IED and
respiratory confounding, cohort generalizability, profile sensitivity, and valid endpoint-specific
inference.

## Acceptance principles retained for current and future summaries

The current replacement results document, and any future rebuild, must:

1. use the analysis and cache versions required by the current README;
2. account for every requested participant in a terminal manifest and reject an interrupted
   `in_progress` run;
3. report endpoint-specific availability and denominators;
4. pass the executable counterexamples and regression tests cited in the issue register;
5. complete blinded raw-data review and document exclusions;
6. use participant-level inference, with per-contact statistics labeled diagnostic;
7. keep all 3D inference disabled until event pairing is repeated inside a calibrated
   time-shift/block null, with an additional matched-contact gate for stage contrasts; and
8. distinguish a null or positive sleep-physiology result from evidence about LC specificity.

V9 materially changes availability but not the construct-validity conclusion: it recovers
subject-level estimates, most consistently SO-linked cardiac acceleration, while 3A is
heterogeneous and 3D remains descriptive. Neither cohort proves human LC tracking.
