# Archived 3A method notes — not the current protocol

> **QUARANTINED DOCUMENT.** Earlier versions called one 3A implementation “validated” and presented
> synthetic control numbers as proof of the complete analysis. That language and those numbers have
> been removed. A helper-level test can expose a bug; it cannot validate raw ECG/iEEG quality,
> staging, anatomy, cohort selection, LC specificity, or a biological conclusion.

Use [METHODS.md](METHODS.md) for the current production contract and
[ISSUE_REGISTER_2026-07.md](ISSUE_REGISTER_2026-07.md) for each defect’s evidence and acceptance
criterion. No legacy cohort output referenced by earlier versions of this page is current.

## Historical diagnostic lessons

The following are retained only as engineering cautions. They are not a declaration that the full
current pipeline or its real-data outputs have been validated.

1. **Resolve returned channels by explicit identity.** A data service may return columns in index
   order rather than request order. ECG and neural channels must be mapped and checked explicitly.
2. **Do not erase time when data are missing.** Deleting nonfinite samples and concatenating the
   remainder changes elapsed time and spectral frequency. The current contract preserves the
   original missing mask and excludes a margin around gaps.
3. **Filtering needs chunk context.** Filtering independent cores creates artificial internal
   endpoints. The current cache pulls overlap and writes only the core.
4. **A coherence null must match the statistic.** A constant circular shift can preserve
   magnitude-squared coherence for a stationary oscillatory pair; it is not automatically a valid
   independence surrogate.
5. **Negative controls must use the same selected window.** Sigma and SWA must be compared in the
   same sigma-defined frequency window rather than after separate peak searches.
6. **Synthetic tests have bounded meaning.** They should demonstrate known failure modes, false
   positive behavior, and planted-signal recovery at the actual analysis geometry. They do not
   replace a corrected real-data rerun or blinded raw review.

## Current interpretation guardrail

An isolated spectral or cross-correlation feature is not proof of an LC rhythm. A defensible report
must pass version/cache lineage checks, missing-data and endpoint accounting, participant-level
inference, negative controls, raw detector QC, and the publication gate in the issue register. Even
then, sigma–cardiac association is an LC-motivated candidate signature unless independently
validated against an LC/NE-sensitive measurement or intervention.
