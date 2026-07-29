# Data inventory — current analysis inputs, not a clinical usability certification

> **SCOPE.** V9 now provides version-gated analysis inventories for the HUP and RESPect inputs used
> here. That proves the pinned files/channels and participant accounting for this run; it does not
> certify clinical suitability, expert staging/anatomy, future service availability, or an
> LC-specific construct.

Before selecting a cohort, consult the current [issue register](ISSUE_REGISTER_2026-07.md),
especially the open staging, anatomy, raw-QC, and construct-validity items. Dataset availability
does not resolve any of those concerns.

## Signal requirements

These requirements describe candidate data, not validated measurements:

| Test | Minimum candidate signals | Important unmet validation |
|---|---|---|
| 3A | Sleep-period iEEG and usable ECG | Expert stage labels or a validated scope, ECG QC, missing-data QC, sufficient continuous duration |
| 3B | Sleep-period iEEG and usable ECG | Expert/raw SO and R-peak validation, artifact and seizure exclusion |
| 3D | Sleep-period iEEG | Contact localization, expert/raw SO and spindle validation, participant-level inference |
| Pupil extension | Simultaneous sleep pupillometry and iEEG | Direct modality availability and construct validation |
| Disease extension | Relevant disease/control labels and comparable physiology | Cohort comparability and an LC/NE-sensitive validation measure |

## Current inputs and historical candidate leads

HUP and RESPect below are the current pinned inputs for this analysis. The remaining archive rows
are leads that must be re-queried before a new availability or exclusivity claim.

| Source | Historical metadata observation | Current status |
|---|---|---|
| iEEG.org HUP `phaseII` | Multi-day epilepsy recordings with iEEG/ECG; a subset has scalp-labeled channels. | **Current analysis-specific inventory:** 25 requested, 24 completed, HUP116 structured skip. The all-25 scalp audit records labels/geometry and derives eight activity-qualified paired C3 records. HUP staging/anatomy/reference remain unvalidated. |
| OpenNeuro `ds003848` (RESPect) | Sleep iEEG with ECG, EMG, EOG, coarse event annotations, and electrode metadata; respiration rows are marked bad. | **Current analysis-specific inventory:** six completed version-pinned caches. Author annotations and pathology metadata are consumed, but N2/N3 are not expert AASM/R&K stages and short/fragmented support limits estimability. |
| Other OpenNeuro iEEG datasets | Earlier survey results suggested several sleep-iEEG candidates but no additional verified sleep-iEEG-plus-ECG cohort. | **Historical survey result only.** Re-run the query before making an exclusivity claim. |
| DANDI, DABI, EBRAINS, and indexed archives | Earlier searches found partial modality matches, awake pupil+iEEG, or subcortical disease recordings rather than a complete match. | **Historical survey result only.** Absence was not established permanently and must be rechecked. |

The current manifests prove raw inputs and complete participant accounting for this analysis.
They do not prove that either cohort is a clinically representative, anatomically validated, or
generally reusable LC-proxy dataset.

## Required inventory acceptance checks

A future inventory may call a cohort usable only after it records, per participant:

1. immutable dataset and file identifiers;
2. actual raw access and readable duration;
3. channel type, units, sampling rate, reference, and `status`;
4. electrode/contact localization evidence rather than numbering inference;
5. missing/acquisition-gap structure;
6. ECG signal quality and detector-review eligibility;
7. stage-source provenance and whether labels are expert-scored or proxy-derived;
8. required endpoint availability, explicit skip/failure reasons, and complete manifest accounting;
9. adequate continuous eligible data for the estimator; and
10. seizure, IED, movement, medication, and respiratory-confound handling where available.

For HUP and RESPect, the analysis-specific inventory is regenerated; the unmet checks above remain
publication limitations. Other dataset names remain discovery leads only.
