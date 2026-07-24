# Historical data-inventory leads — not a current usability certification

> **HISTORICAL / UNVERIFIED SNAPSHOT.** The original inventory was assembled from service and
> dataset metadata in July 2026. Availability, access, channel status, subject counts, recording
> duration, and metadata can change. None of the entries below is certified by this document as
> runnable, analysis-ready, correctly staged, or suitable for an LC-proxy conclusion.

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

## Historical candidate leads

Every row is a lead that must be re-queried and revalidated from current raw metadata.

| Source | Historical metadata observation | Current status |
|---|---|---|
| iEEG.org HUP `phaseII` | Earlier queries found multi-day epilepsy recordings with depth contacts and ECG/EKG in a subset. | **Unverified candidate.** Exact subject/channel counts must be rebuilt. HUP lacks an equivalent expert PSG hypnogram in this workflow, so N2-like/N3-like remains a proxy. |
| OpenNeuro `ds003848` (RESPect) | Earlier BIDS tables indicated sleep iEEG with ECG and EMG; EOG availability and channel roles varied by subject. Some respiration rows were marked bad. | **Unverified candidate.** Read every current `channels.tsv`; exclude bad channels; do not describe rule-based labels as “real” or expert staging. Short or fragmented runs may be non-estimable. |
| Other OpenNeuro iEEG datasets | Earlier survey results suggested several sleep-iEEG candidates but no additional verified sleep-iEEG-plus-ECG cohort. | **Historical survey result only.** Re-run the query before making an exclusivity claim. |
| DANDI, DABI, EBRAINS, and indexed archives | Earlier searches found partial modality matches, awake pupil+iEEG, or subcortical disease recordings rather than a complete match. | **Historical survey result only.** Absence was not established permanently and must be rechecked. |

The former subject-by-subject HUP counts and “two independent runnable cohorts” verdict were removed
because this branch has not regenerated a current, version-gated inventory proving raw access,
eligible channels, estimable endpoints, and complete participant accounting.

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

Until that inventory is regenerated, the dataset names above are discovery leads only.
