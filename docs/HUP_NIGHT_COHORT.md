# Full-night SO coupling across HUP subjects (iEEG.org)

> **WITHDRAWN LEGACY RESULT — DO NOT CITE AS CURRENT LC-PROXY EVIDENCE.** These intervals were
> not verified nights, contacts were not independent participants, SOZ/pathology was not excluded,
> and the estimates do not pass the current v9 provenance, raw-QC, or inference gates. See
> [ISSUE_REGISTER_2026-07.md](ISSUE_REGISTER_2026-07.md).
>
> **CURRENT V9 STATUS:** the neutral cache accounts for all 25 requested participants as 24
> completed plus one structured skip. Endpoint-local estimates are available, but HUP staging and
> anatomy remain unvalidated and 3B/3D inference is disabled. See the v9 report; the legacy values
> below remain withdrawn and do not prove LC tracking.

Extends the single-subject full-night result (HUP165) to additional continuous full nights pulled
from iEEG.org, via `ieeg_pull_night.py` + `hup_night_cohort.py`. Same polarity-robust, cross-validated
method; each subject's value = mean cross-val modulation over its bipolar MTL channels.

## Result (independent full nights)

| Subject | SO→spindle (mean; ch positive) | SO→ripple (mean; ch positive) |
|---|---|---|
| HUP165 | 0.095 (27/28) | 0.051 (28/28) |
| HUP157 | 0.086 (26/28) | 0.088 (27/28) |
| HUP130 | (downloading) | (downloading) |

**Withdrawn interpretation:** these legacy channel-level values were previously described as a
replication, but that conclusion is not supported. The intervals were not verified nights, the
channels are not independent biological replicates, and the analysis omitted the current
participant-level, pathology, provenance, and raw-QC requirements.

## How to read it

- Do not use these values as LC-proxy evidence or as a replication result.
- Statistics over channels are pseudoreplicated because contacts share a participant, reference,
  anatomy, and pathology.
- Blind delta-power staging and failure to exclude SOZ/pathological contacts are unresolved
  confounds, not evidence of robustness.

## Where it fits

It is retained only as a historical artifact. It does not enter the corrected 3A/3B/3D endpoint
summaries.

## Outputs
- `outputs/hup_night_cohort/hup_night_cohort.png` — per-subject (blue) + per-channel (grey) modulation.
- `outputs/hup_night_cohort/cohort_stats.json`.
