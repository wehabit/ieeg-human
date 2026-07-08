# Full-night SO coupling across HUP subjects (iEEG.org)

Extends the single-subject full-night result (HUP165) to additional continuous full nights pulled
from iEEG.org, via `ieeg_pull_night.py` + `hup_night_cohort.py`. Same polarity-robust, cross-validated
method; each subject's value = mean cross-val modulation over its bipolar MTL channels.

## Result (independent full nights)

| Subject | SO→spindle (mean; ch positive) | SO→ripple (mean; ch positive) |
|---|---|---|
| HUP165 | 0.095 (27/28) | 0.051 (28/28) |
| HUP157 | 0.086 (26/28) | 0.088 (27/28) |
| HUP130 | (downloading) | (downloading) |

**Both full nights replicate:** SO→spindle **and** SO→ripple coupling are positive across ~all MTL
channels (26–28 of 28) in each subject independently. This is the within-night, continuous
counterpart to the cross-subject clip replication.

## How to read it (honest framing)
- With only 2–3 subjects, a cohort-level CI over subjects is **not** meaningful; the evidence here is
  the **per-subject, per-channel consistency** (each night shows the effect in nearly every channel).
- Statistics within a subject are over its channels (not independent — shared reference/anatomy).
- Blind NREM staging (delta-power), SOZ not excluded (see the per-channel positivity as the
  pathology-robustness argument: coupling is in ~all channels, not a few).

## Where it fits
- **Cross-subject replication** (population-level) → the three clip cohorts (`nesting_phase_aligned`:
  SO→spindle in all 3, polarity-robust).
- **Within-night continuous** → these full nights (coupling holds across a whole real night, and
  resolves SO→ripple cleanly).

## Outputs
- `outputs/hup_night_cohort/hup_night_cohort.png` — per-subject (blue) + per-channel (grey) modulation.
- `outputs/hup_night_cohort/cohort_stats.json`.
