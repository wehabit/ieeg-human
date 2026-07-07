# Per-patient N3 SO-triggered nesting — does the grand average hide per-patient mess?

**Why.** A grand average (like `nesting_timecourse`) can look convincing while a few patients carry
it. This module breaks the N3 SO→spindle nesting down **per patient** and **per event**, so the
heterogeneity is visible.

**Data.** Normative iEEG Sleep & Wake Atlas (204 Hz, HUP cohort), cached N3 clips, mesiotemporal
channels. At 204 Hz ripples are not resolvable, so this is the SO→**spindle** nesting; SO→ripple
lives in the 1 kHz cohort (`slow_ripple_coupling`, `nesting_timecourse`).

**Method.** Per patient, pool all cached N3 clips; detect slow-oscillation troughs (0.5–1.25 Hz,
≥1 SD, ≥0.6 s apart); average the SO wave and the z-scored spindle-band (11–16 Hz) envelope in a
±1.5 s window. Up-state metric = mean spindle-z in the post-trough window (0–0.75 s).

## Result (honest)

- **49 patients** with ≥20 SO events.
- **Up-state spindle nesting is positive in 42/49 (86%)**, median up-state spindle-z = **0.057**.
- So it is **broadly consistent, not carried by a few** — but it is **genuinely heterogeneous**:
  7/49 patients go the other way, and the median effect is modest. The sorted heatmap shows a clear
  gradient from strong (top) to weak/absent (bottom), and the event-level panel shows the same: many
  SO events have a post-trough spindle increase, many do not.

This is the correct caveat to state alongside the grand-average and the group MI_z (`p = 0.030`):
the effect is real and majority-present, but modest and patient-variable — which is exactly why the
independent 1 kHz replication (`slow_ripple_coupling`, SO→spindle `p = 6×10⁻⁵`) matters.

## Outputs
- `outputs/nesting_by_patient/per_patient_nesting.png` — one small-multiple panel per patient
  (SO wave + spindle envelope), sorted by up-state nesting.
- `outputs/nesting_by_patient/so_trough_locked_spindle_heatmap.png` — SO-trough-locked spindle
  envelope; rows = patients (top) and rows = individual SO events (bottom), columns = time.
- `outputs/nesting_by_patient/per_patient_upstate.csv` — per-patient up-state spindle-z + event count.
