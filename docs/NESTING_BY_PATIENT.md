# Per-patient N3 SO-triggered nesting — does the grand average hide per-patient mess?

> **SEPARATE HISTORICAL STUDY.** These numbers are not outputs of the current LC-proxy v7
> 3A/3B/3D pipelines and must not be used as corrected LC-tracking evidence. Their original
> study-specific validation/inference claims require separate review.

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

## De-circularized validation (`nesting_validation`)

The sorted heatmap above **double-dips** (rows sorted by the same spindle value that is
colour-plotted), so its gradient is *display only, not evidence*. `nesting_validation.py` re-tests
the effect three non-circular ways:

1. **Mean ± CI (no sorting).** Grand up-state spindle-z = **0.058, 95% CI [0.040, 0.077]** — excludes
   zero. The effect is real without any ordering, and visibly small.
2. **Split-half.** Rank patients by up-state spindle on **odd** SO events, display the mean envelope
   from their **even** events. The structure survives: **Spearman(odd, even) r = 0.56, p = 2.7×10⁻⁵**
   → per-patient nesting is a **stable, reproducible trait**, not double-dipping.
3. **Sort by an independent variable (SO trough depth).** Essentially no gradient:
   **Spearman(depth, spindle) r = 0.036** (significant only from n≈15k; negligible size). A deeper
   slow wave does **not** meaningfully predict more spindle.

**Conclusion.** The N3 SO→spindle nesting is **real and reproducible at the patient level**
(CI excludes 0; cross-validates at r=0.56), but **modest** (z≈0.06) and **not a function of
slow-wave size** — reliable in aggregate, mostly stochastic per event. This matches the modest
group MI_z (p=0.030) and the 42/49 heterogeneity, and is why the decisive evidence is the
independent 1 kHz replication (`slow_ripple_coupling`, p=6×10⁻⁵), not this atlas alone.

## Outputs
- `outputs/nesting_by_patient/nesting_validation.png` — the three non-circular views + stats.
- `outputs/nesting_by_patient/validation_stats.json` — split-half r/p, depth r/p, up-state CI.
- `outputs/nesting_by_patient/per_patient_nesting.png` — one small-multiple panel per patient
  (SO wave + spindle envelope), sorted by up-state nesting.
- `outputs/nesting_by_patient/so_trough_locked_spindle_heatmap.png` — SO-trough-locked spindle
  envelope; rows = patients (top) and rows = individual SO events (bottom), columns = time.
- `outputs/nesting_by_patient/per_patient_upstate.csv` — per-patient up-state spindle-z + event count.

## Cross-cohort replication (all three cohorts, stacked)

`nesting_validation_cohorts.py` runs the same non-circular tests on all three cohorts (spindle for
all; ripple for the two ≥1 kHz cohorts; the 204 Hz atlas cannot reach 80–120 Hz). Cohorts are shown
separately (replication), not pooled.

| Cohort | Band | up-state z [95% CI] | >0? | split-half r (p) |
|---|---|---|---|---|
| Atlas 204 Hz (n=49) | spindle | 0.058 [0.041, 0.076] | yes | 0.59 (8e-6) |
| Falach 1 kHz (n=15) | spindle | −0.005 [−0.035, 0.025] | no | 0.58 (0.024) |
| Falach 1 kHz | ripple | −0.003 [−0.013, 0.007] | no | 0.79 (5e-4) |
| Zurich 2 kHz (n=9) | spindle | 0.032 [0.016, 0.048] | yes | 0.67-0.73 (~0.03) |
| Zurich 2 kHz | ripple | −0.001 [−0.006, 0.002] | no | 0.67 (0.050) |

**Key finding.** The fixed post-trough up-state window does **not** uniformly replicate — atlas and
Zurich spindle are >0, but Falach spindle and both ripple bands sit at ~0. Split-half r is high
everywhere (subjects internally consistent), so this is not noise; it is a **polarity/montage
confound**: bipolar SO polarity is not standardized, so "post-trough" maps to different physiological
phases across cohorts.

**This does NOT contradict the coupling result.** M6's SO→spindle/ripple coupling uses **phase-agnostic
Tort MI** (p = 6×10⁻⁵ / 0.0015), which detects SO-phase modulation of amplitude *regardless of where
in the cycle it peaks* — so it survives across cohorts while the fixed-window temporal claim does not.

**Conclusion / what to claim.** Lead with the **phase-amplitude coupling** (robust, phase/polarity
independent). Do **not** headline the "rides the post-trough up-state" temporal narrative — it is
real in the atlas/Zurich spindle but fragile across cohorts. Proper fix: orient SO polarity per
channel (or align to each channel's own preferred SO phase) before the fixed-window average.

## Polarity-robust resolution (`nesting_phase_aligned`)

The fixed-window non-replication above was a **polarity confound**, and fixing it resolves the story.
`nesting_phase_aligned.py` aligns each channel to **its own preferred SO phase**, cross-validated:
half A (interleaved 5 s blocks) gives the preferred phase; half B is measured at that phase. θ* from
A, tested on B → non-circular; per-channel alignment → polarity-independent.

| Cohort | SO→spindle (cross-val modulation [CI]) | SO→ripple |
|---|---|---|
| Atlas 204 Hz (n=51) | 0.047 [0.026, 0.070] ✓>0 | — (204 Hz Nyquist) |
| Falach 1 kHz (n=15) | 0.047 [0.008, 0.097] ✓>0 | 0.011 [−0.008, 0.033] n.s. |
| Zurich 2 kHz (n=9) | 0.024 [0.014, 0.036] ✓>0 | 0.016 [0.010, 0.022] ✓>0 |

**Conclusion.** **SO→spindle coupling replicates across all three independent cohorts** (CI excludes
0 in each) — the earlier Falach "null" was polarity, not biology. **SO→ripple** replicates in Zurich
and is positive-but-underpowered in Falach (n=15) — real but weaker/less consistent than spindle.
This is the defensible, polarity-robust version of the coupling claim; the fixed-window temporal
narrative should still not be headlined.
