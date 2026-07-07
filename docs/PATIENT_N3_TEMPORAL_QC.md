# Patient N3 Temporal QC

This check asks whether the slow-fast nesting is visible patient by patient,
instead of only in a grand average.

## What Was Generated

Script: `analysis/patient_n3_temporal_qc.py`

Output folder: `outputs/patient_n3_temporal_qc/`

- `patients/*_n3_so_triggered_nesting.png`: one figure per patient.
- `all_patients_so_triggered_contact_sheet.png`: all patient averages in one grid.
- `all_patients_kipnis_style_so_heatmap.png`: Kipnis-style patient heatmap.
- `patient_n3_temporal_summary.csv`: clips, channel coverage, duration, and event counts.

## What "Full N3" Means Here

The Pattnaik/Litt atlas is not a continuous full-night dataset. It provides N3
clips, usually up to 10 clips per patient. So "full N3" means all available N3
clips for that patient, across all normative mesiotemporal channels.

## What Each Patient Figure Shows

1. **Raw N3 example**: gray raw LFP, blue slow oscillation, orange spindle envelope.
2. **SO-triggered average**: all detected N3 slow-oscillation troughs aligned at
   time 0; blue is the average slow wave, orange is spindle-envelope timing.
3. **Kipnis-style event heatmap**: each row is one detected slow-oscillation
   cycle. The dark trough at time 0 and red/blue rhythm around it show that the
   slow wave is not only a grand-average artifact.

## Current Coverage

51 of 53 N3/MTL atlas patients produced usable SO-triggered figures.

Two patients (`sub-RID0193`, `sub-RID0294`) were skipped because their N3 EDF
channel names did not overlap the filtered MTL metadata channel names. Treat
these as channel-coverage/QC exclusions, not biological negatives.

## Interpretation

These figures show SO→spindle temporal nesting in N3. They do not show ripples,
because the sleep atlas uses mixed-rate clips that are not appropriate for a
clean ripple analysis across all patients. Ripple nesting is handled in the
separate 1 kHz cohort (`analysis/slow_ripple_coupling.py`).
