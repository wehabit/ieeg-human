# Patient N3 Temporal QC

This check asks whether the slow-fast nesting is visible patient by patient,
instead of only in a grand average.

## What Was Generated

Script: `analysis/patient_n3_temporal_qc.py`

Output folder: `outputs/patient_n3_temporal_qc/`

- `patients/*_n3_so_triggered_nesting.png`: one figure per patient.
- `long_traces/*_n3_long_trace_180s.png`: up to 3 minutes of raw N3 temporal strip per patient.
- `all_patients_so_triggered_contact_sheet.png`: all patient SO-trough-locked averages in one grid.
- `all_patients_so_trough_locked_heatmap.png`: all-patient SO-trough-locked heatmap.
- `patient_n3_temporal_summary.csv`: clips, channel coverage, duration, and event counts.

## Exact Analysis

This QC figure is **not a new p-value test**. It is an event-triggered temporal
visualization that asks: when we align N3 mesiotemporal LFP to detected slow
oscillation troughs, does the repeating slow wave and spindle-envelope timing
remain visible patient by patient?

Steps:

1. Use all available N3 atlas clips for each patient and all normative
   mesiotemporal channels.
2. Detrend each trace and notch 60 Hz line noise.
3. Filter the slow oscillation at **0.5-1.25 Hz**.
4. Detect clear negative SO troughs, at least 0.6 s apart.
5. Extract +/-1.5 s windows around every trough.
6. Plot the SO-filtered LFP and the **11-16 Hz spindle amplitude envelope**
   aligned to time 0.

The repetitive frequency in the heatmap is therefore the **SO band itself:
0.5-1.25 Hz**, or roughly **0.8-2.0 seconds per cycle**. The formal coupling
statistic remains the M2 surrogate-corrected phase-amplitude coupling analysis:
Tort MI_z with circular-shift surrogates, then a paired Wilcoxon signed-rank
test for N3 vs wake across patients (**p = 0.030; N3 > wake in 17/22
patients**).

## What "Full N3" Means Here

The Pattnaik/Litt atlas is not a continuous full-night dataset. It provides N3
clips, usually up to 10 clips per patient. So "full N3" means all available N3
clips for that patient, across all normative mesiotemporal channels.

The long temporal strips use **one representative MTL channel per patient** and
show up to **180 seconds** of N3. These are concatenated 30 s atlas clips, with
clip breaks marked by dotted vertical lines; they should be read as long visual
QC, not as uninterrupted physiology.

## What Each Patient Figure Shows

1. **Raw N3 example**: gray raw LFP, blue slow oscillation, orange spindle envelope.
2. **SO-trough-locked average**: all detected N3 slow-oscillation troughs aligned at
   time 0; blue is the average slow wave, orange is spindle-envelope timing.
3. **SO-trough-locked event heatmap**: each row is one detected slow-oscillation
   cycle. The dark trough at time 0 and red/blue rhythm around it show that the
   slow wave is not only a grand-average artifact.

## What The Long Trace Shows

`long_traces/*_n3_long_trace_180s.png` is the answer to "what does this look
like over minutes?" Each row is 30 s. Gray is the raw LFP, blue is the
0.5-1.25 Hz slow oscillation, and orange is the scaled 11-16 Hz spindle
envelope.

What to infer:

- You can see whether the patient has sustained N3 slow oscillations over
  minutes, not just a hand-picked 8 s example.
- You can see whether spindle power comes in bursts riding on that slow rhythm.
- You can spot noisy clips, discontinuities, or patients with weaker slow waves.

What not to infer:

- The strips are not one continuous 3-minute recording unless no dotted clip
  break appears.
- They are not the formal PAC statistic; they are the time-domain intuition.

## What The Contact Sheet Is For

`all_patients_so_triggered_contact_sheet.png` is a **QC/skepticism plot**. It
does not prove the effect by itself. Its benefit is that it lets you see whether
the group result is broadly present across patients or whether it is being
carried by one or two spectacular examples.

What to infer:

- If many panels show a blue trough at time 0 with flanking slow-wave structure,
  the SO detector is finding a real repeating N3 rhythm across patients.
- If the orange spindle envelope rises at a consistent phase relative to that
  trough/peak, that is the temporal intuition behind SO->spindle nesting.
- If some panels are flat, noisy, or phase-shifted, that is heterogeneity; it is
  useful honesty, not a failure of the package.

What not to infer:

- It is not a clearance measurement.
- It is not a ripple result.
- It is not the main statistical claim; the p-value lives in the M2 PAC analysis.

## Current Coverage

51 of 53 N3/MTL atlas patients produced usable SO-trough-locked figures.

Two patients (`sub-RID0193`, `sub-RID0294`) were skipped because their N3 EDF
channel names did not overlap the filtered MTL metadata channel names. Treat
these as channel-coverage/QC exclusions, not biological negatives.

## Interpretation

These figures show SO→spindle temporal nesting in N3. They do not show ripples,
because the sleep atlas uses mixed-rate clips that are not appropriate for a
clean ripple analysis across all patients. Ripple nesting is handled in the
separate 1 kHz cohort (`analysis/slow_ripple_coupling.py`).
