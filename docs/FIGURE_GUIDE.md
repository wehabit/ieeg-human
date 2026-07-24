# Figure Guide — what each analysis is, in plain English

> **QUARANTINED / NOT AUDITED.** This guide presents withdrawn legacy figures and numerical claims
> as completed findings. It predates the LC-proxy audit and must not be cited or used to select
> current outputs. See `ISSUE_REGISTER_2026-07.md`; the corrected RESPect rerun has no estimable
> 3A/3B endpoint, and the corrected HUP rerun has no estimable 3A/3B or pooled 3D endpoint.

You are testing **one question**: *does the human brain show the Kipnis NREM slow-wave
coordination that the mouse pilot couldn't produce?* Every figure is a piece of that.

The mouse study (Gyuri pointer set) used **single neurons + a vibration stimulus**. The human
atlas has **neither spikes nor a stimulus** — it has **LFP across natural sleep states**. So each
mouse figure either (a) has a human LFP analogue we made, or (b) needs data we don't have yet.

## The Rosetta table: mouse figure type → human module

| Mouse figure (Gyuri) | What it showed | Human analogue | Module | Status |
|---|---|---|---|---|
| raster + PSTH (1–3) | single neurons change firing | — needs spikes | — | not possible here (needs microwires) |
| frequency-specificity (4–5) | response concentrates at 50 Hz | response concentrates in **N3** (state, not freq) | **Slow-power / Spindle-coupling** | ✅ done |
| trial-avg spectrogram (6) | 50 Hz band grows | power spectrum grows in slow band by state | **Slow-power** (`psd_by_state`) | ✅ done |
| 50 Hz-is-pickup controls (7–9) | LFP signal is artifact | wake coupling is artifact → use surrogate **MI_z** | **Spindle-coupling** | ✅ done (built in) |
| spikes-are-real (10–11) | spike QC / ON-OFF summary | — needs spikes | — | not possible here |
| **ripple examples / localization** | **CA1 sharp-wave ripples** | **SO → ripple coupling** | **Ripple-coupling (Phase-1b)** | ✅ done (1 kHz dataset, p=0.0015) |
| unit-by-shank / cell type | where units sit, pyr vs int | — needs spikes | — | not possible here |
| **region processing** (dHPC vs LEC) | regions do opposite things | hippo vs entorhinal vs neocortex | **Region-profile** | ✅ done |
| **theta STA / PLV-by-state** (12–13) | spikes lock to theta wave, more under drive | spindle amp locks to slow-wave phase, more in N3 | **Spindle-coupling + comodulogram** | ✅ done |

## What each human figure actually shows (one line each)

**Slow-power — `slow_power_by_state` / `psd_by_state` / `slow_power_paired_slope`**
The slow wave (0.5–4 Hz) gets stronger going wake→N2→N3. *This is the rhythm existing.*
Rock-solid: 20/21 patients, p=2.9e-6.

**Spindle-coupling — `spindle_coupling_by_state` / `comodulogram` / `tort_phase_amp` / `spindle_coupling_paired_slope`**
That slow wave **organizes** spindles (fast bursts cluster at one phase of the slow wave), and
only in NREM. *This is the coordination — the actual Kipnis mechanism.* Real but noisier:
17/22 patients, p=0.030 (independently replicated in Ripple-coupling, p=6e-5).
- `comodulogram` is the clearest: dark (no coupling) in wake, bright hotspot (coupling) in N3.
- `tort_phase_amp` is the same thing as a bar chart: flat in wake, peaked in N3.

**Phase-preference — the polar plot inside `spindle_coupling_by_state`**
Spindles prefer a specific phase of the slow wave (the up-state), not random. *Coordination is
phase-locked, not coincidence.*

**Region-profile — `mtl_region_profile` (DESCRIPTIVE only)**
Coupling looks strongest in hippocampus, weakest in entorhinal — but **pairwise MI_z differences
are not significant** (all p≥0.13). Treat as a descriptive picture, not a claim.

**Bridge — `mouse_human_bridge`**
The one slide: mouse hit the REM/**theta** column; human hits the NREM/**slow** column — the one
Kipnis *links* to clearance (clearance not measured here). *This is the argument for the grant.*

**result_visuals — `example_traces`**
Raw signal: flat/fast in wake, big slow waves (with a spindle riding one) in N3. *The phenomenon
before any math.*

**Nesting-timecourse — `so_triggered_nesting` / `so_triggered_tf` (temporal / time-domain)**
The same nesting shown in TIME instead of frequency. Align on the slow-oscillation trough and
average: the spindle envelope peaks on the SO up-state (panel B, 10,292 events), and the
time-frequency map shows the slow wave building into the trough with faster bands modulated
around it. *This is the Kipnis Fig 2b / mouse spike-triggered-average view, at the field level.*

**nesting_by_patient — `per_patient_nesting` / `so_trough_locked_spindle_heatmap`**
Sanity check on the grand average: per-patient and per-event N3 SO-triggered spindle nesting.
Positive in 42/49 patients (86%) but heterogeneous (median 0.057) — real and majority-present,
not uniform. See [NESTING_BY_PATIENT.md](NESTING_BY_PATIENT.md).

**HFO-phase — `ripple_so_phase` (Zurich cohort, 2 kHz, event-based)**
Independent 3rd cohort using expert-marked ripples/fast-ripples: what SO phase does each event land
on? Entorhinal/amygdala cluster; hippocampal *ripples* don't — because these are epilepsy HFO
markings and pathological ones couple poorly. *Confirms the mechanism where physiological signal
dominates, and shows SO-coupling tracks physiological (not pathological) events.* Adds the
fast-ripple band. See [HFO_SLOW_PHASE.md](HFO_SLOW_PHASE.md).

## What is still worth generating

- **Ripple-coupling — SO→ripple coupling (Phase-1b).** The one clearly-relevant mouse figure type (ripples) we
  couldn't do at 204 Hz. The Falach/Geva-Sagiv 2024 iEEG sleep dataset is 2 kHz (ripples visible),
  overnight, hippocampus/entorhinal/amygdala, open on Figshare. Adds the ripple half of the human
  coordination story and a true continuous-night view.
- Everything else on the "not possible here" rows needs **single-unit microwire** data
  (Rutishauser/DANDI) — a later, separate Phase-1b if you want the literal spike-level analogue.
