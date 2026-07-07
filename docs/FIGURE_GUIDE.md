# Figure Guide — what each analysis is, in plain English

You are testing **one question**: *does the human brain show the Kipnis NREM slow-wave
coordination that the mouse pilot couldn't produce?* Every figure is a piece of that.

The mouse study (Gyuri pointer set) used **single neurons + a vibration stimulus**. The human
atlas has **neither spikes nor a stimulus** — it has **LFP across natural sleep states**. So each
mouse figure either (a) has a human LFP analogue we made, or (b) needs data we don't have yet.

## The Rosetta table: mouse figure type → human module

| Mouse figure (Gyuri) | What it showed | Human analogue | Module | Status |
|---|---|---|---|---|
| raster + PSTH (1–3) | single neurons change firing | — needs spikes | — | not possible here (needs microwires) |
| frequency-specificity (4–5) | response concentrates at 50 Hz | response concentrates in **N3** (state, not freq) | **M1 / M2** | ✅ done |
| trial-avg spectrogram (6) | 50 Hz band grows | power spectrum grows in slow band by state | **M1** (`psd_by_state`) | ✅ done |
| 50 Hz-is-pickup controls (7–9) | LFP signal is artifact | wake coupling is artifact → use surrogate **MI_z** | **M2** | ✅ done (built in) |
| spikes-are-real (10–11) | spike QC / ON-OFF summary | — needs spikes | — | not possible here |
| **ripple examples / localization** | **CA1 sharp-wave ripples** | **SO → ripple coupling** | **M6 (Phase-1b)** | ✅ done (1 kHz dataset, p=0.0015) |
| unit-by-shank / cell type | where units sit, pyr vs int | — needs spikes | — | not possible here |
| **region processing** (dHPC vs LEC) | regions do opposite things | hippo vs entorhinal vs neocortex | **M4** | ✅ done |
| **theta STA / PLV-by-state** (12–13) | spikes lock to theta wave, more under drive | spindle amp locks to slow-wave phase, more in N3 | **M2 + comodulogram** | ✅ done |

## What each human figure actually shows (one line each)

**M1 — `m1_slow_by_state` / `psd_by_state` / `m1_paired_slope`**
The slow wave (0.5–4 Hz) gets stronger going wake→N2→N3. *This is the rhythm existing.*
Rock-solid: 20/21 patients, p=2.9e-6.

**M2 — `m2_pac_by_state` / `comodulogram` / `tort_phase_amp` / `m2_paired_slope`**
That slow wave **organizes** spindles (fast bursts cluster at one phase of the slow wave), and
only in NREM. *This is the coordination — the actual Kipnis mechanism.* Real but noisier:
17/22 patients, p=0.030 (independently replicated in M6, p=6e-5).
- `comodulogram` is the clearest: dark (no coupling) in wake, bright hotspot (coupling) in N3.
- `tort_phase_amp` is the same thing as a bar chart: flat in wake, peaked in N3.

**M3 — the polar plot inside `m2_pac_by_state`**
Spindles prefer a specific phase of the slow wave (the up-state), not random. *Coordination is
phase-locked, not coincidence.*

**M4 — `m4_region_gradient` (DESCRIPTIVE only)**
Coupling looks strongest in hippocampus, weakest in entorhinal — but **pairwise MI_z differences
are not significant** (all p≥0.13). Treat as a descriptive picture, not a claim.

**M5 — `m5_bridge`**
The one slide: mouse hit the REM/**theta** column; human hits the NREM/**slow** column — the one
Kipnis *links* to clearance (clearance not measured here). *This is the argument for the grant.*

**result_visuals — `example_traces`**
Raw signal: flat/fast in wake, big slow waves (with a spindle riding one) in N3. *The phenomenon
before any math.*

**M7 — `m7_so_triggered_nesting` / `m7_so_triggered_tf` (temporal / time-domain)**
The same nesting shown in TIME instead of frequency. Align on the slow-oscillation trough and
average: the spindle envelope peaks on the SO up-state (panel B, 10,292 events), and the
time-frequency map shows the slow wave building into the trough with faster bands modulated
around it. *This is the Kipnis Fig 2b / mouse spike-triggered-average view, at the field level.*

**M8 — `m8_ripple_so_phase` (Zurich cohort, 2 kHz, event-based)**
Independent 3rd cohort using expert-marked ripples/fast-ripples: what SO phase does each event land
on? Entorhinal/amygdala cluster; hippocampal *ripples* don't — because these are epilepsy HFO
markings and pathological ones couple poorly. *Confirms the mechanism where physiological signal
dominates, and shows SO-coupling tracks physiological (not pathological) events.* Adds the
fast-ripple band. See [M8_ZURICH_RIPPLE_SO_PHASE.md](M8_ZURICH_RIPPLE_SO_PHASE.md).

## What is still worth generating

- **M6 — SO→ripple coupling (Phase-1b).** The one clearly-relevant mouse figure type (ripples) we
  couldn't do at 204 Hz. The Falach/Geva-Sagiv 2024 iEEG sleep dataset is 2 kHz (ripples visible),
  overnight, hippocampus/entorhinal/amygdala, open on Figshare. Adds the ripple half of the human
  coordination story and a true continuous-night view.
- Everything else on the "not possible here" rows needs **single-unit microwire** data
  (Rutishauser/DANDI) — a later, separate Phase-1b if you want the literal spike-level analogue.
