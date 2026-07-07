# Slow-oscillation → ripple (and spindle) coupling  [1 kHz cohort]

**Question:** the mouse study showed CA1 sharp-wave **ripples**; the 204 Hz atlas (Slow-power–Bridge) could
not touch ripples (Nyquist 102 Hz). Do human NREM slow oscillations organize **ripples** (80–120 Hz)
in the medial temporal lobe — and does the spindle coupling from Spindle-coupling replicate on cleaner data?

**Data:** Falach, Geva-Sagiv, Eliashiv et al. 2024 (Figshare, CC-BY-NC), overnight iEEG sleep
segments, **1000 Hz**, SEEG in MTL (amygdala, hippocampus, entorhinal, parahippocampal), with
expert IED annotations and a SOZ flag. NREM-dominant, so this is a **within-NREM** coordination test
(state contrast was Slow-power/Spindle-coupling's job), extended to the ripple band.

## Plain-English answer

Yes on both, and the spindle result strongly replicates.

- **SO → ripple (80–120 Hz): MI_z = 1.14, positive in 13/15 subjects, p = 0.0015.** Human MTL
  ripples are phase-organized by the NREM slow oscillation — the ripple analogue of the mouse
  CA1 figures, and the band the atlas physically could not reach.
- **SO → spindle (11–16 Hz): MI_z = 2.15, 15/15 subjects, p = 6.1e-5.** This **independently
  confirms Spindle-coupling**. The coupling that was real-but-noisy at 204 Hz (p = 0.030) is rock-solid on this
  cleaner, IED-removed, higher-rate dataset.
- Both couplings are present across **all four MTL regions** (hippocampal, entorhinal,
  parahippocampal, amygdala); entorhinal — thin and noisy in Spindle-coupling — is robust here.

## Rigor (why these are physiological, not epileptic, ripples)

- **Bipolar re-referencing** within each electrode shaft (kills reference/volume-conducted artifact).
- **SOZ channels dropped** (`soz_region` flag) — no seizure-onset tissue.
- **IED masking:** every sample within ±0.5 s of an annotated interictal discharge is removed
  before analysis, so pathological fast-ripples around spikes don't drive the coupling.
- **Surrogate-corrected MI_z** (120 circular shifts), same discipline as Spindle-coupling.

## Coverage

129 bipolar MTL channels across 15 subjects: hippocampal 52, amygdala 32, entorhinal 26,
parahippocampal 19. (10 of 25 subjects contributed no channels — all-SOZ shafts or no MTL coverage
after filtering.)

## Where it fits

Ripple-coupling completes the human coordination picture the Gyuri pointer set implied:

| Mouse figure | Human module | Band |
|---|---|---|
| theta STA / PLV | Spindle-coupling + comodulogram | spindle (11–16 Hz) |
| **CA1 ripples** | **Ripple-coupling** | **ripple (80–120 Hz)** |
| dHPC vs LEC | Region-profile | region gradient |

The slow oscillation organizes **both** spindles and ripples in human NREM MTL — the full nesting
that the Kipnis mechanism links to clearance. The 50 Hz mouse drive engaged none of this native
structure; the design implication (Bridge) stands and is now stronger.

## Caveats

- Short NREM-dominant segments (~1–4 min), not full nights → no within-subject state contrast here
  (that is Slow-power/Spindle-coupling). Full-night dynamics (SO→spindle→ripple temporal sequences) still need continuous
  recordings.
- Epilepsy patients; mitigations above reduce but do not eliminate pathological contribution.
- Ripple band fixed at 80–120 Hz; a full ripple-detection (event-based) pipeline is a next step.

## Outputs
- `outputs/slow_ripple_coupling/channel_coupling.csv`
- `outputs/slow_ripple_coupling/ripple_coupling_by_region.{png,svg}`
- `outputs/slow_ripple_coupling/ripple_comodulogram.{png,svg}`
