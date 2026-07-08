# Staresina 2015-style nesting figure (reproduced on HUP165 full night)

Event-based reproduction of the canonical Staresina et al. 2015 (Nat Neurosci) panels on our
continuous full-night data (HUP165, 1024 Hz), for direct visual comparison to the established result.

| Panel | What | Result on our data |
|---|---|---|
| A | SO-trough-triggered spindle amplitude | clear bump — spindles rise on the SO up-state (153k SO troughs) |
| B | preferred SO phase of spindle peaks | clustered in **27/28 channels** (Rayleigh p<0.05) |
| C | spindle-trough-triggered ripple amplitude | bump at spindle trough — ripples rise there (748k troughs) |
| D | preferred spindle phase of ripple peaks | clustered in **16/28 channels (57%)** |

## Two honest points
1. **Reconciles the earlier "spindle→ripple null."** The continuous phase-amplitude method
   (`hierarchical_coupling`) found the spindle→ripple leg ~null; the **event-based** method here (the
   Staresina approach) finds it significant in **16/28 channels** — confirming my caveat that
   event-based detection is more sensitive for the fast leg. So the nesting *is* reproduced:
   SO→spindle strong (27/28), spindle→ripple present in the majority (16/28).
2. **Magnitude is weaker than the original.** Phase clustering is statistically significant (large n)
   but small in effect size (example channels R≈0.03–0.05 vs R≈0.2–0.4 in Staresina 2015). Likely
   causes: blind delta-power NREM staging (not scored PSG), all MTL bipolar channels pooled
   (amygdala/entorhinal, not hand-picked hippocampus), and lenient event thresholds. So: **same
   phenomena, correct directions, but modest magnitude** — consistent with everything else in this
   project (the coupling is real and reproducible but not large).

**Bottom line for a senior:** we reproduce the Staresina 2015 nesting qualitatively on open/continuous
data, at weaker magnitude. This is a *reproduction of established physiology*, not a new finding — the
value is the open, no-patient-access pipeline and the mouse↔human bridge, not the coupling itself.

## Outputs
- `outputs/staresina_style/staresina_style.png` — the 4-panel figure.
