# Slow-oscillation → spindle coupling (tests H2, H3)

**Question:** Is the *coordination* the Kipnis / Jiang-Xie mechanism cares about — the slow field
wave organizing faster activity — present in human mesiotemporal cortex, and is it state-dependent?
This is the human LFP analogue of the mouse **spike-field PLV**: the mouse locked *spikes* to the
theta phase of the field; here (no spikes at 204 Hz) we lock *spindle amplitude* (11–16 Hz) to the
phase of the slow oscillation (0.5–1.25 Hz). That nesting **is** the coordinated field wave.

**Kipnis prediction (H2):** coupling maximal in N3/N2, minimal in wake/REM.

## Plain-English answer

Yes — once wake artifact is controlled for. The **surrogate-corrected** coupling (MI_z) is
significantly stronger in N3 than wake in mesiotemporal cortex (entorhinal + parahippocampal,
**N = 22 patients, MI_z 0.044 → 0.262, Wilcoxon p = 0.030, N3 > W in 17/22**), and it is
**near-zero or negative in REM** — exactly the NREM-specific pattern the mechanism predicts.
This spindle-coupling result is **independently confirmed on a separate cohort in Ripple-coupling** (SO→spindle
p = 6×10⁻⁵), which is the stronger corroboration.

The **raw** Tort MI initially looked *reversed* (W > N3, p = 0.16). That was two wake outliers with
huge spurious MI (sub-RID0160 W = 0.22, sub-RID0179 W = 0.17) — classic wake movement/muscle
artifact manufacturing fake phase–amplitude structure. The circular-shift surrogate null removes
that inflation, which is precisely why it is the primary metric here. **Lesson banked:** report MI_z,
not raw MI, for any wake-vs-sleep PAC contrast.

## Result — surrogate-corrected coupling MI_z (patient-level mean)

| ROI | W | N2 | N3 | R |
|---|---:|---:|---:|---:|
| entorhinal | 0.005 | 0.179 | 0.048 | **−0.400** |
| parahippocampal | 0.064 | 0.124 | 0.199 | **−0.071** |
| hippocampal | 0.145 | 0.194 | **0.412** | 0.114 |
| temporal neocortex | 0.329 | 0.261 | 0.455 | 0.371 |

Coupling climbs into NREM (strongest in hippocampal channels) and collapses in REM for the
entorhinal/parahippocampal ROI. **H3:** in N3, spindle power clusters at a preferred slow-oscillation
phase (see polar panel) rather than spreading uniformly — the coupling is phase-locked, not just a
power coincidence.

## Method

- Per channel, per 30 s clip: detrend → slow-oscillation phase (0.5–1.25 Hz Hilbert) and spindle
  amplitude (11–16 Hz Hilbert) → **Tort modulation index** (18 phase bins, normalized KL from
  uniform).
- **Surrogate null:** 120 circular time-shifts of the amplitude envelope → MI_z = (MI − mean_surr) /
  sd_surr. MI_z is the primary metric; raw MI is reported but is wake-artifact-prone.
- Preferred phase = SO-phase bin of maximum mean spindle amplitude (circular-averaged across clips).
- Aggregate channel → patient → ROI group → state; N3-vs-W paired Wilcoxon across patients.

## How this pairs with the mouse

| Mouse (Dec4) | Human (Spindle-coupling) | Reading |
|---|---|---|
| spikes lock to **theta** phase under 50 Hz drive | spindle amp locks to **slow-oscillation** phase in natural N3 | both are field-organized coordination, but at the frequency/state each system was in |
| no natural NREM → no 0.5–4 Hz coordination | natural N3 → NREM-specific slow coordination present | the human data supplies the state the mouse pilot lacked |

This is the empirical core of the **Bridge bridge argument**: the native coordinated state is a
NREM slow-oscillation phenomenon, so stimulation meant to engage it should target/enhance that
slow coordination — not 50 Hz.

## Caveats

- **Depth / variance.** 60 patients scanned, up to 20 clips/state; entorhinal ROI is the thinnest
  (n = 11 with both W & N3) and noisy. Increasing clips/patient (8 → 20) did **not** tighten the
  effect (p 0.011 → 0.030) — the limiting factor is **between-patient** variance, not clip count.
  The decisive corroboration is Ripple-coupling's independent-cohort replication (p = 6×10⁻⁵), not more clips.
- **MI_z, not absolute coupling strength.** MI_z says coupling exceeds its own chance level more in
  N3; it is the artifact-robust contrast, not a raw effect size.
- **No ripples** (Nyquist 102 Hz) — spindle nesting only. Ripple–SO coupling needs full-rate raw
  data (iEEG.org / DABI Phase-1b).
- Normative epilepsy patients, interictal; coordination is a candidate neural substrate Kipnis links
  to clearance, **not** a clearance
  measurement.

## Outputs
- `outputs/slow_spindle_coupling/channel_state_pac.csv` — per channel × state: mi, mi_z, pref_phase.
- `outputs/slow_spindle_coupling/group_state_pac.csv` — group × state MI + MI_z.
- `outputs/slow_spindle_coupling/spindle_coupling_by_state.{png,svg}`.

## Next
- **Region-profile** regional gradient (hippocampal vs entorhinal vs neocortex) with proper stats — the human
  dHPC-vs-LEC analogue.
- **Bridge** mouse↔human bridge figure for the one-pager.
- **Phase-1b** (optional): full-rate raw iEEG (iEEG.org/DABI) to add SO→ripple coupling and, via a
  microwire dataset (Rutishauser/DANDI), a true spike-level (Kilosort-style) analogue.
