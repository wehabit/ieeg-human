# Human iEEG Test of NREM Hierarchical Nesting (slow oscillation -> spindle -> ripple)

Phase-1 analysis plan. Goal: test whether the coordinated slow-wave field activity that the
Kipnis / Jiang-Xie work links to sleep-dependent CSF perfusion and clearance is present, and
**state-dependent**, in human hippocampus and entorhinal cortex — using open normative iEEG that
requires **no new patient access or committee approval**.

## Why this dataset, why now

The Buzsáki mouse pilot (Dec4, `hpaticmousebuzsakilab`) found that 50 Hz vibrotactile drive
produced **theta (6–10 Hz) spike-field organization** in dorsal hippocampus, but **not** **0.5-4 Hz slow-wave coordination**. See `DEC4_ALL_DHPC_KIPNIS_COORDINATION.md` and
`DEC4_THETA_COUPLING_STATES.md` in that repo. The documented reason: there was **no natural NREM
state** in the mouse pilot — only stimulation epochs. The frequency/state mismatch with the clearance literature was
the central caveat.

Human overnight iEEG closes exactly that gap: it contains **natural N2 / N3 / REM / wake** states,
so we can finally ask whether the slow-band coordination appears where and when the mechanism
predicts.

### Dataset
- **MNI Open iEEG Atlas** (Frauscher et al. 2018) — normative **wake** iEEG, ~1-min bipolar clips,
  100+ epilepsy patients, whole-brain incl. mesiotemporal. https://mni-open-ieegatlas.research.mcgill.ca
- **Multicenter iEEG Sleep Atlas** (LORIS) — normative **sleep** iEEG segmented by stage:
  N2 (1468), N3 (1468), REM (1012); wake (1772). https://ieegatlas.loris.ca
- **Normative iEEG Sleep & Wake Atlases** (Pattnaik & Litt 2024) — processed clips + z-scores,
  CC-BY-NC-SA-4.0, DOI [10.26275/xhte-d11l](https://doi.org/10.26275/xhte-d11l), 4.36 GB.

### Regions of interest
Hippocampus, entorhinal cortex, parahippocampal, amygdala (mesiotemporal), plus a neocortical
comparison set. Mirrors the mouse dHPC vs LEC contrast.

## The nesting prediction, stated as testable hypotheses

- **H1 (state-dependence of coordination).** Slow-wave (0.5–4 Hz) power and coordination in
  hippocampus/EC increase across wake → N2 → N3, and fall in REM. This is the human version of "is
  there slow-band coordination, and is it in the state the mechanism needs."
- **H2 (nesting / phase-amplitude coupling).** The slow-oscillation phase organizes faster events
  (spindles 11–16 Hz; ripples/HFO 80–200 Hz). Coupling strength is maximal in N3/N2, minimal in
  wake/REM. This is the LFP analogue of the mouse **spike-field PLV** — instead of spikes-to-theta,
  it is fast-event-amplitude-to-slow-phase.
- **H3 (phase preference).** Nested events cluster at a preferred slow-oscillation phase (human MTL
  ripples on the up-state). Mirrors the mouse preferred-theta-phase histograms.
- **H4 (regional gradient).** Coordination differs hippocampus vs entorhinal vs neocortex. Mirrors
  the mouse dHPC-vs-LEC 50 Hz gradient figure.

## Analysis modules (each mirrors a mouse-pipeline step)

| Module | Human (LFP, state atlas) | Mouse analogue |
|---|---|---|
| **Slow-power Spectral state map** | Per-channel band power (slow 0.5–4, delta, theta 6–10, spindle 11–16, gamma, HFO 80–200) per state; test wake→N2→N3 slow-band rise | Field theta amplitude by epoch |
| **Spindle-coupling Slow-phase → fast-amp coupling** | Tort modulation index: 0.5–4 Hz phase vs spindle & ripple amplitude, per state | Theta spike-field PLV / PPC |
| **Phase-preference Preferred slow phase** | Phase histogram of nested spindle/ripple events on the SO cycle | Preferred theta-phase histograms |
| **Region-profile Regional gradient** | Hippocampus vs entorhinal vs neocortex coordination | dHPC vs LEC 50 Hz gradient |
| **Bridge Bridge figure** | Overlay mouse (driven theta, 6–10 Hz) vs human (natural NREM slow, 0.5–4 Hz); frame the stim-frequency hypothesis | — |

**Bridge is the payoff for the grant / one-pager.** It states the design hypothesis explicitly: to
engage this substrate, vibrotactile stimulation during NREM should aim to **enhance the native
slow-wave coordination** (slow / Jiang-Xie testing frequencies), not 50 Hz. The mouse 50 Hz result
was target *engagement*; the human atlas defines the *target state*.

## Method notes (carry over the mouse pipeline's rigor)

- Bipolar / region-mean referencing; comb-notch line noise (50 or 60 Hz + harmonics) before banding.
- PAC via Tort MI with surrogate (phase-shuffled) null; report only MI above surrogate CI.
- Bootstrap CIs on all state contrasts; equal-segment resampling so segment-count differences across
  states don't bias the metric (same discipline as the mouse PLV equal-spike resampling).
- Report per-channel and per-patient, not just grand mean.

## Honest caveats (same discipline as the DEC4 docs)

- **Normative epilepsy patients, interictal.** Clips are curated as non-epileptic, but these are not
  healthy brains.
- **LFP only, no spikes.** Coordination here is field-level PAC, not spike-field locking. It is the
  correct level for the field-wave claim but is not identical to the mouse unit analysis.
- **Segmented clips, not continuous overnight.** Within-clip nesting is testable; long-range
  SO→spindle→ripple sequences and full-night dynamics are not (that needs iEEG.org / DABI raw data —
  a Phase-1b extension).
- **No clearance readout.** Coordination is a candidate neural substrate Kipnis links to clearance, not
  a clearance measurement. State this every time.

## Deliverables
- `analysis/` one script per module (Slow-power–Bridge), CLI'd like the mouse scripts.
- `outputs/` figures (png+svg) + csv per module.
- `docs/` a short markdown writeup per module in the DEC4 house style (plain-English answer, method,
  result table, caveats).

## Stack
Python + MNE-Python, numpy/scipy, `tensorpac` (or custom Tort MI), pandas, matplotlib. See
`env/requirements.txt`.

## Sequence
1. **M0 data fetch + loader** — download atlas, parse region/state metadata into a tidy index.
2. **Slow-power spectral state map** — the anchor result; confirms H1 before anything fancy.
3. **Spindle-coupling PAC** — the core coordination metric (H2).
4. **Phase-preference–Region-profile** — phase preference and regional gradient.
5. **Bridge bridge figure** — mouse↔human, feeds the one-pager.
