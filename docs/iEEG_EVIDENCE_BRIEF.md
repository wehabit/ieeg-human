# iEEG Evidence Brief — Human NREM slow-wave coordination in the medial temporal lobe

**One line.** Using only open, public human intracranial EEG — **no new patient access or committee
approval required** — we show that the NREM slow oscillation organizes faster activity (spindles
and ripples) in human medial temporal lobe. This is the human counterpart of the Buzsáki mouse
pilot, and it is a **candidate neural substrate** for the sleep-dependent coordination that Kipnis /
Jiang-Xie link to CSF clearance. **Clearance itself is not measured here.**

**Why it matters for the plan.** The mouse pilot showed *target engagement* (50 Hz vibration
reorganized hippocampal firing, at theta). It could not show the *target state*, because it had no
natural NREM. These human results define that state and can be tested and extended entirely on open
data — with no new patient access or committee approval required.

---

## The evidence (core: two open cohorts; plus a third supporting cohort)

*Slow-power and Spindle-coupling are the same atlas cohort (the slow rhythm, then its coupling). Ripple-coupling is a second,
independent cohort. HFO-phase (below) is a third, supporting cohort.*

### Slow-power — the slow rhythm is present and state-graded
*Normative iEEG Sleep & Wake Atlas (Pattnaik & Litt 2024), HUP cohort, 204 Hz.*
In entorhinal + parahippocampal cortex, relative slow-band (0.5–4 Hz) power rises wake → N3 in
**20 of 21 patients (Wilcoxon p = 2.9×10⁻⁶)**. The same rise appears in hippocampal and
temporal-neocortical channels.
→ *The native slow rhythm exists in human MTL and tracks sleep depth.*

### Spindle-coupling — the slow oscillation organizes spindles (NREM-specific)
*Same atlas; surrogate-corrected phase-amplitude coupling (Tort MI_z).*
Slow-oscillation phase → spindle (11–16 Hz) amplitude coupling is stronger in N3 than wake in
mesiotemporal cortex (**N = 22, MI_z 0.044 → 0.262, p = 0.030, 17/22 patients**), and is
**near-zero or negative in REM** — the NREM-specific pattern predicted. Raw MI is wake-artifact
prone; the surrogate-corrected MI_z is the reported metric.
→ *The slow rhythm coordinates faster events, only in NREM.*

### Ripple-coupling — it also organizes ripples, and independently replicates the spindle result
*Falach / Geva-Sagiv / Eliashiv 2024 (OpenNeuro), 1 kHz, MTL SEEG; SOZ-excluded, IED-masked.*
On a **separate cohort** at higher sampling (ripples now accessible):
- **SO → ripple (80–120 Hz): p = 0.0015** (13/15 subjects).
- **SO → spindle: p = 6×10⁻⁵** (15/15 subjects) — an **independent replication of Spindle-coupling** on cleaner
  data, the decisive corroboration.
→ *The full NREM nesting — slow oscillation organizing both spindles and ripples — is present in
human MTL across independent cohorts.*

*(A third cohort, HFO-phase (Zurich, 2 kHz, expert-marked HFOs), adds fast-ripples and shows the coupling
tracks physiological rather than pathological events — supporting detail, not part of the core
spine.)*

---

## Mouse → Human bridge (figure: `mouse_human_bridge`)

| | Wake | **NREM · slow 0.5–4 Hz** | REM · theta 6–10 Hz |
|---|---|---|---|
| Kipnis Fig 2b (spike-triggered field) | incoherent | slow-wave organizes activity | theta wave |
| **Mouse (Buzsáki Dec4)** | — | *absent* (no natural NREM) | theta present under 50 Hz drive |
| **Human (this work)** | flat coupling | **slow→spindle (Slow-power/Spindle-coupling) + slow→ripple (Ripple-coupling)** | coupling drops |

The mouse hit the REM/theta column; the human data reach the **NREM/slow column** — the coordination
state Kipnis links to clearance. **Working hypothesis:** to engage this substrate, sleep stimulation
should *enhance native slow-wave coordination* (slow / Jiang-Xie frequencies), not 50 Hz.

---

## Honest limits (state these when presenting)

- **No clearance is measured.** We show a candidate *neural substrate*; the CSF/clearance link is
  Kipnis's, and remains a hypothesis in this work.
- **Epilepsy patients, interictal.** Normative filtering (SOZ/resected excluded, IEDs masked)
  mitigates but does not eliminate pathology.
- **Spindle-coupling's effect is modest and between-patient-variable** (p = 0.030); its strength comes from the
  independent Ripple-coupling replication, not from p alone.
- **Region differences (Region-profile) are descriptive** — pairwise MI_z contrasts are not significant.
- Atlas/Falach data are **clips**, not full nights; continuous full-night iEEG (iEEG.org) is set up
  for follow-on within-night analysis but is not part of this brief.

## Reproducibility
Each claim maps to one script and one figure (`analysis/m1…`, `m2…`, `m6…`, `mouse_human_bridge`); final
figures in `outputs/`. Frozen at git commit for this brief.
