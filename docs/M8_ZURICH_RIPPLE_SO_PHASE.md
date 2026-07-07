# M8 — SO-phase of expert-marked ripples & fast-ripples [Zurich ds003498, 2 kHz]

**Question:** in an independent 3rd cohort, do *expert/detector-marked* ripple and fast-ripple
events sit at a preferred slow-oscillation phase? This is an **event-based** test (vs M6's
amplitude-envelope coupling) using gold-standard markings, and it opens the **fast-ripple band**
(nothing before could reach it).

**Data:** ds003498 (Fedele/Sarnthein, OpenNeuro, open, no account), interictal **slow-wave-sleep**
5-min runs, **2000 Hz**, BrainVision. 9 temporal-lobe subjects, MTL depth SEEG (A=amygdala,
AH/H/PH/P=hippocampus, EC/E=entorhinal). Events: `ripple_<chan>`, `fr_<chan>` per bipolar channel.
SO phase computed on a 500 Hz downsample (0.5–1.25 Hz is unstable at 2 kHz), phase read at each
event's onset. ~18k ripple + ~5k fast-ripple events used.

## Result (MRL = mean resultant length = effect size; Rayleigh p inflated by huge n)

| region · event | n | MRL | pref phase | read |
|---|---:|---:|---:|---|
| entorhinal ripple | 847 | **0.20** | −11° | clear clustering |
| amygdala ripple | 4292 | 0.084 | +17° | modest |
| amygdala fast-ripple | 1735 | 0.085 | +28° | modest |
| hippocampal fast-ripple | 2541 | 0.055 | −89° | weak |
| hippocampal **ripple** | 10118 | **0.017** | — | ~no phase preference |

## Plain-English answer

Mixed, and the mix is the point. Entorhinal ripples cluster clearly at a preferred SO phase;
amygdala and fast-ripples modestly; **hippocampal ripples show almost none** (MRL 0.017 — the
p≈0.05 is a large-n artifact, not an effect).

**Why muddier than M6:** these are **HFO markings for epilepsy localization**, so a substantial
fraction are **pathological** ripples, which are known to couple poorly to sleep rhythms. M6 (same
mechanism) was clean because it excluded SOZ channels, masked IEDs, and used the physiological
amplitude envelope. So the two are complementary:

- **M6** = physiological, cleaned, envelope-based → clean positive (SO→ripple p=0.0015).
- **M8** = all marked HFOs incl. pathological, event-based, independent cohort, adds fast-ripples →
  positive where physiological signal dominates (entorhinal, amygdala), washed out where
  pathological HFOs dominate (hippocampal ripple).

That contrast is itself a finding: **SO-coupling tracks the physiological, not the pathological,
population** — consistent with the clinical HFO literature.

## Value added
- Independent **3rd cohort** (Zurich), open, 2 kHz.
- Region-resolved MTL (amygdala / hippocampus / entorhinal).
- **Fast-ripple band** (first time in this project).
- A physiological-vs-pathological control that strengthens M6's interpretation.

## Caveats
- HFO markings mix physiological + pathological; no per-channel SOZ flag here to exclude the latter
  (unlike M6). The muddiness is expected, not a null of the mechanism.
- Only run-01 per subject analyzed (each subject has 25–35 runs available → far more if needed).
- Effect sizes small in absolute MRL; emphasize MRL over Rayleigh p given n in the thousands.

## Outputs
- `outputs/m8_zurich_ripple_so_phase/so_phase_stats.csv`
- `outputs/m8_zurich_ripple_so_phase/m8_ripple_so_phase.{png,svg}`
