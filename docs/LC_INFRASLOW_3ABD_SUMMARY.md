# Testing the locus-coeruleus infraslow hypothesis in human iEEG — summary

**Bottom line: the ~0.02 Hz locus-coeruleus fingerprint is not present**, and there is no N2/N3
difference in any test. Two corrections were made after the fact and are documented rather than
folded in: **3A** looked positive and was **retracted** after a frequency-specificity control, while
**3D** was reported null and is now **positive** — that null was an artifact of not following the
source papers' method. 3B's cited reference turned out to be a review, not a methods paper.

---

## 1. The hypothesis

The locus coeruleus (LC) is proposed to impose a shared **infraslow (~0.02 Hz, ~50 s)** rhythm on
sleep-spindle activity and on heart rate during NREM, coupling the two
(Lecci et al. 2017 *Sci Adv*; Osorio-Forero et al. 2021 *Curr Biol*). Because the LC degenerates
early in neurodegenerative disease, weakened coupling was hypothesised as a candidate biomarker.

None of these signals record the LC directly — they are **proxies**, and that inference rests
entirely on the papers above.

## 2. Dataset

**iEEG.org HUP phase-II** — continuous multi-day intracranial recordings from the Penn epilepsy
monitoring unit.

Discovery step: of 52 candidate datasets, **39 were reachable and every one carries an EKG/ECG
channel** alongside the depth electrodes. This is what made the study possible at all — public iEEG
almost never includes cardiac signals. **25 subjects** have both depth electrodes and EKG.

Analysis used **lateral neocortical contacts** (the highest-numbered contact on each shaft — in SEEG,
contact 1 is deepest/mesial and the highest contact is most lateral). This matters: an initial run on
mesial-temporal contacts was negative, but spindles are thalamo*cortical* and Lecci's effect is
parietal-maximal, so the mesial result did not test the hypothesis. *Caveat: "highest contact =
neocortex" is geometric inference — no electrode localisation files were available.*

Per subject: ~7 h streamed from the highest-delta night, 3 signals derived (1 Hz spindle power,
1 Hz heart rate, slow-oscillation event times), plus 30 s epoch staging features.

### Subject flow

| Stage | n |
|---|---|
| Subjects with depth + EKG (cohort) | **25** |
| Excluded — no usable cortical channels | **2** (HUP116: 4-contact shafts only; HUP141: 2 channels) |
| **Analysed** | **23** |
| Passed staging-quality filter (for N2/N3 contrasts) | **17** |
| Excluded — degenerate stage split | 4 (HUP138 377/2, HUP139 612/4, HUP191 429/12, HUP199 459/45) |
| Excluded — <50 epochs in one stage | 2 (HUP171 33/170, HUP177 70/49) |
| Stage-resolved **3A** possible | **4** only |

## 3. The three tests

| Test | What it asks | Measure | Source |
|---|---|---|---|
| **3A** *(primary)* | Do spindle power and heart rate rise and fall **together every ~50 s**? | Magnitude-squared coherence between the 1 Hz spindle-power series and instantaneous heart rate, read at the pre-specified **0.02 Hz**. Significance via the analytic threshold `1 − α^(1/(K−1))` for K Welch segments. | Lecci 2017; Osorio-Forero 2021 |
| **3B** | Does heart rate shift systematically **around the slow-oscillation trough**? | SO down-state-triggered average HR (SO 0.15–4 Hz, RR at 4 Hz cubic spline, ±5 s), reported as **% above stage-mean HR** and **SO→HR peak latency**. *(Superseded method: modulation depth vs a random-trigger null.)* | **Naji 2019** — *not* the 2022 review |
| **3D** | Does the slow oscillation **organise spindles**? | **Event-locked**: detect SO and spindle events per channel, Rayleigh test on the SO phase at each spindle peak. *(Superseded method: continuous Tort modulation index over all timepoints — this produced a false null.)* | Staresina 2015; Helfrich 2018 |

**Staging caveat.** iEEG has no EOG/EMG, so **true AASM N2/N3 scoring is impossible**. NREM epochs
were split by a 2-component Gaussian mixture on log slow-wave power — high-SWA = "N3-like",
low-SWA = "N2-like". This is the physiological axis the N2/N3 boundary tracks, but it is **not**
scored staging, hence "-like" throughout.

## 4. Results

### 3A — infraslow spindle↔heart coupling: **negative**

Pooled NREM (block length-capped to 55 min so the ~1/K coherence floor is comparable across
subjects): **7 of 23** subjects exceeded their own α=0.05 threshold at 0.02 Hz, against ~1.4
expected. That looked like a strong positive (binomial p = 3×10⁻⁴).

It does not survive the control. Counting exceedances at **every** frequency bin:

| | exceedance |
|---|---|
| 0.02 Hz (pre-specified) | 7/23 = 0.30 |
| 0.05 Hz | **8/23 = 0.35** ← higher |
| median across 0.005–0.25 Hz | 3/23 = **0.13** |
| mean 0.005–0.05 Hz → 0.15–0.25 Hz | 0.21 → 0.07 (broad slope, no peak) |

0.02 Hz ranks **2nd of 63 bins**, and the background exceedance (11.6%) is ~2× the 6% that a
simulation on independent signals predicts. Retested against the *empirical* background rather than
the simulated one, 0.02 Hz gives **p = 0.013** — and it is not the largest bin.

#### Where the 6% and 11.6% come from

These two numbers decide the retraction, so they are worth stating precisely.

- **6% — what chance alone produces.** `validate_3A_false_positive_rate.py` generates **2000 pairs of
  independent synthetic signals** (no coupling by construction) in four spectral shapes — white,
  pink 1/f, brown 1/f², and AR1 ρ=0.99 (strongly autocorrelated, the condition that normally inflates
  coherence) — each 55 min at 1 Hz, pushed through the **identical `msc_block()` code path**. The rate
  at which coherence at 0.02 Hz crosses the analytic threshold is **0.055 / 0.063 / 0.065 / 0.059**.
  So the test delivers ≈ its nominal α = 0.05 *when the two signals share nothing*.
  (The same simulation shows the **band-maximum** test gives **34%** false positives — which is why
  only a pre-specified frequency point is reported anywhere in this work.)

- **11.6% — what the real data produces where no effect is predicted.** For each of the 23 subjects
  the **full** coherence spectrum was computed, and at **every** frequency bin the number of subjects
  exceeding their own threshold was counted. Averaged over the 63 bins spanning 0.005–0.25 Hz, that
  is **11.6%** (≈2.7 of 23 per bin).

The inference: if only chance were operating, the real data would also sit near 6%. It sits at ~2×
that, so genuine shared low-frequency structure exists between sigma power and heart rate — but it is
present at *arbitrary* frequencies. **11.6% is therefore the correct null against which 0.02 Hz must
be judged, and against it the result is neither strong nor the largest bin.** The simulation
established that the *statistic* was sound; it could not establish that the *data* were free of
generic shared structure.

**Interpretation:** real sigma-power and heart-rate series share genuine **broadband low-frequency
structure** (arousals, state changes, drift), strongest at the lowest frequencies and decaying. There
is no ~50 s peak. **The apparent 3A finding is retracted.**

### 3B — SO→heartbeat coupling: **present, but ~25× weaker than published**

**The cited source was wrong.** `10.1073/pnas.2123417119` (Chen, Zhang, Thayer & Mednick 2022) is a
**review** — no heart-rate-burst analysis, no SO-trough method, and its availability statement reads
*"There are no data underlying this work."* The method had been improvised from a one-line
description of it.

The actual methods paper is
[Naji, Krishnan, McDevitt, Bazhenov & Mednick 2019, *J Cogn Neurosci*](https://doi.org/10.1162/jocn_a_01432).
Reimplemented to follow it (`analysis/event_3B_mednick.py`): SO band **0.15–4 Hz**, per-channel
zero-crossing half-wave detection, RR resampled at **4 Hz by piecewise cubic spline**, HR averaged in
a **±5 s** window on the SO down-state trough, reported as **% above that stage's mean HR** plus the
**SO→HR peak latency** — the paper's headline statistic, which the previous version never computed.

| | this work (HUP165) | Naji 2019 (frontal scalp, healthy) |
|---|---|---|
| N2 HR peak | **+0.47%** | **+12.09 ± 1.48%** |
| N3 / SWS HR peak | **+0.48%** | **+3.35 ± 1.01%** |
| SO→HR peak lag | 2.2 s | peak follows the down-state ✓ |
| significance | z = 4.9 (N3), z = −1.8 (N2) | — |

So the coupling **exists** — 9,381 slow oscillations, z = 4.9, with a physiologically sensible ~2.2 s
lag — but at roughly **1/25th** the published magnitude, and Naji's clear **N2 ≫ SWS** pattern does
**not** reproduce (0.47% vs 0.48%).

The most likely explanation is not fixable in code: **Naji recorded frontal scalp EEG (F3/F4), this
uses lateral neocortical iEEG.** Scalp electrodes see large, globally synchronous slow oscillations;
a single lateral intracranial contact sees local ones, and it is the global SOs that plausibly drive
autonomic coupling. Epilepsy patients on anti-seizure medication at ~92 bpm compound it.

### 3D — SO→spindle coupling: **PRESENT** *(corrected)*

**Originally reported null** (N2 MI_z 8.49 vs N3 2.73, p = 0.40; raw Tort MI ~6×10⁻⁵, "negligible").
**That was a method artifact.** The implementation binned *every* timepoint by SO phase into a
continuous modulation index, with no event detection, on a channel-averaged signal — so real events
were averaged together with hours of neither-spindle-nor-SO, driving the estimate toward zero.

Neither source paper does it that way. Helfrich 2018, verbatim: *"we detected SO (0.16–1.25 Hz) and
sleep spindle (12–16 Hz) **events** … phase during the **peak of the detected sleep spindle events**
… Rayleigh z"*. Staresina 2015 uses an *"Event-locked analysis"*.

Re-implemented faithfully (`analysis/event_3D_by_stage.py` — per channel, discrete SO and spindle
events by amplitude threshold, Rayleigh test on the SO phase at each spindle peak):

| | result |
|---|---|
| subjects with ≥1 significant channel | **21 / 23** |
| median fraction of channels significant | **83%** |
| median resultant vector length R | **0.073** |
| **N2-like vs N3-like** (quality-filtered, paired) | median R 0.071 vs 0.069, n = 20, **p = 0.29** — null; fraction significant 0.83 vs 0.58, p = 0.16 |

R ≈ 0.07 agrees with the independent estimate from `master` (R ≈ 0.03–0.05, significant in 27/28
channels) obtained by a separate analysis path — so the two studies are consistent, and the apparent
conflict was entirely the discarded method.

### N2-like vs N3-like: **null in all three tests**

Also null when split by fast (>12 Hz) vs slow (<12 Hz) individual spindle peak (p = 0.27–0.77).

## 5. Conclusions

1. **No evidence for the LC infraslow fingerprint** (~0.02 Hz shared spindle/heart rhythm) in human
   intracranial recordings from 23 epilepsy patients. The one apparently positive result did not
   survive a frequency-specificity control.
2. **SO→spindle coupling (3D) IS present** — 21/23 subjects, median 83% of channels significant,
   R ≈ 0.073 — agreeing with the independent estimate on `master`. This *reverses* an earlier null in
   this document, which was caused by not following the source papers' event-based method.
3. **SO→heartbeat coupling (3B) is present but ~25× weaker** than Naji 2019's frontal-scalp values,
   and their N2 ≫ SWS pattern does not reproduce. Region (lateral iEEG vs frontal scalp) is the most
   likely cause.
4. **No N2-like vs N3-like difference in any test** (3A n/a; 3B 0.47% vs 0.48%; 3D p = 0.29).
5. This is a **negative result for the LC infraslow hypothesis in epilepsy patients**, not a
   refutation of Lecci/Osorio-Forero, who worked in healthy sleepers with scalp EEG and direct LC
   recordings in mouse.

### Method-verification audit

Each test was checked against its source paper. This mattered — two of three deviated:

| Test | Follows source? | Consequence |
|---|---|---|
| **3A** | Partially — Lecci uses 4-s epochs, 10–15 Hz sigma, and **cross-correlation**; this used 1-s bins, individual spindle peak ±1 Hz, and **coherence**. But the decisive test (spectrum of the sigma-power time course, Lecci Fig 1C) was run their way. | negative stands |
| **3B** | ❌ No — cited a **review**; method improvised. Now reimplemented per Naji 2019. | verdict changed |
| **3D** | ❌ No — continuous MI over all timepoints instead of event-locked circular statistics. Now reimplemented per Staresina/Helfrich. | **verdict reversed** |

## 6. Limitations

- **Epilepsy patients on anti-seizure medication**, many tachycardic during NREM (~92 bpm in HUP165),
  which may blunt autonomic modulation.
- **Staging is not scored** — no EOG/EMG (see §3).
- **Slow oscillations were detected on the channel average**, which attenuates them when contacts are
  not synchronous. This plausibly suppressed both 3B and 3D and should be redone per channel.
- **No respiration/SpO₂ channel** exists, so apnoea-driven heart-rate swings cannot be excluded.
- **Electrode locations are inferred from contact numbering**, not localisation files.
- **Stage-resolved 3A rests on only 4 subjects** — N2 and N3 alternate faster than a 0.02 Hz
  coherence estimate needs.

## 7. What is nonetheless established

- **25 public iEEG subjects with simultaneous depth electrodes and EKG** were identified — the
  enabling resource ([DATA_INVENTORY_LC_INFRASLOW.md](DATA_INVENTORY_LC_INFRASLOW.md)).
- A **calibrated pipeline**: the analytic coherence threshold delivers 5.5–6.5% false positives on
  independent signals, while the band-maximum test gives **34%** — so only a pre-specified frequency
  point is defensible ([validate_3A_false_positive_rate.py](../analysis/validate_3A_false_positive_rate.py)).
- Documented methodological traps ([DEC_3A_METHOD_NOTES.md](DEC_3A_METHOD_NOTES.md)), notably that
  shift/phase surrogates are **invalid** for oscillatory coherence, and that band-passing before
  computing coherence destroys the estimate.

## 8. Next step

**OpenNeuro ds003848** — the one public iEEG sleep dataset with a verified ECG channel (n≈6,
1 h continuous, 2048 Hz), recorded with **EOG and EMG**, so real sleep staging is possible, on
cortical grids at an independent site. If the effect is absent there too, this becomes a clean,
well-controlled negative.

---
---

## 9. References

### The tests

**3A — infraslow spindle↔heart coupling**
- Lecci S, Fernandez LMJ, Weber FD, Cardis R, Chatton J-Y, Born J, Lüthi A (2017). *Coordinated
  infraslow neural and cardiac oscillations mark fragility and offline periods in mammalian sleep.*
  **Science Advances** 3:e1602026. [10.1126/sciadv.1602026](https://doi.org/10.1126/sciadv.1602026)
  — defines the 0.02 Hz sigma-power oscillation and its cardiac counterpart; source of the ~50 s
  periodicity, the S2 > SWS claim, and the fast-spindle-peak band.
- Osorio-Forero A, Cardis R, Vantomme G, Guillaume-Gentil A, Katsioudi G, Devenoges C, Fernandez LMJ,
  Lüthi A (2021). *Noradrenergic circuit control of non-REM sleep substates.* **Current Biology**
  31:5009–5023. [10.1016/j.cub.2021.09.041](https://doi.org/10.1016/j.cub.2021.09.041)
  — shows the LC drives the infraslow spindle-clustering rhythm and coordinates infraslow heart rate.

**3B — slow-oscillation ↔ heartbeat coupling**
- Chen P-C, Zhang J, Thayer JF, Mednick SC (2022). *Understanding the roles of central and autonomic
  activity during sleep in the improvement of working memory and episodic memory.* **PNAS**
  119(44):e2123417119. [10.1073/pnas.2123417119](https://doi.org/10.1073/pnas.2123417119)
  — the event-based autonomic–central coupling framework this test follows.

**3D — slow-oscillation → spindle coupling**
- Staresina BP, Bergmann TO, Bonnefond M, van der Meij R, Jensen O, Deuker L, Elger CE, Axmacher N,
  Fell J (2015). *Hierarchical nesting of slow oscillations, spindles and ripples in the human
  hippocampus.* **Nature Neuroscience** 18:1679–1686. [10.1038/nn.4119](https://doi.org/10.1038/nn.4119)
- Helfrich RF, Mander BA, Jagust WJ, Knight RT, Walker MP (2018). *Old brains come uncoupled in
  sleep: slow wave–spindle synchrony, brain atrophy and forgetting.* **Neuron** 97:221–230.
  [10.1016/j.neuron.2017.11.020](https://doi.org/10.1016/j.neuron.2017.11.020)
- Winer JR, Mander BA, Helfrich RF, Maass A, Harrison TM, Baker SL, Knight RT, Jagust WJ, Walker MP
  (2019). *Sleep as a potential biomarker of tau and β-amyloid burden in the human brain.*
  **J Neurosci** 39:6315–6324.
  [10.1523/JNEUROSCI.0503-19.2019](https://doi.org/10.1523/JNEUROSCI.0503-19.2019)
  — the disease link motivating the whole study: impaired SO–spindle coupling predicts medial
  temporal tau.

### Why the LC is the target
- Jacobs HIL et al. (2021). *In vivo and neuropathology data support locus coeruleus integrity as
  indicator of Alzheimer's disease pathology and cognitive decline.* **Sci Transl Med** 13:eabj2511.
  [10.1126/scitranslmed.abj2511](https://doi.org/10.1126/scitranslmed.abj2511)

### Methods
- Tort ABL, Komorowski R, Eichenbaum H, Kopell N (2010). *Measuring phase-amplitude coupling between
  neuronal oscillations of different frequencies.* **J Neurophysiol** 104:1195–1210. — the modulation
  index used in 3D.
- Welch PD (1967). *The use of FFT for the estimation of power spectra.* **IEEE Trans Audio
  Electroacoust** 15:70–73. — the segmented spectral estimator underlying 3A.
- Halliday DM, Rosenberg JR, Amjad AM, Breeze P, Conway BA, Farmer SF (1995). *A framework for the
  analysis of mixed time series/point process data.* **Prog Biophys Mol Biol** 64:237–278. — source
  of the analytic coherence confidence limit `1 − α^(1/(K−1))`.
- Pan J, Tompkins WJ (1985). *A real-time QRS detection algorithm.* **IEEE Trans Biomed Eng**
  32:230–236. — R-peak detection (used before switching to NeuroKit2).
- Makowski D, Pham T, Lau ZJ, et al. (2021). *NeuroKit2: a Python toolbox for neurophysiological
  signal processing.* **Behavior Research Methods** 53:1689–1696. — the validated R-peak detection
  actually used.

### Data
- Pattnaik AR, Litt B et al. — Normative iEEG Sleep & Wake Atlas, Pennsieve
  [10.26275/xhte-d11l](https://doi.org/10.26275/xhte-d11l)
- iEEG.org (IEEG Portal) — HUP phase-II continuous recordings; the cohort used here.
- OpenNeuro **ds003848** — the replication target (iEEG sleep + verified ECG + EOG/EMG).

*PDFs of the primary references are in the companion literature repository under
`papers/lc_infraslow_sleep_coupling/`.*


---

*Scripts: `analysis/cohort_stages_3ABD.py` (3A/3B/3D per stage), `analysis/frequency_specificity_3A.py`
(the control that retracted 3A), `analysis/summarize_cohort_stages.py`,
`analysis/validate_3A_false_positive_rate.py`. Per-subject JSON in `outputs/cohort_stages_3ABD/` and
`outputs/freq_specificity_3A/`.*