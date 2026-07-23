# Testing the locus-coeruleus infraslow hypothesis in human iEEG — summary

**Bottom line: the ~0.02 Hz locus-coeruleus fingerprint is not present**, and there is no N2/N3
difference in any test.

> **Method correction (2026-07).** An audit found real bugs in the first version of this analysis.
> The code was rebuilt, every fix ships with a test that fails on the old code and passes on the new
> (`analysis/test_*.py`), and the corrected cohort numbers live in
> `outputs/corrected_3AB/COHORT_SUMMARY.txt`. The negative conclusion is unchanged in direction but is
> now **stronger and correctly specified**:
>
> 1. **3A NaN handling** deleted missing seconds and spliced the survivors, compressing the time axis
>    so a true 0.02 Hz rhythm could read up to 0.029 Hz (median 8.7%, max 30.8% of seconds dropped).
>    Replaced with a gap-aware coherence estimator (`analysis/spectral_gapped.py`), calibrated to a
>    5% false-positive rate at the real analysis geometry (K≈107, 34 bouts).
> 2. **3A now follows Lecci** instead of testing one hard-coded 0.0195 Hz coherence bin. It fits each
>    subject's own infraslow spectral peak over all NREM bouts (Lecci Fig 1G) and runs the
>    **cross-correlation** with heart rate (Lecci Fig 6) — the paper's actual coupling statistic,
>    which had never been computed here (`analysis/lecci_faithful_3A.py`).
> 3. **3B's surrogate null was stage-blind** — it drew surrogate triggers from the whole night while
>    scoring "% above *this stage's* mean HR", so the null measured the stage-vs-night HR offset, not
>    SO-locking. On zero-effect synthetic data that alone produced z = −21.8 / +24.0. Fixed to a
>    stage-matched null and run on the **whole cohort (n=23)** rather than one subject
>    (`analysis/event_3B_cached.py`). The old single-subject "z = 4.9" was an artifact of this bug.
>
> The Results sections below have been updated to the corrected numbers. Two further corrections were
> made earlier and remain: **3D** was reported null and is now **positive** (an artifact of not using
> the source papers' event-based method), and 3B's cited reference turned out to be a review, not a
> methods paper.

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
| **3A** *(primary)* | Do spindle power and heart rate rise and fall **together every ~50 s**? | **Step 1:** per-subject infraslow Morlet spectrum of the sigma-power time course over all NREM bouts ≥120 s, Gaussian peak fit — do peaks cluster at ~0.019 Hz (Lecci Fig 1G)? **Step 2:** gap-aware coherence (read at 0.02 Hz and each subject's own peak) **and** cross-correlation with heart rate (Lecci Fig 6). *(Superseded: a single hard-coded 0.0195 Hz coherence bin on a delete-and-spliced 55-min block.)* | Lecci 2017; Osorio-Forero 2021 |
| **3B** | Does heart rate shift systematically **around the slow-oscillation trough**? | SO down-state-triggered average HR (SO 0.15–4 Hz, RR at 4 Hz cubic spline, ±5 s), % above stage-mean HR + SO→HR latency, against a **stage-matched** random-trigger null. *(Superseded nulls: whole-night triggers; modulation-depth vs random.)* | **Naji 2019** — *not* the 2022 review |
| **3D** | Does the slow oscillation **organise spindles**? | **Event-locked**: detect SO and spindle events per channel, Rayleigh test on the SO phase at each spindle peak. *(Superseded method: continuous Tort modulation index over all timepoints — this produced a false null.)* | Staresina 2015; Helfrich 2018 |

**Staging caveat.** iEEG has no EOG/EMG, so **true AASM N2/N3 scoring is impossible**. NREM epochs
were split by a 2-component Gaussian mixture on log slow-wave power — high-SWA = "N3-like",
low-SWA = "N2-like". This is the physiological axis the N2/N3 boundary tracks, but it is **not**
scored staging, hence "-like" throughout.

## 4. Results

### 3A — infraslow spindle↔heart coupling: **negative** (Lecci's own two-step method)

Rebuilt to follow Lecci rather than test one hard-coded coherence bin (`lecci_faithful_3A.py`, n=23),
pooling **all** NREM bouts ≥120 s (coherence K median **59**, versus 7–23 on the original single
55-min block) with the gap-aware estimator.

**Step 1 — is there a ~50 s sigma rhythm at all? (Lecci Fig 1G)** Per subject, a duration-weighted
Morlet spectrum of the sigma-power time course over every NREM bout ≥120 s, with a Gaussian peak fit.
Individual subjects *do* show infraslow bumps (some very strong — HUP212 prominence 173, HUP157 22.9),
but at **scattered** frequencies. The decisive cohort question is whether the 23 fitted peaks cluster
at ~0.019 Hz like Lecci's (a real shared rhythm → simulated peak SD ≈ 0.0008 Hz) or scatter like 1/f
noise (SD ≈ 0.013 Hz):

| | value | expected if… |
|---|---|---|
| observed peak SD | **0.0137 Hz** | real rhythm ≈ 0.0008; noise ≈ 0.013 |
| this cohort's own scale-free surrogates (peak SD) | 0.0139 Hz | — |
| bootstrap: real peaks tighter than surrogates? | **p = 0.49** | — |
| fraction of peaks in 0.015–0.025 Hz | 0.30 | true-in-band ≈ 0.94; noise ≈ 0.25 |
| sigma more prominent than the SWA control (Lecci's own control) | Wilcoxon **p = 0.11** | — |

The peaks are indistinguishable from where 1/f noise happens to bump. **There is no shared ~50 s sigma
rhythm at the cohort level.** (The per-subject scale-free significance control is deliberately reported
but is weak — a single ~2 h night gives a noisy prominence estimate; on synthetic data even a planted
rhythm scores p ≈ 0.21. Peak *clustering* across subjects, above, is the powered test.)

**Step 2 — does heart rate track it? (Lecci Fig 6)** Two independent controls agree it does not.

- **Frequency-specificity** (the recomputed version of the control that retracted the original 3A):
  0.02 Hz is significant in 6/23 subjects, but ranks only **10th of 63 bins** against a 16.5%
  background (binomial p = 0.17). Higher K here lowers each subject's threshold, so more bins cross —
  but 0.02 Hz is even less special than before.
- **Cross-correlation** — *Lecci's actual coupling statistic, which had never been run in this work*:
  z-scored 120 s intervals, HR as source wave, averaged within then across subjects. Group |r| =
  **0.036** (lag −1 s); per-subject median |r| = 0.056, with lags scattered (IQR −10 to −1 s) and no
  consistent sign (t on peak r vs 0: **p = 0.23**). A positive result would be a clear correlogram
  peak at a consistent lag; this is flat.

**Interpretation:** with Lecci's own method, better powered than the original, sigma power and heart
rate do not share a ~50 s rhythm in these patients — neither signal reliably carries the oscillation,
and they do not track each other. This *confirms and strengthens* the original negative, which had
reached the right conclusion through a mis-specified single-bin test.

<details><summary>Original single-bin analysis (superseded, retained for the record)</summary>

The first version tested one hard-coded bin at 0.0195 Hz on a single 55-min block: 7/23 subjects
exceeded threshold (binomial p = 3×10⁻⁴ against a calibrated 6% false-positive rate), but a
frequency-specificity control put 0.02 Hz only **2nd of 63 bins** (0.05 Hz scored 8/23) against an
**11.6%** empirical background — retested against which 0.02 Hz gave p = 0.013 and was not the largest
bin. Right conclusion, but the test was not Lecci's and used the delete-and-splice coherence that the
gap-aware estimator replaces. `frequency_specificity_3A.py` / `validate_3A_false_positive_rate.py`.
</details>

### 3B — SO→heartbeat coupling: **weak and heterogeneous** (n=23, stage-matched null)

**Two things were fixed since the single-subject version.** (1) The cited source was wrong:
`10.1073/pnas.2123417119` (Chen, Zhang, Thayer & Mednick 2022) is a **review** with no such analysis;
the actual methods paper is
[Naji, Krishnan, McDevitt, Bazhenov & Mednick 2019, *J Cogn Neurosci*](https://doi.org/10.1162/jocn_a_01432)
(SO band 0.15–4 Hz, per-channel zero-crossing half-waves, RR at 4 Hz cubic spline, HR averaged ±5 s on
the down-state trough, reported as % above that stage's mean HR plus SO→HR latency). (2) The **surrogate
null was stage-blind** — random triggers were drawn from the whole night while the statistic was scored
against the *stage* mean HR, so the null measured the stage-vs-night HR offset rather than SO-locking.
On zero-effect synthetic data this alone produced z = −21.8 (N2) / +24.0 (N3), and it fully explains
the old single-subject "z = 4.9". Fixed to a stage-matched null and run on the whole cohort
(`event_3B_cached.py`).

*(Naji scored with Rechtschaffen & Kales, whose deep-sleep stage is **SWS** = R&K stages 3+4; AASM
merged these into **N3**, so **SWS = N3** — the same stage, different scoring vocabulary.)*

| | corrected cohort (n=23) | Naji 2019 (frontal scalp, healthy) |
|---|---|---|
| N2 (Stage 2) HR peak | **+1.8%** (SD 8.2) → ~7× weaker | **+12.09 ± 1.48%** |
| N3 (= SWS) HR peak | **+3.5%** (SD 12.7) → **matches** | **+3.35 ± 1.01%** |
| significant per subject (z > 1.96) | **7/23** (N2), **6/20** (N3) | — |
| cohort z-test (t on z vs 0) | N2 **p = 0.014**, N3 **p = 0.020** | — |
| N2 ≫ N3 (SWS)? (Naji's headline pattern) | paired **p = 0.57 — not reproduced** | 3.6× |
| SO→HR lag | ~1–2 s (follows the down-state ✓) | follows the down-state ✓ |

So SO→heartbeat coupling is **real but weak and present in only a minority of subjects** — the cohort
z-test is significantly nonzero, but the median subject shows almost nothing (median z ≈ 0.5–0.7) and
a handful of subjects (HUP139 z=7.8, HUP150 z=10.0, HUP151 z=7.6, HUP143 z=6.7) carry the effect.
Notably the **N3 (= SWS) magnitude (+3.5%) lands right on Naji's SWS value (+3.35%)**, while N2 is far
weaker and Naji's N2 ≫ N3 pattern is absent.

The most likely reason the pattern does not fully match is not fixable in code: **Naji recorded frontal
scalp EEG (F3/F4), this uses lateral neocortical iEEG.** Scalp electrodes see large, globally
synchronous slow oscillations; a single lateral intracranial contact sees local ones, and it is the
global SOs that plausibly drive autonomic coupling. Epilepsy patients on anti-seizure medication at
~92 bpm compound it. *(The earlier claim "present but ~25× weaker" was HUP165 alone, before the null
fix; it is superseded by these cohort numbers.)*

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
   intracranial recordings from 23 epilepsy patients, now on Lecci's own two-step method. The
   sigma-power infraslow peaks **scatter like 1/f noise** rather than clustering at 0.019 Hz (peak SD
   0.014 Hz ≈ this cohort's surrogates), and heart rate does not track them — the **cross-correlation
   is null** (median |r| 0.06, p = 0.23) and 0.02 Hz ranks 10th of 63 bins.
2. **SO→spindle coupling (3D) IS present** — 21/23 subjects, median 83% of channels significant,
   R ≈ 0.073 — agreeing with the independent estimate on `master`. This *reverses* an earlier null in
   this document, which was caused by not following the source papers' event-based method.
3. **SO→heartbeat coupling (3B) is weak and heterogeneous** — significant in only 7/23 (N2) and 6/20
   (N3) subjects, though the cohort z-test is nonzero (p ≈ 0.01–0.02). N3/SWS magnitude (+3.5%) matches
   Naji 2019, but N2 (+1.8%) is ~7× weaker and Naji's N2 ≫ SWS pattern does not reproduce (p = 0.57).
   Region (lateral iEEG vs frontal scalp) is the most likely cause. *(The prior "~25× weaker" figure
   was one subject before the stage-matched-null fix.)*
4. **No N2-like vs N3-like difference in any test** (3A: peaks scatter in both; 3B: N2 vs N3 paired
   p = 0.57; 3D p = 0.29).
5. This is a **negative result for the LC infraslow hypothesis in epilepsy patients**, not a
   refutation of Lecci/Osorio-Forero, who worked in healthy sleepers with scalp EEG and direct LC
   recordings in mouse.

### Method-verification audit

Each test was checked against its source paper, and the analysis was rebuilt where it deviated (see
the Method-correction note at the top). Every fix ships with a test in `analysis/test_*.py`.

| Test | Follows source? | Consequence |
|---|---|---|
| **3A** | ✅ Now yes — per-subject infraslow spectral peak fit over all NREM bouts (Lecci Fig 1G) **and** cross-correlation with heart rate (Lecci Fig 6), on a gap-aware coherence estimator. The original single-bin coherence test used delete-and-splice NaN handling, now replaced. | negative **strengthened** (better powered, Lecci's own statistic) |
| **3B** | ✅ Now yes — Naji 2019 method with a **stage-matched surrogate null** (the earlier whole-night null was invalid), run on the full cohort rather than n=1. | verdict refined: weak/heterogeneous, not "25× weaker" |
| **3D** | ✅ Now yes — event-locked circular statistics per Staresina/Helfrich, replacing a continuous MI over all timepoints. | **verdict reversed** (null → present) |

## 6. Limitations

- **Epilepsy patients on anti-seizure medication**, many tachycardic during NREM (~92 bpm in HUP165),
  which may blunt autonomic modulation.
- **Staging is not scored** — no EOG/EMG in HUP iEEG, so NREM is the GMM-on-slow-wave-power proxy and
  cannot exclude REM or quiet wake. That residual wake/REM contamination is exactly the broadband
  low-frequency structure that could smear a real infraslow rhythm, so it is the single biggest threat
  to the 3A negative. **OpenNeuro ds003848 (EMG + EOG) is the planned test of this** — see §8.
- **The corrected 3B uses per-channel SO detection** (the cached series stores per-channel trough
  times); the older channel-average version, which attenuated non-synchronous SOs, is superseded. 3D
  on `master` should still be redone per channel.
- **No respiration/SpO₂ channel** in HUP, so apnoea-driven heart-rate swings cannot be excluded here
  (ds003848 *does* carry thoracic/abdominal belts, though marked bad).
- **Electrode locations are inferred from contact numbering**, not localisation files.
- **Stage-resolved 3A is not the primary analysis** — the corrected 3A pools all NREM (N2 and N3
  alternate faster than an infraslow spectral estimate needs); the N2/N3 contrast rests on the 3B/3D
  event-based tests.

## 7. What is nonetheless established

- **25 public iEEG subjects with simultaneous depth electrodes and EKG** were identified — the
  enabling resource ([DATA_INVENTORY_LC_INFRASLOW.md](DATA_INVENTORY_LC_INFRASLOW.md)).
- A **calibrated pipeline**: the analytic coherence threshold delivers 5.5–6.5% false positives on
  independent signals, while the band-maximum test gives **34%** — so only a pre-specified frequency
  point is defensible ([validate_3A_false_positive_rate.py](../analysis/validate_3A_false_positive_rate.py)).
- Documented methodological traps ([DEC_3A_METHOD_NOTES.md](DEC_3A_METHOD_NOTES.md)), notably that
  shift/phase surrogates are **invalid** for oscillatory coherence, and that band-passing before
  computing coherence destroys the estimate.

## 8. Independent replication — OpenNeuro ds003848 (done)

**OpenNeuro ds003848** (Utrecht RESPect, n=6, 1 h continuous @ 2048 Hz, ECG + EMG + EOG) was run to
test the biggest caveat above: that the 3A negative might be a staging artifact, since HUP NREM is a
GMM proxy that cannot exclude REM/wake. With **real EMG/EOG staging** that removes 6–41% of epochs as
wake/REM, both tests remain consistent with the HUP negative — no infraslow peak clustering at
0.019 Hz (fitted peaks sit at ~0.04 Hz, 0/5 in band), no HR coupling, and a null 3B. **So the negative
is not merely a staging artifact.**

The caveat cuts the other way too: the 1 h recordings, fragmented by real staging, leave very few NREM
bouts ≥120 s (coherence K median **7** vs 59 in HUP; only 2/6 subjects estimable), so ds003848 is a
**directional** replication, not a powered one. A definitive test needs longer continuous iEEG sleep
with cardiac + EOG/EMG, which no current public dataset provides. Full write-up:
**[DS003848_REPLICATION.md](DS003848_REPLICATION.md)**; numbers in
`outputs/corrected_3AB/DS003848_SUMMARY.txt`.

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