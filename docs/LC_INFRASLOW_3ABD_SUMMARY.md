# Testing the locus-coeruleus infraslow hypothesis in human iEEG — summary

**Bottom line: negative.** Across 23 analysed subjects we find no evidence that spindle activity and
heart rate share the ~0.02 Hz (~50 s) locus-coeruleus fingerprint, and no N2/N3 difference in any of
the three tests. One apparently positive result (3A) did not survive a frequency-specificity control
and is retracted below.

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

| Test | What it asks | Measure |
|---|---|---|
| **3A** *(primary)* | Do spindle power and heart rate rise and fall **together every ~50 s**? | Magnitude-squared coherence between the 1 Hz spindle-power series and instantaneous heart rate, read at the pre-specified **0.02 Hz**. Significance via the analytic threshold `1 − α^(1/(K−1))` for K Welch segments. |
| **3B** | Does heart rate shift systematically **around the slow-oscillation trough**? | SO-trough-triggered average heart rate; modulation depth vs a random-trigger null (z). |
| **3D** | Does the slow oscillation **organise spindles** (context/QC)? | Tort modulation index between SO phase and spindle amplitude, surrogate-corrected. |

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

**Interpretation:** real sigma-power and heart-rate series share genuine **broadband low-frequency
structure** (arousals, state changes, drift), strongest at the lowest frequencies and decaying. There
is no ~50 s peak. **The apparent 3A finding is retracted.**

### 3B — SO→heartbeat coupling: **null**

N2-like z = 0.70 vs N3-like z = 1.89 (paired, n = 17, Wilcoxon **p = 0.11**). Direction favours
N3-like but does not reach significance.

### 3D — SO→spindle coupling: **null, and effect size negligible**

N2-like MI_z = 8.49 vs N3-like MI_z = 2.73 (paired, n = 17, **p = 0.40**).

More importantly, **raw Tort MI is ~6×10⁻⁵** (median; MI is normalised 0–1) — negligible coupling.
Large MI_z values (up to 212) are an artifact of huge sample counts shrinking the surrogate standard
deviation. **Effect size, not z, is the number to quote.**

### N2-like vs N3-like: **null in all three tests**

Also null when split by fast (>12 Hz) vs slow (<12 Hz) individual spindle peak (p = 0.27–0.77).

## 5. Conclusions

1. **No evidence for the LC infraslow fingerprint** (~0.02 Hz shared spindle/heart rhythm) in human
   intracranial recordings from 23 epilepsy patients.
2. **No N2/N3 difference** in infraslow coupling, SO→heartbeat coupling, or SO→spindle coupling.
3. **SO→spindle coupling is negligible** on lateral neocortical contacts in this pipeline
   (raw MI ~6×10⁻⁵) — note this contrasts with the robust SO→spindle coupling this repo reports in
   *mesial temporal* contacts, and is likely partly a methodological artifact (see limitations).
4. This is a **negative result in epilepsy patients**, not a refutation of Lecci/Osorio-Forero, who
   worked in healthy sleepers with scalp EEG and (in mouse) direct LC recordings.

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

*Scripts: `analysis/cohort_stages_3ABD.py` (3A/3B/3D per stage), `analysis/frequency_specificity_3A.py`
(the control that retracted 3A), `analysis/summarize_cohort_stages.py`,
`analysis/validate_3A_false_positive_rate.py`. Per-subject JSON in `outputs/cohort_stages_3ABD/` and
`outputs/freq_specificity_3A/`.*
