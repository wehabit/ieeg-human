# The `lc-infraslow-3ABD` study — what we did and how

## The overarching goal: is heart rate a readable proxy for the locus coeruleus?

The locus coeruleus (LC) degenerates early in Alzheimer's, so a non-invasive LC readout would be a
valuable biomarker. The proposed readout is indirect: the LC is claimed to impose a shared
**~0.02 Hz (~50 s) infraslow rhythm** on both sleep-spindle (sigma) power and heart rate during NREM.
If that holds, heart rate reads out the LC.

We tested this in **human intracranial EEG** — 25 Penn epilepsy patients (iEEG.org HUP `phaseII`) who
carry simultaneous depth electrodes *and* an EKG channel, which public iEEG almost never has. Three
questions, each motivated by a specific paper.

---

## 3A — Do spindle power and heart rate oscillate together every ~50 s?

**Motivation:** [Lecci 2017, *Sci Adv*](https://doi.org/10.1126/sciadv.1602026) — the 0.02 Hz
sigma/cardiac oscillation; [Osorio-Forero 2021, *Curr Biol*](https://doi.org/10.1016/j.cub.2021.09.041)
— the LC drives it. This is the load-bearing test for "read the LC off the heartbeat."

**How we implemented the answer** (Lecci's two-step method):

1. **Build the two time series.** Per subject, stream ~7 h from the highest-delta night. From the
   lateral neocortical contacts: notch → band-pass at the individual fast-spindle peak (±1 Hz) →
   Hilbert envelope → **sigma-power series at 1 Hz**, IED-masked. From the EKG: NeuroKit2 R-peaks →
   RR intervals → **instantaneous heart rate**.
2. **Step 1 — is a ~50 s rhythm even present?** A duration-weighted Morlet spectrum of the sigma-power
   series over every NREM bout ≥120 s, with a Gaussian peak fit — one fitted infraslow peak per
   subject. Cohort test: do the peaks **cluster** at ~0.019 Hz (a real shared rhythm) or **scatter**
   like 1/f noise? Compared against each subject's own scale-free surrogates.
3. **Step 2 — does heart rate track it?** Magnitude-squared coherence between sigma power and heart
   rate (read at 0.02 Hz and at each subject's own peak), plus the **cross-correlation** (z-scored
   120 s intervals, HR as source wave) — the lag-resolved statistic that is Lecci's actual coupling
   measure.

**Answer: no.** The peaks scatter like noise (peak SD 0.014 Hz ≈ this cohort's surrogates 0.014;
0/… clustering at 0.019), sigma is not more prominent than the SWA control (p = 0.11), coherence at
0.02 Hz ranks 10th of 63 bins, and the cross-correlation is flat (group |r| = 0.036, p = 0.23).
**Heart rate does not carry the LC's proposed fingerprint.**

---

## 3B — Does heart rate shift systematically around the slow-oscillation trough?

**Motivation:** [Naji 2019, *J Cogn Neurosci*](https://doi.org/10.1162/jocn_a_01432) — timing between
cortical slow oscillations and heart-rate bursts.

**How we implemented the answer:** Detect slow oscillations per channel (0.15–4 Hz, zero-crossing
negative half-waves with duration + amplitude criteria). Resample RR to 4 Hz by cubic spline. Average
heart rate in a ±5 s window locked to each SO down-state trough, expressed as **% above that stage's
mean HR**, with the SO→HR latency. Significance against a **stage-matched** random-trigger null
(surrogate triggers drawn from the same sleep stage as the real ones). Run per stage (N2 vs N3),
across all 23 subjects.

**Answer: weak and heterogeneous.** Significant in only 7/23 (N2) and 6/20 (N3); the N3 magnitude
(+3.5%) matches Naji's SWS but N2 (+1.8%) is ~7× weaker, and Naji's N2 ≫ SWS pattern does not
reproduce (paired p = 0.57).

---

## 3D — Does the slow oscillation organise spindles?

**Motivation:** [Staresina 2015, *Nat Neurosci*](https://doi.org/10.1038/nn.4119) ·
[Helfrich 2018, *Neuron*](https://doi.org/10.1016/j.neuron.2017.11.020) (SO→spindle nesting);
clinically [Winer 2019, *J Neurosci*](https://doi.org/10.1523/JNEUROSCI.0503-19.2019) (this coupling
predicts medial-temporal tau) and [Jacobs 2021, *Sci Transl Med*](https://doi.org/10.1126/scitranslmed.abj2511)
(LC integrity ↔ AD). This is the disease link that makes the study matter.

**How we implemented the answer:** Event-locked circular statistics — detect SO and spindle events per
channel, take the SO phase at each spindle peak, Rayleigh test per channel.

**Answer: present.** 21/23 subjects, median 83% of channels significant, R ≈ 0.073 — agreeing with the
independent estimate on the `master` branch. No N2/N3 difference.

---

## Staging (common to all three)

iEEG has no EOG/EMG, so NREM and the N2/N3 split were derived from a Gaussian mixture on slow-wave
power — a physiological proxy, labelled "-like" throughout, not scored AASM staging.

---

## Independent replication — OpenNeuro ds003848

To test whether the 3A negative was merely a staging artifact (residual wake/REM smearing a real
rhythm), we ran the **same** 3A/3B pipeline on ds003848 (Utrecht RESPect, n=6) — the one public iEEG
sleep dataset with **ECG + EMG + EOG**, which let us stage with **real REM/wake exclusion**. Staging
was built from EMG tone + EOG movement + iEEG slow-wave power, validated to separate Wake/REM/N2/N3 as
the physiology predicts.

**Answer:** consistent with the HUP negative — excluding REM/wake did **not** rescue an infraslow
rhythm, and 3B was null. Directional, not definitive: the 1-hour recordings, fragmented by real
staging, are underpowered for 3A (coherence K median 7 vs 59 in HUP). Full write-up:
[DS003848_REPLICATION.md](DS003848_REPLICATION.md).

---

## Overall

Across a 23-patient cohort and a properly-staged 6-patient replication, the LC's proposed ~50 s
spindle↔heart fingerprint (3A) is **not detectable** in human iEEG, and the SO→heartbeat coupling (3B)
is weak — so **heart rate is not a validated LC proxy in this population**. The cortical SO→spindle
coupling (3D) that links to tau *is* present. This is a negative **in epilepsy patients**, not a
refutation of Lecci/Osorio-Forero, who worked in healthy sleepers with scalp EEG and direct LC
recordings in mouse.

*Full method/results write-up: [LC_INFRASLOW_3ABD_SUMMARY.md](LC_INFRASLOW_3ABD_SUMMARY.md). Corrected
cohort numbers: `outputs/corrected_3AB/COHORT_SUMMARY.txt` (HUP) and `DS003848_SUMMARY.txt`
(replication).*
