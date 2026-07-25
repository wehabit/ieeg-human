# Key Literature

> **SCOPE NOTE.** Literature descriptions may remain useful, but historical repository-output
> numbers on this page are not v8 3A/3B/3D results or evidence of LC tracking. See the
> [v8 QC-sensitivity report](QC_SENSITIVITY_RESULTS_2026-07.md) and
> [issue register](ISSUE_REGISTER_2026-07.md).

The primary papers do **not** provide a 70%, 74%, 75%, 80%, or 90% recording-coverage
qualification rule. Lecci uses artifact-free NREM bouts ≥120 s; Naji uses uninterrupted 3-minute
bins; Staresina/Helfrich use artifact-free NREM plus event morphology/amplitude criteria. Their
75th percentiles concern event amplitude, not usable-recording coverage.

Fact-checked reference notes for the human iEEG NREM-coordination work. Each dataset/method
entry below has been verified against the primary source (paper, dataset record, or patent) and,
where relevant, against this repo's own analysis outputs. Verification status is marked per section.

---

## Part 1 — Datasets & methods (verified)

### 1. Geva-Sagiv, Mankin, Eliashiv, … Nir & Fried (2023) — *Nature Neuroscience* 26:1100–1110
*doi 10.1038/s41593-023-01324-5 · PMID 37264156 · PMC10244181*

In 18 UCLA epilepsy patients implanted with Behnke–Fried macro–micro depth electrodes, the authors
ran a within-participant, two-night, counterbalanced experiment (undisturbed sleep vs sleep with
real-time closed-loop stimulation) to ask whether coordination between hippocampal ripples,
thalamocortical spindles, and cortical slow waves causally supports memory consolidation. MTL iEEG
was streamed at 2 kHz, band-passed 0.5–4 Hz, and thresholded adaptively (80 µV default, re-estimated
every 400 s) so that brief 50 ms bursts delivered to prefrontal white matter landed on the
slow-oscillation up-state — achieving a mean 241.3 ms delay to the MTL peak versus 373.3 ms in a
mixed-phase control arm that received identical stimulation without phase-locking. Synchronized
stimulation improved overnight recognition memory in 6/6 prefrontal-white-matter participants
(one-sided binomial ≈ 0.016), driven by fewer false alarms rather than more hits, with pairing
accuracy unaffected; it also raised 9–16 Hz spindle power brain-wide across 565 contacts (figure
caption P = 1.4×10⁻³⁹; the Results text prints P < 10⁻³⁰ — an internal inconsistency, so cite the
figure), lifted the fraction of non-MTL single units phase-locked to MTL slow waves from 34.0% to
50.0%, and — critically — increased ripple↔neocortical-slow-wave co-occurrence (P = 8.2×10⁻⁴)
without changing ripple rate at all. That co-occurrence increase tracked memory gain at ρ = 0.8,
P = 0.007 (n = 30 pairs, 8 patients), while the stricter triple co-occurrence did not reach
significance (ρ = 0.7, P = 0.2, n = 12 pairs). Mixed-phase stimulation produced nothing, and
sometimes degraded performance.

**Why it matters here.** This is our intervention already executed in humans, with memory as the
readout. It establishes that (a) timing is the entire effect — same energy, wrong phase, no benefit —
and (b) the mechanism enhances coordination, not event count. Detection code is public
(github.com/mgevasagiv/); Fried, Nir & Geva-Sagiv hold **US Patent 12,654,013** on the method
(granted 2026-06-16; assignees UC Regents + Ramot/Tel Aviv University).

### 2. Pattnaik & Ong, … Litt — Normative iEEG Sleep & Wake Atlas (Pennsieve dataset 414, 2024; preprint 2025)
*Dataset DOI 10.26275/xhte-d11l (v1, Oct 2024, 4.36 GB, 3,164 files) · Preprint medRxiv 2025.09.10.25335533 (Sep 2025)*

This is a data resource, not a peer-reviewed finding. The Pennsieve dataset supplies processed iEEG
clips at 204 Hz with z-scores computed against separate wake and sleep normative atlases plus
deidentified electrode metadata, and the accompanying manuscript — "Normative intracranial EEG
highlights epileptic abnormalities across wakefulness and sleep" — remains a medRxiv preprint, not
peer-reviewed. The normative reference is built from 106 subjects and tested on 30, summarizing
spectral power and coherence across six frequency bands. Importantly, that 106-subject normative
baseline is largely the MNI Open iEEG Atlas (Frauscher et al. 2018) supplemented with normal HUP
channels — so it is **not** statistically independent of the MNI atlas. Released under
CC-BY-NC-SA-4.0.

**Why it matters here.** Our largest normative resource and the source of the wake→N2→N3 slow-power
result. Different analyses draw different subsets of this atlas: N = 21 for slow-power by state, N = 22
for slow→spindle coupling, n = 49 for per-patient nesting (53 N3/MTL patients total). Two constraints:
cite it as a dataset + preprint, and note it is the **only non-commercial-licensed dataset** we use.
At 204 Hz it cannot reach ripples (Nyquist 102 Hz) — spindles only.

### 3. Falach, Geva-Sagiv, Eliashiv, … Fried & Nir (2024) — *Scientific Data* 11:1354
*doi 10.1038/s41597-024-04187-y · dataset figshare 26131978*

A curated benchmark for interictal epileptiform discharge (IED) detection, not a sleep-physiology
study. It pools depth-electrode iEEG from 25 drug-resistant epilepsy patients across two centers
(Tel Aviv Sourasky 9, UCLA 16), 857 channels, acquired at 2 kHz and released downsampled to 1 kHz,
with 852 IEDs manually annotated by two expert neurologists working on non-overlapping patient sets.
Sleep was staged by AASM criteria from scalp EEG where available and by a validated automatic
algorithm for the 10 patients lacking it. The authors trained a LightGBM detector on
kurtosis/entropy/Teager–Kaiser features, reaching 94.3% accuracy, 94.4% precision and 94.3%
sensitivity under five-fold cross-validation, and 91.2%/86.7% under leave-one-out generalization. The
crucial practical fact — verified directly against our downloaded files — is that the release is
~76 minutes in total (76.7 min measured locally), a median of ~3.0 minutes per patient, excerpted
from overnight recordings rather than being overnight data. The **dataset is CC BY 4.0** on figshare;
the CC-BY-NC-ND licence applies only to the article text.

**Why it matters here.** Three minutes per patient is why our Falach SO→ripple estimate came back at
0.011 [−0.008, 0.033] — underpowered by design, never a biological null, and the docs now say so. The
dataset is deliberately enriched for epileptiform activity, but its expert IED labels are a gift: use
them to mask rather than trusting an automatic detector. Commercial use is permitted.

### 4. Fedele, Burnos, Boran, Krayenbühl, Hilfiker, Grunwald & Sarnthein (2017) — *Scientific Reports* 7:13836
*doi 10.1038/s41598-017-13064-1 · OpenNeuro ds003498 (CC0)*

The Zurich group recorded interictal intracranial EEG during slow-wave sleep in 20 epilepsy patients
and asked whether high-frequency oscillations mark epileptogenic tissue well enough to predict
surgical outcome in the individual patient — finding that resection of HFO-generating tissue
predicted seizure freedom. For our purposes the value is the expert/detector-validated event
markings: each bipolar channel carries labelled ripple (80–250 Hz) and fast-ripple (250–500 Hz)
events in BIDS `events.tsv` files, recorded at 2 kHz. Adam Li converted the release to BIDS and used
it to validate mne-hfo. Licensed CC0 (public domain) — the most permissive dataset we hold.

**Why it matters here.** Our only cohort with gold-standard, independently marked ripple events —
which is what makes the ripple-triggered slow-wave average non-circular (z > 2 in 54/63 channels,
median z = 7.92). Two caveats: the markings mix pathological with physiological HFOs (which is why
hippocampal ripples failed to phase-cluster while entorhinal ones did, MRL = 0.20), and we hold only
9 of the 20 subjects locally.

### 5. Demuru, van Blooijs, Zweiphenning, … Zijlmans & the RESPect group (2022) — *Neuroinformatics* 20:727–736
*doi 10.1007/s12021-022-09567-6 · OpenNeuro ds003848 v1.0.3*

[The primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC9440951/) presents a workflow for
organizing clinical intraoperative and long-term iEEG in BIDS; it is a dataset/methods paper, not
validation of an LC proxy or expert AASM/R&K N2/N3 scoring. In the pinned public RESPect snapshot,
the sleep-run `events.tsv` sidecars contain coarse sleep/transition, artifact, seizure, and curated
SWS/REM-selection intervals, while `electrodes.tsv` contains atlas and pathology fields.

**Why it matters here.** Those sidecars are analysis inputs, not optional descriptions. The legacy
RESPect cache ignored both and therefore could misclassify author REM, retain author-marked
disturbances, and include pathological/non-cortical contacts. V8 consumes them
conservatively, while leaving author-unknown sleep unclassified and treating parietal/frontal iEEG
ROIs as motivated scalp-source adaptations rather than direct LC measurements.

---

## Cross-cutting note

None of the Pattnaik, Falach, or Zurich holdings is a continuous full night: Pattnaik is clips,
Falach is ~3 min/patient, and Zurich is 5-min runs. The RESPect sleep runs used here are
approximately one hour and may provide too little strictly annotated continuous data for a cohort
endpoint. The historical HUP analyses used longer continuous high-delta candidate intervals, but
those intervals have not been verified as lights-off-to-wake nights. V8 recovers subject-level
endpoint-local estimates but retains unvalidated staging/anatomy and disabled 3B/3D inference;
availability is not a biological result.

**License summary:** Falach = CC BY 4.0 (commercial OK) · Zurich ds003498 = CC0 (public domain) ·
Pattnaik atlas = CC-BY-NC-SA-4.0 (non-commercial). Only one of the three is non-commercial.

---

## Part 2 — Conceptual literature (NOT YET VERIFIED)

> The references below are placeholders. Each needs its citation (year, journal, DOI) and the
> specific claim we rely on verified against the primary source before it is finalized here.

- Kipnis — (glymphatic / meningeal-lymphatic clearance; sleep-dependent coordination link)
- Hauglund —
- Uji / Tamaki —
- Dagum —
- Murdock —
- Miao —
- Staresina —
- Mander —
- Zhang —
- Mlinarič —
- Azadian —
