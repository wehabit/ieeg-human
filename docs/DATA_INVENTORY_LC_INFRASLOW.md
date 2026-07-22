# Data inventory — which datasets can run the LC-infraslow coupling tests

Deliverable #1 of the LC-infraslow brief: for every candidate source, does it carry the signals the
five tests need, and which test can it actually drive? Surveyed 2026-07-21 by live query
(iEEG.org authenticated scan; OpenNeuro GraphQL; DANDI API over all 876 dandisets; DABI project
listing; EBrains + openlists index). **Rule: a signal is "yes" only if verified from real metadata
(channel tables / participants / API), never inferred from a paper's topic.**

The five tests and their data requirements:
- **3A** RR ↔ spindle-envelope infraslow (~0.02 Hz) coherence — needs iEEG **sleep + ECG**
- **3B** SO ↔ heartbeat event coupling — needs iEEG **sleep + ECG**
- **3C** pupil vs spindle/HR — needs **sleep pupillometry + iEEG**
- **3D** SO–spindle PAC (QC/context) — needs iEEG **sleep** only
- **3E** neurodegenerative vs control comparison — needs **disease labels**

## Bottom line

| Test | Runnable on public data? | Where |
|---|---|---|
| 3A (primary) | ✅ **Yes — two independent cohorts** | iEEG.org HUP `phaseII` (n≥25) **+** OpenNeuro **ds003848** (n=6) |
| 3B | ✅ Yes — same two cohorts | HUP `phaseII` + ds003848 |
| 3D | ✅ Yes — abundant | HUP nights + atlas + Falach + Zurich (ds003498) + several pediatric NREM sets |
| 3C pupil | ❌ **No** | pupil+iEEG exists but is **awake task only** (EBrains Kucewicz; DANDI 000623/001613) — no sleep pupillometry+iEEG anywhere |
| 3E disease | ❌ **No (usable)** | every iEEG **sleep** dataset is epilepsy-only; DABI has Parkinson's invasive data but it is **subcortical STN-LFP, not cortical iEEG, not co-recorded with sleep spindles** — cannot be joined |

**So the core coupling story (3A/3B/3D) is runnable now, and even replicable across two independent
cohorts. The disease test (3E) and pupil test (3C) have no public home — they require the user's own
labeled/instrumented cohort.**

## Primary cohort — iEEG.org HUP `phaseII` (Penn epilepsy monitoring)

Live scan reached 39 of 52 candidate `HUP{n}_phaseII` datasets; **all 39 carry an EKG/ECG channel**
(labels `EKG1/EKG2`, or `ECG1/ECG2` on HUP160/212; older subjects use `EEG EKG 01-Ref` naming).
**25 confirmed with both mesial-temporal depth AND EKG** — the usable set for 3A/3B (firm lower
bound: HUP142 and older grid subjects have depths under non-anatomical naming that the regex missed,
so true n may be higher). Multi-day continuous → contain natural NREM, long enough to resolve the
~50 s rhythm across many cycles.

Usable subjects (all epilepsy-only; no hypnogram → stage blind from delta, as `ieeg_pull_night.py` does; no respiration/SpO₂ → apnea cannot be regressed out):

| Subject | sfreq | length | Subject | sfreq | length |
|---|---|---|---|---|---|
| HUP116 | 500 | 96 h | HUP178 | 512 | 186 h |
| HUP130 | 1024 | 126 h | HUP182 | 512 | 260 h |
| HUP133 | 512 | 187 h | HUP185 | 512 | 233 h |
| HUP138 | 1024 | 172 h | HUP187 | 512 | 177 h |
| HUP139 | 1024 | 259 h | HUP191 | 512 | 216 h |
| HUP141 | 512 | 146 h | HUP199 | 512 | 148 h |
| HUP143 | 1024 | 239 h | HUP205 | 512 | 198 h |
| HUP150 | 512 | 113 h | HUP211 | 512 | 210 h |
| HUP151 | 512 | 187 h | HUP212 | 1024 | 204 h |
| HUP157 | 1024 | 191 h | HUP160 | 1024 | 336 h |
| HUP165 | 1024 | 459 h | HUP171 | 512 | 286 h |
| HUP172 | 512 | 240 h | HUP173 | 256 | 237 h |
| HUP177 | 512 | 192 h | | | |

Already pulled (MTL only, no EKG yet): HUP165, HUP157, HUP130 → re-pull with EKG added.

## Independent replication cohort — OpenNeuro ds003848

The **single OpenNeuro iEEG sleep dataset with a verified ECG channel**. Utrecht RESPect long-term
iEEG; n=6 long-term patients (a **mix of ECoG grid and SEEG depth** — e.g. RESP0521 = 64 ECoG,
RESP0749 = 67 SEEG, RESP0800 = 110 SEEG), 1 h continuous `task-Sleep` runs @ **2048 Hz**, **50 Hz**
line. Short per subject but continuous — enough cycles for infraslow coherence — and, crucially,
**an independent cohort/site**, mirroring how this repo already validated SO–spindle coupling across
three cohorts. Value: cross-cohort replication of 3A/3B, and — unlike HUP — **EMG + EOG enable real
REM/wake exclusion and staging** rather than the GMM-on-slow-wave proxy.

**Verified against the raw `channels.tsv` (2026-07), correcting two earlier claims here:**
- **Respiration belts DO exist** — `thor+` (thoracic) and `abdo+` (abdominal) are present (marked
  `status: bad`, units `Adim.`), so apnoea screening is possible in principle. The earlier "no
  respiration belt" was wrong.
- **Trust `channels.tsv`, not the `ieeg.json` sidecar counts.** For RESP0521 the JSON reads
  `EOGChannelCount: 0`, but `channels.tsv` lists `orb+` as an EOG channel. EOG presence must be checked
  per subject from `channels.tsv` — it looks present on the ECoG subject and absent on the SEEG
  subjects' JSON, so read all six before relying on it. ECG (`ECG+`) and EMG (`emg1+`/`emg2+`) are
  present across the subjects checked.

## The negatives (surveyed, not assumed)

- **OpenNeuro** (72 iEEG datasets): only ds003848 has iEEG-sleep+ECG. Zero have neurodegenerative
  labels — all epilepsy etiologies (HS, FCD, tumor, cavernoma, Sturge-Weber). Many good 3D sets.
- **DANDI** (876 dandisets; 58 human, ~20 human-iEEG): none combine sleep+ECG. AJILE12 (000055) =
  55 days ECoG spanning sleep but no ECG, no staging → weak 3D only. Pupil+iEEG exists (000623,
  001613) but awake movie/rest only. No neurodegenerative iEEG.
- **DABI** (222 projects, mostly 401-access-gated): no verifiable iEEG-sleep+ECG. Fried/UCLA sleep
  depth set (U01NS108930) has **0 files deposited**; Zelmann/MGH SEEG "natural sleep" is stim-evoked;
  PD/Mayberg sleep sets are **subcortical LFP**, not cortical iEEG. Holds Parkinson's invasive
  cohorts but not in a form joinable to the sleep-spindle tests.
- **EBrains / openlists index**: no iEEG-sleep+ECG. ds003498 `ieeg.json` shows `ECGChannelCount=0`;
  Pennsieve-414/MNI atlas is bipolar-iEEG-only (ECG=false). iEEG+pupil (EBrains GKNT-T3X, Kucewicz
  n=10) is a **waking memory task**. Neurodegenerative invasive data (figshare 31254136) = STN-DBS
  LFP + MEG, not iEEG, not sleep.

## Implication for scope

Run **3A/3B on HUP `phaseII` (n up to 25) with ds003848 as an independent replication**, and **3D**
as QC. Write **3C and 3E as feasibility/"what's needed"** — they cannot be run on any public data and
require the user's own instrumented (pupil) or labeled (neurodegenerative) cohort.
