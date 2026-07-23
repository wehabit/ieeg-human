# ds003848 replication — LC-infraslow tests in a properly-staged iEEG sleep cohort

**Bottom line: consistent with the HUP negative, and it removes the biggest caveat — but this cohort
is underpowered for 3A, so it is a directional replication, not a definitive one.**

The HUP `phaseII` result rests on a staging *proxy* (a Gaussian mixture on slow-wave power, because
iEEG has no EOG/EMG). The obvious objection: maybe the ~50 s rhythm is real but smeared by wake/REM
epochs that the proxy could not exclude. OpenNeuro **ds003848** (Utrecht RESPect long-term iEEG) is
the one public iEEG sleep dataset that carries **ECG + EMG + EOG**, so NREM can be scored with real
REM/wake exclusion. This is the test of that objection.

## Cohort and staging

Six patients, ~1 h continuous `task-[Ss]leep` @ 2048 Hz, 50 Hz line; **3 ECoG grid + 3 SEEG depth**.
Every subject verified (from raw `channels.tsv`) to carry iEEG, ECG, EMG, EOG and (bad) respiration
belts. Pipeline: `analysis/stage_ds003848.py` (MNE BrainVision reader → derived series identical to
the HUP cache → EMG/EOG staging) then the **same** corrected code
(`lecci_faithful_3A.py`, `event_3B_cached.py`) via a real-stage adapter.

**Staging is rule-based and physiologically validated.** True AASM scoring needs scalp EEG, which
this dataset lacks; what EMG + EOG buy is REM/Wake exclusion. Per 30 s epoch: submental EMG RMS,
EOG movement variance, iEEG slow-wave power, robust-z within subject. Wake = high EMG; REM = atonia +
phasic eye movement + low SWA; NREM split N2/N3 by GMM on SWA. On RESP0521 the labels separate exactly
as the physiology predicts:

| stage | n | EMG z | EOG z | SWA z |
|---|---|---|---|---|
| Wake | 15 | **+1.28** | +0.31 | +0.25 |
| REM | 3 | −0.55 | **+1.37** | −0.54 |
| N2 | 32 | −0.16 | −0.18 | −0.73 |
| N3 | 48 | −0.02 | +0.08 | **+0.51** |

Per-subject stage counts (W / R / N2 / N3, of ~120–128 epochs):
0521 15/3/32/48 · 0699 35/17/25/25 · 0724 19/1/6/75 · 0749 12/3/22/22 · 0779 25/5/39/41 · 0800 19/1/46/42.
Wake+REM is **6–41%** of epochs — real contamination that the HUP proxy could not remove, and here is
removed. Heart rates are also healthier than HUP (e.g. 70 bpm vs HUP's ~92 bpm tachycardia).

## 3A — underpowered here, but what signal exists does not support Lecci

The catch is duration. HUP streams ~7 h/subject; ds003848 is **1 h**, and real staging fragments NREM
into short runs, leaving very few bouts ≥120 s and a low Welch segment count:

| subject | bouts ≥120 s | coherence K | own peak (Hz) |
|---|---|---|---|
| RESP0521 | 5 | — | 0.035 |
| RESP0699 | 3 | — | 0.041 |
| RESP0724 | 4 | — | 0.052 |
| RESP0749 | **0** | — | — |
| RESP0779 | 3 | (low) | 0.027 |
| RESP0800 | 7 | (low) | 0.044 |

- **Step 1 (peaks).** 5/6 subjects yield a fitted peak; they sit at **0.027–0.052 Hz (mean 0.040)**,
  **0/5 inside 0.015–0.025 Hz**, and significantly *above* Lecci's 0.019 Hz (t p = 0.008). The peaks
  are somewhat tighter than this cohort's scale-free surrogates (SD 0.0093 vs 0.0137) but the bootstrap
  is not significant (p = 0.21), and they cluster near **0.04 Hz, not 0.019** — so there is no Lecci
  fingerprint. Sigma is not reliably more prominent than the SWA control (p = 0.31).
- **Step 2 (HR coupling).** Only **2/6** subjects have enough contiguous NREM for a coherence estimate
  at all (K median **7**, versus 59 in HUP) — 0/2 significant. The cross-correlation is flat (group
  |r| = 0.066; per-subject signed peak r vs 0: p = 0.09, and the trend is *negative* if anything).

So 3A here is **too underpowered to stand alone**, but every direction it points agrees with HUP: no
infraslow peak at 0.019 Hz, no heart-rate coupling. Crucially, **excluding REM and wake did not rescue
a rhythm** — the HUP negative is not merely a staging artifact.

## 3B — null (better powered than 3A, event-based)

3B needs no contiguity, so it uses thousands of SOs per stage and is the more trustworthy test on this
cohort. It is **null**:

| | ds003848 (n=6) | HUP (n=23) | Naji 2019 |
|---|---|---|---|
| N2 HR peak | +0.65% (0/5 sig; t p=0.88) | +1.8% (7/23 sig; t p=0.014) | +12.09% |
| N3 HR peak | +0.99% (0/6 sig; t p=0.28) | +3.5% (6/20 sig; t p=0.020) | +3.35% |
| N2 ≫ N3? | no (p=0.63) | no (p=0.57) | 3.6× |

No subject reaches z > 1.96 in either stage, and the cohort z-tests are null. Note the contrast with
HUP, whose 3B cohort z-test *was* significant (driven by ~7/23 subjects): with only 6 subjects here
this could be low power, or it could mean part of HUP's weak 3B signal came from wake/arousal
contamination that real staging removes. **n = 6 cannot distinguish these**, so this is flagged, not
concluded.

## What this replication does and does not establish

- **Does:** the HUP 3A negative survives genuine REM/wake exclusion — the strongest single objection to
  it (staging contamination) does not rescue an infraslow LC rhythm. Direction of every test agrees
  with HUP. Staging from EMG/EOG works and is validated.
- **Does not:** provide a *powered* independent 3A test. 1 h/subject, fragmented by real staging, gives
  median K = 7 and only 2/6 usable coherence estimates. A definitive replication needs longer
  continuous sleep with cardiac + EOG/EMG, which no current public iEEG dataset provides.

## Honest limitations

- **1 h recordings** → few NREM bouts ≥120 s → 3A badly underpowered (RESP0749 unusable for 3A).
- **Rule-based staging without scalp EEG** — validated to separate the stages, but not AASM-scored;
  the N2/N3 split is still SWA-driven, so it inherits that axis.
- **n = 6**, mixed ECoG/SEEG montage.
- Respiration belts exist but are marked `bad`, so apnoea screening was not attempted.

## Reproduce

```
.venv/bin/python analysis/stage_ds003848.py              # download → stage → cache (deletes raw .eeg)
.venv/bin/python analysis/run_ds003848_replication.py    # corrected 3A + 3B on real-staged NREM
.venv/bin/python analysis/summarize_corrected_3AB.py \
    --a-dir outputs/ds003848_3A --b-dir outputs/ds003848_3B --label "ds003848 replication"
```
Full numbers: `outputs/corrected_3AB/DS003848_SUMMARY.txt`. Cache: `data/derived/ds003848/` (gitignored).
