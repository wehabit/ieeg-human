# Methods — 3A, 3B, 3D

## Datasets

**Two datasets, three scopes:** each test was first built on a single worked subject (**HUP165**,
n=1), then run on the full **iEEG.org HUP phaseII** cohort (n=23), then replicated on the independent
**Utrecht RESPect** cohort (n=6). HUP165 is one subject *within* HUP phaseII (`HUP165_phaseII`), not a
separate dataset.

| | HUP phaseII (primary) | Utrecht RESPect (replication) |
|---|---|---|
| Source | iEEG.org (Penn epilepsy monitoring unit) | OpenNeuro **ds003848** (UMC Utrecht) |
| n analysed | 23 (of 25 with depth + EKG); developed on 1 worked subject (HUP165) | 6 |
| Recording | continuous multi-day; ~7 h/subject streamed from the highest-delta night | 1 h continuous `task-[Ss]leep` |
| Sampling | 256–1024 Hz (per subject) | 2048 Hz, 50 Hz line |
| iEEG electrodes | **SEEG depth**; lateral neocortical contacts (see below) | **3 ECoG grid + 3 SEEG depth** |
| Cardiac | `EKG1`/`EKG2` | `ECG` |
| EMG / EOG | none | **EMG + EOG** (+ thoracic/abdominal respiration belts) |
| Localisation | none — depth ordering only (contact 1 = deepest/mesial, high = lateral) | **MNI coordinates + Destrieux atlas labels** |

## Electrode selection

- **HUP:** the highest-numbered contact on each depth shaft (shafts with ≥6 contacts), up to 6 channels
  — the most **lateral neocortical** end of each SEEG trajectory (contact 1 is deepest/mesial). Chosen
  because spindles are thalamo*cortical* and Lecci's rhythm is parietal-maximal. No anatomical labels
  exist, so region beyond mesial-vs-lateral depth is not recoverable.
- **RESPect:** all good iEEG channels (54–93/subject), classified by role from `channels.tsv` and by
  Destrieux atlas region from `electrodes.tsv`.

## Common preprocessing

- Line notch (60 Hz HUP / 50 Hz RESPect + harmonics).
- **Fast-spindle peak (FSP):** per subject, whiten the PSD (remove the log-log 1/f trend), take the
  local peak in 10.5–16 Hz (fallback 13 Hz).
- **Sigma power:** Hilbert envelope of the band-passed signal (FSP ± 1 Hz, or fixed 10–15 Hz), squared,
  IED-masked, averaged across contacts, binned to **1 Hz**.
- **Heart rate:** R-peaks (NeuroKit2) → RR intervals, gated to physiological 0.33–1.5 s → instantaneous
  HR resampled to a 1 Hz grid (3A) or **4 Hz piecewise cubic spline** (3B).
- **Staging (30 s epochs):**
  - *HUP* (no EOG/EMG): NREM by a 2-component Gaussian mixture on delta ratio; N2/N3 split by a
    2-component GMM on log slow-wave (0.5–4 Hz) power. A proxy — labelled "N2-like/N3-like".
  - *RESPect* (has EOG/EMG): rule-based. Per epoch, robust-z of submental **EMG** RMS (10–100 Hz),
    **EOG** movement variance (0.3–6 Hz), and slow-wave power. **Wake** = EMG z > 1.0; **REM** = EMG
    z < −0.3 & SWA z < 0 & EOG z > 0.5; remaining = **NREM**, split N2/N3 by GMM on log SWA. Real
    REM/wake exclusion.

---

## 3A — Do spindle power and heart rate oscillate together every ~50 s?

**Source:** Lecci 2017 (*Sci Adv*); Osorio-Forero 2021 (*Curr Biol*). **Datasets:** HUP (n=23) + RESPect (n=6).

**Method (per subject, on pooled NREM):**
1. Two 1 Hz series — sigma power (FSP ± 1 Hz, IED-masked) and instantaneous heart rate.
2. **Step 1 — is there a ~50 s rhythm?** Morlet spectrum (3-cycle) of the sigma-power time course over
   **every NREM bout ≥ 120 s**, 0.004–0.12 Hz at 0.001 Hz resolution, duration-weighted across bouts,
   normalised, **Gaussian-fit** to locate that subject's own infraslow peak. Cohort test: do the peaks
   **cluster** near 0.019 Hz, versus each subject's scale-free (1/f) surrogates?
3. **Step 2 — does HR track it?** Magnitude-squared coherence (Welch, nperseg 256, 0.005 Hz drift
   high-pass, gap-aware) read at 0.02 Hz **and** at the subject's own fitted peak; significance via the
   analytic threshold `1 − α^(1/(K−1))`. Plus **cross-correlation** of z-scored 120 s intervals (HR as
   source wave, ±60 s lag).

| Step | Lecci 2017 | Our implementation |
|---|---|---|
| Signal | **parietal scalp** EEG, 4-s epochs, sigma 10–15 Hz | **lateral neocortical iEEG**, 1-s bins, FSP ± 1 Hz |
| Rhythm | Morlet spectrum + Gaussian peak fit per bout | same (Morlet + Gaussian fit) |
| Coupling | cross-correlation | cross-correlation **+ gap-aware coherence** |
| Staging | R&K scored | GMM proxy (HUP) / EMG-EOG (RESPect) |

---

## 3B — Does heart rate shift around the slow-oscillation trough?

**Source:** Naji 2019 (*J Cogn Neurosci*). **Datasets:** HUP (n=23) + RESPect (n=6).

**Method (per channel, per stage):**
- **SO detection:** zero-phase band-pass **0.15–4 Hz**; negative half-waves with **duration 0.3–1.0 s**
  and negative-peak amplitude **and** peak-to-peak ≥ **75th percentile** within channel.
- **HR:** R-peaks → RR → **4 Hz cubic spline**.
- **Statistic:** mean HR in a ±5 s window on the SO down-state trough; effect = **peak of the
  post-trough curve as % above that stage's mean HR**; plus SO→HR peak latency. Significance = a
  **stage-matched random-trigger null** (200 surrogates drawn from the same stage), z-scored.

| Step | Naji 2019 | Our implementation |
|---|---|---|
| EEG band | zero-phase 0.15–4 Hz | 0.15–4 Hz ✓ |
| SO detection | Dang-Vu 2008 criteria (p2p amp, up-state amp, down/up-state duration) on **F3/F4 scalp** | negative half-waves, dur 0.3–1.0 s, amp & p2p ≥ **75th pct**, on **lateral iEEG** (percentile substitutes for µV) |
| RR → HR | 4 Hz cubic spline | identical ✓ |
| Statistic | peak of mean HR curve, % above stage mean; SO→HR interval | identical ✓ |
| Significance | **none** (descriptive, mean ± SEM); real endpoint = ΔT ↔ behaviour (Pearson) | **stage-matched null + z**; no behaviour available |
| Region | frontal scalp (global SOs) | lateral neocortical iEEG (local SOs) |

**3B follow-ups (RESPect only, from Destrieux labels):** SO-globality gradient (consensus across
channels), region contrast (autonomic-adjacent vs posterior), and an N2 K-complex proxy (isolated
frontal N2 events). See [3B_REGION_GLOBALITY.md](3B_REGION_GLOBALITY.md).

---

## 3D — Does the slow oscillation organise spindles?

**Source:** Staresina 2015 (*Nat Neurosci*); Helfrich 2018 (*Neuron*). **Dataset:** HUP (n=23).

**Method (per channel, per stage — event-based):**
- **SO phase:** band-pass **0.5–1.25 Hz**, Hilbert phase.
- **SO events:** troughs by `find_peaks(−SO, height = SD, min spacing 0.8 s)`.
- **Spindle events:** Hilbert envelope of FSP ± 1 Hz, z-scored; peaks at **z ≥ 1.5, min spacing 0.3 s**.
- **Coupling:** take the SO phase at each spindle peak → **Rayleigh test** on that phase distribution →
  resultant vector length R (≥20 events per channel per stage).

| Step | Staresina/Helfrich | Our implementation |
|---|---|---|
| Approach | **event-locked** — detect SO and spindle events, SO phase at spindle peak, Rayleigh | same |
| SO band | 0.16–1.25 Hz (Helfrich) | 0.5–1.25 Hz |
| Spindle band | 12–16 Hz | FSP ± 1 Hz (individual) |

---

*Code: `analysis/lecci_faithful_3A.py` (3A), `analysis/event_3B_mednick.py` + `event_3B_cached.py`
(3B), `analysis/event_3D_by_stage.py` (3D); derived-series cache built by `cache_lc_series.py` (HUP)
and `stage_ds003848.py` (RESPect). Results in `outputs/`.*
