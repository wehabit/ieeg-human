# Locus-coeruleus infraslow coupling in human iEEG — **result: negative**

> **Standalone branch.** This branch (`lc-infraslow-3ABD`) holds one self-contained study and is
> **not intended to be merged**. `master` holds a different study — NREM hierarchical nesting
> (slow oscillation → spindle → ripple) in human MTL, after
> [Staresina 2015](https://doi.org/10.1038/nn.4119) — with its own README; some shared
> infrastructure in `analysis/` belongs to it.

The **locus coeruleus (LC)** is proposed to coordinate the cardinal events of NREM sleep — **slow
oscillations (SO), sleep spindles, and heart rate** — and it degenerates early in neurodegenerative
disease, so a measurable coupling among these would be a candidate biomarker. This study asks whether
that coordination is detectable in human intracranial EEG, through **three questions** linking the slow
oscillation, spindles, and heart rate (RR), each split by lighter (**N2**) vs deeper (**N3**) NREM.

**Answer:** the proposed LC infraslow spindle↔heart rhythm (**3A**) is **absent**, and SO→heart-rate
coupling (**3B**) is **weak and non-replicating**; only the cortical SO→spindle coupling (**3D**) is
**present**. No N2-vs-N3 difference in any test.

**📄 Methods: [docs/METHODS.md](docs/METHODS.md)** · **Full write-up:
[docs/LC_INFRASLOW_3ABD_SUMMARY.md](docs/LC_INFRASLOW_3ABD_SUMMARY.md)** — hypothesis, datasets, tests,
subject flow, results, limitations, references.

---

## The three tests

| | Question (each also tested N2 vs N3) | Result | Source paper |
|---|---|---|---|
| **3A** | Do **spindle (sigma) power and heart rate share a ~0.02 Hz (~50 s) infraslow rhythm** during NREM — the proposed LC signature? | **negative** (HUP n=23, K median 59). Step 1: sigma-power peaks **scatter like 1/f noise** (SD 0.014 Hz ≈ this cohort's own surrogates 0.014; a true shared rhythm would give ~0.001) and are no more prominent than the SWA control (p=0.11). Step 2: HR does **not** track them — cross-correlation null (median \|r\|=0.06, p=0.23), and 0.02 Hz ranks only **10th of 63** bins in the frequency-specificity control. **RESPect (n=6):** consistent — fitted peaks at ~0.04 Hz (0/5 in the infraslow band), coherence estimable in only 2/6 and non-significant; but underpowered (1 h nights → K median 7). | [Lecci 2017, *Sci Adv*](https://doi.org/10.1126/sciadv.1602026) · [Osorio-Forero 2021, *Curr Biol*](https://doi.org/10.1016/j.cub.2021.09.041) |
| **3B** | Does **heart rate rise at a fixed latency after the slow-oscillation down-state** (SO→HR coupling)? | **weak and heterogeneous** (n=23, stage-matched null). Significant in only 7/23 (N2) and 6/20 (N3); cohort z-test nonzero (t p=0.014 / 0.020). N3 HR peak **+3.5% ≈ Naji's +3.35%**; N2 **+1.8% is ~7× weaker** than Naji's +12.1%, and Naji's **N2 ≫ N3 pattern does not reproduce** (paired p=0.57). SO→HR lag ~1–2 s. **Null** in the RESPect replication. | [Naji & Mednick 2019, *J Cogn Neurosci*](https://doi.org/10.1162/jocn_a_01432) |
| **3D** | Are **sleep spindles phase-locked to the slow oscillation** (SO→spindle coupling)? | **present** — 21/23 subjects, median **83%** of channels significant, R = 0.073. No N2/N3 difference (p = 0.29). | [Staresina 2015, *Nat Neurosci*](https://doi.org/10.1038/nn.4119) · [Helfrich 2018, *Neuron*](https://doi.org/10.1016/j.neuron.2017.11.020) |

Why it mattered clinically: [Winer 2019, *J Neurosci*](https://doi.org/10.1523/JNEUROSCI.0503-19.2019)
(impaired SO–spindle coupling predicts medial-temporal tau) and
[Jacobs 2021, *Sci Transl Med*](https://doi.org/10.1126/scitranslmed.abj2511)
(LC integrity indexes Alzheimer's pathology and cognitive decline).

### 3A — no infraslow spindle↔heart rhythm

Two-step method after Lecci, pooling all NREM (coherence K median **59**) with a gap-aware estimator
([`analysis/lecci_faithful_3A.py`](analysis/lecci_faithful_3A.py)):

- **Step 1 — is there a ~50 s sigma rhythm?** Per subject, a duration-weighted Morlet spectrum of the
  sigma-power time course over all NREM bouts ≥120 s, with a Gaussian peak fit (Lecci Fig 1G). Cohort
  question: do the 23 peaks **cluster** at ~0.019 Hz (a real shared rhythm → simulated SD ≈ 0.001 Hz)
  or **scatter** like 1/f noise? Observed peak SD = **0.0137 Hz**, essentially identical to this
  cohort's own scale-free surrogates (0.0139 Hz); bootstrap finds the real peaks **not** tighter
  (p = 0.49). Sigma is **not** more prominent than the SWA control (Wilcoxon p = 0.11). Individual
  subjects have infraslow bumps, but at scattered frequencies — no shared ~50 s rhythm.
- **Step 2 — does heart rate track it?** Two controls agree it does not: 0.02 Hz ranks **10th of 63**
  bins in the frequency-specificity test (background 16.5%, p = 0.17), and the **cross-correlation**
  (Lecci's coupling statistic) is null — group \|r\| = 0.036, per-subject median 0.056, lags scattered
  (IQR −10 to −1 s), t vs 0 p = 0.23.

**Replication — RESPect (n=6, real EMG/EOG staging):** consistent with the HUP negative. Fitted infraslow
peaks sit at ~0.04 Hz (mean 0.040, **0/5 inside 0.015–0.025 Hz**), heart-rate coherence is estimable in
only 2/6 subjects and significant in none, and the cross-correlation is flat. Excluding REM/wake did
**not** rescue an infraslow rhythm — so the negative is not a staging artifact. But the 1 h recordings
are underpowered for 3A (coherence K median 7 vs 59), so this is a directional, not definitive,
replication. See [docs/DS003848_REPLICATION.md](docs/DS003848_REPLICATION.md).

### 3B — SO→heartbeat coupling is weak and heterogeneous

Method after Naji 2019 ([`analysis/event_3B_mednick.py`](analysis/event_3B_mednick.py) + `event_3B_cached.py`):

- **SO detection** — per channel, band-pass 0.15–4 Hz, negative half-waves with duration 0.3–1.0 s and
  amplitude & peak-to-peak ≥ 75th percentile within channel.
- **Heart rate** — R-peaks → RR → resampled to 4 Hz by piecewise cubic spline.
- **Statistic** — mean HR in a ±5 s window on the SO down-state trough; effect = peak of the post-trough
  curve as **% above that stage's mean HR**, plus the SO→HR latency; significance against a
  **stage-matched** random-trigger null (200 surrogates), z-scored.

**HUP (n=23):**

| | N2 | N3 | Naji 2019 (scalp) |
|---|---|---|---|
| HR peak (% above stage mean) | **+1.8%** | **+3.5%** | +12.09% (N2) / +3.35% (N3) |
| significant (z > 1.96) | 7/23 | 6/20 | — |
| cohort z-test (t vs 0) | p = 0.014 | p = 0.020 | — |
| SO→HR lag | ~1–2 s | ~1–2 s | follows down-state ✓ |

So the coupling is significant at the group level but carried by only ~30% of patients (median subject
near zero). **N3 matches Naji's +3.35%; N2 is ~7× weaker** than Naji's +12.09%, and Naji's **N2 ≫ N3
pattern does not reproduce** (paired p = 0.57).

**Replication — RESPect (n=6):** **null** — N2 +0.6% (0/5 significant, t p = 0.88), N3 +1.0% (0/6,
t p = 0.28), N2 ≫ N3 not reproduced.

**Why lateral iEEG likely misses it (follow-ups on RESPect):** the coupling lives in large slow waves,
and a lateral depth contact sees mostly local ones. Restricting to globally-synchronous SOs strengthened
the effect in N3 (global > local in 5/5 subjects) but stayed ~+2%; cortical region (frontal/cingulate/
insula vs posterior) did **not** matter; and isolating N2 K-complexes was too rare to test. Nothing
approached scalp magnitudes. See [docs/3B_METHOD_COMPARISON.md](docs/3B_METHOD_COMPARISON.md) and
[docs/3B_REGION_GLOBALITY.md](docs/3B_REGION_GLOBALITY.md).

### 3D — SO→spindle coupling is present (agrees with `master`)

Per channel: discrete SO and spindle events, Rayleigh test on the SO phase at each spindle peak
([`analysis/event_3D_by_stage.py`](analysis/event_3D_by_stage.py)), after Staresina 2015 / Helfrich 2018:

| | value |
|---|---|
| subjects with ≥1 significant channel | **21 / 23** |
| median fraction of channels significant | **83%** |
| median resultant vector length R | **0.073** |

R ≈ 0.07 matches `master`'s independent estimate (R ≈ 0.03–0.05, 27/28 channels) from a separate
analysis path. **N2 vs N3 null**: median R 0.071 vs 0.069 (n = 20, Wilcoxon p = 0.29).

**Subject flow:** 25 subjects with depth + EKG → 2 excluded (no usable cortical channels) → **23
analysed** → 17 passed the staging-quality filter. No N2-vs-N3 difference in any test.

## Datasets

**Two datasets, examined at three scopes** (n=1 worked subject → n=23 cohort → n=6 replication):

| Dataset | Rate | Coverage | n (scope) |
|---|---|---|---|
| **iEEG.org HUP `phaseII`** — continuous multi-day recordings, Penn epilepsy monitoring unit | 256–1024 Hz | lateral neocortical contacts (highest per shaft) + **`EKG1`/`EKG2`**; no electrode localisation | **1** — worked single subject (HUP165), method developed here |
| — same dataset, full cohort | 256–1024 Hz | as above | **25** with depth + EKG → **23** analysed |
| **OpenNeuro ds003848** — Utrecht **RESPect** long-term iEEG (independent replication) | 2048 Hz (50 Hz line) | 3 ECoG grid + 3 SEEG depth; ECG + **EMG + EOG** + respiration belts; **Destrieux atlas** labels + MNI coords | **6** |

*HUP165 is one subject of iEEG.org HUP `phaseII` (dataset name `HUP165_phaseII`), not a separate
dataset — each test was first built/visualised on it, then run on the full cohort.*

The enabling find: of 52 candidate datasets, 39 were reachable and **every one carries an EKG
channel** alongside the depth electrodes — public iEEG almost never includes cardiac signals. Full
survey of iEEG.org, OpenNeuro, DANDI, DABI and EBrains, with a test-by-test feasibility verdict, in
[docs/DATA_INVENTORY_LC_INFRASLOW.md](docs/DATA_INVENTORY_LC_INFRASLOW.md).

**On the replication (ds003848 / RESPect):** staged with real REM/wake exclusion
(`analysis/stage_ds003848.py`), then the same 3A/3B. Consistent with the HUP negative — excluding
REM/wake does **not** rescue an infraslow rhythm — but the 1 h recordings are underpowered for 3A
(coherence K median 7 vs 59). Directional, not definitive. See
[docs/DS003848_REPLICATION.md](docs/DS003848_REPLICATION.md).

## Scripts

| Script | What it does |
|---|---|
| `probe_physio_channels.py` | discovers which iEEG.org datasets carry EKG/ECG |
| `ekg_quality_check.py` | verifies the clinical EKG is usable for R-peak detection |
| `infraslow_rr_sigma_coherence.py` | the 3A pipeline on one subject (+ overlay/coherence figure) |
| `results_3A_tutorial_style.py` | single-subject 3A result drawn in the teaching-figure style |
| `build_3A_results_page.py` | builds the self-contained HTML results page |
| `cohort_3A_cortical.py` | cortical-channel selection, night finding, first cohort pass |
| **`cohort_stages_3ABD.py`** | 3A/3B/3D per sleep stage, one streaming pass per subject |
| `summarize_cohort_stages.py` | cohort aggregation with staging-quality filters + paired tests |
| `validate_3A_false_positive_rate.py` | calibrates the coherence test on independent signals |
| **`frequency_specificity_3A.py`** | the frequency-specificity control (exceedance at every bin) |
| `tutorial_signal_walkthrough.py` | synthetic teaching figures explaining 3A/3B/3D |
| **`spectral_gapped.py`** | **gap-aware coherence/PSD**; `test_spectral_gapped.py`, `test_coherence_calibration.py` |
| **`cache_lc_series.py`** | streams each night once → cached derived series (`data/derived/lc_infraslow/`) |
| **`lecci_faithful_3A.py`** | **3A** — Lecci per-subject peak fit + cross-correlation; `test_lecci_faithful.py` |
| **`event_3B_cached.py`** | **3B** — Naji method, stage-matched null, whole cohort; `test_3B_null.py` |
| `summarize_corrected_3AB.py` | cohort summary → `outputs/corrected_3AB/COHORT_SUMMARY.txt` |

*(Other scripts in `analysis/` — `slow_power_by_state`, `slow_ripple_coupling`, `hfo_slow_phase`,
`atlas.py`, … — belong to the separate study on `master`.)*

## Outputs

```
outputs/corrected_3AB/               cohort summary (3A peak fit + cross-correlation, 3B stage-matched)
outputs/lecci_faithful_3A/           per-subject 3A (peak fit, coherence, cross-correlation)
outputs/event_3B_cached/             per-subject 3B (stage-matched null)
outputs/cohort_stages_3ABD/          per-subject 3A/3B/3D by stage (JSON) + cohort CSV
outputs/freq_specificity_3A/         full coherence spectra (frequency-specificity control)
outputs/results_3A_tutorial_style/   single-subject result figures (png/svg/json)
outputs/ekg_quality_check/           EKG quality verification
outputs/signal_tutorial/             synthetic method-explainer figures
reports/HUP165_phaseII_3A_results.html   worked single-subject report (self-contained HTML)
docs/LC_INFRASLOW_3ABD_SUMMARY.md    the write-up
docs/DATA_INVENTORY_LC_INFRASLOW.md  what public data can and cannot support these tests
docs/DEC_3A_METHOD_NOTES.md          method notes / traps
docs/HUP165_N2_N3_RESULTS.md         full N2/N3 breakdown for the worked subject
```

## Method notes worth reading before reusing any of this

[docs/DEC_3A_METHOD_NOTES.md](docs/DEC_3A_METHOD_NOTES.md) documents methodological points that matter
for these signals:

- `ds.get_data()` returns channels in **ascending index order, not the order requested** — remap before
  assigning EKG vs brain channels
- **do not band-pass before computing coherence** — it destroys the estimate (a known-coupled control
  pair scores 0.09 instead of 0.98)
- **shift and phase-randomisation surrogates are invalid** for oscillatory coupling — a shifted 50 s
  rhythm is still coherent with itself; use the analytic threshold `1 − α^(1/(K−1))`
- the **band-maximum** coherence test has a **34%** false-positive rate; only a pre-specified
  frequency point is defensible
- report **effect size, not z** for phase-amplitude coupling (Tort MI_z is inflated by sample size)
- coherence is biased by window length (floor ≈ 1/K), so any contrast between conditions must use
  **length-matched** windows

## Limitations

Epilepsy patients on anti-seizure medication, many tachycardic in NREM · staging is **not** scored
AASM (no EOG/EMG in iEEG; N2/N3 approximated by a Gaussian mixture on slow-wave power, hence
"-like") · slow oscillations detected on the channel average, which likely suppressed 3B/3D and
should be redone per channel · no respiration channel, so apnoea cannot be excluded · electrode
locations inferred from contact numbering, not localisation files · stage-resolved 3A rests on only
4 subjects.

This is a **negative result in epilepsy patients**, not a refutation of Lecci/Osorio-Forero, who
worked in healthy sleepers with scalp EEG and direct LC recordings in mouse.

## Stack

Python + numpy/scipy · **NeuroKit2** (R-peak detection) · scikit-learn (staging) · matplotlib ·
pandas · `ieeg` (iEEG.org client). Raw data and credentials are gitignored; see
[docs/IEEG_ORG_SETUP.md](docs/IEEG_ORG_SETUP.md) for portal access.
