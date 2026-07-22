# Locus-coeruleus infraslow coupling in human iEEG — **result: negative**

> **Standalone branch.** This branch (`lc-infraslow-3ABD`) holds one self-contained study and is
> **not intended to be merged**. `master` holds a different study — NREM hierarchical nesting
> (slow oscillation → spindle → ripple) in human MTL, after
> [Staresina 2015](https://doi.org/10.1038/nn.4119) — with its own README; some shared
> infrastructure in `analysis/` belongs to it.

Does the **locus coeruleus** leave a shared **~0.02 Hz (~50 s) rhythm** in sleep spindles and heart
rate during NREM ([Lecci 2017](https://doi.org/10.1126/sciadv.1602026); [Osorio-Forero 2021](https://doi.org/10.1016/j.cub.2021.09.041)), and does that coupling differ between lighter and
deeper NREM? The LC degenerates early in neurodegenerative disease, so weakened coupling was a
candidate biomarker.

**Answer: no** for the infraslow (~50 s) rhythm — and no N2/N3 difference in any test.
But two of the three tests were **corrected after the fact**: 3D reverses from null to clearly
**positive** once measured the way the source papers do, and 3B's cited source turned out to be a
review rather than a methods paper. Both corrections are documented below rather than quietly folded in.

**📄 Full write-up: [docs/LC_INFRASLOW_3ABD_SUMMARY.md](docs/LC_INFRASLOW_3ABD_SUMMARY.md)** —
hypothesis, dataset, tests, subject flow, results, limitations, references.

---

## The three tests

| | Question | Result | Source paper |
|---|---|---|---|
| **3A** | Do spindle power and heart rate rise and fall together every ~50 s? | **negative** — retracted after a frequency-specificity control | [Lecci 2017, *Sci Adv*](https://doi.org/10.1126/sciadv.1602026) · [Osorio-Forero 2021, *Curr Biol*](https://doi.org/10.1016/j.cub.2021.09.041) |
| **3B** | Does heart rate shift around the slow-oscillation trough? | ⚠️ **present but ~25× weaker than published** — HR peak +0.5% above stage mean vs Naji's +12.1% (Stage 2); SO→HR lag 2.2 s | [Naji & Mednick 2019, *J Cogn Neurosci*](https://doi.org/10.1162/jocn_a_01432) — *the actual methods paper; the 2022 PNAS reference is a review with no such analysis* |
| **3D** | Does the slow oscillation organise spindles? | ✅ **PRESENT** — 21/23 subjects, median **83%** of channels significant, R = 0.073. *(Corrected: an earlier continuous-MI version reported this as null; that was a method artifact.)* No N2/N3 difference (p = 0.29). | [Staresina 2015, *Nat Neurosci*](https://doi.org/10.1038/nn.4119) · [Helfrich 2018, *Neuron*](https://doi.org/10.1016/j.neuron.2017.11.020) |

Why it mattered clinically: [Winer 2019, *J Neurosci*](https://doi.org/10.1523/JNEUROSCI.0503-19.2019)
(impaired SO–spindle coupling predicts medial-temporal tau) and
[Jacobs 2021, *Sci Transl Med*](https://doi.org/10.1126/scitranslmed.abj2511)
(LC integrity indexes Alzheimer's pathology and cognitive decline).

### 3D, corrected: SO→spindle coupling **is** present — and agrees with `master`

An earlier version of 3D reported this as null (raw modulation index ~6×10⁻⁵). **That was wrong, and
the cause was my implementation, not the data.** It binned *every* timepoint by slow-oscillation
phase into a continuous Tort modulation index, with **no event detection at all**, on a
**channel-averaged** signal. Most of a night is neither spindle nor slow oscillation, so that
averages the real events together with hours of nothing and drives the estimate toward zero.

Neither [Staresina 2015](https://doi.org/10.1038/nn.4119) nor
[Helfrich 2018](https://doi.org/10.1016/j.neuron.2017.11.020) does it that way. Helfrich, verbatim:
*"we detected SO (0.16–1.25 Hz) and sleep spindle (12–16 Hz) **events** … phase during the **peak of
the detected sleep spindle events** … Rayleigh z"*. Re-implemented that way — per channel, discrete
SO and spindle events, Rayleigh test on the SO phase at each spindle peak
(`analysis/event_3D_by_stage.py`) — the result reverses:

| | event-based (correct) | old continuous MI |
|---|---|---|
| subjects with ≥1 significant channel | **21 / 23** | — |
| median fraction of channels significant | **83%** | — |
| median resultant vector length R | **0.073** | raw MI 6×10⁻⁵ ("negligible") |
| verdict | **coupling present** | "null" |

R ≈ 0.07 lands on `master`'s independent estimate (R ≈ 0.03–0.05, 27/28 channels) from a completely
separate analysis path — so the two studies **agree**, and the apparent conflict was entirely an
artifact of the discarded method.

**N2-like vs N3-like remains null**: median R 0.071 vs 0.069 (n = 20, Wilcoxon p = 0.29); fraction of
channels significant 0.83 vs 0.58 (p = 0.16). Direction favours N2 but does not reach significance.

**Subject flow:** 25 subjects with depth electrodes + EKG → 2 excluded (no usable cortical channels)
→ **23 analysed** → 17 passed the staging-quality filter. **No N2-like vs N3-like difference in any
test**, including when split by fast vs slow individual spindle peak.

**Why 3A was retracted.** 7/23 subjects cleared threshold at 0.02 Hz, which looked like a strong
positive. But counting exceedances at *every* frequency gave a background of **11.6%** — against the
**6%** that a simulation on independent signals predicts — with 0.02 Hz only **2nd of 63 bins**
(0.05 Hz scored higher). Spindle power and heart rate do share genuine broadband low-frequency
structure, but there is no ~50 s peak.

## Dataset

| Dataset | Rate | Coverage | n |
|---|---|---|---|
| **iEEG.org HUP `phaseII`** — continuous multi-day recordings, Penn epilepsy monitoring unit | 256–1024 Hz | lateral neocortical contacts + **`EKG1`/`EKG2`** | **25** with depth + EKG |

The enabling find: of 52 candidate datasets, 39 were reachable and **every one carries an EKG
channel** alongside the depth electrodes — public iEEG almost never includes cardiac signals. Full
survey of iEEG.org, OpenNeuro, DANDI, DABI and EBrains, with a test-by-test feasibility verdict, in
[docs/DATA_INVENTORY_LC_INFRASLOW.md](docs/DATA_INVENTORY_LC_INFRASLOW.md).

**Replication target (not run):** OpenNeuro **ds003848** — the one public iEEG sleep dataset with a
verified ECG channel, plus **EOG/EMG** so real sleep staging is possible.

## Scripts

| Script | What it does |
|---|---|
| `probe_physio_channels.py` | discovers which iEEG.org datasets carry EKG/ECG |
| `ekg_quality_check.py` | verifies the clinical EKG is usable for R-peak detection |
| `infraslow_rr_sigma_coherence.py` | the 3A pipeline on one subject (+ overlay/coherence figure) |
| `results_3A_tutorial_style.py` | single-subject 3A result drawn in the teaching-figure style |
| `build_3A_results_page.py` | builds the self-contained HTML results page |
| `cohort_3A_cortical.py` | cortical-channel selection, night finding, first cohort pass |
| **`cohort_stages_3ABD.py`** | **main analysis** — 3A/3B/3D per sleep stage, one streaming pass per subject |
| `summarize_cohort_stages.py` | cohort aggregation with staging-quality filters + paired tests |
| `validate_3A_false_positive_rate.py` | calibrates the coherence test on independent signals |
| **`frequency_specificity_3A.py`** | **the control that retracted 3A** |
| `tutorial_signal_walkthrough.py` | synthetic teaching figures explaining 3A/3B/3D |

*(Other scripts in `analysis/` — `slow_power_by_state`, `slow_ripple_coupling`, `hfo_slow_phase`,
`atlas.py`, … — belong to the separate study on `master`.)*

## Outputs

```
outputs/cohort_stages_3ABD/          per-subject 3A/3B/3D by stage (JSON) + cohort CSV
outputs/freq_specificity_3A/         full coherence spectra — the retraction evidence
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

[docs/DEC_3A_METHOD_NOTES.md](docs/DEC_3A_METHOD_NOTES.md) documents pitfalls that each
independently invalidated an earlier version of these results:

- `ds.get_data()` returns channels in **ascending index order, not the order requested** — so the
  heartbeat detector was briefly running on a brain channel
- **band-passing before computing coherence** destroys the estimate: a known-coupled control pair
  scored 0.09 instead of 0.98
- **shift and phase-randomisation surrogates are invalid** for oscillatory coupling — a shifted 50 s
  rhythm is still coherent with itself; use the analytic threshold `1 − α^(1/(K−1))`
- the **band-maximum** coherence test has a **34%** false-positive rate; only a pre-specified
  frequency point is defensible
- report **effect size, not z** — Tort MI_z reached 212 where the raw MI was ~6×10⁻⁵ (negligible)
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
