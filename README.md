# Locus-coeruleus infraslow coupling in human iEEG — **result: negative**

> **Standalone branch.** This branch (`lc-infraslow-3ABD`) holds one self-contained study and is
> **not intended to be merged**. `master` holds a different study (Kipnis NREM slow-oscillation
> nesting in human MTL) with its own README; some shared infrastructure in `analysis/` belongs to it.

Does the **locus coeruleus** leave a shared **~0.02 Hz (~50 s) rhythm** in sleep spindles and heart
rate during NREM ([Lecci 2017](https://doi.org/10.1126/sciadv.1602026); [Osorio-Forero 2021](https://doi.org/10.1016/j.cub.2021.09.041)), and does that coupling differ between lighter and
deeper NREM? The LC degenerates early in neurodegenerative disease, so weakened coupling was a
candidate biomarker.

**Answer: no.** No infraslow spindle–heart coupling, and no N2/N3 difference, in 23 subjects.

**📄 Full write-up: [docs/LC_INFRASLOW_3ABD_SUMMARY.md](docs/LC_INFRASLOW_3ABD_SUMMARY.md)** —
hypothesis, dataset, tests, subject flow, results, limitations, references.

---

## The three tests

| | Question | Result | Source paper |
|---|---|---|---|
| **3A** | Do spindle power and heart rate rise and fall together every ~50 s? | **negative** — retracted after a frequency-specificity control | [Lecci 2017, *Sci Adv*](https://doi.org/10.1126/sciadv.1602026) · [Osorio-Forero 2021, *Curr Biol*](https://doi.org/10.1016/j.cub.2021.09.041) |
| **3B** | Does heart rate shift around the slow-oscillation trough? | **null** (p = 0.11) | [Chen & Mednick 2022, *PNAS*](https://doi.org/10.1073/pnas.2123417119) |
| **3D** | Does the slow oscillation organise spindles? | **null** (p = 0.40), raw effect ~6×10⁻⁵ | [Staresina 2015, *Nat Neurosci*](https://doi.org/10.1038/nn.4119) · [Helfrich 2018, *Neuron*](https://doi.org/10.1016/j.neuron.2017.11.020) |

Why it mattered clinically: [Winer 2019, *J Neurosci*](https://doi.org/10.1523/JNEUROSCI.0503-19.2019)
(impaired SO–spindle coupling predicts medial-temporal tau) and
[Jacobs 2021, *Sci Transl Med*](https://doi.org/10.1126/scitranslmed.abj2511)
(LC integrity indexes Alzheimer's pathology and cognitive decline).

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

## Honest limitations

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
