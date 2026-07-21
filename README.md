# iEEG — human intracranial tests of sleep-coordination mechanisms

This repository holds **two independent studies** on open human intracranial EEG, both requiring
**no new patient access or committee approval**.

| | Study | Result |
|---|---|---|
| **A** | **Locus-coeruleus infraslow coupling** — do spindles and heart rate share a ~50 s rhythm? *(this branch)* | **negative** |
| **B** | **Kipnis NREM nesting** — does the slow oscillation organise spindles and ripples in human MTL? *(`master`)* | positive |

---

# Study A — the locus-coeruleus infraslow test  ·  **result: negative**

Does the **locus coeruleus** leave a shared **~0.02 Hz (~50 s) rhythm** in sleep spindles and heart
rate during NREM (Lecci 2017; Osorio-Forero 2021), and does that coupling differ between lighter and
deeper NREM? Motivation: the LC degenerates early in neurodegenerative disease, so weakened coupling
was a candidate biomarker.

**📄 Full write-up: [docs/LC_INFRASLOW_3ABD_SUMMARY.md](docs/LC_INFRASLOW_3ABD_SUMMARY.md)** —
hypothesis, dataset, tests, subject flow, results, limitations, references.

## The three tests

| | Question | Result | Source |
|---|---|---|---|
| **3A** | Do spindle power and heart rate rise and fall together every ~50 s? | **negative** — retracted after a frequency-specificity control | Lecci 2017; Osorio-Forero 2021 |
| **3B** | Does heart rate shift around the slow-oscillation trough? | **null** (p = 0.11) | Chen & Mednick 2022 |
| **3D** | Does the slow oscillation organise spindles? | **null** (p = 0.40), raw effect ~6×10⁻⁵ | Staresina 2015; Helfrich 2018 |

**Subject flow:** 25 subjects with depth electrodes + EKG → 2 excluded (no usable cortical channels)
→ **23 analysed** → 17 passed the staging-quality filter. **No N2-like vs N3-like difference in any
test.**

**Why 3A was retracted.** 7/23 subjects cleared threshold at 0.02 Hz — but counting exceedances at
*every* frequency gave a background of **11.6%** (vs the **6%** a simulation on independent signals
predicts), with 0.02 Hz only **2nd of 63 bins**. Spindle power and heart rate do share genuine
broadband low-frequency structure, but there is no ~50 s peak.

## Dataset

| Dataset | Rate | Coverage | n |
|---|---|---|---|
| **iEEG.org HUP `phaseII`** — continuous multi-day recordings, Penn epilepsy monitoring unit | 256–1024 Hz | lateral neocortical contacts + **`EKG1`/`EKG2`** | **25** with depth + EKG |

The enabling find: of 52 candidate datasets, 39 were reachable and **every one carries an EKG
channel** alongside the depth electrodes — public iEEG almost never includes cardiac signals. Full
survey of iEEG.org, OpenNeuro, DANDI, DABI and EBrains in
[docs/DATA_INVENTORY_LC_INFRASLOW.md](docs/DATA_INVENTORY_LC_INFRASLOW.md).

**Replication target (not yet run):** OpenNeuro **ds003848** — the one public iEEG sleep dataset with
a verified ECG channel, plus **EOG/EMG** so real sleep staging is possible.

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

## Outputs

```
outputs/cohort_stages_3ABD/      per-subject 3A/3B/3D by stage (JSON) + cohort CSV
outputs/freq_specificity_3A/     full coherence spectra -- the retraction evidence
outputs/results_3A_tutorial_style/   single-subject result figures (png/svg/json)
outputs/ekg_quality_check/       EKG quality verification
outputs/signal_tutorial/         synthetic method-explainer figures
reports/HUP165_phaseII_3A_results.html   worked single-subject report (self-contained)
```

## Method notes worth reading before reusing this

[docs/DEC_3A_METHOD_NOTES.md](docs/DEC_3A_METHOD_NOTES.md) documents traps that each independently
invalidated earlier results:

- `ds.get_data()` returns channels in **ascending index order, not the order requested**
- **band-passing before computing coherence** destroys the estimate (a known-coupled pair scored 0.09 instead of 0.98)
- **shift/phase surrogates are invalid** for oscillatory coupling — use the analytic threshold
- the **band-maximum** coherence test has a **34%** false-positive rate; only a pre-specified frequency point is defensible
- report **effect size, not z** — Tort MI_z reached 212 where the raw MI was ~6×10⁻⁵

---

# Study B — Kipnis NREM nesting  ·  result: positive

*(the repository's original study; primary line of work on `master`)*

Does the NREM slow oscillation organize faster activity (spindles, ripples) in **human medial
temporal lobe** — the candidate neural substrate that Kipnis / Jiang-Xie link to sleep-dependent CSF
clearance? (Clearance itself is not measured.) Sister project to the Buzsáki mouse pilot
(`hpaticmousebuzsakilab`).

**Start here:** [docs/iEEG_EVIDENCE_BRIEF.md](docs/iEEG_EVIDENCE_BRIEF.md) (2-page spine) ·
[docs/FIGURE_GUIDE.md](docs/FIGURE_GUIDE.md) · [docs/PLAN.md](docs/PLAN.md)

## Datasets

| # | Dataset | Rate | Coverage |
|---|---|---|---|
| 1 | **Normative iEEG Sleep & Wake Atlas** — Pattnaik & Litt 2024, Pennsieve [10.26275/xhte-d11l](https://doi.org/10.26275/xhte-d11l) (30 s clips, W/N2/N3/R) | 204 Hz | mesiotemporal + neocortex |
| 2 | **iEEG sleep + IED/HFO annotations** — Falach/Geva-Sagiv/Eliashiv 2024, [Figshare 26131978](https://figshare.com/articles/dataset/26131978) (~3 min NREM/subject) | 1 kHz | MTL SEEG |
| 3 | **Zurich HFO sleep atlas** — Fedele/Sarnthein, OpenNeuro [ds003498](https://openneuro.org/datasets/ds003498) | 2 kHz | MTL SEEG |

Fetchers: dataset 1 via `analysis/atlas.py` (Pennsieve API); datasets 2–3 via public download.
iEEG.org access setup in [docs/IEEG_ORG_SETUP.md](docs/IEEG_ORG_SETUP.md).

## Analysis modules (one script → one result → one figure)

| Script | What it shows | Dataset |
|---|---|---|
| `slow_power_by_state` | slow-band (0.5–4 Hz) power rises wake→N3 (p = 2.9×10⁻⁶) | 1 |
| `slow_spindle_coupling` | slow-oscillation → spindle coupling, NREM-specific (p = 0.030) | 1 |
| `mtl_region_profile` | slow-power & coupling by MTL region (descriptive; pairwise n.s.) | 1 |
| `slow_ripple_coupling` | slow-oscillation → ripple (p = 0.0015) **and** → spindle (p = 6×10⁻⁵) | 2 |
| `nesting_timecourse` | SO-triggered averages: spindle/ripple ride the SO up-state | 2 |
| `hfo_slow_phase` | SO-phase of expert-marked ripples/fast-ripples | 3 |
| `mouse_human_bridge` | the mouse-theta vs human-slow bridge figure | — |
| `make_result_visuals` | explanatory raw-signal / spectrum / slopegraph figures | 1 |
| `hup165_night_coupling`, `hup_night_cohort`, `fullnight_strength`, `staresina_style` | full-night continuous analyses (iEEG.org) | 4 |

---

## Layout

```
analysis/   Study A: cohort_stages_3ABD, frequency_specificity_3A, infraslow_rr_sigma_coherence, ...
            Study B: slow_power_by_state, slow_spindle_coupling, slow_ripple_coupling, ...
            fetchers: atlas.py (Pennsieve), ieeg_portal.py + ieeg_pull_night.py (iEEG.org)
outputs/    figures (png/svg) + summary csv/json, one folder per module
reports/    self-contained HTML result pages
docs/       write-ups, data inventory, method notes, evidence brief, figure guide
data/       raw data + credentials — GITIGNORED (reproduced by the fetchers)
env/        requirements.txt
```

Python + numpy/scipy, MNE-Python, **NeuroKit2** (R-peak detection), scikit-learn (staging),
matplotlib, pandas, `ieeg` (iEEG.org client).
