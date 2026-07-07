# iEEG — Human test of the Kipnis sleep-coordination substrate

Does the NREM slow oscillation organize faster activity (spindles, ripples) in **human medial
temporal lobe** — the candidate neural substrate that Kipnis / Jiang-Xie link to sleep-dependent
CSF clearance? We test this on **open, public** human intracranial EEG, with **no new patient
access or committee approval required**. (Clearance itself is not measured here.)

Sister project to the Buzsáki mouse pilot (`hpaticmousebuzsakilab`): the mouse showed 50 Hz
vibrotactile drive organized dorsal-hippocampal firing at **theta**, but not at the Kipnis
**0.5–4 Hz slow-wave** band — because it had no natural NREM. This project supplies that state.

**Start here:** [docs/iEEG_EVIDENCE_BRIEF.md](docs/iEEG_EVIDENCE_BRIEF.md) (2-page spine) ·
[docs/FIGURE_GUIDE.md](docs/FIGURE_GUIDE.md) (plain-English map of every figure) ·
[docs/PLAN.md](docs/PLAN.md).

## Datasets used (all open; raw data gitignored)

| # | Dataset | Rate | Coverage | Used by |
|---|---|---|---|---|
| 1 | **Normative iEEG Sleep & Wake Atlas** — Pattnaik & Litt 2024, Pennsieve DOI [10.26275/xhte-d11l](https://doi.org/10.26275/xhte-d11l) (HUP cohort; 30 s clips, W/N2/N3/R) | 204 Hz | mesiotemporal + neocortex | `slow_power_by_state`, `slow_spindle_coupling`, `mtl_region_profile` |
| 2 | **iEEG sleep + IED/HFO annotations** — Falach/Geva-Sagiv/Eliashiv 2024, [Figshare 26131978](https://figshare.com/articles/dataset/26131978) (~3 min NREM/subject) | 1 kHz | MTL SEEG (amygdala, hippocampus, entorhinal) | `slow_ripple_coupling`, `nesting_timecourse` |
| 3 | **Zurich HFO sleep atlas** — Fedele/Sarnthein, OpenNeuro [ds003498](https://openneuro.org/datasets/ds003498) (5-min SWS runs × 25–35/subject; expert HFO markings) | 2 kHz | MTL SEEG | `hfo_slow_phase` |
| 4 | **iEEG.org continuous** — e.g. `HUP165_phaseII` (~19 days continuous; one full NREM night pulled to `data/ieeg_portal/HUP165_night1/`) | 1024 Hz | MTL depth (LA/LB/LC/LH) | follow-on within-night analysis |

Fetchers: dataset 1 via `analysis/atlas.py` (Pennsieve API); datasets 2–3 via public download;
dataset 4 via `analysis/ieeg_portal.py` + `analysis/ieeg_pull_night.py` (see
[docs/IEEG_ORG_SETUP.md](docs/IEEG_ORG_SETUP.md)).

## Analysis modules (one script → one result → one figure)

| Script | What it shows | Dataset |
|---|---|---|
| `slow_power_by_state` | slow-band (0.5–4 Hz) power rises wake→N3 (p = 2.9×10⁻⁶) | 1 |
| `slow_spindle_coupling` | slow-oscillation → spindle coupling, NREM-specific (p = 0.030) | 1 |
| `mtl_region_profile` | slow-power & coupling by MTL region (descriptive; pairwise n.s.) | 1 |
| `slow_ripple_coupling` | slow-oscillation → ripple (p = 0.0015) **and** → spindle (p = 6×10⁻⁵, replicates coupling) | 2 |
| `nesting_timecourse` | time-domain SO-triggered averages: spindle/ripple ride the SO up-state | 2 |
| `hfo_slow_phase` | SO-phase of expert-marked ripples/fast-ripples (independent cohort) | 3 |
| `mouse_human_bridge` | the mouse-theta vs human-slow bridge figure | — |
| `make_result_visuals` | explanatory raw-signal / spectrum / slopegraph figures | 1 |

## Layout
```
analysis/       one script per module (above) + data fetchers (atlas, ieeg_portal, ieeg_pull_night)
outputs/        figures (png/svg) + summary csv, one folder per module; all_figures/ = consolidated
docs/           evidence brief, figure guide, per-module writeups, plan, iEEG.org setup
data/           downloaded raw data + credentials — GITIGNORED (reproduced by the fetchers)
env/            requirements.txt
```
