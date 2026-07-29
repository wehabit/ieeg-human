# Participant-level visual guide

This file is generated from `figure_values.json`; do not edit its numerical
claims by hand. The figures were recomputed with analysis
`2026-07-qc-sensitivity-v9` and locked profile `overlap11_endpoint_local`.
They contain deidentified participant-derived measurements, not synthetic
teaching signals. Examples are selected for endpoint availability, not for a
favorable biological direction.

## Paper map

| Question | Main literature source | What is tested here |
|---|---|---|
| 3A | Lecci et al. 2017 | Infraslow 10–15 Hz sigma spectrum and its relationship to heart rate |
| 3B | Naji et al. 2019; Dang-Vu et al. for SO morphology | Heart-rate/RR response following an SO down-state trough |
| 3D | Staresina et al. 2015 and Helfrich et al. 2018 | Spindle timing relative to SO phase |

Naji is not the source method for 3A. None of these endpoints directly
measures human LC firing or norepinephrine.

## 3A — RESP0699

`3A_RESP0699_all_three_measurements` shows an accepted sigma-spectrum peak at
0.01686 Hz, fixed-target coherence of
0.679 against a nominal analytic threshold
of 0.632, and paper-direction
cross-correlation `r = 0.303` at a best lag of
0 s. Availability passes for all three
measurements, but the lag does not reproduce the approximately +5 s timing
expected from Lecci.

### Coherence-threshold explanation

RESP0699 contributed `K = 4` accepted
256-second Welch windows with
128-second overlap. Their starts were
990, 1118, 1246, 1650 s.
For pointwise `alpha = 0.05`,
`1 - alpha^(1/(K-1)) = 0.63159685`.
The evaluated Fourier bin was
0.01953125 Hz, where coherence was
0.679199; its nominal uncorrected
pointwise probability was 0.033.
An illustrative three-participant Bonferroni sensitivity gives a threshold of
0.745,
which this participant does not exceed. This is not a preregistered cohort
decision rule and is not evidence of LC specificity.

## 3B — HUP160

| State | SO events | Local HR change | Mean contact HR-peak time |
|---|---:|---:|---:|
| N2-like | 635 | +0.606% | 2.04 s |
| N3-like | 902 | +0.145% | 2.50 s |
| Pooled NREM | 1,878 | +0.557% | 2.88 s |

These measurements are descriptive. Event-locking p/z values remain disabled
because whole-stage shifts do not preserve local nonstationary cardiac trends
and clustered SO timing.

### RESP0699 check against Naji

RESP0699 has 0 scalp EEG channels and no
F3/F4 labels. Its nearest available comparison is therefore a
right-frontal ECoG adaptation using 5 contacts,
510 s of stable support, and
167 SOs. The HR peak was
+1.414% above the whole-stage mean,
the local pre/post change was
+0.455%, and the mean
contact-specific peak lag was 2.35 s.
Naji reported an SWS group mean of
3.35% and a Stage-2 group
mean of 12.09%.
This is qualitative timing context, not an F3/F4 replication, group
comparison, behavioral replication, or LC validation.

## 3D — HUP172

Qualified contacts are LD12 (764 pairs), LJ12 (766 pairs), LM11 (789 pairs). Their equal-contact participant vector
has `R = 0.193` and preferred phase
48.5°. The 80% valid-NREM
contact-support rule shown in the plot is repository-defined, not
paper-derived. The vector is descriptive only; inference is disabled because
nearest-event selection within a finite SO-centered window can create apparent
phase structure under independence.

## Reproduction

```bash
.venv/bin/python visualization/make_participant_result_visuals.py
```

PNG files are convenient for slides; SVG files are preferable for editing.
Exact source and result hashes are recorded in `RUN_MANIFEST.json`.
