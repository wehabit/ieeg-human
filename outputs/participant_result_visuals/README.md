# Participant-level visual guide

These plots are recomputed from the final v8 neutral caches with the locked
`overlap11_endpoint_local` profile. They contain real deidentified participant-derived data, not
synthetic teaching signals. The examples were selected because all relevant measurements were
available, not because they gave the most favorable biological result.

## Paper map

| Question | Main literature source | What is tested here |
|---|---|---|
| 3A | Lecci et al. 2017 | Infraslow 10–15 Hz sigma spectrum and its relationship to heart rate |
| 3B | Naji et al. 2019; Dang-Vu et al. for SO morphology | Heart-rate/RR response following an SO down-state trough |
| 3D | Staresina et al. 2015 and Helfrich et al. 2018 | Spindle timing relative to SO phase |

Naji is not the source method for 3A. None of these endpoints directly measures human LC firing or
norepinephrine.

## 3A — RESP0699

`3A_RESP0699_all_three_measurements` shows:

1. an actual NREM excerpt of sigma power and heart rate;
2. an accepted sigma-spectrum peak at 0.01686 Hz;
3. fixed-0.02-Hz coherence of 0.679, above this participant's analytic threshold of 0.632; and
4. a positive paper-direction cross-correlation of `r = 0.303`, whose best lag is 0 s rather than
   the approximately +5 s Lecci expectation.

The participant passes availability for all three measurements, but the timing does not fully
match the paper prediction.

### Saved 3A threshold explanation

`3A_RESP0699_coherence_threshold_explained` records how the nominal coherence threshold was
calculated. RESP0699 contributed four accepted 256-second Welch windows, with 128-second overlap:

| Window | Recording interval |
|---:|---:|
| 1 | 990–1,246 s |
| 2 | 1,118–1,374 s |
| 3 | 1,246–1,502 s |
| 4 | 1,650–1,906 s |

For `K = 4` windows and pointwise `alpha = 0.05`, the saved analysis uses

`1 - 0.05^(1 / (4 - 1)) = 0.63159685`.

The actual evaluated Fourier bin is 0.01953125 Hz, where coherence is 0.679199. The corresponding
nominal uncorrected pointwise probability is approximately 0.033. This threshold was not learned
from the participant, is unrelated to F3/F4, and is not corrected across participants. An
illustrative Bonferroni sensitivity across the three estimable RESP participants gives a threshold
of approximately 0.745, which this participant does not exceed. This sensitivity is not a new
preregistered decision rule.

## 3B — HUP160

`3B_HUP160_event_locked_heart_rate` shows SO-trough-aligned heart-rate curves:

| State | SO events | Local HR change | Mean contact HR-peak time |
|---|---:|---:|---:|
| N2-like | 635 | +0.606% | 2.04 s |
| N3-like | 902 | +0.145% | 2.50 s |
| Pooled NREM | 1,878 | +0.557% | 2.88 s |

These are descriptive measurements. The event-locking p/z fields remain disabled because the
current whole-stage shift does not preserve local nonstationary cardiac trends and clustered SO
timing.

### Saved RESP0699 check against Naji

`3B_RESP0699_naji_frontal_ecog_check` records the nearest available Naji-style comparison for the
same RESP0699 participant used in the 3A example.

The source metadata declare 64 right-hemisphere ECoG contacts, zero scalp EEG channels, and no F3
or F4. Consequently, this is a right-frontal ECoG adaptation, not an F3/F4 replication.

| Measurement | RESP0699 result | Naji reference |
|---|---:|---:|
| Stable SWS/N3-like support | 510 s | uninterrupted 3-min bins |
| Event-supported contacts | 5 | scalp F3 and F4 |
| SO events | 167 | not directly comparable |
| HR peak above stage mean | +1.414% | SWS: +3.35 ± 1.01% |
| Local pre/post HR change | +0.455% | not the paper's reported baseline |
| Mean contact HR-peak time | 2.35 s | Stage 2: approximately 2.10–2.23 s |
| N2 endpoint | unavailable | Stage 2: +12.09 ± 1.48% |

RESP0699 agrees only qualitatively: heart rate rises after the SO trough on an approximately
2–3-second timescale. It does not match Naji's published SWS group magnitude, cannot test the
Stage-2-versus-SWS contrast, has no valid event-locking p/z, and lacks the paper's behavioral
processing-speed endpoint. The source paper is
[Naji et al. 2019](https://escholarship.org/uc/item/5393b9zk).

The complete written interpretation and evidence map are saved in
`RESP0699_THRESHOLD_AND_NAJI_CHECK.md`.

## 3D — HUP172

`3D_HUP172_SO_spindle_phase` shows:

1. which contacts meet the base support requirements;
2. the separately normalized phase distribution for each qualified contact; and
3. equal-contact complex-vector aggregation.

Qualified contacts are LD12 (802 pairs), LJ12 (799), and LM11 (820). Their equal-contact
participant vector has `R = 0.192` and preferred phase 49.7°. The base 80% valid-NREM contact rule
shown in panel A is a repository sensitivity rule, not a paper-derived biological cutoff.

The vector is descriptive only. Inference is disabled because nearest-event selection inside a
finite SO-centered window can induce apparent phase structure under independence.

## Reproduction

```bash
.venv/bin/python visualization/make_participant_result_visuals.py
```

Exact plotted values and source paths are recorded in `figure_values.json`. PNG files are convenient
for slides; SVG files are preferable when editing or exporting at larger sizes.
