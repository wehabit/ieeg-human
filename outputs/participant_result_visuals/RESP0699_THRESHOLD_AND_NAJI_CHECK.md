# RESP0699 threshold and Naji-style check

This file is generated from `figure_values.json`; do not edit its numerical
claims by hand. It accompanies the RESP0699 3A threshold and 3B Naji-check
figures produced by analysis `2026-07-qc-sensitivity-v9` under locked
profile `overlap11_endpoint_local`.

## 1. The 3A coherence threshold

The nominal pointwise threshold is
`Ccrit = 1 - alpha^(1 / (K - 1))`. RESP0699 had
`K = 4` accepted
256-second windows,
`alpha = 0.05`, and therefore
`Ccrit = 0.63159685`. Window starts were
990, 1,118, 1,246, 1,650 s, with 128 s overlap.

The closest Fourier bin to 0.02 Hz was
0.01953125 Hz. Coherence there was
0.679199, with nominal uncorrected
pointwise probability 0.033. A
three-participant Bonferroni sensitivity gives
0.745,
which the observed value does not exceed. These are nominal analytic
sensitivities—not a participant-learned threshold, cohort result, F3/F4
measurement, or proof of LC/norepinephrine activity.

## 2. Whether RESP0699 can replicate Naji et al.

It cannot provide an exact F3/F4 replication. Source metadata show
0 scalp EEG channels, no F3, and no F4.
The available adaptation uses 5 right-frontal ECoG
contacts (C01, C17, C25, C33, C41) during 510 s of
author-SWS-selected N3-like data. Across 167 SOs:

- HR peak above the whole-stage mean:
  +1.414%
- local pre/post HR change:
  +0.455%
- mean contact-specific HR-peak time:
  2.35 s
- participant-average curve peak:
  3.00 s

Naji reported group means of
3.35 ±
1.01% in SWS and
12.09 ±
1.48% in Stage 2. RESP0699 lacks the
paper's bilateral referenced scalp montage, expert Stage-2 endpoint,
Pan–Tompkins-plus-visual ECG validation, and behavioral texture-discrimination
measure. Its percentile-based intracranial SO detector is also an adaptation.
Event-locking p/z inference remains disabled.

The correct conclusion is: RESP0699 shows a small SO-following cardiac
acceleration in an intracranial adaptation. It does not establish an F3/F4
Naji replication or validate an LC proxy.

Primary paper: [Naji et al. 2019](https://escholarship.org/uc/item/5393b9zk).
