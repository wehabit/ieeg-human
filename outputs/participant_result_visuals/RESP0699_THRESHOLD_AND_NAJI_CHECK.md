# RESP0699 threshold and Naji-style check

This document accompanies:

- `3A_RESP0699_coherence_threshold_explained.png` and `.svg`
- `3B_RESP0699_naji_frontal_ecog_check.png` and `.svg`

The figures are recomputed from the final v8 neutral cache under the locked
`overlap11_endpoint_local` profile. They contain deidentified participant-derived measurements,
not synthetic signals.

## 1. The 3A coherence threshold

The reported 0.63159685 value is a nominal pointwise 5% analytic threshold for
magnitude-squared coherence:

`Ccrit = 1 - alpha^(1 / (K - 1))`.

RESP0699 had `K = 4` accepted 256-second Welch windows and `alpha = 0.05`, so:

`Ccrit = 1 - 0.05^(1 / 3) = 0.63159685`.

The accepted windows began at 990, 1,118, 1,246, and 1,650 seconds. Adjacent windows in the long
NREM bout overlap by 128 seconds. Five of the 1,500 valid NREM seconds were short-gap fills.

The closest Fourier bin to the prespecified 0.02-Hz target is 0.01953125 Hz. RESP0699's
coherence at that bin is 0.679199, which exceeds 0.631597. Under the same nominal independent-window
null, the uncorrected pointwise probability is approximately:

`(1 - 0.679199)^3 = 0.033`.

This means only that the fixed-frequency point exceeds the repository's nominal individual
cutoff. It is not a threshold learned from this participant, not an F3/F4 measurement, not a
cohort-level result, and not proof of LC or norepinephrine activity.

The threshold is not adjusted for choosing a positive participant from the three RESP participants
with estimable coherence. As an illustrative sensitivity, dividing alpha by three gives a cutoff
of approximately 0.745, which 0.679 does not exceed. That calculation illustrates selection and
multiplicity risk; it is not a new preregistered endpoint.

Implementation evidence:

- `analysis/spectral_gapped.py`, function `analytic_msc_threshold`
- `analysis/run_qc_grid.py`, function `analyse_3a`
- `outputs/qc_grid/staging_window_support_v1/respect_qc_grid.json`

## 2. Whether RESP0699 can replicate Naji et al.

It cannot provide an exact F3/F4 replication.

The original recording metadata declare:

- 64 right-hemisphere ECoG contacts
- zero scalp EEG channels
- no F3 or F4 labels
- one ECG channel

Naji et al. used bilateral scalp F3/A2 and F4/A1, visually scored Stage 2 and SWS in uninterrupted
three-minute bins, Pan–Tompkins R-peak detection with visual confirmation, a 4-Hz RR series, and
slow oscillations detected after 0.15–4-Hz filtering.

The nearest current adaptation uses five right-frontal ECoG contacts—C01, C17, C25, C33, and
C41—during 510 seconds of stable author-SWS-selected N3-like data. Across 167 slow oscillations:

- HR peak above the whole-stage mean: +1.414%
- HR peak above the local pre-event mean: +0.455%
- mean of contact-specific HR-peak times: 2.35 seconds
- peak of the participant-average curve: 3.00 seconds

Naji reported an SWS group HR increase of +3.35 ± 1.01% and a Stage-2 increase of
+12.09 ± 1.48%. Their reported Stage-2 SO–HR delays were approximately 2.10–2.23 seconds across
visits.

Therefore, RESP0699 has partial qualitative agreement—the cardiac response points upward and
occurs on an approximately 2–3-second timescale—but it is not a direct replication:

- the montage is right-frontal ECoG rather than bilateral F3/F4 scalp EEG;
- expert Stage-2 scoring and a usable N2 endpoint are absent;
- the intracranial SO amplitude rule is percentile-based rather than the paper's scalp criterion;
- ECG peaks use NeuroKit and were not visually validated;
- the comparable SWS magnitude is smaller than the published group mean;
- valid event-locking p/z inference is disabled; and
- no texture-discrimination outcome is available.

The correct conclusion is:

> RESP0699 shows a small SO-following cardiac acceleration in an intracranial adaptation. It does
> not establish an F3/F4 Naji replication or validate an LC proxy.

Primary paper:

- Naji M, Krishnan GP, McDevitt EA, Bazhenov M, Mednick SC. 2019.
  https://escholarship.org/uc/item/5393b9zk

Local evidence:

- `data/ds003848_raw/sub-RESP0699_ses-1_task-sleep_run-030608_ieeg.json`
- `data/ds003848_raw/sub-RESP0699_ses-1_task-sleep_run-030608_channels.tsv`
- `data/ds003848_raw/sub-RESP0699_ses-1_task-sleep_run-030608_events.tsv`
- `outputs/qc_grid/staging_window_support_v1/respect_qc_grid.json`
