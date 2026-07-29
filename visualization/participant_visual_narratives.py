"""Deterministic Markdown narratives for participant-level teaching figures.

The prose in ``outputs/participant_result_visuals`` contains numerical claims.
Keeping those files as hand-edited companions allowed stale values to be
re-hashed after a result rebuild.  These renderers make every reported number a
direct function of ``figure_values.json`` instead.
"""
from __future__ import annotations


def _stage_row(stage: str, value: dict) -> str:
    labels = {
        "N2": "N2-like",
        "N3": "N3-like",
        "NREM": "Pooled NREM",
    }
    return (
        f"| {labels[stage]} | {int(value['n_so_total']):,} | "
        f"{float(value['local_change_pct']):+.3f}% | "
        f"{float(value['mean_channel_peak_lag_s']):.2f} s |"
    )


def render_visual_readme(values: dict) -> str:
    """Render the complete visual guide from one figure-values payload."""
    figures = values["figures"]
    result_3a = figures["3A"]
    threshold = figures["3A_threshold_explainer"]
    result_3b = figures["3B"]
    naji = figures["3B_RESP0699_Naji_check"]
    result_3d = figures["3D"]
    stages = result_3b["stages"]
    event_counts = result_3d.get("paired_events_by_contact", {})
    if event_counts:
        contact_text = ", ".join(
            f"{contact} ({int(event_counts[contact]):,} pairs)"
            for contact in result_3d["qualified_contacts"]
        )
    else:
        contact_text = ", ".join(result_3d["qualified_contacts"])

    return f"""# Participant-level visual guide

This file is generated from `figure_values.json`; do not edit its numerical
claims by hand. The figures were recomputed with analysis
`{values['analysis_version']}` and locked profile `{values['profile_id']}`.
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
{float(result_3a['peak_hz']):.5f} Hz, fixed-target coherence of
{float(result_3a['coherence_0p02']):.3f} against a nominal analytic threshold
of {float(result_3a['coherence_threshold']):.3f}, and paper-direction
cross-correlation `r = {float(result_3a['xcorr_r']):.3f}` at a best lag of
{float(result_3a['xcorr_lag_s']):.0f} s. Availability passes for all three
measurements, but the lag does not reproduce the approximately +5 s timing
expected from Lecci.

### Coherence-threshold explanation

RESP0699 contributed `K = {int(threshold['K'])}` accepted
{int(threshold['welch_window_seconds'])}-second Welch windows with
{int(threshold['welch_overlap_seconds'])}-second overlap. Their starts were
{", ".join(str(int(value)) for value in threshold['welch_window_starts_s'])} s.
For pointwise `alpha = {float(threshold['alpha']):.2f}`,
`1 - alpha^(1/(K-1)) = {float(threshold['analytic_threshold']):.8f}`.
The evaluated Fourier bin was
{float(threshold['evaluated_frequency_hz']):.8f} Hz, where coherence was
{float(threshold['observed_coherence']):.6f}; its nominal uncorrected
pointwise probability was {float(threshold['nominal_pointwise_p']):.3f}.
An illustrative three-participant Bonferroni sensitivity gives a threshold of
{float(threshold['three_participant_bonferroni_threshold_sensitivity']):.3f},
which this participant does not exceed. This is not a preregistered cohort
decision rule and is not evidence of LC specificity.

## 3B — HUP160

| State | SO events | Local HR change | Mean contact HR-peak time |
|---|---:|---:|---:|
{_stage_row('N2', stages['N2'])}
{_stage_row('N3', stages['N3'])}
{_stage_row('NREM', stages['NREM'])}

These measurements are descriptive. Event-locking p/z values remain disabled
because whole-stage shifts do not preserve local nonstationary cardiac trends
and clustered SO timing.

### RESP0699 check against Naji

RESP0699 has {int(naji['scalp_eeg_channel_count'])} scalp EEG channels and no
F3/F4 labels. Its nearest available comparison is therefore a
right-frontal ECoG adaptation using {int(naji['n_channels'])} contacts,
{int(naji['stable_seconds'])} s of stable support, and
{int(naji['n_so_total'])} SOs. The HR peak was
{float(naji['pct_above_stage_mean']):+.3f}% above the whole-stage mean,
the local pre/post change was
{float(naji['event_locked_local_change_pct']):+.3f}%, and the mean
contact-specific peak lag was {float(naji['mean_contact_peak_lag_s']):.2f} s.
Naji reported an SWS group mean of
{float(naji['naji_reference']['sws_hr_peak_pct']):.2f}% and a Stage-2 group
mean of {float(naji['naji_reference']['stage2_hr_peak_pct']):.2f}%.
This is qualitative timing context, not an F3/F4 replication, group
comparison, behavioral replication, or LC validation.

## 3D — HUP172

Qualified contacts are {contact_text}. Their equal-contact participant vector
has `R = {float(result_3d['participant_R']):.3f}` and preferred phase
{float(result_3d['participant_preferred_phase_deg']):.1f}°. The 80% valid-NREM
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
"""


def render_resp0699_check(values: dict) -> str:
    """Render the detailed RESP0699 threshold/Naji companion document."""
    figures = values["figures"]
    threshold = figures["3A_threshold_explainer"]
    naji = figures["3B_RESP0699_Naji_check"]
    contacts = ", ".join(naji["event_contact_ids"])
    starts = ", ".join(
        f"{int(value):,}" for value in threshold["welch_window_starts_s"])
    reference = naji["naji_reference"]

    return f"""# RESP0699 threshold and Naji-style check

This file is generated from `figure_values.json`; do not edit its numerical
claims by hand. It accompanies the RESP0699 3A threshold and 3B Naji-check
figures produced by analysis `{values['analysis_version']}` under locked
profile `{values['profile_id']}`.

## 1. The 3A coherence threshold

The nominal pointwise threshold is
`Ccrit = 1 - alpha^(1 / (K - 1))`. RESP0699 had
`K = {int(threshold['K'])}` accepted
{int(threshold['welch_window_seconds'])}-second windows,
`alpha = {float(threshold['alpha']):.2f}`, and therefore
`Ccrit = {float(threshold['analytic_threshold']):.8f}`. Window starts were
{starts} s, with {int(threshold['welch_overlap_seconds'])} s overlap.

The closest Fourier bin to 0.02 Hz was
{float(threshold['evaluated_frequency_hz']):.8f} Hz. Coherence there was
{float(threshold['observed_coherence']):.6f}, with nominal uncorrected
pointwise probability {float(threshold['nominal_pointwise_p']):.3f}. A
three-participant Bonferroni sensitivity gives
{float(threshold['three_participant_bonferroni_threshold_sensitivity']):.3f},
which the observed value does not exceed. These are nominal analytic
sensitivities—not a participant-learned threshold, cohort result, F3/F4
measurement, or proof of LC/norepinephrine activity.

## 2. Whether RESP0699 can replicate Naji et al.

It cannot provide an exact F3/F4 replication. Source metadata show
{int(naji['scalp_eeg_channel_count'])} scalp EEG channels, no F3, and no F4.
The available adaptation uses {int(naji['n_channels'])} right-frontal ECoG
contacts ({contacts}) during {int(naji['stable_seconds'])} s of
{naji['stage']} data. Across {int(naji['n_so_total'])} SOs:

- HR peak above the whole-stage mean:
  {float(naji['pct_above_stage_mean']):+.3f}%
- local pre/post HR change:
  {float(naji['event_locked_local_change_pct']):+.3f}%
- mean contact-specific HR-peak time:
  {float(naji['mean_contact_peak_lag_s']):.2f} s
- participant-average curve peak:
  {float(naji['participant_curve_peak_lag_s']):.2f} s

Naji reported group means of
{float(reference['sws_hr_peak_pct']):.2f} ±
{float(reference['sws_plus_minus']):.2f}% in SWS and
{float(reference['stage2_hr_peak_pct']):.2f} ±
{float(reference['stage2_plus_minus']):.2f}% in Stage 2. RESP0699 lacks the
paper's bilateral referenced scalp montage, expert Stage-2 endpoint,
Pan–Tompkins-plus-visual ECG validation, and behavioral texture-discrimination
measure. Its percentile-based intracranial SO detector is also an adaptation.
Event-locking p/z inference remains disabled.

The correct conclusion is: RESP0699 shows a small SO-following cardiac
acceleration in an intracranial adaptation. It does not establish an F3/F4
Naji replication or validate an LC proxy.

Primary paper: [Naji et al. 2019](https://escholarship.org/uc/item/5393b9zk).
"""
