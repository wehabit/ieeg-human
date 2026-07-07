# Slow-power by state — anchor result (tests H1)

**Question:** In human mesiotemporal cortex, does slow-wave (0.5–4 Hz) power rise from wake into
deep NREM (W → N2 → N3), as the Kipnis / Jiang-Xie sleep-coordination mechanism predicts? This is
the natural-sleep contrast the Buzsáki mouse pilot could not run (it had only 50 Hz drive, no
scored NREM), and it is the entorhinal analogue of the mouse LEC arm.

**Data:** Normative iEEG Sleep & Wake Atlas (Pattnaik & Litt 2024, Pennsieve DOI 10.26275/xhte-d11l),
HUP cohort, bipolar, 30 s clips, 204 Hz. No new patient access or committee approval required.

## Plain-English answer (full cohort)

Yes. In entorhinal + parahippocampal channels, relative slow-band power rose from wake to N3 in
**20 of 21 patients** (mean 0.586 → 0.775; Wilcoxon **p = 2.9×10⁻⁶**, N = 21). The same W → N2 → N3
rise appears in hippocampal and temporal-neocortical channels.

This is the human slow-wave state signature the mouse pilot lacked. It sets up the central contrast:
mouse 50 Hz drive produced **theta** organization; natural human NREM produces **slow-band**
organization.

## Result (relative slow 0.5–4 Hz, patient-level mean; full cohort)

| ROI | W | N2 | N3 | R | N (W/N3) |
|---|---:|---:|---:|---:|---:|
| entorhinal | 0.621 | 0.664 | 0.739 | 0.626 | 10/10 |
| parahippocampal | 0.593 | 0.664 | 0.804 | 0.707 | 26/20 |
| hippocampal | 0.620 | 0.678 | 0.754 | 0.708 | 48/33 |
| temporal neocortex | 0.574 | 0.678 | 0.835 | 0.691 | 56/43 |

## Method

- Per channel, per 30 s clip: detrend → 60 Hz notch → Welch PSD (4 s windows, 50% overlap) →
  band power normalized to 0.5–45 Hz total. Averaged across clips → channel; channel → patient →
  ROI group; patient-level bootstrap 95% CI.
- Normative filter: SOZ and resected channels dropped.
- H1 test: entorhinal+parahippocampal, wake vs N3, paired Wilcoxon across patients (N = 21 with
  both states). (A 12-patient pilot gave the same direction at p = 0.078; the full cohort above is
  the reported result.)

## Caveats

- **Relative power, not coordination.** Slow-power shows the slow rhythm is *present and state-graded*.
  Whether it *coordinates* faster events (the field-wave claim) is Spindle-coupling / Ripple-coupling.
- **Nyquist 102 Hz** → no ripple/HFO band in this dataset; spindle coupling is the accessible
  nesting signature here (ripples are added in Ripple-coupling on a 1 kHz dataset).
- Normative epilepsy patients, interictal; not healthy brains. **No clearance readout** — slow-wave
  coordination is a *candidate neural substrate* that Kipnis links to clearance, not a measurement
  of clearance.

## Outputs
- `outputs/slow_power_by_state/channel_state_bandpower.csv` — per channel × state, all bands.
- `outputs/slow_power_by_state/group_state_slow.csv` — group × state slow power + CI.
- `outputs/slow_power_by_state/slow_power_by_state.{png,svg}`.
