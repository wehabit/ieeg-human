# Branch summary — current v9 evidence

This branch tests three distinct pieces of human sleep physiology. They are motivated by LC/NE
work, but none measures the locus coeruleus (LC) or norepinephrine directly.

| Label | Question | Current interpretation |
|---|---|---|
| 3A | Do infraslow sigma-power dynamics covary with cardiac dynamics? | Individual candidate features exist, but the cohort pattern is heterogeneous; not an LC biomarker |
| 3B | Does RR/heart rate change around cortical slow-oscillation troughs? | A small descriptive acceleration is the most consistent observation; inference is disabled and the effect is not LC-specific |
| 3D | At what SO phase do spindle events occur? | Generic descriptive SO–spindle nesting/QC; inference is disabled and this is not an LC proxy |

## Current v9 availability

Availability means that the endpoint could be computed under the locked
`overlap11_endpoint_local` profile. It does not mean that its hypothesis was supported.

| Cohort | Requested/analyzed | 3A spectrum/coherence/xcorr | 3B N2/N3/pooled | 3D descriptive |
|---|---:|---:|---:|---:|
| HUP | 25/24; HUP116 structured skip | 20/14/16 | 6/5/15 | 7 |
| RESPect | 6/6 | 3/3/3 | 0/1/2 | outside scope |

RESPect does not reach the default five-participant reporting gate. HUP reaches that availability
gate, but the tracked artifacts do not contain a valid cohort inferential test for 3A and
explicitly disable 3B and 3D inference.

## What changed during the final rebuild

- Endpoint-local 3A now uses shared finite EEG support for spectra and a separate EEG-plus-HR mask
  for cardiac endpoints. HUP availability changed from 19/14/15 to 20/14/16 because HUP138 gained
  a valid spectrum and cross-correlation, while coherence remained unavailable.
- In the paired scalp–iEEG comparison, this correction removed HUP187's accepted iEEG spectral
  peak. The both-accepted paired peak sample changed from three to two (HUP160 and HUP212).
  Paired coherence, cross-correlation, and all paired 3B estimates were unchanged.
- The flat export no longer double-counts iEEG 3B rows for F3-present participants. The separate
  role-pair table intentionally retains one comparator per role.
- HUP138 is the only participant with geometry-matched F3 and F4 labels, but both streamed signals
  are exactly flat. An exact bilateral Naji-style participant check is therefore unavailable for
  an objective signal-quality reason.

## What remains unresolved

- HUP N2-like/N3-like values are unvalidated algorithmic proxies, not expert stages.
- HUP contacts need coordinate, tissue, region, pathology/SOZ, and clinical review.
- RESPect and HUP iEEG sensors are adaptations of scalp sources, not validated homologues.
- ECG, SO, spindle, artifact, and missing-boundary detectors need blinded raw-data validation.
- HUP intervals are sparse-probe high-delta candidates, not verified sleep onset.
- 3B needs a local-trend/dependence-preserving null; 3D needs a pairing-aware complete-train null.
- No current signal validates LC specificity. That requires an independent LC/NE-sensitive
  measurement or intervention.

The authoritative numerical interpretation is
[QC_SENSITIVITY_RESULTS_2026-07.md](QC_SENSITIVITY_RESULTS_2026-07.md). Detailed methods are in
[METHODS.md](METHODS.md), the simultaneous sensor comparison is in
[PAIRED_SCALP_IEEG_RESULTS_2026-07.md](PAIRED_SCALP_IEEG_RESULTS_2026-07.md), and every confirmed
defect/open concern is tracked in [ISSUE_REGISTER_2026-07.md](ISSUE_REGISTER_2026-07.md).
