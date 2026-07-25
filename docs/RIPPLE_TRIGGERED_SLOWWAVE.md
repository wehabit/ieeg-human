# Ripple-triggered slow-wave average in human MTL

> **SEPARATE HISTORICAL STUDY.** These numbers are not outputs of the current LC-proxy v8
> 3A/3B/3D pipelines and must not be used as corrected LC-tracking evidence. Their original
> study-specific validation/inference claims require separate review.

The direct human parallel to the mouse **spike-triggered LFP** (Kipnis Fig 1c): trigger on each
expert-marked ripple and average the surrounding **slow-wave (0.5–4 Hz) LFP**. If ripples ride a
consistent slow-wave field, the average is a clear slow wave; if not, it is flat.

- **Mouse:** spike → surrounding LFP wave.
- **Human (this figure):** ripple/HFO event → surrounding slow-wave LFP.

**Data:** Zurich ds003498 (2 kHz), **expert-marked ripples** (independent of the slow wave → the
average is non-circular). Per bipolar MTL channel with ≥30 ripples; SO computed at 500 Hz.

**Polarity handled correctly:** bipolar SO polarity varies by channel, so the analysis is **per
channel** (no cross-channel averaging of signed waves). Significance is the peak-to-peak amplitude of
the ripple-triggered slow-wave average vs a **shuffled-ripple-time surrogate** (200 shuffles) — an
amplitude test, so polarity-free.

## Result
- **63 channels** (≥30 ripples). Ripple-triggered slow-wave modulation **z > 2 in 54/63 (86%)**,
  **median z = 7.92**. In the large majority of channels, ripples sit on a consistent slow wave.

This is the **time-domain, waveform form** of the coupling shown elsewhere: `hfo_slow_phase` tested
the SO *phase* of the same marked ripples (Rayleigh); the phase-aligned analysis tested amplitude
modulation; this shows the actual slow wave the ripples ride. Together: **human MTL slow oscillations
organize ripple activity**, shown three complementary ways.

## Caveats
- Zurich HFO markings include some **pathological** ripples; the surrogate test shows the slow-wave
  locking is real, but pathological and physiological ripples are not separated here.
- Per-channel (channels not independent within subject); 8/9 subjects contribute.

## Outputs
- `outputs/ripple_triggered_slowwave/ripple_triggered_slowwave.png` — example per-channel
  ripple-triggered slow waves (top by z) + z distribution.
- `outputs/ripple_triggered_slowwave/channel_z.csv` — per-channel z + ripple count.
