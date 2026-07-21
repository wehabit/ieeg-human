# 3A method notes — infraslow RR ↔ sigma coherence, and three ways to get it wrong

Validated implementation: `analysis/infraslow_rr_sigma_coherence.py`.
Every choice below was forced by a control test, not by preference.

## The validated method

1. **EKG → R-peaks** with NeuroKit2, run in **chunks** (~300 s). Its artifact correction misbehaves
   on multi-million-sample arrays; per-chunk it is reliable. Verified: 91.5 bpm, 0.02% RR dropped.
2. **RR → instantaneous HR** on a 4 Hz grid; drop only physiologically impossible RR
   (<0.33 s / >1.5 s). Do **not** additionally reject beats for deviating from a local median —
   that discards the very HR variability being measured (an early version flagged 61% of beats).
3. **MTL/cortex → sigma envelope**: notch, band-pass 11–16 Hz, Hilbert amplitude, per-channel
   median-normalise, then average channels; resample to the same 4 Hz grid.
4. **High-pass both at 0.005 Hz** — removes drift only, keeps the whole infraslow band.
5. **Magnitude-squared coherence (Welch)** on those broadband series; read the value **at 0.02 Hz**
   and as the 0.01–0.03 Hz max.
6. **Significance: analytic MSC threshold** `crit = 1 − α^(1/(K−1))` for K Welch segments;
   Bonferroni over band bins for the band-max. No surrogates.

## Three errors that each independently invalidated results

**1. Channel columns come back re-ordered.** `ds.get_data(...)` returns columns in **ascending
channel-index order**, *not* the order requested. `EKG1` is index 0 in HUP165, so it returns as
column 0. Reading "the last column" as EKG meant running heartbeat detection on an LFP channel
(and the sigma average silently included the EKG). Symptom: beat counts of ~34 bpm instead of ~91.
Fix: remap columns via `argsort` after the pull.

**2. Do not band-pass before computing coherence.** Filtering both signals to 0.01–0.03 Hz and then
running a frequency-resolved coherence is redundant and destroys the estimate.
*Control:* two signals that literally share a 0.02 Hz rhythm scored **0.087** (i.e. "nothing") with
the pre-filter, and **0.978** without it. Band-pass is for the display overlay only.

**3. Shift/phase surrogates are invalid for oscillatory coupling.** MSC is invariant to a constant
phase offset, so a circularly-shifted 50 s rhythm is *still perfectly coherent* with the original.
Phase-randomisation fails identically: it assigns a random but **globally constant** phase, which
coherence reads as perfect consistency across segments.
*Control:* a genuinely coupled pair gave observed 0.978 vs surrogate 95th-percentile 0.982 → p=0.31,
i.e. the surrogate declared real coupling non-significant. Use the analytic threshold instead.

## Control results the method must reproduce

| Control | @0.02 Hz | crit (bin) | verdict |
|---|---|---|---|
| shared 0.02 Hz rhythm + independent noise | 0.978 | 0.095 | SIGNIFICANT ✅ |
| two independent noise series | 0.020 | 0.095 | not significant ✅ |

## Interpretation guardrail

A single bin crossing the **uncorrected** per-bin threshold inside an otherwise flat coherence
spectrum is **not** evidence of coupling — check (a) the band-corrected threshold and (b) whether a
~0.02 Hz peak exists in each signal's own power spectrum at all. If neither signal has an infraslow
spectral peak, there is no rhythm to couple and a marginal crossing is a threshold artifact.
