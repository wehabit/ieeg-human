# ds003848 RESPect — current v9 sensitivity

> **CURRENT V9 RESULT:** the complete six-subject neutral cache and sensitivity grids pass exact
> cache-code, input-identity, acquisition, calibration, profile, result-file, and runtime checks. The
> endpoint-local sensitivity recovers individual 3A and 3B estimates, but no endpoint reaches the
> default cohort minimum of five. This is limited availability, not positive or negative evidence
> for an LC proxy.

## Current v9 result

- Cache: 6 completed, 0 skipped, 0 failed. Each subject uses seven files pinned to OpenNeuro
  snapshot 1.0.3 by S3 version ID, byte size, and SHA-256.
- Endpoint-local 3A: spectrum 3/6, cross-correlation 3/6, and fixed-0.02-Hz coherence 3/6.
  RESP0699 and RESP0724 have accepted peaks at 0.0169 and 0.0221 Hz; only RESP0699's pointwise
  coherence exceeds its analytic threshold. The Lecci-direction maxima occur at 0, 0, and 5 s.
- Endpoint-local 3B: N2 0/6, N3 1/6, and exploratory pooled NREM 2/6. RESP0699 N3 has a
  descriptive local HR change of +0.455% at a mean-channel latency of 2.35 s.
- Historical `audit80`: 3A remains 0/6 because several locally introduced gates stack; 3B is
  0/1/2 because the current endpoint-local design no longer lets a global power gate suppress
  endpoint-local cardiac estimates.
- None reaches the default cohort minimum of five. No 3B p/z value is produced.

The fixed 80% gate was neither a paper requirement nor preregistered. RESP0521 had approximately
79.2% historical aggregate coverage and RESP0699 approximately 79.7%, illustrating why a binary
79-versus-80 declaration is not scientifically justified. The locked one-axis grid shows which
support choices change availability; the endpoint-local profile was designed without using the
biological result direction. See
[QC_SENSITIVITY_RESULTS_2026-07.md](QC_SENSITIVITY_RESULTS_2026-07.md).

The audit found three additional reasons the old values cannot be interpreted:

- the legacy cache ignored the author `events.tsv`, including NREM/REM/SWS selections,
  transition, artifact, and seizure intervals; direct overlap showed that the proxy called many
  author-REM epochs N2/N3;
- it ignored `electrodes.tsv`, allowing documented SOZ/resected/edge and other non-normative
  contacts into endpoints; and
- 3B did not enforce Naji's uninterrupted 3-minute stable-stage rule.

V9 makes author annotations primary, leaves author-unknown sleep unclassified,
excludes annotated disturbances, filters documented pathological/non-cortical contacts, and
intersects motivated parietal/frontal Destrieux ROIs. These corrections can sharply reduce the
number of estimable participants. They do not provide expert AASM/R&K N2/N3 scoring, validate the
iEEG ROIs as homologues of the cited scalp sensors, remove possible postictal effects, or measure
LC directly.

Everything from “Historical cohort and staging” through the old numerical sections is retained only
to document what was withdrawn; it is not the v9 result.

The separate v9 HUP regeneration has 24 completed caches and one explicit skip, with many
endpoint-local subject estimates but unresolved staging/anatomy/inference limitations. Historical
HUP comparisons below remain withdrawn.

The HUP `phaseII` result rests on a staging *proxy* (a Gaussian mixture on slow-wave power, because
iEEG has no EOG/EMG). The historical RESPect analysis attempted to test whether wake/REM
contamination explained that result. It did not establish real REM/wake exclusion because it failed
to consume the dataset's author annotations.

## Historical cohort and staging

Six patients, ~1 h continuous `task-[Ss]leep` @ 2048 Hz, 50 Hz line; **3 ECoG grid + 3 SEEG depth**.
Every subject verified (from raw `channels.tsv`) to carry iEEG, ECG, EMG, EOG and (bad) respiration
belts. Pipeline: `analysis/stage_ds003848.py` (MNE BrainVision reader → neutral
cache plus EMG/EOG staging), followed by the authoritative profile-materialized
3A/3B/3D path in `analysis/run_qc_grid.py`.

**Staging is rule-based and not validated against expert labels.** True AASM scoring needs scalp EEG, which
this dataset lacks; what EMG + EOG buy is REM/Wake exclusion. Per 30 s epoch: submental EMG RMS,
EOG movement variance, iEEG slow-wave power, robust-z within subject. Wake = high EMG; REM = atonia +
phasic eye movement + low SWA; NREM split N2/N3 by GMM on SWA. On RESP0521 the labels separate exactly
as the physiology predicts:

| stage | n | EMG z | EOG z | SWA z |
|---|---|---|---|---|
| Wake | 15 | **+1.28** | +0.31 | +0.25 |
| REM | 3 | −0.55 | **+1.37** | −0.54 |
| N2 | 32 | −0.16 | −0.18 | −0.73 |
| N3 | 48 | −0.02 | +0.08 | **+0.51** |

Per-subject stage counts (W / R / N2 / N3, of ~120–128 epochs):
0521 15/3/32/48 · 0699 35/17/25/25 · 0724 19/1/6/75 · 0749 12/3/22/22 · 0779 25/5/39/41 · 0800 19/1/46/42.
Wake+REM is **6–41%** of epochs — real contamination that the HUP proxy could not remove, and here is
removed. Heart rates are also healthier than HUP (e.g. 70 bpm vs HUP's ~92 bpm tachycardia).

## Historical 3A — withdrawn

The catch is duration. HUP streams ~7 h/subject; ds003848 is **1 h**, and the stage proxy fragments NREM
into short runs, leaving very few bouts ≥120 s and a low Welch segment count:

| subject | bouts ≥120 s | coherence K | own peak (Hz) |
|---|---|---|---|
| RESP0521 | 5 | — | 0.035 |
| RESP0699 | 3 | — | 0.041 |
| RESP0724 | 4 | — | 0.052 |
| RESP0749 | **0** | — | — |
| RESP0779 | 3 | (low) | 0.027 |
| RESP0800 | 7 | (low) | 0.044 |

- **Step 1 (peaks).** 5/6 subjects yield a fitted peak; they sit at **0.027–0.052 Hz (mean 0.040)**,
  **0/5 inside 0.015–0.025 Hz**, and significantly *above* Lecci's 0.019 Hz (t p = 0.008). The peaks
  are somewhat tighter than this cohort's scale-free surrogates (SD 0.0093 vs 0.0137) but the bootstrap
  is not significant (p = 0.21), and they cluster near **0.04 Hz, not 0.019** — so there is no Lecci
  fingerprint. Sigma is not reliably more prominent than the SWA control (p = 0.31).
- **Step 2 (HR coupling).** Only **2/6** subjects have enough contiguous NREM for a coherence estimate
  at all (K median **7**, versus 59 in HUP) — 0/2 significant. The cross-correlation is flat (group
  |r| = 0.066; per-subject signed peak r vs 0: p = 0.09, and the trend is *negative* if anything).

The superseded analysis was previously interpreted as agreeing with a “negative” HUP result and as
showing that REM/wake exclusion did not rescue a rhythm. That interpretation is withdrawn: v9
recovers individual 3A estimates but does not provide a powered RESPect cohort result.

## Historical 3B — withdrawn

The superseded 3B analysis was previously described as a more trustworthy “null” because it used
thousands of SOs. Its values are retained below only as withdrawn provenance:

| | ds003848 (n=6) | HUP (n=23) | Naji 2019 |
|---|---|---|---|
| N2 HR peak | +0.65% (0/5 sig; t p=0.88) | +1.8% (7/23 sig; t p=0.014) | +12.09% |
| N3 HR peak | +0.99% (0/6 sig; t p=0.28) | +3.5% (6/20 sig; t p=0.020) | +3.35% |
| N2 ≫ N3? | no (p=0.63) | no (p=0.57) | 3.6× |

The legacy write-up contrasted these values with a nominally significant HUP result and speculated
about power or wake/arousal contamination. V9 recovers descriptive 3B estimates but disables
inference, so that legacy numerical contrast and its p-value interpretation remain withdrawn.

## What this replication does and does not establish

- **Does not establish:** that the HUP result survives genuine REM/wake exclusion. That claim
  requires a powered cohort and validation against expert staging; RESPect remains below the
  default cohort minimum.
- **Does not:** provide a *powered* independent 3A test. The legacy one-hour proxy analysis reported
  median K = 7 and only 2/6 usable coherence estimates. A definitive replication needs longer
  continuous sleep with cardiac + EOG/EMG, which no current public iEEG dataset provides.

## Honest limitations

- **1 h recordings** → few NREM bouts ≥120 s → 3A badly underpowered (RESP0749 unusable for 3A).
- **Rule-based staging without scalp EEG** — not validated against expert labels and not AASM-scored;
  the N2/N3 split is still SWA-driven, so it inherits that axis.
- **Strict current annotations** — explicit NREM/REM exists for only part of the cohort; unknown
  sleep remains unclassified and can make current endpoints unavailable.
- **Source transfer** — conservative Destrieux parietal/frontal iEEG intersections are not validated
  equivalents of Lecci's C3/parietal or Naji's F3/F4 scalp sensors and do not measure LC.
- **Seizure/postictal confounding** — exact author-marked seizures are excluded, but longer
  postictal autonomic/sleep effects need a clinically justified sensitivity exclusion.
- **n = 6**, mixed ECoG/SEEG montage.
- Respiration belts exist but are marked `bad`, so apnoea screening was not attempted.

## Reproduce

```bash
.venv/bin/python analysis/stage_ds003848.py --force
.venv/bin/python analysis/calibrate_staging_windows.py \
  --update-profile-pin

for grid in coverage_oat_v1 staging_window_support_v1 \
            auxiliary_window_support_v1 event_count_oat_v1; do
  .venv/bin/python analysis/run_qc_grid.py \
    --cache-dir data/derived/ds003848 --grid "$grid"
done
```

The withdrawn `run_ds003848_replication.py` and standalone endpoint writers
exit without producing results. Full grid artifacts are written under the
gitignored `outputs/qc_grid/`; the cache is
`data/derived/ds003848/` (also gitignored).
