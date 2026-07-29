# 3B method comparison — current descriptive adaptation vs Naji 2019

> The method comparison below describes the current RR-domain adaptation. Only the numerical
> section at the end is legacy/withdrawn. Whole-stage shifts are diagnostic only, no 3B
> event-locking p/z claim is available, and the direct legacy entry point is hard-stopped.

**Purpose.** A step-by-step comparison of how we implemented the slow-oscillation → heart-rate
coupling test against the method it follows. Written because the coupling was cited to the wrong
paper and described as "validated in our own data" — this documents what we actually did and how it
differs from the source.

## Which paper is the source

The empirical methods paper is **Naji, Krishnan, McDevitt, Bazhenov & Mednick (2019)**, *"Timing
between Cortical Slow Oscillations and Heart Rate Bursts during Sleep Predicts Temporal Processing
Speed, but Not Offline Consolidation"*, J Cogn Neurosci 31(10):1484–1490
([DOI](https://doi.org/10.1162/jocn_a_01432) · [lab PDF](http://sleepandcognitionlab.org/wp-content/uploads/2019/07/Timing-between-Cortical-Slow-Oscillations-and-Heart-Rate-Bursts-during-Sleep-Predicts-Temporal-Processing-Speed-but-Not-Offline-Consolidation.pdf)).
It is a scalp PSG/ECG and behavioral study, not a direct LC measurement.

The **PNAS 2022** paper (Chen, Zhang, Thayer & Mednick, *"Understanding the roles of central and
autonomic activity during sleep…"*, 119(44)) is a **review** — its Data Availability statement reads
verbatim *"There are no data underlying this work."* It proposes a framework and reviews others'
findings; it does not itself demonstrate the SO→HR relationship. Mednick is a coauthor on Naji 2019,
so "Mednick and colleagues" is fine when pointed at Naji 2019.

**There is no public code repository for Naji 2019** (no GitHub/OSF). The comparison below is against
the paper's Methods (which cite Dang-Vu et al. 2008 for the SO-detection criteria). To obtain their
actual code, contact the Mednick lab (UC Irvine Sleep & Cognition Lab).

## Our code

- [`analysis/event_3b_estimators.py`](../analysis/event_3b_estimators.py) —
  the single reusable RR-domain estimator implementation.
- [`analysis/event_3B_mednick.py`](../analysis/event_3B_mednick.py) —
  withdrawn command-line entry point and compatibility re-exports only; it
  does not maintain a second estimator.
- [`analysis/event_3B_cached.py`](../analysis/event_3B_cached.py) —
  importable profile-materialized estimator; its standalone writer is withdrawn.
- [`analysis/run_qc_grid.py`](../analysis/run_qc_grid.py) —
  the authoritative HUP and RESPect 3A/3B/3D runner.
- [`analysis/cohort_stages_3ABD.py`](../analysis/cohort_stages_3ABD.py) —
  a hard-stopped compatibility shim; its superseded whole-night-null
  estimator has been removed.
- [`analysis/test_3B_null.py`](../analysis/test_3B_null.py) — regression test for the null fix.
- [`analysis/scalp_f3_f4_exploratory.py`](../analysis/scalp_f3_f4_exploratory.py) —
  a separate HUP138-only bilateral scalp sensitivity. It is not imported into
  paired/group reporting and cannot run until a current terminal sidecar proves
  both F3 and F4.
- [`analysis/scalp_f3_f4_artifact_validation.py`](../analysis/scalp_f3_f4_artifact_validation.py) —
  clean-checkout validation of any published bilateral bytes; publication mode
  additionally requires and recomputes from the ignored private sidecar.

---

## A. Slow-oscillation detection

**Naji 2019:** *"EEG signals were filtered (zero-phase bandpass, .15–4 Hz). SO were detected from the
F3 and F4 electrodes based on a set of criteria for peak-to-peak amplitude, up-state amplitude, and
duration of down- and up-states (Dang-Vu et al. 2008)."*
The analyzed scalp derivations were F3/A2 and F4/A1.

| Aspect | Naji 2019 (Dang-Vu criteria) | Current 3B candidate detector | Match |
|---|---|---|---|
| Filter band | zero-phase bandpass **0.15–4 Hz** | notch + `sosfiltfilt` **0.15–4 Hz** | ⚠️ same SO band, additional notch |
| Where detected | **2 frontal scalp derivations**, F3/A2 and F4/A1 | conservative frontal iEEG ROI, per contact | ❌ modality/region/reference |
| Wave definition | full SO event (down-state **and** up-state) | complete negative/down and positive/up half-waves | ✅ structural definition aligned |
| Duration gate | duration of down- and up-states is named, but Naji does not report numerical limits | negative-to-positive 0.3–1.5 s, then positive-to-negative <1.0 s | ❌ implementation-specific; exact paper match cannot be claimed |
| Amplitude gate | peak-to-peak and up-state amplitude criteria are named, but Naji does not report their numerical cutoffs; the cited Dang-Vu scalp method uses µV criteria | up-state amp **≥ 75th pct** AND peak-to-peak **≥ 75th pct**, within channel | ❌ relative iEEG adaptation |
| Up-state amplitude | criterion named but not numerically specified in Naji | positive up-state and peak-to-peak amplitudes are percentile-gated separately | ⚠️ unverified adaptation |
| Staging | R&K visual scoring; continuous, undisturbed 3-min **Stage 2 or SWS** bins | GMM N2-like/N3-like proxy (HUP) or annotation-constrained labels (ds003848) | ❌ not R&K Stage 2/SWS |

**Why the deviations, and what they do:**

- **Sensor/reference transfer (F3/A2 and F4/A1 scalp → frontal iEEG)** is a major difference.
  The withdrawn numerical gap cannot be assigned uniquely to modality, region, reference, spatial
  scale, or channel aggregation.
- **Paper amplitude criteria → 75th-percentile threshold** is an explicit analysis choice, not a
  paper-derived rule. Scalp µV criteria do not transfer directly to depth-contact scale/reference,
  but other iEEG calibration strategies are possible. Both up-state and peak-to-peak amplitudes
  must clear their within-channel percentile; retention is therefore relative and can be smaller
  than 25%.
- **Complete half-wave timing is enforced in the current estimator**, but Naji does not publish the
  numerical duration limits needed to call it an exact match. Code-level morphology alignment also
  does not validate iEEG polarity or scalp-event homology.

---

## B. The statistics

This is the biggest conceptual difference and the one most relevant to the "validated" claim.

**What Naji actually reports and tests.** Naji reports the HR peak as **% above the stage mean** —
**12.09 ± 1.48%** in Stage 2 and **3.35 ± 1.01%** in SWS — as a descriptive measurement
(mean ± SEM across subjects). The paper does not report a null test for whether the event-locked
curve itself differs from chance. Its inferential claim is a **Pearson correlation between SO–HR
timing (ΔT ≈ 2.1–2.2 s) and texture-discrimination speed**. Timing was not associated with offline
performance improvement, and cross-visit correlations were poor.

**What the current descriptive estimator computes.** `subject_so_triggered`:

| Aspect | Naji 2019 | Current v9 adaptation | Match |
|---|---|---|---|
| R peaks | Pan–Tompkins detector followed by visual confirmation | NeuroKit automated detector; no visual validation | ❌ detector/validation transfer |
| RR series | RR → **4 Hz piecewise cubic spline** | **4 Hz PCHIP** inside continuous beat runs, followed by explicit long-gap re-masking | ⚠️ intentionally different |
| Window | 10 s centered on the trough (−5 to +5 s) | ±5 s (`HALF_WIN=5`) | ✅ |
| Effect statistic | HR peak after SO, reported as % above average Stage 2/SWS HR; conversion/aggregation order is not fully specified | minimum of the participant-average RR curve, converted once to HR and % above the stage baseline | ⚠️ related quantity, implementation order is repository-specific |
| SO→HR timing (ΔT) | trough → HR peak | `peak_lag_s` | ✅ same |
| **Significance test** | **none** for the descriptive coupling curve | no exposed event-locking p/z; a shared stage shift is retained as a diagnostic only | ✅ no inferential coupling claim |
| Aggregation | average the F3 and F4 **HR-maximum times** for ΔT; the paper does not fully specify the two-derivation magnitude aggregation | average contact RR curves before one participant magnitude; average contact-specific timings | ⚠️ timing concept aligned, magnitude order and sensors are not |
| Behavioural endpoint | **ΔT ↔ TDT speed** (Pearson) — their headline | **not tested** (no behaviour in iEEG) | ❌ we cannot do their actual claim |

**Three things this means:**

1. **The current output is descriptive.** The repository does not expose the
   old random-trigger z or p value as evidence. Its shared-shift statistic is a
   diagnostic because a whole-stage shift does not necessarily preserve local
   cardiac trends and clustered SO timing.
2. **The magnitude has a related definition, but is not an apples-to-apples
   replication.** Both use a post-trough curve extremum relative to a stage
   baseline, but interpolation, staging, SO criteria, sensor location,
   reference, and channel aggregation differ. The withdrawn numerical gap
   cannot be attributed uniquely to scalp-versus-iEEG region or spatial scale.
3. **We did not test Naji's actual finding.** Their result is a visit-specific Pearson association
   between SO→HR timing and perceptual speed. We have no behavioural task, so we cannot replicate
   that — we only measured the coupling itself.

### Bounded HUP138 F3/F4 scalp-only sensitivity

HUP138 is the only audited HUP participant with unique, geometry-matched F3 and
F4 labels. Its acquisition plan now includes both channels, but no current
terminal sidecar or bilateral result is implied by the plan alone. Once the
private sidecar is regenerated, the separate runner applies the locked 3B gates
to F3 and F4 independently and calls `subject_so_triggered` once with both event
sets in the RR domain (`minimum_channels=2`). It does **not** average two
unilateral scalar results: magnitude comes from the mean F3/F4 RR curve, while
`peak_lag_s` is the mean of the two electrode-specific RR-minimum/HR-maximum
lags.

N2 and N3 remain unavailable whenever there is no uninterrupted 180-s proxy
stage support. Pooled NREM is labeled a non-Naji sensitivity. All z/p fields
remain null; the shared stage shift is diagnostic only. The artifact is one
participant, scalp-only, and cannot support a paired iEEG, group, equivalence,
behavioural, or LC claim. Unknown online references, automated ECG without
visual adjudication, PCHIP rather than the paper's spline, repository-specific
duration and percentile gates, and trough-only candidate artifact QC also
prevent an exact-method claim.

The runner writes deterministic JSON/CSV/PNG/SVG plus a byte-exact terminal
manifest under `outputs/scalp_f3_f4_exploratory/`. It hard-stops on dirty source
or nonterminal input lineage. Offline validation requires current source and
profile hashes plus well-formed, manifest-consistent producing-runtime metadata,
without opportunistically reading private data. Publication validation also
requires that producing runtime and recomputes the endpoint and figure bytes
from the pinned private inputs. Source tests do not generate those artifacts:

```bash
.venv/bin/python analysis/test_scalp_f3_f4_exploratory.py
.venv/bin/python analysis/scalp_f3_f4_artifact_validation.py --artifact-validation offline
```

---

## What the legacy pipeline reported (withdrawn context)

| Cohort | N2 HR peak | N3 HR peak | Significant / n | N2 ≫ N3? |
|---|---|---|---|---|
| **HUP iEEG (n=23)** | +1.8% (t p=0.014) | +3.5% (t p=0.020) | 7/23, 6/20 | no (p=0.57) |
| **ds003848 legacy proxy (n=6; withdrawn)** | +0.6% (t p=0.88) | +1.0% (t p=0.28) | 0/5, 0/6 | no (p=0.63) |
| **Naji 2019 (scalp, healthy)** | **+12.09%** | **+3.35%** | — | yes (3.6×) |

The legacy write-up called HUP coupling weak/heterogeneous, RESPect null, and
3A absent. Those biological verdicts are withdrawn. The historical v8
endpoint-local grid recovered descriptive 3B estimates for 6/5/15 HUP
participants and 0/1/2 RESPect participants in N2/N3/pooled NREM,
respectively; inference remained disabled and RESPect remained below the
default cohort minimum. These are not regenerated v9 results. See
[QC_SENSITIVITY_RESULTS_2026-07.md](QC_SENSITIVITY_RESULTS_2026-07.md) and
[ISSUE_REGISTER_2026-07.md](ISSUE_REGISTER_2026-07.md).
