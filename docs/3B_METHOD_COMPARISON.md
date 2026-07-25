# WITHDRAWN 3B comparison — superseded implementation vs Naji 2019

> **LEGACY / WITHDRAWN NUMBERS.** This comparison predates the RR-domain v8 descriptive
> estimator. Whole-stage shifts are now diagnostic only and no 3B event-locking p/z claim is
> available. The direct legacy entry point is hard-stopped. V8 recovers subject-level N2-like,
> N3-like, and pooled estimates, but inference remains explicitly disabled.

**Purpose.** A step-by-step comparison of how we implemented the slow-oscillation → heart-rate
coupling test against the method it follows. Written because the coupling was cited to the wrong
paper and described as "validated in our own data" — this documents what we actually did and how it
differs from the source.

## Which paper is the source

The empirical methods paper is **Naji, Krishnan, McDevitt, Bazhenov & Mednick (2019)**, *"Timing
between Cortical Slow Oscillations and Heart Rate Bursts during Sleep Predicts Temporal Processing
Speed, but Not Offline Consolidation"*, J Cogn Neurosci 31(10):1484–1490
([DOI](https://doi.org/10.1162/jocn_a_01432) · [lab PDF](http://sleepandcognitionlab.org/wp-content/uploads/2019/07/Timing-between-Cortical-Slow-Oscillations-and-Heart-Rate-Bursts-during-Sleep-Predicts-Temporal-Processing-Speed-but-Not-Offline-Consolidation.pdf)).

The **PNAS 2022** paper (Chen, Zhang, Thayer & Mednick, *"Understanding the roles of central and
autonomic activity during sleep…"*, 119(44)) is a **review** — its Data Availability statement reads
verbatim *"There are no data underlying this work."* It proposes a framework and reviews others'
findings; it does not itself demonstrate the SO→HR relationship. Mednick is a coauthor on Naji 2019,
so "Mednick and colleagues" is fine when pointed at Naji 2019.

**There is no public code repository for Naji 2019** (no GitHub/OSF). The comparison below is against
the paper's Methods (which cite Dang-Vu et al. 2008 for the SO-detection criteria). To obtain their
actual code, contact the Mednick lab (UC Irvine Sleep & Cognition Lab).

## Our code

- [`analysis/event_3B_mednick.py`](../analysis/event_3B_mednick.py) — retained helper functions only;
  its direct command-line analysis is withdrawn,
  `rr_to_hr_4hz()`, `so_triggered()` (the SO-triggered HR average + stage-matched null).
- [`analysis/event_3B_cached.py`](../analysis/event_3B_cached.py) — cohort runner (HUP + ds003848).
- [`analysis/cohort_stages_3ABD.py`](../analysis/cohort_stages_3ABD.py) `test_3B` — the earlier,
  superseded version (whole-night null; now fixed).
- [`analysis/test_3B_null.py`](../analysis/test_3B_null.py) — regression test for the null fix.

---

## A. Slow-oscillation detection

**Naji 2019:** *"EEG signals were filtered (zero-phase bandpass, .15–4 Hz). SO were detected from the
F3 and F4 electrodes based on a set of criteria for peak-to-peak amplitude, up-state amplitude, and
duration of down- and up-states (Dang-Vu et al. 2008)."*

| Aspect | Naji 2019 (Dang-Vu criteria) | Our `detect_so_halfwaves` | Match |
|---|---|---|---|
| Filter band | zero-phase bandpass **0.15–4 Hz** | notch + `sosfiltfilt` **0.15–4 Hz** | ✅ identical |
| Where detected | **2 frontal scalp** electrodes (F3, F4) | **every lateral neocortical iEEG contact** (54–93/subject), per channel | ❌ modality/region |
| Wave definition | full SO event (down-state **and** up-state) | **negative half-wave** only (zero-crossing → argmin trough) | ⚠️ simplified |
| Duration gate | duration of **down- and up-states** | **down-state only**, 0.3–1.0 s | ⚠️ up-state gate dropped |
| Amplitude gate | **absolute µV** (peak-to-peak + up-state amplitude) | negative-peak amp **≥ 75th pct** AND peak-to-peak **≥ 75th pct**, within channel | ❌ relative, not absolute |
| Up-state amplitude | explicit criterion | not gated separately; p2p uses max of next ≤1 s as a proxy up-state | ⚠️ approximated |
| Staging | R&K visual scoring, 3-min stable N2 / N3 bins (Naji's paper labels N3 "SWS") | GMM-on-slow-wave proxy (HUP) or EMG/EOG rule-based (ds003848) | ❌ not R&K |

**Why the deviations, and what they do:**

- **Region (F3/F4 scalp → lateral iEEG)** is the dominant difference and is unavoidable — a
  depth-electrode patient has no scalp montage. Scalp F3/F4 see **large, globally synchronous** SOs;
  a single lateral intracranial contact sees **local** ones. This is the main reason our HR effect is
  ~7× smaller than Naji's (+1.8% vs +12%).
- **Absolute µV → 75th-percentile threshold** is forced: intracranial amplitudes do not map to scalp
  µV, so a fixed Dang-Vu µV bar is meaningless on depth contacts. Consequence: our detector always
  keeps the top 25% of waves per channel, so SO *count* is set by the percentile rather than by an
  absolute amplitude.
- **Half-wave-only + dropped up-state duration gate** makes ours slightly more permissive than
  Dang-Vu's full-wave detector.

---

## B. The statistics

This is the biggest conceptual difference and the one most relevant to the "validated" claim.

**What Naji actually tests.** Naji reports the HR peak as **% above the stage mean** — **12.09 ± 1.48%**
(N2), **3.35 ± 1.01%** (N3) — as a **descriptive** measurement (mean ± SEM across subjects). They
**do not test the coupling against a null**; the coupling is taken as given (their Fig 1). Their actual
*statistical* claim is a **Pearson correlation between the SO→HR timing (ΔT ≈ 2.1–2.2 s) and
texture-discrimination speed** — i.e. the *timing* predicts behaviour, and even that was **negative for
memory consolidation** (hence the paper's title).

**What we test.** `so_triggered`:

| Aspect | Naji 2019 | Our `so_triggered` | Match |
|---|---|---|---|
| HR series | RR → **4 Hz cubic spline** | identical (`rr_to_hr_4hz`) | ✅ |
| Window | 10 s around trough; HR in **±5 s** post-trough | ±5 s (`HALF_WIN=5`) | ✅ |
| Effect statistic | **peak** of mean HR curve, % above stage mean | `curve[post].max()`, % above stage mean | ✅ same statistic |
| SO→HR timing (ΔT) | trough → HR peak | `peak_lag_s` | ✅ same |
| **Significance test** | **none** (descriptive, mean ± SEM) | **stage-matched random-trigger null**, 200 surrogates, **z-score** | ❌ we added a null |
| Aggregation | average F3/F4, then across subjects | average z across all channels, then **cohort t-test** on per-subject z | ⚠️ many correlated channels |
| Behavioural endpoint | **ΔT ↔ TDT speed** (Pearson) — their headline | **not tested** (no behaviour in iEEG) | ❌ we cannot do their actual claim |

**Three things this means:**

1. **We are stricter than Naji.** They never asked "is the coupling above chance?" — they reported the
   % descriptively. We built a stage-matched surrogate null and z-scored it. So our "significant in
   7/23" is *our* bar; by Naji's own descriptive standard you would simply report the +1.8% / +3.5%.
   The added null is more rigorous, but "not significant" partly reflects a test Naji never applied.
2. **The magnitude comparison is fair.** Both take the peak of the averaged HR curve as % above stage
   mean, so +1.8% (ours, Stage 2) vs +12.09% (theirs) is apples-to-apples — the ~7× gap is real and
   driven by region (local iEEG vs global frontal scalp). Because both the observed effect and the
   surrogates take `max()` over the post-trough window, that max's upward bias is symmetric and
   cancels in the z (so the z is fair, even though the raw +1.8%, like Naji's peak, is a max-based
   number).
3. **We did not test Naji's actual finding.** Their result is that SO→HR *timing predicts perceptual
   speed*. We have no behavioural task, so we cannot replicate that — we only measured the coupling
   itself.

---

## What the legacy pipeline reported (withdrawn context)

| Cohort | N2 HR peak | N3 HR peak | Significant / n | N2 ≫ N3? |
|---|---|---|---|---|
| **HUP iEEG (n=23)** | +1.8% (t p=0.014) | +3.5% (t p=0.020) | 7/23, 6/20 | no (p=0.57) |
| **ds003848 legacy proxy (n=6; withdrawn)** | +0.6% (t p=0.88) | +1.0% (t p=0.28) | 0/5, 0/6 | no (p=0.63) |
| **Naji 2019 (scalp, healthy)** | **+12.09%** | **+3.35%** | — | yes (3.6×) |

The legacy write-up called HUP coupling weak/heterogeneous, RESPect null, and 3A absent. Those
biological verdicts are withdrawn. In the corrected runs, neither cohort has an estimable 3B
endpoint and neither has an estimable 3A endpoint, which does not establish absence. See
[BRANCH_SUMMARY.md](BRANCH_SUMMARY.md) and
[LC_INFRASLOW_3ABD_SUMMARY.md](LC_INFRASLOW_3ABD_SUMMARY.md).
