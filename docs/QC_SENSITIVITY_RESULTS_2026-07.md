# LC-proxy QC sensitivity results — v8

## Bottom line

The v8 rerun changes the **availability** conclusion, but it does not establish a human
locus-coeruleus (LC) biomarker.

This document describes the exact hash-pinned v8 baseline. Subsequent raw-channel selector and
numerical flat-line fixes change the cache-builder digest, so the current code intentionally
rejects these caches as reusable current-build artifacts. A complete HUP cache/grid rebuild is
still required to measure how those fixes change the cohort-wide results.

- The former fixed-80% pipeline stacked several locally introduced coverage and count gates. It
  made every HUP endpoint unavailable even though the cited papers do not specify a 70%, 75%, 80%,
  or 90% recording-coverage cutoff.
- A locked, outcome-blind endpoint-local sensitivity profile recovers subject-level estimates:
  HUP has 19/14/15 usable 3A spectrum/coherence/cross-correlation records, 6/5/15 usable 3B
  N2-like/N3-like/pooled records, and 7 descriptive 3D records. RESPect has 3 usable 3A records,
  one N3 3B record, and two exploratory pooled-NREM 3B records.
- RESPect still does not reach the default cohort minimum of five. HUP reaches the nominal
  availability minimum, but the saved artifacts contain no valid cohort inferential test for 3A,
  explicitly disable 3B inference, and explicitly disable 3D inference.
- The most consistent recovered observation is a small, descriptive cardiac acceleration about
  2–3 seconds after cortical slow-oscillation troughs. That is compatible with downstream
  CNS–autonomic timing; it is not LC-specific.
- 3A is heterogeneous and does not reproduce a clean group-level approximately five-second
  HR-leading-sigma pattern. 3D preferred phases are dispersed and are descriptive only.

Current contract:

```text
analysis_version = 2026-07-qc-sensitivity-v8
cache_schema_version = 2026-07-neutral-per-contact-gap-aware-source-pin-v8
```

## What SWA means

**SWA** means **slow-wave activity**, here broadband EEG/iEEG power from approximately
**0.5–4 Hz**. It is commonly associated with deeper NREM sleep. It is not:

- a discrete slow oscillation (SO);
- the same thing as SO–spindle coupling; or
- a direct or validated LC measurement.

## Why 80% is not a paper-derived qualification rule

Lecci used consecutive artifact-free 30-second NREM epochs, bouts lasting at least 120 seconds,
and the first 210 minutes from sleep onset. Naji used visually scored uninterrupted three-minute
sleep bins and visually checked R peaks. Staresina and Helfrich used artifact-free NREM and
event-amplitude/duration rules. None of those methods specifies a 70%, 74%, 75%, 80%, or 90%
whole-record coverage cutoff, a minimum of three iEEG contacts, 30 SOs per contact, 1,200 valid
seconds per contact, 200 paired events, or a five-participant availability gate.

The “75%” in the SO/spindle literature is an **event-amplitude percentile**, not a percentage of
the recording that must survive QC. The former 80% rules entered this repository after the
biological analyses and were neither prospective nor preregistered.

It is therefore not defensible to say that 79% is scientifically invalid while 80% is valid.
It is also not defensible to select 70% because it produces a preferred result. v8 handles this by:

1. retaining reversible numerator/denominator and observation masks in neutral caches;
2. defining profiles before inspecting the direction of the LC-proxy results;
3. changing one coverage/count axis at a time; and
4. reporting endpoint availability and estimates across the locked grid.

Primary-method evidence:

- [Lecci et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5298853/) defines human SWA
  0.5–4 Hz, sigma 10–15 Hz, visually scored artifact-free NREM, bouts ≥120 s, and the first
  210 minutes after sleep onset. Its human observations do not directly measure LC.
- [Naji et al. 2019](https://escholarship.org/uc/item/5393b9zk) uses visually scored
  uninterrupted 3-minute bins, F3/F4, a 4-Hz RR tachogram, and visually checked R peaks; its main
  inference relates SO–HR timing to behavior, which this repository does not measure.
- [Dang-Vu et al. 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2567508/) supplies the cited scalp
  SO morphology and absolute voltage rules; the latter cannot be transferred literally to iEEG.
- [Staresina et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4625581/) and
  [Helfrich et al. 2018](https://helfrich-lab.com/uploads/papers/2018a_Helfrich_Neuron.pdf)
  motivate the fixed SO/spindle bands and event-amplitude/duration choices, not a recording-
  coverage percentage or an LC measure.

## Profiles and calibration

| Profile | Purpose | Key behavior |
|---|---|---|
| `audit80` | Historical post-audit reference | Fixed contact sets plus the former stacked 80% coverage rules; not preregistered and not paper-derived |
| `overlap11_endpoint_local` | Outcome-blind exact-support sensitivity | Uses an overlap-connected component and median polish; resolves power, staging, cardiac, and event support at their own endpoints; keeps paper method constants fixed |

`overlap11_endpoint_local` is not a “70% profile.” Relative to `audit80`, it changes the
aggregation architecture and several support/count gates at once. The large availability change
between those two profiles must not be attributed to one percentage.

The “11” refers only to iEEG staging spectra: at least 11 complete 4-second Hann/Welch windows,
with 2-second overlap, out of 14 possible windows in a 30-second epoch. An outcome-blind
reconstruction calibration used 42,075 records from all six RESPect participants. The first exact
support stratum meeting the locked engineering targets was 11 windows; the first pooled
at-or-above stratum passed at 9. At exactly 11 windows:

- `n = 4,430` reconstruction records from all six participants;
- SWA Spearman correlation = 0.982;
- delta-ratio Spearman correlation = 0.972;
- median/p95 absolute log-SWA error = 0.079/0.293; and
- median/p95 absolute delta-ratio error = 0.0107/0.0455.

These are engineering measurement-stability targets, not values supplied by a paper. The
calibration does not validate the artifact detector, sleep stages, or LC construct. It used
RESPect only; applying 11 windows to HUP is explicitly an **unvalidated transport sensitivity**.
EMG/EOG still require 14 complete windows in the endpoint-local base profile because the iEEG
calibration does not validate auxiliary channels.

The calibration artifact is byte-pinned and records the exact calibration-code hashes, runtime,
cache manifest, and per-cache hashes. Its current SHA-256 is
`c6d4f1deccd3e0c7a03875185b2378116b1e33aab78323db91a39dc0b828af1f`.

## Endpoint availability

Counts below mean “a numerical subject-level estimate was available under this profile.” They do
not mean that the hypothesis was supported.

| Cohort/profile | 3A spectrum | 3A coherence | 3A cross-correlation | 3B N2 | 3B N3 | 3B pooled NREM | 3D descriptive |
|---|---:|---:|---:|---:|---:|---:|---:|
| RESPect `audit80` | 0 | 0 | 0 | 0 | 1 | 2 | N/A |
| RESPect `overlap11_endpoint_local` | 3 | 3 | 3 | 0 | 1 | 2 | N/A |
| HUP `audit80` | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| HUP `overlap11_endpoint_local` | 19 | 14 | 15 | 6 | 5 | 15 | 7 |

HUP requested 25 participants. Twenty-four were analyzed and HUP116 was an explicit structured
skip because it had no cortical-contact candidates. All 24 HUP records had enough neutral fields
to evaluate the 3D support rules even when the descriptive effect failed those rules.

## Coverage sensitivity: 70% versus 75% versus 80%

The coverage grid changes one axis at a time around `audit80`.

- For HUP, changing any single 80% coverage threshold to 70%, 75%, or 90% recovers **zero**
  endpoints because other fixed-80/count blockers remain. This proves that the old
  disqualification was caused by a stack, not a single scientifically meaningful cliff.
- For RESPect, a 70% or 75% **power-bin contact fraction** admits RESP0521 for all three 3A
  endpoints; 80% and 90% do not.
- A 70% or 75% **aggregate-power coverage** admits spectra for RESP0521 and RESP0699 and
  coherence for RESP0521; 80% and 90% admit none.
- A 25% per-second clean-sample rule admits RESP0521, whereas 50%, 75%, and 90% do not.
- The other one-at-a-time coverage changes do not alter RESPect endpoint counts.

Thus 70% and 75% are often identical, but not universally. For example, HUP150 has approximately
72.4% historical aggregate coverage and passes a 70% threshold but not 75%; HUP177 has
approximately 74.6% and behaves the same way. Neither fact alone establishes endpoint validity:
the endpoint-specific observation graph, stage support, cardiac support, contacts, and event counts
still matter.

## Staging and event-count sensitivity

RESPect availability is invariant across 1–14 valid iEEG Welch windows:
3A = 3/3/3, 3B = 0/1/2, and 3D is outside scope. HUP is more sensitive:

| Valid iEEG windows | HUP 3A spectrum/coh/xcorr | HUP 3B N2/N3/pooled | HUP 3D |
|---:|---:|---:|---:|
| 1 | 20/16/16 | 5/5/17 | 6 |
| 4 | 20/16/16 | 7/7/17 | 6 |
| 7 | 18/14/14 | 5/5/15 | 5 |
| 9 | 19/14/15 | 8/8/15 | 7 |
| 11 | 19/14/15 | 6/5/15 | 7 |
| 12 | 17/13/14 | 2/2/15 | 7 |
| 14 | 11/7/11 | 3/2/11 | 7 |

The HUP stage-specific counts are nonmonotonic because the support rule changes the distribution
used to fit proxy stages; it is not merely removing a nested subset. This is model-specification
sensitivity and another reason not to treat HUP N2-like/N3-like labels as expert stages.

The auxiliary 1–14-window grid is an expected endpoint-level negative control:

- HUP has no EMG/EOG input, so auxiliary support is not applicable.
- In RESPect, the threshold changes finite EMG/EOG intermediates and changes some proxy labels,
  but author-constrained final labels remain identical. v8 now serializes finite auxiliary counts,
  proxy-label counts, final-label counts, and their disagreement so this invariance is auditable.

For 3B, HUP N2/N3 counts generally remain 5–6 and pooled remains 15 across the locked
SO-count/contact/amplitude-percentile grid. RESPect N3 disappears at 50 SOs per contact or at the
90th amplitude percentile. These are fragile subject-availability changes, not inferential
results.

For 3D, the largest availability driver is minimum contacts:

| Minimum contacts | HUP descriptive records |
|---:|---:|
| 1 | 13 |
| 2 | 10 |
| 3 (base) | 7 |
| 4 | 5 |

Raising events per contact to 100 or total paired events to 400 reduces 7 to 6. Changing the
valid-NREM duration from 120 through 1,200 seconds does not change the base count. No setting
enables inference.

## What the recovered estimates suggest

### 3A — infraslow sigma/cardiac candidate

RESPect has three subject spectra. Two have accepted algorithmic peaks, at 0.0169 and 0.0221 Hz,
bracketing Lecci's approximately 0.019-Hz report. Only one of three fixed-0.02-Hz coherence values
exceeds its individual analytic threshold. The Lecci-direction cross-correlation maxima are at
0, 0, and 5 seconds, with median `r = 0.303`, but `n = 3` is below the cohort minimum.

HUP has 19 spectra and 9 accepted peaks. Their median is 0.0166 Hz, but their range is broad
(0.0104–0.0325 Hz). Four of 14 fixed-0.02-Hz coherence values exceed their individual analytic
thresholds. Fourteen of 15 Lecci-direction correlations are positive, but the median is only
`r = 0.068`; lags span 0–15 seconds and are concentrated at zero rather than around five seconds.

**Interpretation:** there are individual infraslow sigma features, but the current artifacts do
not contain a valid cohort test and do not show a clean, robust replication of the proposed
human timing signature. Even a replication would be LC-motivated, not LC-specific.

### 3B — SO-linked cardiac timing

HUP descriptive local HR changes are:

| State | Available participants | Positive direction | Median local HR change | Median mean-channel latency |
|---|---:|---:|---:|---:|
| N2-like | 6 | 6/6 | +0.402% | about 2.7 s |
| N3-like | 5 | 5/5 | +0.529% | 2.5 s |
| Pooled NREM, exploratory | 15 | 14/15 | +0.245% | 2.0 s |

RESPect contributes one N3 estimate (`+0.455%`, 2.35 s) and two pooled estimates
(`+0.891%` and `+1.613%`, approximately 1.9 and 2.7 s).

Every saved 3B estimate has `p_upper = null`, `z = null`, and an explicit
`inference_status` saying that whole-stage shifts do not preserve local nonstationary HR trends or
event-density clustering. Low values in fields explicitly named `*_diagnostic` are not valid
p-values.

**Interpretation:** the direction and approximately 2–3-second timing are the clearest recovered
physiological convergence. They support describing SO-linked autonomic timing, not an inferential
Naji replication, behavioral replication, or LC proxy validation.

### 3D — SO–spindle nesting

Seven HUP participants pass the descriptive base support rules. Participant vector lengths range
from 0.050 to 0.292, while preferred phases range from about −43° to 161°. A direct, unvalidated
equal-participant arithmetic mean of the stored vectors has `R = 0.087` and phase approximately
0.08°, illustrating weak common direction; it is not a stored inferential result.

Every record sets `inference_enabled = false` and `inferential_p_value = null`. The nearest-event
pairing step is known to induce common phase under independence unless a complete-train
time-shift/block null repeats detection, pairing, and aggregation.

**Interpretation:** 3D is a generic sleep-coordination/QC measure with dispersed descriptive
phase, not an LC proxy and not a valid positive or null cohort result.

## Did the fixes change the results?

Yes, materially:

- **Before v8:** the documentation said all HUP and nearly all RESPect endpoints were unavailable.
- **After v8:** endpoint-local materialization recovers many HUP estimates and three RESPect 3A
  estimates. The previous all-unavailable conclusion was partly an engineering artifact.
- **What did not change:** there is still no proof of LC specificity, no valid 3B or 3D inference,
  no expert validation of HUP stages/anatomy, and no clean cohort-level 3A replication claim.

The defensible conclusion is therefore not “the data are all disqualified,” and it is not “the LC
proxy is proven.” It is:

> Subject-level LC-motivated downstream physiology is measurable in these data after endpoint-local
> QC. SO-linked cardiac acceleration is descriptively consistent, while the infraslow sigma and
> SO–spindle results are heterogeneous. LC specificity and the main inferential claims remain
> unvalidated.

## Newly confirmed implementation/release defects and fixes

| ID | Concern | Evidence | Resolution |
|---|---|---|---|
| A81 | The 80%/count stack was treated as paper-derived and globally disqualified endpoints | Paper-method comparison and one-axis grid | Preserve `audit80` only as historical; add locked endpoint-local and one-axis sensitivities |
| A82 | “Pooled NREM” matched only the literal `NREM` label, excluding N2/N3 and their transitions | Synthetic N2→N3→NREM regression | Use the union of NREM/N2/N3 while preserving continuity |
| A83 | Configured minimum finite RR support was not enforced | Adversarial profile regression | Enforce whole-record and stage-local finite-RR requirements |
| A84 | Skip caches could require pickled `None`, and JSON could contain nonstandard NaN/Infinity | `allow_pickle=False` and strict-JSON regressions | Store typed scalar/string fields and fail strict serialization |
| A85 | Sigma and SWA could be aggregated from different contact/time support | Joint-support regression | Use one joint observation graph/contact set |
| A86 | The overlap component selector could choose too few contacts despite another valid component | Disconnected-component regression | Select only components meeting minimum contacts; deterministic tie break |
| A87 | The 50%-clean per-second power decision was irreversible | 25% versus 50% reconstruction regression | Cache per-contact power numerators and clean-sample denominators |
| A88 | A brief gap erased an entire staging/auxiliary epoch | Planted-gap Welch-window regressions | Cache complete 4-second window powers/support and select support offline |
| A89 | 3B used symmetric 0.3–1.0-second half-waves rather than the cited Dang-Vu timing | Direct method comparison | Use negative half 0.3–1.5 s and following positive half <1 s |
| A90 | Cached HUP data could not reconstruct the 3D detector | Required-field regression | Cache 20-Hz RMS, SO phase, valid masks, and pre-threshold SO candidates |
| A91 | Finite placeholder arrays made missing 3D acquisition look like 100% support | All-unmeasured contact regression | Use sample-weighted measured acquisition and event-valid fractions |
| A92 | The cache-source digest omitted `spectral_gapped.py` | Digest membership audit | Include it in the cache-producer hash |
| A93 | Calibration trusted incomplete lineage and was not immutably pinned | Tampered-cache/artifact tests | Verify manifest/cache bytes and embedded identity; pin artifact and purpose-specific source hashes |
| A94 | Median polish could return a finite nonconverged estimate | Forced-low-iteration regression | Fail closed and suppress downstream stages/endpoints |
| A95 | New profile/overlap/3D/calibration regressions were absent from CI | Workflow audit | Add every new suite to CI |
| A96 | Grid outputs did not fully validate pipeline/runtime/exact hash keys | Adversarial manifest tests | Validate and record pipeline, runtime, manifest, cache, profile, and source hashes |
| A97 | 3A peak/control diagnostics could remain populated when profile support failed | Below-support regression | Withhold spectra, peaks, controls, coherence, and cross-correlation |
| A98 | Transient portal failures could be retried without preserving the failed-run trail | Real HUP retry audit | Add a finalizer that verifies retry bytes and records recovery provenance |
| A99 | Release validation demanded obsolete v7 directories masquerade as current | Audit-script failure on quarantined outputs | Validate v8 QC artifacts and require explicit legacy markers for old directories |
| A100 | An invariant auxiliary grid looked indistinguishable from an inactive code path | Real-cache intermediate audit | Serialize finite auxiliary, proxy-stage, final-stage, and disagreement diagnostics |
| A101 | The real-cache coherence test read obsolete raw stage aggregates and tested no v8 subject | `--require-real-cache` failed while endpoint-local records existed | Materialize the exact v8 3A profile; HUP133/HUP139 real gap geometries passed calibrated FPR checks on the frozen run; report later stale cache lineage as pending by default while keeping publication mode fail-closed |

## Remaining real blockers

1. None of 3A/3B/3D is a validated human LC measurement.
2. HUP proxy stages need expert PSG validation; HUP contacts need coordinates, tissue/ROI,
   pathology/SOZ, and clinical review.
3. RESPect/HUP sensors are not validated homologues of the cited scalp sources.
4. ECG, SO, spindle, artifact, and gap-boundary detectors need blinded real-data validation.
5. HUP timing is a sparse-probe high-delta interval, not verified sleep onset.
6. The Lecci signal construction remains an approximation to the FieldTrip/Morlet implementation.
7. 3B needs a local-trend/dependence-preserving null; its behavioral endpoint is absent.
8. 3D needs a complete-train pairing-aware null and matched-contact stage contrast.
9. Intracranial polarity/reference is not calibrated to scalp downstate definitions.
10. Respiration/apnea, medication, arousal, seizures/postictal physiology, and other common drivers
    remain incompletely controlled.

## Reproducible artifacts

- `analysis/qc_profiles_v1.json` — immutable profile definitions and grids
- `outputs/qc_calibration/staging_window_calibration.json` — outcome-blind calibration
- `outputs/qc_grid_public/coverage_oat_v1/` — every 70/75/80/90 profile,
  cohort count, and participant endpoint-availability record
- `outputs/qc_grid_public/staging_window_support_v1/` — compact 1–14
  iEEG-window evidence
- `outputs/qc_grid_public/auxiliary_window_support_v1/` — compact 1–14
  EMG/EOG-window control evidence
- `outputs/qc_grid_public/event_count_oat_v1/` — compact
  event/contact/cohort-count sensitivity evidence
- `outputs/qc_grid_public/locked/overlap11_endpoint_local__hup.json` — complete
  frozen locked-profile subject records required by the paired exact-match
  analysis

The full numeric grids are deterministic ignored products under
`outputs/qc_grid/`; regenerate them before the compact evidence with
`analysis/compact_qc_grid_artifacts.py`. See
[ARTIFACT_POLICY.md](ARTIFACT_POLICY.md). All committed JSON is strict, all
cache/result bytes are hash-pinned, and the old v7 endpoint folders remain
explicitly quarantined as legacy outputs.
