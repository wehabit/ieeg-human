# WITHDRAWN 3B locality/globality follow-up — superseded estimators

> **LEGACY / WITHDRAWN NUMBERS.** This follow-up predates the RR-domain v8 descriptive
> path and unvalidated-anatomy audit. Its whole-stage-shift inference is invalid under local
> nonstationarity, its command-line entry points are hard-stopped, and its numerical results are
> historical only. V8 recovers descriptive 3B estimates but not this follow-up or a valid
> event-locking p/z claim, so none of the old locality/globality conclusions is current.

**Withdrawn historical bottom line: partial, mixed, and none reached scalp magnitudes.** Three predictions,
three answers: (1) coupling is stronger for **global** slow waves — directional support in N3 (global
> local in 5/5 subjects, marginal at n=5); (2) it is stronger on **autonomic-adjacent cortex** — **not**
supported; (3) it is carried by **N2 K-complexes** — **underpowered and unrecovered** (isolated frontal
N2 events are too rare in iEEG). Across all three, even the best iEEG condition stays far below Naji's
scalp values, so this neither rescues nor refutes the coupling — it explains *why* lateral iEEG misses it.

**Two axes, kept separate (an earlier draft conflated them).** "Globally-synchronous" (spatial: how many
contacts share a wave) and "K-complex" (a stage/type: the large *isolated* slow wave that defines **N2**)
are different things. K-complexes are an **N2** hallmark tied to micro-arousals and phasic autonomic
bursts (de Zambotti 2016) — that is what Naji's **N2 > N3** effect points to. Dense, mutually-synchronous
slow-wave *trains* are an **N3** property. Test 1 (globality) turns out to measure the N3 axis; Test 3
(K-complexes) targets the N2 axis directly.

## Why we ran this

The legacy lateral-iEEG 3B was described as weak/null
([3B_METHOD_COMPARISON.md](3B_METHOD_COMPARISON.md)). The SO→HR coupling
Naji measured is a property of large slow waves, and a lateral neocortical contact sees mostly *local*
slow waves. Three falsifiable predictions follow, all testable on the Utrecht RESPect cohort (OpenNeuro
ds003848) from the cache (all-channel SO troughs) plus its Destrieux atlas labels — no re-streaming. Code:
[`analysis/region_global_3B.py`](../analysis/region_global_3B.py) (Tests 1–2) and
[`analysis/kcomplex_n2_3B.py`](../analysis/kcomplex_n2_3B.py) (Test 3). Both are withdrawn because
their stage-shift null does not preserve local HR trends or SO event-density clustering.

## Test 1 — globality gradient

Per-channel SO troughs are clustered into consensus events; each event's **globality** = fraction of
channels participating within ±150 ms. Events are binned local (<15%) / regional (15–40%) /
global (≥40%) and run through the same SO-triggered HR test.

**N3 — the predicted rising gradient appears:**

| globality | pooled HR peak (% above stage mean) | pooled z |
|---|---|---|
| local | +0.13% | +0.22 |
| regional | +0.97% | +0.97 |
| **global** | **+1.95%** | +0.86 |

Paired local→global: **global > local in 5/5 subjects**, mean diff +1.81% (paired t p = 0.10,
Wilcoxon p = 0.06 — the floor for n = 5). So the direction is consistent and sizeable but only
marginally significant at this n.

**Stage 2 (N2) — flat/absent:** pooled local +0.32% → regional +0.24% → global +0.37%, and global
events are rare in N2 (0–40 per subject vs 12–131 in N3). No gradient.

**Important — this metric measures the N3 axis, not K-complexes.** A cross-channel "global" event needs
many waves co-occurring, which happens in the **dense slow-wave trains of N3**, not among the **sparse,
isolated K-complexes of N2** (hence 0–40 global events in N2 vs 12–131 in N3). So Test 1 confirms that
*N3 (slow-wave-sleep) synchrony* couples to HR — it does **not** test Naji's N2/K-complex effect. That needs a
K-complex-specific detector (Test 3).

## Test 2 — region contrast (autonomic-adjacent vs posterior/lateral)

Contacts grouped by Destrieux label into **autonomic-adjacent** (frontal / cingulate / insula) vs
**posterior/lateral** (temporal / parietal / occipital); per-channel SO-triggered test averaged within
each group.

| stage | autonomic-adjacent (pooled z) | posterior/lateral (pooled z) |
|---|---|---|
| N2 | −0.50 | +0.75 |
| N3 | +0.48 | +1.05 |

**Not supported — the opposite, if anything.** Autonomic-adjacent > posterior in only 2/6 subjects
(N3, paired p = 0.30). The one subject with real **insula** coverage (RESP0749, 8 insula contacts —
the key cortical autonomic hub) shows autonomic z ≈ 0. So the coupling is not focal to autonomic
cortex; it is diffuse, present about equally wherever a global SO reaches.

## Test 3 — N2 K-complexes (Naji's actual effect)

Naji's headline is the **N2** effect (+12% vs +3.35% in N3), and N2's signature autonomic event is the
**K-complex** — a large, *isolated* slow wave on frontal/midline cortex. Test 1 could not see this (its
consensus metric captures N3 trains). So Test 3 targets it directly
([`kcomplex_n2_3B.py`](../analysis/kcomplex_n2_3B.py)): a **K-complex proxy** = a slow-wave event that is
(1) in **N2**, (2) on **frontal/cingulate** contacts, and (3) **isolated** (no other frontal event within
±2.5 s), versus three controls (N2 frontal *train*, N2 *posterior* isolated, N3 frontal isolated).

*Caveat: the cache holds trough times only, so this is an "isolated large N2 frontal slow-wave event"
proxy — it does not verify biphasic K-complex morphology (that needs re-downloading the raw BrainVision).*

| condition | pooled z | pooled % above N2 mean | usable subjects |
|---|---|---|---|
| **N2 frontal-isolated (K-complex proxy)** | **+0.09** | **+0.88%** | 3 / 6 |
| N2 frontal-train | −1.49 | −0.17% | 5 |
| N2 posterior-isolated | −0.84 | −0.37% | 2 |
| N3 frontal-isolated | +1.08 | +1.38% | 2 |

**Underpowered, and Naji's N2 effect is not recovered.** Isolated K-complex-like events are **rare** in
iEEG — only 1–12% of frontal events, so just 3/6 subjects reached the 30-event minimum. The K-complex
proxy's z ≈ 0. But the *direction* among N2 conditions is as predicted: isolated frontal (+0.88%) >
train (−0.17%) > posterior (−0.37%) — isolated frontal N2 events couple more than N2 trains or posterior
events, just far too weakly and sparsely to reach significance. The most likely reason: K-complexes are
frontal-**midline** (Fz/Cz) events, and these "frontal" iEEG contacts are *lateral* frontal cortex, not
midline — so even a K-complex-targeted test on depth electrodes barely captures them.

## Interpretation

- **What our lateral 3B missed was SO globality/scale, not lobe.** Restricting to globally-coordinated
  SOs recovers a modest, direction-consistent HR modulation in N3 (5/5 subjects); restricting to
  "autonomic" cortex does not; and isolating N2 K-complexes is too sparse to test on iEEG. A global slow
  wave is global everywhere, so any participating contact captures it via cross-channel *consensus* — you
  don't isolate the coupling by picking a frontal contact.
- **Both N2 and N3 mechanisms stay far below scalp.** Global N3 SOs peak at ~+2%, the N2 K-complex proxy
  at ~+0.9%, versus Naji's +12% (N2) / +3.35% (N3) on frontal scalp. A depth-electrode consensus cannot
  reconstruct what scalp spatial-averaging integrates — and the N2/K-complex effect additionally needs
  frontal-midline coverage that these clinical montages lack.
- **Net for the "validated" claim:** unchanged. iEEG neither validates nor refutes the SO→HR coupling;
  these analyses show *why* (locality/scale dilutes it, and midline K-complexes are barely sampled) and
  recover only weak, direction-consistent traces, nothing approaching the published effect.

## Caveats

- **n = 6**, split by stage and globality/region — exploratory and directional, not powered.
- **Clinical coverage**, not autonomic-targeted: RESP0724 is almost entirely temporal; only RESP0749
  has insula contacts. Per-subject channel counts are in the JSON.
- Consensus globality is bounded by how many channels each subject has (35–93); a subject with fewer
  channels reaches "global" less often.
- "Autonomic-adjacent" uses cortical Destrieux labels; the deepest autonomic hubs (amygdala, brainstem)
  are not sampled by these montages.
- **K-complex proxy is trough-time-based**, not morphology-verified, and "frontal" here is lateral
  frontal cortex, not the frontal-midline (Fz/Cz) where K-complexes are maximal — a morphology-verified,
  midline-targeted test would need re-downloading the raw data.

Numbers: `outputs/region_global_3B/*.json`, `outputs/kcomplex_n2_3B/*.json`. Companion:
[3B_METHOD_COMPARISON.md](3B_METHOD_COMPARISON.md), [BRANCH_SUMMARY.md](BRANCH_SUMMARY.md).
