# 3B — is the null a locality artifact? SO globality and region tests (ds003848)

**Bottom line: partial support, honestly mixed.** The prediction that the coupling is carried by
**global** slow waves gets directional support in SWS (global SOs modulate HR ~15× more than local
ones, 5/5 subjects, though only marginally significant at n=5). The prediction that it is stronger on
**autonomic-adjacent cortex** (frontal/cingulate/insula) is **not** supported — those contacts are no
better than posterior/lateral ones. So it looks like SO *globality*, not cortical *location*, is what
our lateral single-channel 3B missed — but even global iEEG SOs stay far below Naji's scalp magnitudes,
so this neither fully rescues nor cleanly refutes the coupling.

## Why we ran this

Our lateral-iEEG 3B was weak/null ([3B_METHOD_COMPARISON.md](3B_METHOD_COMPARISON.md)). The coupling
Naji 2019 measured lives in large, global, frontally-maximal (K-complex-type) slow waves near the
central autonomic network; a lateral neocortical contact sees mostly *local* slow waves far from that
network. Two falsifiable predictions follow, both testable on ds003848 from the cache (all-channel SO
troughs) plus its Destrieux atlas labels — no re-streaming. Code:
[`analysis/region_global_3B.py`](../analysis/region_global_3B.py), reusing the corrected
`so_triggered` (stage-matched null) unchanged.

## Test 1 — globality gradient

Per-channel SO troughs are clustered into consensus events; each event's **globality** = fraction of
channels participating within ±150 ms. Events are binned local (<15%) / regional (15–40%) /
global (≥40%) and run through the same SO-triggered HR test.

**SWS (N3) — the predicted rising gradient appears:**

| globality | pooled HR peak (% above stage mean) | pooled z |
|---|---|---|
| local | +0.13% | +0.22 |
| regional | +0.97% | +0.97 |
| **global** | **+1.95%** | +0.86 |

Paired local→global: **global > local in 5/5 subjects**, mean diff +1.81% (paired t p = 0.10,
Wilcoxon p = 0.06 — the floor for n = 5). So the direction is consistent and sizeable but only
marginally significant at this n.

**Stage 2 (N2) — flat/absent:** pooled local +0.32% → regional +0.24% → global +0.37%, and global
events are rare in N2 (only one subject had a populated global bin). No gradient.

That N3-but-not-N2 pattern is physiologically coherent: globally-coordinated slow waves dominate SWS,
so that is where a globality effect should show.

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

## Interpretation

- **What our lateral 3B missed was SO globality, not lobe.** Restricting to globally-coordinated SOs
  recovers a modest, direction-consistent HR modulation in SWS (5/5 subjects); restricting to
  "autonomic" cortex does not. This makes sense — a global slow wave is global everywhere, so any
  participating contact captures it, and you isolate it by cross-channel *consensus*, not by picking a
  frontal contact.
- **But it does not reach scalp magnitudes.** Even global iEEG SOs peak at ~+2% (N3), versus Naji's
  +12% (S2) / +3.35% (SWS) on frontal scalp. A consensus of depth contacts still cannot reconstruct
  what scalp spatial-averaging integrates, so a fully powered iEEG replication would need dense frontal
  coverage and far more than 6 subjects.
- **Net for the "validated" claim:** unchanged. iEEG neither validates nor refutes the SO→HR coupling;
  this analysis shows *why* (locality dilutes it) and recovers a weak globality-dependent trace in SWS,
  but nothing approaching the published effect.

## Caveats

- **n = 6**, split by stage and globality/region — exploratory and directional, not powered.
- **Clinical coverage**, not autonomic-targeted: RESP0724 is almost entirely temporal; only RESP0749
  has insula contacts. Per-subject channel counts are in the JSON.
- Consensus globality is bounded by how many channels each subject has (35–93); a subject with fewer
  channels reaches "global" less often.
- "Autonomic-adjacent" uses cortical Destrieux labels; the deepest autonomic hubs (amygdala, brainstem)
  are not sampled by these montages.

Numbers: `outputs/region_global_3B/*.json`. Companion:
[3B_METHOD_COMPARISON.md](3B_METHOD_COMPARISON.md), [BRANCH_SUMMARY.md](BRANCH_SUMMARY.md).
