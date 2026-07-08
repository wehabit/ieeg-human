# Hierarchical coupling test (SO→spindle→ripple) — the missing spindle→ripple leg

The Staresina 2015 hierarchy is SO ⊃ spindle ⊃ ripple. Earlier analyses here tested only the two
SO→ legs; this adds the **spindle→ripple** leg (polarity-robust, cross-validated) on the
ripple-capable cohorts (the 204 Hz atlas cannot reach ripples).

| Cohort | SO→spindle | spindle→ripple | SO→ripple |
|---|---|---|---|
| Falach 1 kHz (n=15) | 0.047 [0.008,0.097] ✓ | 0.002 [−0.003,0.006] **n.s.** | 0.011 [−0.008,0.033] n.s. |
| Zurich 2 kHz (n=9) | 0.024 [0.014,0.036] ✓ | 0.004 [0.001,0.006] ✓ (tiny) | 0.016 [0.010,0.022] ✓ |

**Verdict: the full hierarchy is NOT supported by these data.** The spindle→ripple leg is null in
Falach and only marginally/tinily positive in Zurich. Therefore we drop the word "hierarchical" and
claim only what holds: **SO→spindle coupling** (robust, 3 cohorts, polarity-robust) and **SO→ripple**
(Zurich + HUP full nights).

**Caveat (method sensitivity):** spindle→ripple nesting is conventionally shown event-based (ripple
*rate* locked to spindle phase; Staresina 2015), which is more sensitive than the continuous
phase-amplitude method used here. Absence here is "not detected with this method," not proof of
absence. An event-based spindle→ripple test is the proper follow-up if the leg matters.

## Outputs
- `outputs/hierarchical_coupling/hierarchical_coupling.png` — 2 cohorts × 3 legs, aligned curves + CI.
