"""Aggregate the staged cohort with the quality filters applied, and report honestly.

Filters (declared before looking at the result):
  * a stage must hold >= MIN_EP epochs, and
  * the minority stage must be >= MIN_BAL of NREM epochs
    -- a 612/4 "split" is a GMM outlier cluster, not a sleep stage.
Pooled 3A is additionally tested as a binomial: under the null, 5% of subjects should exceed
their own alpha=0.05 coherence threshold.
"""
import glob, json, os
import numpy as np
from scipy import stats

OUT = "/Users/paris/Documents/iEEG/outputs/cohort_stages_3ABD"
MIN_EP = 50
MIN_BAL = 0.15
FAST_FSP = 12.0

rows = []
for f in sorted(glob.glob(os.path.join(OUT, "*.json"))):
    d = json.load(open(f))
    if d.get("status") != "ok":
        continue
    n2, n3 = d.get("n_N2", 0), d.get("n_N3", 0)
    tot = max(n2 + n3, 1)
    rows.append(dict(
        sub=d["subject"].replace("_phaseII", ""), fsp=d.get("fast_spindle_peak_hz"),
        n2=n2, n3=n3, bal=min(n2, n3) / tot, sep=d.get("gmm_separation"),
        d_n2=d["N2"].get("3D_mi_z"), d_n3=d["N3"].get("3D_mi_z"),
        b_n2=d["N2"].get("3B_z"), b_n3=d["N3"].get("3B_z"),
        a_n2=(d["N2"].get("3A") or {}).get("at"), a_n3=(d["N3"].get("3A") or {}).get("at"),
        pool_at=(d.get("NREM_pooled", {}).get("3A") or {}).get("at"),
        pool_crit=(d.get("NREM_pooled", {}).get("3A") or {}).get("crit"),
        pool_K=(d.get("NREM_pooled", {}).get("3A") or {}).get("K"),
    ))

print(f"subjects with results: {len(rows)}\n")
print(f"{'sub':8s} {'FSP':>6s} {'N2':>5s} {'N3':>5s} {'bal':>5s} | {'3D N2':>8s} {'3D N3':>8s} | "
      f"{'3B N2':>6s} {'3B N3':>6s} | {'3Apool':>7s} {'crit':>6s} {'':3s}")
for r in rows:
    ok = (r["n2"] >= MIN_EP and r["n3"] >= MIN_EP and r["bal"] >= MIN_BAL)
    sig = ""
    if r["pool_at"] is not None and r["pool_crit"] is not None:
        sig = "*" if r["pool_at"] > r["pool_crit"] else ""
    fmt = lambda v, w=8, p=1: (f"{v:{w}.{p}f}" if isinstance(v, (int, float)) else f"{'n/a':>{w}s}")
    print(f"{r['sub']:8s} {fmt(r['fsp'],6,2)} {r['n2']:5d} {r['n3']:5d} {r['bal']:5.2f} | "
          f"{fmt(r['d_n2'])} {fmt(r['d_n3'])} | {fmt(r['b_n2'],6)} {fmt(r['b_n3'],6)} | "
          f"{fmt(r['pool_at'],7,3)} {fmt(r['pool_crit'],6,3)} {sig:3s}  {'' if ok else '<- excluded'}")

use = [r for r in rows if r["n2"] >= MIN_EP and r["n3"] >= MIN_EP and r["bal"] >= MIN_BAL]
print(f"\npassing staging-quality filter: {len(use)}/{len(rows)} "
      f"(>= {MIN_EP} epochs per stage and minority stage >= {MIN_BAL:.0%})")


def paired(tag, key2, key3, subset, higher):
    a = np.array([r[key2] for r in subset], float)
    b = np.array([r[key3] for r in subset], float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5:
        print(f"  {tag:22s} n={int(m.sum())} -- too few"); return
    p = stats.wilcoxon(a[m], b[m])[1]
    print(f"  {tag:22s} N2 median {np.median(a[m]):8.3f} | N3 median {np.median(b[m]):8.3f} | "
          f"n={int(m.sum())} p={p:.4f} | N2>N3 in {int((a[m]>b[m]).sum())}/{int(m.sum())}"
          f"   [{higher}]")


print("\n=== paired N2 vs N3, quality-filtered ===")
paired("3D SO->spindle MI_z", "d_n2", "d_n3", use, "Lecci predicts N2 > N3")
paired("3B SO->heartbeat z", "b_n2", "b_n3", use, "")
paired("3A coherence@0.02", "a_n2", "a_n3", use, "")

print("\n=== pooled-NREM 3A (all subjects, length-capped so K is comparable) ===")
pa = [(r["pool_at"], r["pool_crit"]) for r in rows if r["pool_at"] is not None and r["pool_crit"] is not None]
hits = sum(1 for a, c in pa if a > c)
if pa:
    bp = stats.binomtest(hits, len(pa), 0.05, alternative="greater").pvalue
    print(f"  {hits}/{len(pa)} subjects exceed their own alpha=0.05 threshold "
          f"(expected {0.05*len(pa):.1f} by chance) -> binomial p={bp:.4f}")
    print(f"  median coherence {np.median([a for a,_ in pa]):.3f}")

print("\n=== split by spindle peak (Lecci's effect is at the FAST spindle peak) ===")
for lab, sel in (("fast FSP > 12 Hz", [r for r in use if (r["fsp"] or 0) > FAST_FSP]),
                 ("slow FSP < 12 Hz", [r for r in use if (r["fsp"] or 99) <= FAST_FSP])):
    print(f"-- {lab}: n={len(sel)}")
    if len(sel) >= 5:
        paired("   3D SO->spindle MI_z", "d_n2", "d_n3", sel, "")
        paired("   3B SO->heartbeat z", "b_n2", "b_n3", sel, "")
