"""Cohort summary for the EVENT-BASED 3D (SO->spindle), N2-like vs N3-like.

Same staging-quality filter as summarize_cohort_stages.py so the two are comparable:
a stage needs >= MIN_EVENTS spindle events and >= MIN_CH channels tested.
"""
import glob, json, os
import numpy as np
from scipy import stats

OUT = "/Users/paris/Documents/iEEG/outputs/event_3D_by_stage"
MIN_EVENTS = 200      # spindle events per stage
MIN_CH = 3            # channels tested per stage

rows = []
for f in sorted(glob.glob(os.path.join(OUT, "*.json"))):
    d = json.load(open(f))
    if d.get("status") != "ok":
        continue
    r = dict(sub=d["subject"].replace("_phaseII", ""), fsp=d.get("fsp"),
             n2_ep=d.get("n_N2"), n3_ep=d.get("n_N3"))
    for st in ("N2", "N3", "ALL"):
        s = d.get(st) or {}
        r[f"{st}_sig"] = s.get("n_channels_significant")
        r[f"{st}_tested"] = s.get("n_channels_tested")
        r[f"{st}_R"] = s.get("median_R")
        r[f"{st}_ev"] = s.get("n_spindle_events")
        r[f"{st}_phase"] = s.get("mean_preferred_phase_deg")
    rows.append(r)

print(f"subjects with results: {len(rows)}\n")
hdr = f"{'sub':9s}{'FSP':>6s} | {'N2 sig':>7s}{'N2 R':>8s}{'N2 ev':>8s} | {'N3 sig':>7s}{'N3 R':>8s}{'N3 ev':>8s} | {'ALL sig':>8s}{'ALL R':>8s}"
print(hdr); print("-" * len(hdr))
for r in rows:
    f2 = lambda v, p=3: (f"{v:.{p}f}" if isinstance(v, float) else "  n/a")
    ok = (r["N2_ev"] or 0) >= MIN_EVENTS and (r["N3_ev"] or 0) >= MIN_EVENTS
    print(f"{r['sub']:9s}{r['fsp'] or 0:6.2f} | {str(r['N2_sig'])+'/'+str(r['N2_tested']):>7s}"
          f"{f2(r['N2_R']):>8s}{r['N2_ev'] or 0:8d} | {str(r['N3_sig'])+'/'+str(r['N3_tested']):>7s}"
          f"{f2(r['N3_R']):>8s}{r['N3_ev'] or 0:8d} | {str(r['ALL_sig'])+'/'+str(r['ALL_tested']):>8s}"
          f"{f2(r['ALL_R']):>8s}" + ("" if ok else "   <- excluded"))

use = [r for r in rows if (r["N2_ev"] or 0) >= MIN_EVENTS and (r["N3_ev"] or 0) >= MIN_EVENTS
       and (r["N2_tested"] or 0) >= MIN_CH and (r["N3_tested"] or 0) >= MIN_CH]
print(f"\npassing filter (>= {MIN_EVENTS} spindle events and >= {MIN_CH} channels per stage): {len(use)}/{len(rows)}")

# Is SO->spindle coupling present at all? (pooled NREM, per subject)
allr = [r for r in rows if r["ALL_R"] is not None]
sig_frac = [r["ALL_sig"] / r["ALL_tested"] for r in allr if r["ALL_tested"]]
n_any = sum(1 for r in allr if (r["ALL_sig"] or 0) >= 1)
print("\n=== Is SO->spindle coupling PRESENT? (pooled NREM) ===")
print(f"  subjects with >=1 significant channel : {n_any}/{len(allr)}")
print(f"  median fraction of channels significant: {np.median(sig_frac):.2f}")
print(f"  median resultant vector length R       : {np.median([r['ALL_R'] for r in allr]):.3f}")

# N2 vs N3
print("\n=== N2-like vs N3-like (quality-filtered, paired) ===")
for tag, key in (("median R", "R"), ("fraction of channels significant", "sig")):
    if key == "R":
        a = np.array([r["N2_R"] for r in use], float); b = np.array([r["N3_R"] for r in use], float)
    else:
        a = np.array([r["N2_sig"] / r["N2_tested"] for r in use], float)
        b = np.array([r["N3_sig"] / r["N3_tested"] for r in use], float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() >= 5:
        p = stats.wilcoxon(a[m], b[m])[1]
        print(f"  {tag:34s} N2 {np.median(a[m]):.3f} | N3 {np.median(b[m]):.3f} | "
              f"n={int(m.sum())} Wilcoxon p={p:.4f} | N2>N3 in {int((a[m]>b[m]).sum())}/{int(m.sum())}")
