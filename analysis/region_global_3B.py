"""Does SO->heartbeat coupling depend on SO globality and cortical region?

Motivation. Our lateral-iEEG 3B was weak/null, but the coupling Naji 2019 measured lives in LARGE,
GLOBAL, frontally-maximal (K-complex-type) slow waves near the central autonomic network (insula,
cingulate, medial PFC). A single lateral neocortical contact sees mostly LOCAL slow waves far from
that network, so a null there neither validates nor refutes the relationship. This script tests the
two falsifiable predictions of that account, on ds003848 (which has all-channel SO troughs cached and
Destrieux atlas labels):

  TEST 1 -- globality gradient. Cluster per-channel SO troughs into consensus events; each event's
    globality = how many channels participate within +/-DELTA. Bin events local -> regional -> global
    and run the SAME so_triggered (stage-matched null) on each bin. Prediction: HR peak and z RISE
    with globality. A rising gradient => the earlier null was a locality artifact; a flat one => the
    coupling is genuinely absent.

  TEST 2 -- region contrast. Group contacts by Destrieux label into autonomic-adjacent
    (frontal/cingulate/insula) vs posterior/lateral (temporal/parietal/occipital); run per-channel
    so_triggered within each group and average. Prediction: stronger in the autonomic-adjacent group.

Both run from the cache -- no re-streaming. n=6 and coverage is clinical (some subjects lack frontal/
insular contacts), so this is exploratory/directional, reported with per-subject coverage.

    .venv/bin/python analysis/region_global_3B.py
"""
import csv, glob, json, os
import numpy as np

import lecci_faithful_3A as L
from event_3B_mednick import so_triggered, FS_RR
from stage_ds003848 import OUT as DS_CACHE
from spectral_gapped import fill_short_gaps
from cohort_stages_3ABD import EPOCH

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ELEC = os.path.join(ROOT, "data", "ds003848_raw", "electrodes")
OUT = os.path.join(ROOT, "outputs", "region_global_3B")

DELTA_S = 0.15        # co-occurrence tolerance for "same" SO across channels
GRID_S = 0.05         # consensus raster resolution
GLOBALITY_BINS = [(0.0, 0.15, "local"), (0.15, 0.40, "regional"), (0.40, 1.01, "global")]
AUTONOMIC = {"frontal", "cingulate", "insula"}
POSTERIOR = {"temporal", "parietal", "occipital"}


def region_of(label):
    l = (label or "").lower()
    if "front" in l:
        return "frontal"
    if "cingul" in l:
        return "cingulate"
    if "insula" in l or "ins_" in l or "insular" in l:
        return "insula"
    if "temporal" in l or "hippocamp" in l or "amygd" in l or "fusi" in l:
        return "temporal"
    if "pariet" in l or "precuneus" in l or "supramarg" in l or "angular" in l:
        return "parietal"
    if "occipital" in l or "calcarine" in l or "lingual" in l or "cuneus" in l:
        return "occipital"
    return "other"


def load_regions(subject):
    fp = os.path.join(ELEC, f"{subject}_electrodes.tsv")
    if not os.path.exists(fp):
        return {}
    out = {}
    for r in csv.DictReader(open(fp), delimiter="\t"):
        if r.get("name"):
            out[r["name"]] = region_of(r.get("Destrieux_label_text") or r.get("Destrieux_label") or "")
    return out


def stage_pool_and_mean(hr, lab, stage):
    eps = np.where(lab == stage)[0]
    if len(eps) < 10:
        return None, None, None
    m = np.zeros(len(hr), bool)
    for e in eps:
        m[int(e * EPOCH * FS_RR):int((e + 1) * EPOCH * FS_RR)] = True
    m &= np.isfinite(hr)
    if m.sum() < 400:
        return None, None, None
    return float(hr[m].mean()), np.where(m)[0], set(eps.tolist())


def consensus_events(so_by_ch, total_s, n_ch):
    """Cluster per-channel troughs into consensus events with a globality count.

    Each channel contributes at most 1 to a location (its trough train is dilated by +/-DELTA and
    binarised); participation summed across channels = globality. Returns (event_times_s, globality)."""
    T = int(total_s / GRID_S) + 1
    box = np.ones(int(2 * DELTA_S / GRID_S) + 1)
    partic = np.zeros(T)
    for ch, tr in so_by_ch.items():
        if not len(tr):
            continue
        idx = (np.asarray(tr) / GRID_S).astype(int)
        idx = idx[(idx >= 0) & (idx < T)]
        raster = np.zeros(T)
        raster[idx] = 1.0
        partic += (np.convolve(raster, box, "same") > 0).astype(float)
    # local maxima of participation = consensus events
    from scipy.signal import find_peaks
    pk, _ = find_peaks(partic, height=1, distance=int(DELTA_S / GRID_S))
    return pk * GRID_S, partic[pk] / max(n_ch, 1)


def analyse(subject):
    L.set_cache(DS_CACHE)
    d = L.load(subject)
    if d is None:
        return None
    regions = load_regions(subject)
    hr, _, _ = fill_short_gaps(d["hr_4"], FS_RR, max_gap_s=5.0)
    lab, nrem, _ = L.stages_for(d)
    ctx = [str(c) for c in d["cortical_chans"]]
    so_by_ch = {c: (d[f"so_t_{c}"] if f"so_t_{c}" in d.files else np.array([])) for c in ctx}
    total_s = len(d["hr_4"]) / FS_RR
    rng = np.random.RandomState(0)

    rec = dict(subject=subject, n_channels=len(ctx),
               region_counts={g: sum(1 for c in ctx if regions.get(c) == g)
                              for g in ("frontal", "cingulate", "insula", "temporal",
                                        "parietal", "occipital", "other")})

    ev_t, glob = consensus_events(so_by_ch, total_s, len(ctx))

    for stage in ("N2", "N3"):
        stage_mean, pool, keep = stage_pool_and_mean(hr, lab, stage)
        if pool is None:
            rec[stage] = None
            continue
        out = dict(stage_mean_hr=stage_mean)

        # --- TEST 1: globality gradient (consensus events in this stage) ---
        in_stage = np.array([int(t // EPOCH) in keep for t in ev_t])
        grad = []
        for lo, hi, name in GLOBALITY_BINS:
            sel = in_stage & (glob >= lo) & (glob < hi)
            tt = ev_t[sel]
            r = so_triggered(hr, tt, stage_mean, pool, rng=rng) if len(tt) >= 30 else None
            grad.append(dict(bin=name, glob_range=[lo, hi], n_events=int(sel.sum()),
                             pct=(None if not r else r["pct_above_stage_mean"]),
                             z=(None if not r else r["z"]),
                             lag=(None if not r else r["peak_lag_s"])))
        out["globality_gradient"] = grad

        # --- TEST 2: region contrast (per-channel within group, averaged) ---
        for grp_name, grp in (("autonomic", AUTONOMIC), ("posterior", POSTERIOR)):
            chans = [c for c in ctx if regions.get(c) in grp]
            zs, pcts, nso = [], [], 0
            for c in chans:
                tt = [float(t) for t in so_by_ch[c] if int(t // EPOCH) in keep]
                r = so_triggered(hr, tt, stage_mean, pool, rng=rng)
                if r:
                    zs.append(r["z"]); pcts.append(r["pct_above_stage_mean"]); nso += r["n_so"]
            out[f"region_{grp_name}"] = dict(
                n_channels=len(chans), n_channels_usable=len(zs), n_so=nso,
                z=(float(np.mean(zs)) if zs else None),
                pct=(float(np.mean(pcts)) if pcts else None))
        rec[stage] = out
    return rec


def main():
    os.makedirs(OUT, exist_ok=True)
    from stage_ds003848 import SUBJECTS
    rows = []
    for s in SUBJECTS:
        try:
            r = analyse(s)
        except Exception as e:
            import traceback
            print(f"[{s}] ERROR {type(e).__name__}: {e}", flush=True); traceback.print_exc(); continue
        if r is None:
            print(f"[{s}] no cache", flush=True); continue
        json.dump(r, open(os.path.join(OUT, f"{s}.json"), "w"))
        rows.append(r)

    # ---- report ----
    print("\n" + "=" * 92)
    print("TEST 1 -- GLOBALITY GRADIENT (does SO->HR coupling rise from local to global SOs?)")
    print("=" * 92)
    for stage in ("N2", "N3"):
        print(f"\n[{stage}]  per-subject % above stage-mean HR (z) by SO globality bin")
        print(f"{'subject':16s} {'local':>18s} {'regional':>18s} {'global':>18s}")
        agg = {b[2]: {"pct": [], "z": []} for b in GLOBALITY_BINS}
        for r in rows:
            st = r.get(stage)
            if not st:
                print(f"{r['subject']:16s} {'(stage n/a)':>18s}"); continue
            cells = []
            for g in st["globality_gradient"]:
                if g["pct"] is None:
                    cells.append(f"{'-':>18s}")
                else:
                    cells.append(f"{g['pct']:+6.2f}% (z{g['z']:+.1f}) ".rjust(18))
                    agg[g["bin"]]["pct"].append(g["pct"]); agg[g["bin"]]["z"].append(g["z"])
            print(f"{r['subject']:16s} " + "".join(cells))
        gg = " -> ".join(f"{b[2]} {np.mean(agg[b[2]]['pct']):+.2f}% (z{np.mean(agg[b[2]]['z']):+.2f}, n={len(agg[b[2]]['pct'])})"
                         if agg[b[2]]['pct'] else f"{b[2]} n/a" for b in GLOBALITY_BINS)
        print(f"  POOLED: {gg}")

    print("\n" + "=" * 92)
    print("TEST 2 -- REGION CONTRAST (autonomic-adjacent front/cing/insula vs posterior/lateral)")
    print("=" * 92)
    for stage in ("N2", "N3"):
        print(f"\n[{stage}]  % above stage-mean HR (z) | usable channels")
        print(f"{'subject':16s} {'AUTONOMIC(f/c/i)':>26s} {'POSTERIOR(t/p/o)':>26s}")
        a_z, p_z = [], []
        for r in rows:
            st = r.get(stage)
            if not st:
                print(f"{r['subject']:16s} {'(stage n/a)':>26s}"); continue
            a, p = st["region_autonomic"], st["region_posterior"]
            af = f"{a['pct']:+.2f}% z{a['z']:+.1f} [{a['n_channels_usable']}ch]" if a['z'] is not None else f"no ch [{a['n_channels']}]"
            pf = f"{p['pct']:+.2f}% z{p['z']:+.1f} [{p['n_channels_usable']}ch]" if p['z'] is not None else f"no ch [{p['n_channels']}]"
            print(f"{r['subject']:16s} {af:>26s} {pf:>26s}")
            if a['z'] is not None: a_z.append(a['z'])
            if p['z'] is not None: p_z.append(p['z'])
        if a_z and p_z:
            from scipy import stats
            print(f"  POOLED z: autonomic mean {np.mean(a_z):+.2f} (n={len(a_z)}) | "
                  f"posterior mean {np.mean(p_z):+.2f} (n={len(p_z)})")
    print(f"\nper-subject JSON -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
