"""Does the SO->heartbeat coupling live in N2 K-complexes? (the test our globality analysis missed)

Naji 2019's headline effect is in STAGE 2 (+12% vs +3.35% in SWS), and Stage 2 is the home of the
K-complex -- a large, ISOLATED, frontally-maximal slow wave tied to micro-arousals and phasic
autonomic bursts (de Zambotti 2016). Our earlier globality test accidentally measured N3 instead:
its cross-channel-consensus metric fires on the DENSE slow-wave trains of SWS, and N2 has almost no
"global" events by that metric (0-40 per subject vs 12-131 in N3). So it never tested the K-complex
mechanism.

This does. On the Utrecht RESPect cohort (OpenNeuro ds003848 -- real EMG/EOG staging + Destrieux atlas
labels), it isolates a K-complex proxy and tests its HR coupling against three controls:

  K-complex proxy = a slow-wave event that is
     (1) in N2,                                   [K-complexes are an N2 hallmark]
     (2) on FRONTAL / CINGULATE (midline) contacts, [K-complex generators]
     (3) ISOLATED -- no other frontal event within +/-2.5 s.   [not a SWS train]
  Controls:  N2 frontal TRAIN (non-isolated) · N2 POSTERIOR isolated · N3 frontal isolated.

Prediction (if Naji's mechanism is present and visible in iEEG): the N2 frontal-isolated proxy shows
the strongest SO->HR coupling.

CAVEAT: the cache stores SO trough TIMES only, not raw morphology, so this is an "isolated large
slow-wave event in N2 on frontal contacts" proxy -- it does NOT verify biphasic K-complex morphology
(that would need re-downloading the raw BrainVision). All cached troughs already passed the
top-25%-amplitude SO criterion, so they are the large waves. Exploratory; n<=6, sparse frontal
coverage in some subjects.

    .venv/bin/python analysis/kcomplex_n2_3B.py
"""
import csv, json, os
import numpy as np

import lecci_faithful_3A as L
from event_3B_mednick import so_triggered, FS_RR
from stage_ds003848 import OUT as DS_CACHE, SUBJECTS
from region_global_3B import load_regions, ELEC, stage_pool_and_mean
from spectral_gapped import fill_short_gaps
from cohort_stages_3ABD import EPOCH

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs", "kcomplex_n2_3B")

MERGE_S = 0.15        # merge near-simultaneous troughs across channels into one event
ISO_S = 2.5          # isolation: no other same-group event within this window => K-complex-like
FRONTAL_GEN = {"frontal", "cingulate"}          # K-complex generators (midline frontal-central)
POSTERIOR = {"temporal", "parietal", "occipital"}


def merged_events(so_by_ch, chans):
    """Pool troughs from `chans`, merge those within MERGE_S into single events (median time)."""
    t = np.concatenate([np.asarray(so_by_ch[c], float) for c in chans if len(so_by_ch.get(c, []))]) \
        if chans else np.array([])
    if len(t) == 0:
        return np.array([])
    t.sort()
    events, cur = [], [t[0]]
    for x in t[1:]:
        if x - cur[-1] <= MERGE_S:
            cur.append(x)
        else:
            events.append(np.median(cur)); cur = [x]
    events.append(np.median(cur))
    return np.array(events)


def isolated(events, iso_s=ISO_S):
    """Boolean mask: event has no neighbour within +/-iso_s (K-complex-like isolation)."""
    if len(events) == 0:
        return np.array([], bool)
    prev = np.concatenate(([np.inf], np.diff(events)))
    nxt = np.concatenate((np.diff(events), [np.inf]))
    return (prev > iso_s) & (nxt > iso_s)


def run_condition(hr, events, keep_epochs, stage_mean, pool, rng):
    tt = np.array([t for t in events if int(t // EPOCH) in keep_epochs])
    if len(tt) < 30:
        return dict(n=int(len(tt)), pct=None, z=None, lag=None)
    r = so_triggered(hr, tt, stage_mean, pool, rng=rng)
    if not r:
        return dict(n=int(len(tt)), pct=None, z=None, lag=None)
    return dict(n=r["n_so"], pct=r["pct_above_stage_mean"], z=r["z"], lag=r["peak_lag_s"])


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

    front = [c for c in ctx if regions.get(c) in FRONTAL_GEN]
    post = [c for c in ctx if regions.get(c) in POSTERIOR]
    ev_front = merged_events(so_by_ch, front)
    ev_post = merged_events(so_by_ch, post)
    iso_front = isolated(ev_front)
    iso_post = isolated(ev_post)

    rec = dict(subject=subject, n_frontal_ch=len(front), n_posterior_ch=len(post),
               n_frontal_events=int(len(ev_front)), frontal_isolated_frac=float(iso_front.mean()) if len(iso_front) else None)
    rng = np.random.RandomState(0)

    for stage in ("N2", "N3"):
        stage_mean, pool, keep = stage_pool_and_mean(hr, lab, stage)
        if pool is None:
            rec[stage] = None
            continue
        conds = {}
        conds["frontal_isolated"] = run_condition(hr, ev_front[iso_front], keep, stage_mean, pool, rng)
        conds["frontal_train"] = run_condition(hr, ev_front[~iso_front] if len(iso_front) else np.array([]),
                                               keep, stage_mean, pool, rng)
        conds["posterior_isolated"] = run_condition(hr, ev_post[iso_post] if len(iso_post) else np.array([]),
                                                    keep, stage_mean, pool, rng)
        rec[stage] = dict(stage_mean_hr=stage_mean, **conds)
    return rec


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for s in SUBJECTS:
        try:
            r = analyse(s)
        except Exception as e:
            import traceback
            print(f"[{s}] ERROR {type(e).__name__}: {e}"); traceback.print_exc(); continue
        if r is None:
            print(f"[{s}] no cache"); continue
        json.dump(r, open(os.path.join(OUT, f"{s}.json"), "w"))
        rows.append(r)

    def cell(c):
        return f"{c['pct']:+.2f}% z{c['z']:+.1f} (n={c['n']})" if c and c.get("z") is not None \
            else (f"n={c['n']}<30" if c else "-")

    print("\n" + "=" * 100)
    print("K-COMPLEX TEST -- SO->HR coupling for N2 frontal-ISOLATED events (K-complex proxy) vs controls")
    print("=" * 100)
    print(f"\nfrontal/cingulate generator contacts, isolated = no neighbour within +/-{ISO_S}s\n")
    print(f"{'subject':15s} {'front ch':>8s} {'N2 front-ISO (K-cplx)':>24s} {'N2 front-train':>20s} "
          f"{'N2 post-ISO':>18s} {'N3 front-ISO':>18s}")
    agg = {"kc": [], "train": [], "post": [], "n3iso": []}
    for r in rows:
        n2 = r.get("N2"); n3 = r.get("N3")
        kc = cell(n2["frontal_isolated"]) if n2 else "(N2 n/a)"
        tr = cell(n2["frontal_train"]) if n2 else "-"
        po = cell(n2["posterior_isolated"]) if n2 else "-"
        n3iso = cell(n3["frontal_isolated"]) if n3 else "(N3 n/a)"
        print(f"{r['subject']:15s} {r['n_frontal_ch']:8d} {kc:>24s} {tr:>20s} {po:>18s} {n3iso:>18s}")
        if n2 and n2["frontal_isolated"].get("z") is not None: agg["kc"].append(n2["frontal_isolated"])
        if n2 and n2["frontal_train"].get("z") is not None: agg["train"].append(n2["frontal_train"])
        if n2 and n2["posterior_isolated"].get("z") is not None: agg["post"].append(n2["posterior_isolated"])
        if n3 and n3["frontal_isolated"].get("z") is not None: agg["n3iso"].append(n3["frontal_isolated"])

    print("\nPOOLED (mean z, mean % above stage-mean HR):")
    for key, name in (("kc", "N2 frontal-ISOLATED (K-complex proxy)"), ("train", "N2 frontal-train"),
                      ("post", "N2 posterior-isolated"), ("n3iso", "N3 frontal-isolated")):
        a = agg[key]
        if a:
            print(f"  {name:38s} z={np.mean([x['z'] for x in a]):+.2f}  "
                  f"%={np.mean([x['pct'] for x in a]):+.2f}  (n={len(a)} subjects)")
        else:
            print(f"  {name:38s} no usable subjects")
    print(f"\nPrediction: if Naji's N2/K-complex mechanism is visible in iEEG, the K-complex proxy "
          f"is strongest.\nper-subject JSON -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
