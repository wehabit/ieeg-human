"""3B (Naji 2019 SO -> heart-rate coupling) across the WHOLE cohort, from the cached series.

Two problems with the previous 3B are fixed here:

  1. SCOPE. `event_3B_mednick.py` re-streamed each night and had only ever been run on HUP165, yet
     the summary reported "SO->heartbeat coupling is present but ~25x weaker than published" as a
     cohort-level finding. It was n=1. This runs every cached subject.
  2. THE NULL. Surrogate triggers were drawn from the whole night while the statistic was expressed
     as "% above THIS STAGE's mean HR", so the null measured the stage-vs-night mean-HR offset
     rather than SO-locking. On zero-effect synthetic data that produced z = -21.8 / +24.0 in the
     two stages (test_3B_null.py). Surrogates are now drawn from the same stage as the real
     triggers, via the fixed `so_triggered`.

Method is otherwise Naji 2019 as implemented in event_3B_mednick: SO 0.15-4 Hz per-channel
zero-crossing half-waves, RR at 4 Hz by cubic spline, HR averaged +/-5 s on the down-state trough,
reported as % above that stage's mean HR plus the SO->HR peak latency.

    .venv/bin/python analysis/event_3B_cached.py
"""
import argparse, json, os
import numpy as np

from event_3B_mednick import so_triggered, FS_RR
from lecci_faithful_3A import load, CACHE
from cohort_stages_3ABD import stage_epochs, EPOCH
from cohort_3A_cortical import COHORT
from spectral_gapped import fill_short_gaps

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs", "event_3B_cached")

NAJI = {"N2": 12.09, "N3": 3.35}          # % above stage mean HR, frontal scalp, healthy sleepers


def analyse(subject):
    d = load(subject)
    if d is None:
        return None
    hr, _, _ = fill_short_gaps(d["hr_4"], FS_RR, max_gap_s=5.0)
    if np.isfinite(hr).sum() < 1000:
        return dict(subject=subject, status="skip", reason="no usable RR")
    ep = dict(dr=d["ep_dr"], swa=d["ep_swa"], clean=d["ep_clean"])
    lab, nrem, sep = stage_epochs(ep)
    ctx = [str(c) for c in d["cortical_chans"]]
    rec = dict(subject=subject, status="ok", n_N2=int((lab == "N2").sum()),
               n_N3=int((lab == "N3").sum()), n_channels=len(ctx),
               gmm_separation=(None if sep is None else float(sep)),
               whole_night_mean_hr=float(np.nanmean(hr)))
    rng = np.random.RandomState(0)
    for stage in ("N2", "N3"):
        eps = np.where(lab == stage)[0]
        if len(eps) < 20:
            rec[stage] = None; continue
        m = np.zeros(len(hr), bool)
        for e in eps:
            m[int(e * EPOCH * FS_RR):int((e + 1) * EPOCH * FS_RR)] = True
        m &= np.isfinite(hr)
        if m.sum() < 400:
            rec[stage] = None; continue
        stage_mean = float(hr[m].mean())
        pool = np.where(m)[0]
        keep = set(eps.tolist())
        per_ch = []
        for c in ctx:
            key = f"so_t_{c}"
            if key not in d.files:
                continue
            tt = [float(x) for x in d[key] if int(x // EPOCH) in keep]
            r = so_triggered(hr, tt, stage_mean, pool, rng=rng)
            if r:
                r.pop("curve"); r.pop("lag_s"); r["ch"] = c
                per_ch.append(r)
        if not per_ch:
            rec[stage] = None; continue
        rec[stage] = dict(
            stage_mean_hr=stage_mean, n_channels=len(per_ch),
            n_so_total=int(sum(x["n_so"] for x in per_ch)),
            pct_above_stage_mean=float(np.mean([x["pct_above_stage_mean"] for x in per_ch])),
            peak_lag_s=float(np.mean([x["peak_lag_s"] for x in per_ch])),
            z=float(np.mean([x["z"] for x in per_ch])),
            null_mean_pct=float(np.mean([x["null_mean_pct"] for x in per_ch])),
            per_channel=per_ch)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default=",".join(map(str, COHORT)))
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for n in [int(x) for x in a.subjects.split(",") if x.strip()]:
        name = f"HUP{n}_phaseII"
        try:
            r = analyse(name)
        except Exception as e:
            print(f"[{name}] ERROR {type(e).__name__}: {e}", flush=True); continue
        if r is None:
            print(f"[{name}] no cache", flush=True); continue
        if r.get("status") != "ok":
            print(f"[{name}] {r.get('reason')}", flush=True); continue
        json.dump(r, open(os.path.join(OUT, f"{name}.json"), "w"))
        rows.append(r)
        parts = []
        for st in ("N2", "N3"):
            v = r.get(st)
            parts.append(f"{st} n/a" if not v else
                         f"{st} {v['pct_above_stage_mean']:+.2f}% (null {v['null_mean_pct']:+.2f}) "
                         f"z={v['z']:+.2f} lag={v['peak_lag_s']:.1f}s nSO={v['n_so_total']}")
        print(f"[{name}] " + " | ".join(parts), flush=True)
    print(f"\n{len(rows)} subjects -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
