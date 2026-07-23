"""Run the corrected 3A (Lecci-faithful) and 3B (stage-matched null) on the ds003848 cohort.

ds003848 is the independent replication set: an OpenNeuro iEEG sleep dataset (Utrecht RESPect) with
ECG + EMG + EOG, so NREM is REM/wake-excluded via real staging (`stage_ds003848.py`) rather than the
GMM-on-slow-wave proxy used for HUP. This reuses the exact corrected analysis code, just pointed at
the ds003848 cache and real stage labels.

    .venv/bin/python analysis/stage_ds003848.py           # build the cache first
    .venv/bin/python analysis/run_ds003848_replication.py
    .venv/bin/python analysis/summarize_corrected_3AB.py \
        --a-dir outputs/ds003848_3A --b-dir outputs/ds003848_3B --label "ds003848 replication"
"""
import json, os
import lecci_faithful_3A as L
import event_3B_cached as B
from stage_ds003848 import SUBJECTS, OUT as DS_CACHE

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_A = os.path.join(ROOT, "outputs", "ds003848_3A")
OUT_B = os.path.join(ROOT, "outputs", "ds003848_3B")


def main():
    L.set_cache(DS_CACHE)          # both L.analyse and B.analyse resolve subjects via L.load
    os.makedirs(OUT_A, exist_ok=True)
    os.makedirs(OUT_B, exist_ok=True)
    for s in SUBJECTS:
        if not os.path.exists(os.path.join(DS_CACHE, f"{s}.npz")):
            print(f"[{s}] no cache -- run stage_ds003848.py first", flush=True); continue
        # 3A
        try:
            r = L.analyse(s, n_sur=200)
        except Exception as e:
            print(f"[{s}] 3A ERROR {type(e).__name__}: {e}", flush=True); r = None
        if r and r.get("status") == "ok":
            json.dump(r, open(os.path.join(OUT_A, f"{s}.json"), "w"))
            pk = r["peak"]; co = r.get("coherence") or {}; xc = r.get("xcorr") or {}
            print(f"[{s}] 3A bouts={r['n_bouts']} peak="
                  f"{'none' if not pk.get('peak_hz') else format(pk['peak_hz'],'.4f')} "
                  f"coh@own={co.get('at_own_peak')} xcorr r={xc.get('peak_r')}", flush=True)
        else:
            print(f"[{s}] 3A skip: {r.get('reason') if r else 'no result'}", flush=True)
        # 3B
        try:
            rb = B.analyse(s)
        except Exception as e:
            print(f"[{s}] 3B ERROR {type(e).__name__}: {e}", flush=True); rb = None
        if rb and rb.get("status") == "ok":
            json.dump(rb, open(os.path.join(OUT_B, f"{s}.json"), "w"))
            for st in ("N2", "N3"):
                v = rb.get(st)
                if v:
                    print(f"[{s}] 3B {st}: {v['pct_above_stage_mean']:+.2f}% z={v['z']:+.2f} "
                          f"nSO={v['n_so_total']}", flush=True)
        else:
            print(f"[{s}] 3B skip: {rb.get('reason') if rb else 'no result'}", flush=True)
    print(f"\n3A -> {OUT_A}\n3B -> {OUT_B}", flush=True)


if __name__ == "__main__":
    main()
