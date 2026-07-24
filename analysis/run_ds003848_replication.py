"""Run the corrected 3A approximation and 3B stage-matched analysis on ds003848.

ds003848 is the independent replication set: an OpenNeuro iEEG sleep dataset (Utrecht RESPect) with
ECG + EMG + EOG. Author-provided sleep/NREM/REM/SWS and transition annotations define the primary
coarse states; the unvalidated multimodal proxy is sensitivity-only within author-unknown sleep.
These are not expert AASM N2/N3 labels. This reuses the corrected analysis code, pointed at the
ds003848 cache.

    .venv/bin/python analysis/stage_ds003848.py           # build the cache first
    .venv/bin/python analysis/run_ds003848_replication.py [--force]
    .venv/bin/python analysis/summarize_corrected_3AB.py \
        --a-dir outputs/ds003848_3A --b-dir outputs/ds003848_3B --label "ds003848 replication"
"""
import argparse, json, os
import lecci_faithful_3A as L
import event_3B_cached as B
from stage_ds003848 import SUBJECTS, OUT as DS_CACHE
from pipeline_version import (ANALYSIS_VERSION, CACHE_SCHEMA_VERSION, atomic_json_dump,
                              cache_code_sha256, file_sha256, source_tree_sha256, start_run_manifest,
                              validated_complete_run_exists, write_run_manifest)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_A = os.path.join(ROOT, "outputs", "ds003848_3A")
OUT_B = os.path.join(ROOT, "outputs", "ds003848_3B")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--force", action="store_true",
        help="replace results that are not two exact reusable completed endpoint runs")
    args = parser.parse_args()
    L.set_cache(DS_CACHE)          # both L.analyse and B.analyse resolve subjects via L.load
    os.makedirs(OUT_A, exist_ok=True)
    os.makedirs(OUT_B, exist_ok=True)
    config_a = dict(band="fixed", smooth_4s=True, n_sur=200,
                    analysis_version=ANALYSIS_VERSION,
                    cache_schema_version=CACHE_SCHEMA_VERSION,
                    cache_code_sha256=cache_code_sha256(ROOT))
    config_b = dict(tachogram_domain="rr", n_surrogates=1000,
                    analysis_version=ANALYSIS_VERSION,
                    cache_schema_version=CACHE_SCHEMA_VERSION,
                    cache_code_sha256=cache_code_sha256(ROOT),
                    null_method="shared circular shift in eligible stage-time",
                    contact_qc=B.CONTACT_QC_METHOD,
                    minimum_event_channels=B.MIN_EVENT_CHANNELS,
                    stable_stage_minimum_s=int(B.MIN_STABLE_STAGE_EPOCHS * B.EPOCH))
    requested = list(SUBJECTS)
    existing_a = any(
        os.path.exists(os.path.join(OUT_A, f"{subject}.json"))
        for subject in requested)
    existing_b = any(
        os.path.exists(os.path.join(OUT_B, f"{subject}.json"))
        for subject in requested)
    if not args.force and (existing_a or existing_b):
        expected_cache_inputs = {
            subject: L.cache_lineage_entry(L.cache_lineage(subject))
            for subject in requested
        }
        reusable_a = (
            existing_a
            and validated_complete_run_exists(
                OUT_A, pipeline="ds003848_lecci_3A", requested=requested,
                config={**config_a, "cache_inputs": expected_cache_inputs},
                suffix=".json", require_current_source_tree=True)
        )
        reusable_b = (
            existing_b
            and validated_complete_run_exists(
                OUT_B, pipeline="ds003848_event_3B", requested=requested,
                config={**config_b, "cache_inputs": expected_cache_inputs},
                suffix=".json", require_current_source_tree=True)
        )
        if reusable_a and reusable_b:
            print(
                f"validated existing complete ds003848 3A/3B runs "
                f"({len(requested)} subjects)",
                flush=True)
            return
        raise RuntimeError(
            "only one ds003848 endpoint run is reusable; rerun both with --force")
    run_id_a = start_run_manifest(
        OUT_A, pipeline="ds003848_lecci_3A", requested=requested, config=config_a)
    run_id_b = start_run_manifest(
        OUT_B, pipeline="ds003848_event_3B", requested=requested, config=config_b)
    tree_digest = source_tree_sha256(ROOT)
    completed_a, completed_b, skipped_a, skipped_b, failed_a, failed_b = [], [], [], [], [], []
    cache_inputs_a, cache_inputs_b = {}, {}
    for s in SUBJECTS:
        if not os.path.exists(os.path.join(DS_CACHE, f"{s}.npz")):
            error = "no cache -- run stage_ds003848.py first"
            print(f"[{s}] {error}", flush=True)
            failed_a.append(dict(subject=s, error=error))
            failed_b.append(dict(subject=s, error=error))
            continue
        # 3A
        try:
            r = L.analyse(s, band="fixed", smooth_4s=True, n_sur=200)
        except Exception as e:
            print(f"[{s}] 3A ERROR {type(e).__name__}: {e}", flush=True)
            failed_a.append(dict(subject=s, error=f"{type(e).__name__}: {e}"))
            r = None
        if r is None:
            error = "analysis returned no record (cache is unavailable or not status=ok)"
            print(f"[{s}] 3A ERROR {error}", flush=True)
            if not any(value.get("subject") == s for value in failed_a):
                failed_a.append(dict(subject=s, error=error))
        else:
            r["run_id"] = run_id_a
            r["source_tree_sha256"] = tree_digest
            atomic_json_dump(r, os.path.join(OUT_A, f"{s}.json"))
            cache_inputs_a[s] = L.cache_lineage_entry(r)
        if r and r.get("status") == "ok":
            completed_a.append(s)
            pk = r["peak"]; co = r.get("coherence") or {}; xc = r.get("xcorr") or {}
            print(f"[{s}] 3A bouts={r['n_bouts']} peak="
                  f"{'none' if not pk.get('peak_hz') else format(pk['peak_hz'],'.4f')} "
                  f"coh@own={co.get('at_own_peak')} xcorr r={xc.get('peak_r')}", flush=True)
        elif r is not None:
            print(f"[{s}] 3A skip: {r.get('reason') if r else 'no result'}", flush=True)
            skipped_a.append(dict(subject=s, reason=r.get("reason", "no result")))
        # 3B
        try:
            rb = B.analyse(s)
        except Exception as e:
            print(f"[{s}] 3B ERROR {type(e).__name__}: {e}", flush=True)
            failed_b.append(dict(subject=s, error=f"{type(e).__name__}: {e}"))
            rb = None
        if rb is None:
            error = "analysis returned no record (cache is unavailable or not status=ok)"
            print(f"[{s}] 3B ERROR {error}", flush=True)
            if not any(value.get("subject") == s for value in failed_b):
                failed_b.append(dict(subject=s, error=error))
        else:
            rb["run_id"] = run_id_b
            rb["source_tree_sha256"] = tree_digest
            atomic_json_dump(rb, os.path.join(OUT_B, f"{s}.json"))
            cache_inputs_b[s] = L.cache_lineage_entry(rb)
        if rb and rb.get("status") == "ok":
            completed_b.append(s)
            for st in ("N2", "N3"):
                v = rb.get(st)
                if v:
                    print(f"[{s}] 3B {st}: raw {v['pct_above_stage_mean']:+.2f}% "
                          f"local Δ={v['event_locked_local_change_pct']:+.2f}% "
                          f"(inference disabled) nSO={v['n_so_total']}", flush=True)
        elif rb is not None:
            print(f"[{s}] 3B skip: {rb.get('reason') if rb else 'no result'}", flush=True)
            skipped_b.append(dict(subject=s, reason=rb.get("reason", "no result")))
    write_run_manifest(
        OUT_A, pipeline="ds003848_lecci_3A", requested=requested,
        completed=completed_a, skipped=skipped_a, failed=failed_a,
        config={**config_a, "cache_inputs": cache_inputs_a},
        run_id=run_id_a,
        result_files_sha256={
            subject: file_sha256(os.path.join(OUT_A, f"{subject}.json"))
            for subject in completed_a + [value["subject"] for value in skipped_a]
        })
    write_run_manifest(
        OUT_B, pipeline="ds003848_event_3B", requested=requested,
        completed=completed_b, skipped=skipped_b, failed=failed_b,
        config={**config_b, "cache_inputs": cache_inputs_b},
        run_id=run_id_b,
        result_files_sha256={
            subject: file_sha256(os.path.join(OUT_B, f"{subject}.json"))
            for subject in completed_b + [value["subject"] for value in skipped_b]
        })
    print(f"\n3A -> {OUT_A}\n3B -> {OUT_B}", flush=True)
    if failed_a or failed_b:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
