"""
Pull one full night of NREM from a continuous iEEG.org dataset (default HUP165_phaseII).

No sleep staging exists in the portal, so we stage blind from the iEEG's own delta power:
  A. NIGHT FINDER  - probe 1 channel every 15 min over the first SCAN_H hours; the longest
     sustained high-delta stretch = a night.
  B. STAGE         - pull that 1 channel across the night; score each 30 s epoch NREM
     (relative delta 0.5-4 Hz dominant) vs not.
  C. PULL NREM     - download the MTL channel set for the contiguous NREM blocks only, saved
     block-by-block (memory-safe) to data/ieeg_portal/<out>/.

Output: per-block .npz (data[samples x ch], ch_names, sfreq, start_sec) + index.csv + meta.json.
"""
import argparse, json, os, re, time
import numpy as np, pandas as pd
from scipy import signal
from ieeg.auth import Session

MAX_REQ_S = 600       # iEEG.org caps get_data at ~600 s per request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CRED = os.path.join(ROOT, "data", "ieeg_secret", "credentials.json")
EPOCH = 30.0          # s, staging epoch
DELTA = (0.5, 4.0); DENOM = (0.5, 25.0)


def sess():
    c = json.load(open(CRED))
    return Session(c["username"], c["password"])


def get(ds, idx, start_s, dur_s, tries=3):
    """get_data with usec; retries transient portal errors. Caps at MAX_REQ_S."""
    dur_s = min(dur_s, MAX_REQ_S)
    for k in range(tries):
        try:
            return ds.get_data(int(start_s * 1e6), int(dur_s * 1e6), idx)
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(1.5)


def delta_ratio(x, sf):
    m = np.isfinite(x)
    if m.mean() < 0.7:                    # too much gap -> unusable
        return np.nan
    if not m.all():                        # interpolate short gaps
        x = x.copy(); idx = np.arange(len(x))
        x[~m] = np.interp(idx[~m], idx[m], x[m])
    if np.std(x) < 1e-9:
        return np.nan
    f, p = signal.welch(x, sf, nperseg=int(min(4 * sf, len(x))))
    num = np.trapezoid(p[(f >= DELTA[0]) & (f < DELTA[1])], f[(f >= DELTA[0]) & (f < DELTA[1])])
    den = np.trapezoid(p[(f >= DENOM[0]) & (f < DENOM[1])], f[(f >= DENOM[0]) & (f < DENOM[1])])
    return float(num / den) if den > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="HUP165_phaseII")
    ap.add_argument("--scan-h", type=float, default=30.0)
    ap.add_argument("--scan-step-min", type=float, default=15.0)
    ap.add_argument("--night-h", type=float, default=9.0)
    ap.add_argument("--mtl-regex", default=r"^L[ABCH]\d+$")   # left mesial-temporal shafts
    ap.add_argument("--max-contact", type=int, default=8)
    ap.add_argument("--nrem-pct", type=float, default=55.0)   # keep epochs above this delta pctile
    ap.add_argument("--block-min", type=float, default=10.0)  # max download block length (<=10 min portal cap)
    ap.add_argument("--night-start-h", type=float, default=None)  # skip the scan if known
    ap.add_argument("--out", default="HUP165_night1")
    a = ap.parse_args()
    outdir = os.path.join(ROOT, "data", "ieeg_portal", a.out); os.makedirs(outdir, exist_ok=True)

    s = sess(); ds = s.open_dataset(a.dataset)
    labels = ds.get_channel_labels()
    sf = ds.get_time_series_details(labels[0]).sample_rate
    total_h = ds.get_time_series_details(labels[0]).duration / 3.6e9
    lab_idx = {l: i for i, l in enumerate(labels)}

    def contact(l):
        m = re.match(r"[A-Za-z]+(\d+)$", l)
        return int(m.group(1)) if m else 999
    mtl = [l for l in labels if re.match(a.mtl_regex, l) and contact(l) <= a.max_contact]
    mtl_idx = [lab_idx[l] for l in mtl]
    stage_ch = next((l for l in ["LB4", "LC4", "LB6", "LH4", "LC6"] if l in lab_idx), mtl[0])
    print(f"[night] {a.dataset}: {sf:.0f} Hz, {total_h:.0f} h. MTL={len(mtl)} ch {mtl[:8]}...; "
          f"stage on {stage_ch}", flush=True)

    # ---- A. night finder (skippable) ----
    if a.night_start_h is not None:
        best_t = a.night_start_h * 3600
        print(f"[night] using provided night start = {a.night_start_h:.1f} h (scan skipped)", flush=True)
    else:
        step = a.scan_step_min * 60
        probes = int(a.scan_h * 3600 / step)
        tr = []
        for k in range(probes):
            t = k * step
            try:
                x = get(ds, [lab_idx[stage_ch]], t, 6.0)[:, 0]
                tr.append((t, delta_ratio(x, sf)))
            except Exception:
                tr.append((t, np.nan))
        tr = pd.DataFrame(tr, columns=["t", "dr"]).dropna()
        tr.to_csv(os.path.join(outdir, "scan_delta.csv"), index=False)
        win = int(a.night_h * 3600 / step)
        best_t, best_m = 0, -1
        for i in range(len(tr) - win):
            m = tr.iloc[i:i + win].dr.mean()
            if m > best_m:
                best_m, best_t = m, tr.t.iloc[i]
        print(f"[night] chosen night start = {best_t/3600:.1f} h (mean delta-ratio {best_m:.2f})", flush=True)

    # ---- B. stage the night (1 channel, full window, 30 s epochs) ----
    n_ep = int(a.night_h * 3600 / EPOCH)
    drs = np.full(n_ep, np.nan)
    chunk = MAX_REQ_S                      # 10-min chunks (portal cap), epoch offline
    n_chunks = int(np.ceil(a.night_h * 3600 / chunk))
    for ci in range(n_chunks):
        cstart = best_t + ci * chunk
        cdur = min(chunk, best_t + a.night_h * 3600 - cstart)
        try:
            x = get(ds, [lab_idx[stage_ch]], cstart, cdur)[:, 0]
        except Exception:
            continue
        for k in range(int(cdur / EPOCH)):
            gi = int(ci * chunk / EPOCH) + k
            if gi < n_ep:
                seg = x[int(k * EPOCH * sf):int((k + 1) * EPOCH * sf)]
                drs[gi] = delta_ratio(seg, sf)
    valid = drs[np.isfinite(drs)]
    thr = np.percentile(valid, a.nrem_pct) if len(valid) else np.nan
    nrem = np.isfinite(drs) & (drs >= thr)
    print(f"[night] staged {n_ep} epochs; NREM threshold(delta_ratio)={thr:.2f}; "
          f"NREM epochs={nrem.sum()} ({100*nrem.mean():.0f}%), ~{nrem.sum()*EPOCH/3600:.1f} h", flush=True)
    pd.DataFrame({"epoch": range(n_ep), "start_s": best_t + np.arange(n_ep) * EPOCH,
                  "delta_ratio": drs, "nrem": nrem}).to_csv(os.path.join(outdir, "hypnogram.csv"), index=False)

    # ---- C. contiguous NREM blocks -> pull MTL channels, save per block ----
    blocks = []
    e = 0
    while e < n_ep:
        if nrem[e]:
            j = e
            while j < n_ep and nrem[j] and (j - e) * EPOCH < a.block_min * 60:
                j += 1
            blocks.append((best_t + e * EPOCH, (j - e) * EPOCH)); e = j
        else:
            e += 1
    print(f"[night] {len(blocks)} NREM blocks to download on {len(mtl)} MTL channels", flush=True)
    idxrows = []
    for bi, (t0, dur) in enumerate(blocks):
        try:
            data = get(ds, mtl_idx, t0, dur).astype(np.float32)   # samples x ch
        except Exception as ex:
            print(f"  block {bi} FAILED: {ex}", flush=True); continue
        fn = f"block_{bi:04d}.npz"
        np.savez_compressed(os.path.join(outdir, fn), data=data, ch_names=mtl, sfreq=sf, start_sec=t0)
        idxrows.append({"block": bi, "file": fn, "start_sec": t0, "dur_sec": dur,
                        "n_samples": data.shape[0]})
        if bi % 5 == 0:
            print(f"  saved block {bi+1}/{len(blocks)} ({data.shape}) t={t0/3600:.2f}h", flush=True)
    pd.DataFrame(idxrows).to_csv(os.path.join(outdir, "index.csv"), index=False)
    json.dump({"dataset": a.dataset, "sfreq": sf, "night_start_h": best_t / 3600,
               "night_h": a.night_h, "mtl_channels": mtl, "stage_ch": stage_ch,
               "nrem_threshold": float(thr), "n_blocks": len(idxrows),
               "nrem_hours": float(nrem.sum() * EPOCH / 3600)},
              open(os.path.join(outdir, "meta.json"), "w"), indent=2)
    tot = sum(r["n_samples"] for r in idxrows) / sf / 3600
    print(f"[night] DONE: {len(idxrows)} blocks, {tot:.2f} h NREM on {len(mtl)} ch -> {outdir}", flush=True)
    s.close_dataset(ds)


if __name__ == "__main__":
    main()
