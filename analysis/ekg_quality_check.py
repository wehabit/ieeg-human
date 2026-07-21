"""
EKG-quality check on one HUP night (default HUP165) — de-risks tests 3A/3B before scaling.

Pulls EKG1/EKG2 from iEEG.org for the LONGEST already-staged NREM block of a pulled night,
runs a self-contained Pan-Tompkins R-peak detector (no neurokit2 dependency), and reports
whether the clinical-rig EKG is clean enough to build an RR series:
  - HR range / median, RR distribution
  - artifact fraction (physiologically implausible RR: <0.33 s [>180 bpm] or >1.5 s [<40 bpm])
  - EKG1 vs EKG2 vs bipolar(EKG1-EKG2) beat-count agreement
Figure: EKG trace + detected R-peaks (20 s) + RR tachogram + HR histogram.

Usage: .venv/bin/python analysis/ekg_quality_check.py [--night HUP165_night1] [--dataset HUP165_phaseII]
"""
import argparse, csv, json, os, time
import numpy as np
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ieeg.auth import Session

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CRED = os.path.join(ROOT, "data", "ieeg_secret", "credentials.json")
MAX_REQ_S = 600
LINE_HZ = 60.0   # Penn (US) mains


def sess():
    c = json.load(open(CRED))
    return Session(c["username"], c["password"])


def get(ds, idx, start_s, dur_s, tries=3):
    dur_s = min(dur_s, MAX_REQ_S)
    for k in range(tries):
        try:
            return ds.get_data(int(start_s * 1e6), int(dur_s * 1e6), idx)
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(1.5)


def pan_tompkins(x, sf):
    """Return R-peak sample indices. Classic PT: bandpass 5-15, derivative, square,
    moving-window integrate, adaptive peak pick."""
    x = np.nan_to_num(x.astype(float))
    # notch mains + harmonic
    for f0 in (LINE_HZ, 2 * LINE_HZ):
        if f0 < sf / 2:
            b, a = signal.iirnotch(f0, 30, sf)
            x = signal.filtfilt(b, a, x)
    # bandpass 5-15 Hz (QRS energy)
    b, a = signal.butter(2, [5 / (sf / 2), 15 / (sf / 2)], btype="band")
    xf = signal.filtfilt(b, a, x)
    deriv = np.gradient(xf)
    sq = deriv ** 2
    win = max(1, int(0.15 * sf))
    integ = np.convolve(sq, np.ones(win) / win, mode="same")
    # adaptive threshold; enforce 200 ms refractory (max ~300 bpm ceiling for detection)
    thr = 0.3 * np.mean(integ) + 0.2 * np.median(integ)
    peaks, _ = signal.find_peaks(integ, height=thr, distance=int(0.30 * sf))
    # refine each peak to local max of bandpassed signal within +/-100 ms
    r = []
    w = int(0.10 * sf)
    for p in peaks:
        lo, hi = max(0, p - w), min(len(xf), p + w)
        r.append(lo + int(np.argmax(np.abs(xf[lo:hi]))))
    return np.array(sorted(set(r))), integ, xf


def beats_and_rr(x, sf):
    r, _, _ = pan_tompkins(x, sf)
    rr = np.diff(r) / sf
    return r, rr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--night", default="HUP165_night1")
    ap.add_argument("--dataset", default="HUP165_phaseII")
    a = ap.parse_args()

    nd = os.path.join(ROOT, "data", "ieeg_portal", a.night)
    rows = list(csv.DictReader(open(os.path.join(nd, "index.csv"))))
    blk = max(rows, key=lambda r: float(r["dur_sec"]))
    start_s, dur_s = float(blk["start_sec"]), float(blk["dur_sec"])
    print(f"[ekg] {a.dataset}: longest NREM block = {dur_s:.0f}s @ t={start_s:.0f}s (block {blk['block']})", flush=True)

    s = sess(); ds = s.open_dataset(a.dataset)
    labels = ds.get_channel_labels()
    ekg = [l for l in labels if l.upper().replace(" ", "").startswith(("EKG", "ECG")) or "EKG" in l.upper() or "ECG" in l.upper()]
    print(f"[ekg] EKG/ECG channels: {ekg}", flush=True)
    if not ekg:
        raise SystemExit("No EKG channel found in dataset labels.")
    idx = [labels.index(l) for l in ekg[:2]]
    sf = ds.get_time_series_details(labels[idx[0]]).sample_rate

    data = get(ds, idx, start_s, min(dur_s, MAX_REQ_S))   # samples x 2
    print(f"[ekg] pulled {data.shape} @ {sf:.0f} Hz", flush=True)

    chans = {ekg[0]: data[:, 0]}
    if data.shape[1] > 1:
        chans[ekg[1]] = data[:, 1]
        chans[f"{ekg[0]}-{ekg[1]}"] = data[:, 0] - data[:, 1]

    summary = {}
    best = None
    for name, x in chans.items():
        r, rr = beats_and_rr(x, sf)
        if len(rr) < 3:
            summary[name] = dict(n_beats=len(r), note="too few beats")
            continue
        bad = ((rr < 0.33) | (rr > 1.5)).mean()
        hr = 60.0 / rr
        summary[name] = dict(
            n_beats=int(len(r)),
            hr_median=float(np.median(hr)),
            hr_iqr=[float(np.percentile(hr, 25)), float(np.percentile(hr, 75))],
            rr_median_s=float(np.median(rr)),
            artifact_frac=float(bad),
            mean_beats_per_min=float(len(r) / (dur_s / 60.0)),
        )
        if best is None or bad < summary[best]["artifact_frac"]:
            best = name
    print("[ekg] per-channel summary:")
    print(json.dumps(summary, indent=2))
    print(f"[ekg] BEST channel = {best} (lowest artifact fraction)", flush=True)

    # ---- figure ----
    xb = chans[best]
    r, integ, xf = pan_tompkins(xb, sf)
    rr = np.diff(r) / sf
    hr = 60.0 / rr
    t = np.arange(len(xb)) / sf

    fig, ax = plt.subplots(3, 1, figsize=(11, 9))
    seg = (t >= 20) & (t < 40)                       # 20 s window
    ax[0].plot(t[seg], xf[seg], lw=0.7, color="0.2")
    rin = r[(r / sf >= 20) & (r / sf < 40)]
    ax[0].plot(rin / sf, xf[rin], "rv", ms=6)
    ax[0].set(title=f"{a.dataset}  EKG={best}  (bandpassed 5-15 Hz, 20 s)", xlabel="s", ylabel="a.u.")
    ax[1].plot(r[1:] / sf, rr, ".-", ms=3, lw=0.6)
    ax[1].axhspan(0.33, 1.5, color="g", alpha=0.06)
    ax[1].set(title=f"RR tachogram  ({len(r)} beats, {len(r)/(dur_s/60):.0f} bpm avg)", xlabel="s", ylabel="RR (s)")
    ax[2].hist(hr, bins=40, color="0.4")
    ax[2].set(title=f"HR distribution  median={np.median(hr):.0f} bpm, artifact={((rr<0.33)|(rr>1.5)).mean()*100:.1f}%",
              xlabel="bpm", ylabel="count")
    fig.tight_layout()
    outdir = os.path.join(ROOT, "outputs", "ekg_quality_check")
    os.makedirs(outdir, exist_ok=True)
    fp = os.path.join(outdir, f"{a.dataset}_ekg_quality.png")
    fig.savefig(fp, dpi=130)
    json.dump({"dataset": a.dataset, "block": blk, "sfreq": sf, "best_channel": best, "summary": summary},
              open(os.path.join(outdir, f"{a.dataset}_ekg_quality.json"), "w"), indent=2)
    print(f"[ekg] figure -> {fp}", flush=True)


if __name__ == "__main__":
    main()
