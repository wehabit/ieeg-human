"""
LEGACY 3A — retained only for shared acquisition/filter helpers.

The executable analysis is quarantined because it averages raw contacts before power, uses a
single-bin coherence endpoint, and lacks corrected cache/version controls. Use
``cache_lc_series.py`` followed by ``lecci_faithful_3A.py`` for the current approximation.

The primary LC-fingerprint measure (Lecci 2017 / Osorio-Forero 2021): during NREM, the ~0.02 Hz
(~50 s) rhythm in sigma (spindle) power co-varies with the ~0.02 Hz rhythm in heart rate.

METHOD (validated against positive/negative controls -- see docs/DEC_3A_METHOD_NOTES.md):
  EKG  -> NeuroKit2 R-peaks (chunked) -> RR -> drop impossible RR -> 4 Hz instantaneous-HR grid
  MTL  -> notch -> sigma (11-16 Hz) -> Hilbert envelope -> mean across channels -> 4 Hz grid
  both -> high-pass 0.005 Hz (remove drift only; NOT band-passed to the infraslow band)
  coherence = magnitude-squared coherence (Welch) read at 0.02 Hz and as the 0.01-0.03 Hz max
  significance = ANALYTIC MSC threshold  crit = 1 - alpha^(1/(K-1))  for K Welch segments
                 (Bonferroni over band bins for the band-max)

Three errors this implementation deliberately avoids:
  1. ds.get_data returns columns in ASCENDING channel-index order, not the requested order.
  2. Band-passing to 0.01-0.03 Hz *before* coherence destroys the estimate (a known-coupled
     control pair scored 0.09 instead of 0.98).
  3. Circular-shift / phase-randomised surrogates are INVALID for oscillatory coupling: a shifted
     50 s rhythm is still a 50 s rhythm and MSC ignores a constant phase offset.

Outputs: outputs/infraslow_rr_sigma_coherence/<ds>_3A.{json,png,svg}

Usage: .venv/bin/python analysis/infraslow_rr_sigma_coherence.py \
         --dataset HUP165_phaseII --night HUP165_night1 --win-start-s 41280 --win-min 140
"""
import argparse, json, os, time
import numpy as np
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import neurokit2 as nk
from ieeg.auth import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CRED = os.path.join(ROOT, "data", "ieeg_secret", "credentials.json")
MAX_REQ_S = 600
LINE_HZ = 60.0
SIGMA = (11.0, 16.0)
INFRA = (0.01, 0.03)      # infraslow band of interest
F_TARGET = 0.02           # the ~50 s LC rhythm
FS_IS = 4.0               # analysis grid (Hz)
DRIFT_HP = 0.005          # high-pass to remove drift (keeps the whole infraslow band)
EDGE_TRIM_S = 200.0       # drop high-pass transient
EPOCH = 30.0
NREM_DR = 0.90            # absolute delta-ratio threshold for NREM (matches night-staging)
ALPHA = 0.05
IEEG_CONNECT_TIMEOUT_S = 10.0
IEEG_READ_TIMEOUT_S = 90.0


def configure_http_session(http, timeout=(IEEG_CONNECT_TIMEOUT_S, IEEG_READ_TIMEOUT_S)):
    """Bound portal requests and retry transient metadata failures.

    ``ieeg==1.6`` does not supply a Requests timeout, so one stalled response can otherwise block
    an entire cohort regeneration indefinitely. Data POSTs are already retried explicitly by
    :func:`get`; the adapter retry policy therefore remains limited to idempotent GET metadata
    requests.
    """
    original_request = http.request

    def request_with_timeout(method, url, **kwargs):
        kwargs.setdefault("timeout", timeout)
        return original_request(method, url, **kwargs)

    http.request = request_with_timeout
    retry = Retry(
        total=2,
        connect=2,
        read=2,
        status=2,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(("GET",)),
    )
    adapter = HTTPAdapter(max_retries=retry)
    http.mount("https://", adapter)
    http.mount("http://", adapter)
    return http


def sess():
    c = json.load(open(CRED))
    session = Session(c["username"], c["password"])
    configure_http_session(session.api.http)
    return session


def get(ds, idx, start_s, dur_s, tries=3):
    dur_s = min(dur_s, MAX_REQ_S)
    for k in range(tries):
        try:
            return ds.get_data(int(start_s * 1e6), int(dur_s * 1e6), idx)
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(1.5)


def pull_continuous(ds, idx, start_s, dur_s):
    """Pull [start_s, start_s+dur_s) for channel indices idx across 600 s chunks.
    ds.get_data returns columns in ASCENDING channel-index order, so remap to requested order."""
    out, t = [], 0.0
    while t < dur_s:
        d = min(MAX_REQ_S, dur_s - t)
        out.append(get(ds, idx, start_s + t, d))
        t += d
    data = np.concatenate(out, axis=0)
    order = np.argsort(idx)
    inv = np.empty(len(idx), dtype=int); inv[order] = np.arange(len(idx))
    return data[:, inv]


def notch(x, sf):
    y = x.astype(float)
    for f0 in (LINE_HZ, 2 * LINE_HZ):
        if f0 < sf / 2:
            b, a = signal.iirnotch(f0, 30, sf)
            y = signal.filtfilt(b, a, y)
    return y


def r_peaks(x, sf, chunk_s=300):
    """NeuroKit2 R-peaks, run in CHUNKS (its artifact correction misbehaves on
    multi-million-sample arrays; per-chunk it is reliable)."""
    sf = int(sf); step = int(chunk_s * sf); out = []
    for st in range(0, len(x), step):
        seg = np.nan_to_num(x[st:st + step].astype(float))
        if len(seg) < sf * 5:
            continue
        cleaned = nk.ecg_clean(seg, sampling_rate=sf, method="neurokit")
        _, info = nk.ecg_peaks(cleaned, sampling_rate=sf, method="neurokit", correct_artifacts=True)
        out.append(np.asarray(info["ECG_R_Peaks"], dtype=int) + st)
    return np.concatenate(out) if out else np.array([], dtype=int)


def drift_highpass(x):
    sos = signal.butter(3, DRIFT_HP, btype="high", fs=FS_IS, output="sos")
    return signal.sosfiltfilt(sos, x)


def bandpass_infra(x):
    """ONLY for the visual overlay -- never fed to the coherence estimate."""
    sos = signal.butter(3, [INFRA[0], INFRA[1]], btype="band", fs=FS_IS, output="sos")
    return signal.sosfiltfilt(sos, x)


def hr_series(x_ekg, sf, t_grid):
    """EKG -> instantaneous HR on t_grid; returns (broadband detrended, infraslow, QC)."""
    r = r_peaks(x_ekg, sf)
    bt = r / sf
    rr = np.diff(bt)
    bad = (rr < 0.33) | (rr > 1.5)      # physiologically impossible only
    hr = np.interp(t_grid, bt[1:][~bad], 60.0 / rr[~bad])
    hr_bb = drift_highpass(hr - hr.mean())
    return hr_bb, bandpass_infra(hr_bb), float(bad.mean()), int(len(r)), float(np.median(60.0 / rr[~bad]))


def sigma_series(x_mtl, sf, t_grid):
    """MTL channels -> mean sigma Hilbert envelope on t_grid; (broadband detrended, infraslow)."""
    envs = []
    for j in range(x_mtl.shape[1]):
        y = notch(np.nan_to_num(x_mtl[:, j].astype(float)), sf)
        b, a = signal.butter(4, [SIGMA[0] / (sf / 2), SIGMA[1] / (sf / 2)], btype="band")
        env = np.abs(signal.hilbert(signal.filtfilt(b, a, y)))
        envs.append(env / (np.median(env) + 1e-12))
    env = np.mean(envs, axis=0)
    env_g = np.interp(t_grid, np.arange(len(env)) / sf, env)
    sig_bb = drift_highpass(env_g - env_g.mean())
    return sig_bb, bandpass_infra(sig_bb), env_g


def coherence_analytic(a, b):
    """Broadband MSC + analytic significance threshold (no surrogates)."""
    nperseg = int(min(2048, (len(a) // 4) // 2 * 2))
    f, cxy = signal.coherence(a, b, fs=FS_IS, nperseg=nperseg, noverlap=nperseg // 2)
    K = int((len(a) - nperseg // 2) // (nperseg // 2))          # ~independent Welch segments
    band = (f >= INFRA[0]) & (f <= INFRA[1]); nb = int(band.sum())
    at_t = float(cxy[np.argmin(np.abs(f - F_TARGET))])
    bmax = float(cxy[band].max()); fpk = float(f[band][np.argmax(cxy[band])])
    crit_bin = 1 - ALPHA ** (1 / (K - 1))
    crit_band = 1 - (ALPHA / nb) ** (1 / (K - 1))
    return dict(f=f, cxy=cxy, K=K, nbins=nb, at_target=at_t, band_max=bmax,
                band_peak_hz=fpk, crit_bin=crit_bin, crit_band=crit_band,
                sig_at_target=bool(at_t > crit_bin), sig_band=bool(bmax > crit_band), nperseg=nperseg)


def nrem_fraction(x_stage, sf, dur):
    n_ep = int(dur / EPOCH); dr = np.full(n_ep, np.nan)
    for e in range(n_ep):
        seg = x_stage[int(e * EPOCH * sf):int((e + 1) * EPOCH * sf)]
        if len(seg) < sf:
            continue
        f, p = signal.welch(seg, sf, nperseg=int(min(4 * sf, len(seg))))
        num = np.trapezoid(p[(f >= 0.5) & (f < 4)], f[(f >= 0.5) & (f < 4)])
        den = np.trapezoid(p[(f >= 0.5) & (f < 25)], f[(f >= 0.5) & (f < 25)])
        dr[e] = num / den if den > 0 else np.nan
    return float(np.nanmean(dr >= NREM_DR)), float(np.nanmedian(dr))


def main():
    raise SystemExit(
        "LEGACY/WITHDRAWN 3A entry point: use cache_lc_series.py followed by "
        "lecci_faithful_3A.py. Shared helper functions remain importable.")
    raise SystemExit(
        "LEGACY 3A QUARANTINED: use cache_lc_series.py and lecci_faithful_3A.py. "
        "Only helper functions from this module remain supported.")
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="HUP165_phaseII")
    ap.add_argument("--night", default="HUP165_night1")
    ap.add_argument("--win-start-s", type=float, default=None)
    ap.add_argument("--win-min", type=float, default=140.0)
    ap.add_argument("--mtl-chans", default="LB2,LB4,LC2,LC4,LH2,LH4")
    ap.add_argument("--ekg", default="EKG1")
    a = ap.parse_args()

    meta = json.load(open(os.path.join(ROOT, "data", "ieeg_portal", a.night, "meta.json")))
    win_start = a.win_start_s if a.win_start_s is not None else meta["night_start_h"] * 3600.0
    win_s = a.win_min * 60.0

    s = sess(); ds = s.open_dataset(a.dataset)
    labels = ds.get_channel_labels()
    sf = ds.get_time_series_details(labels[0]).sample_rate
    want = [c for c in a.mtl_chans.split(",") if c in labels]
    ekg = a.ekg if a.ekg in labels else next(l for l in labels if l.upper().startswith(("EKG", "ECG")))
    idx = [labels.index(c) for c in want] + [labels.index(ekg)]
    print(f"[3A] {a.dataset} @ {sf:.0f} Hz | MTL={want} | EKG={ekg} | "
          f"{win_start/3600:.2f}-{(win_start+win_s)/3600:.2f} h ({a.win_min:.0f} min)", flush=True)

    data = pull_continuous(ds, idx, win_start, win_s)
    x_mtl, x_ekg = data[:, :len(want)], data[:, len(want)]
    dur = len(x_ekg) / sf
    t_grid = np.arange(0, dur, 1.0 / FS_IS)

    hr_bb, hr_is, frac_bad, nbeats, med_hr = hr_series(x_ekg, sf, t_grid)
    sig_bb, sig_is, env_g = sigma_series(x_mtl, sf, t_grid)
    print(f"[3A] {nbeats} beats, median HR {med_hr:.1f} bpm, {frac_bad*100:.2f}% RR dropped", flush=True)

    k = int(EDGE_TRIM_S * FS_IS)
    hr_bb, sig_bb, hr_is, sig_is, t_t = hr_bb[k:-k], sig_bb[k:-k], hr_is[k:-k], sig_is[k:-k], t_grid[k:-k]

    C = coherence_analytic(hr_bb, sig_bb)

    # lag from infraslow cross-correlation (+ => sigma follows HR)
    az = (hr_is - hr_is.mean()) / hr_is.std(); bz = (sig_is - sig_is.mean()) / sig_is.std()
    xc = signal.correlate(bz, az, mode="full") / len(az)
    lags = signal.correlation_lags(len(bz), len(az), mode="full") / FS_IS
    m = np.abs(lags) <= 100; i_pk = int(np.argmax(xc[m]))

    nrem_frac, dr_med = nrem_fraction(x_mtl[:, 0], sf, dur)

    out = dict(dataset=a.dataset, window_h=[win_start / 3600, (win_start + win_s) / 3600],
               dur_min=dur / 60, mtl_chans=want, ekg=ekg, sf=sf,
               nbeats=nbeats, median_hr_bpm=med_hr, rr_dropped_frac=frac_bad,
               nrem_frac=nrem_frac, delta_ratio_median=dr_med,
               coherence_at_0p02Hz=C["at_target"], crit_bin=C["crit_bin"],
               significant_at_0p02Hz=C["sig_at_target"],
               coherence_band_max=C["band_max"], band_peak_hz=C["band_peak_hz"],
               crit_band=C["crit_band"], significant_band=C["sig_band"],
               welch_segments_K=C["K"], nperseg=C["nperseg"],
               xcorr_peak=float(xc[m][i_pk]), xcorr_lag_s=float(lags[m][i_pk]))
    outdir = os.path.join(ROOT, "outputs", "infraslow_rr_sigma_coherence")
    os.makedirs(outdir, exist_ok=True)
    json.dump(out, open(os.path.join(outdir, f"{a.dataset}_3A.json"), "w"), indent=2)
    print("[3A] RESULT:\n" + json.dumps(out, indent=2), flush=True)

    # ---- figure ----
    f, cxy = C["f"], C["cxy"]
    fig, ax = plt.subplots(3, 1, figsize=(11.5, 10))
    ax[0].plot(t_t / 60, bz, color="#1b6ca8", lw=1.0, label="sigma-power infraslow (11-16 Hz env)")
    ax[0].plot(t_t / 60, az, color="#c0392b", lw=1.0, label="heart-rate infraslow")
    ax[0].set(title=f"{a.dataset}  infraslow overlay (0.01-0.03 Hz, display only)  |  "
                    f"lag={out['xcorr_lag_s']:+.0f}s  r={out['xcorr_peak']:.2f}", ylabel="z", xlabel="min")
    ax[0].legend(fontsize=8, loc="upper right")

    sig_txt = "SIGNIFICANT" if C["sig_at_target"] else "not significant"
    ax[1].plot(f * 1000, cxy, color="0.2", lw=1.5)
    ax[1].axvspan(INFRA[0] * 1000, INFRA[1] * 1000, color="g", alpha=0.10, label="infraslow band")
    ax[1].axhline(C["crit_bin"], color="r", ls="--", lw=1, label=f"analytic crit @bin ({C['crit_bin']:.3f})")
    ax[1].axhline(C["crit_band"], color="darkorange", ls=":", lw=1, label=f"crit band-corrected ({C['crit_band']:.3f})")
    ax[1].plot(F_TARGET * 1000, C["at_target"], "v", color="#c0392b", ms=10)
    ax[1].set(title=f"MSC (RR x sigma)  |  @0.02 Hz = {C['at_target']:.3f} -> {sig_txt}  "
                    f"(K={C['K']} segments)", xlabel="mHz", ylabel="coherence", xlim=[0, 60])
    ax[1].legend(fontsize=8)

    for sgl, lab, col in ((sig_bb, "sigma envelope", "#1b6ca8"), (hr_bb, "heart rate", "#c0392b")):
        fp, pp = signal.welch(sgl, fs=FS_IS, nperseg=C["nperseg"])
        ax[2].semilogy(fp * 1000, pp / pp.max(), color=col, lw=1.3, label=lab)
    ax[2].axvspan(INFRA[0] * 1000, INFRA[1] * 1000, color="g", alpha=0.10)
    ax[2].set(title="power spectra (normalised) — is there a ~0.02 Hz rhythm in each signal at all?",
              xlabel="mHz", ylabel="norm. power", xlim=[0, 60])
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, f"{a.dataset}_3A.{ext}"), dpi=130)
    print(f"[3A] figures -> {outdir}", flush=True)


if __name__ == "__main__":
    raise SystemExit(
        "LEGACY/WITHDRAWN 3A entry point: use cache_lc_series.py followed by "
        "lecci_faithful_3A.py. Shared helper functions remain importable.")
