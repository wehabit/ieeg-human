"""OpenNeuro ds003848 (Utrecht RESPect long-term iEEG) -> staged derived-series cache.

The independent replication cohort for 3A/3B. Unlike HUP, every subject carries EMG + EOG, so NREM
can be scored with real REM/wake exclusion instead of the GMM-on-slow-wave proxy -- the single
biggest threat to the HUP 3A negative (residual wake/REM is exactly the broadband low-frequency
structure that could smear a real infraslow rhythm).

Six patients, ~1 h continuous `task-[Ss]leep` runs @ 2048 Hz, 50 Hz line. Verified against the raw
channels.tsv (2026-07): all six have iEEG (3 ECoG grid + 3 SEEG depth), ECG, EMG, EOG, and (bad)
respiration belts. Channel ROLE is taken from channels.tsv by row index, which BIDS guarantees
matches the data-column order, so naming variants (ECG+/ecg1+, emg+/EMG2, orb+/Orb+) don't matter.

Pipeline per subject: download the BrainVision triple if absent -> stream in chunks with MNE ->
derive the SAME series the HUP cache stores (so lecci_faithful_3A / event_3B_cached consume it
unchanged) PLUS EMG/EOG staging features and a scored `stage_lab` -> save data/derived/ds003848/.
The raw .eeg (~3.9 GB) is optionally deleted after caching so disk stays bounded.

STAGING (rule-based, honest about its limits). True AASM scoring needs scalp EEG, which this dataset
lacks. What EMG + EOG buy is REM/Wake exclusion, and that is what this does:
  * per 30 s epoch: submental EMG RMS (notch + 10-100 Hz), EOG movement variance (0.3-6 Hz),
    slow-wave power (0.5-4 Hz, from the iEEG mean), all robust-z-scored within subject
  * Wake  = high muscle tone (EMG z > 1.0)
  * REM   = muscle atonia (EMG z < -0.3) + low SWA + phasic eye movement (EOG z > 0.5)
  * NREM  = everything else with adequate delta; split N2/N3 by a 2-component GMM on log SWA
The N2/N3 split is still SWA-driven (as in Lecci's S2/SWS axis), but now on epochs from which Wake
and REM have been removed -- which the HUP proxy could not do.

    .venv/bin/python analysis/stage_ds003848.py [--subjects sub-RESP0521,...] [--keep-raw]
"""
import argparse, os, time, urllib.request
import numpy as np
import csv, io
from scipy import signal, interpolate
import neurokit2 as nk
import mne

from infraslow_rr_sigma_coherence import ROOT
from cohort_stages_3ABD import (fsp_from, band_sos, delta_ratio, EPOCH, SWA_BAND, CHUNK_S)
from cache_lc_series import detect_so_halfwaves, SIGMA_FIXED, SWA_BAND_L, SO_BAND_NAJI, FS_RR
from results_3A_tutorial_style import ied_clean_mask

mne.set_log_level("ERROR")

BASE = "https://openneuro.org/crn/datasets/ds003848/snapshots/1.0.3/files"
RAW = os.path.join(ROOT, "data", "ds003848_raw")
OUT = os.path.join(ROOT, "data", "derived", "ds003848")

# subject -> (session, sleep-run basename without extension). One representative sleep run each.
SUBJECTS = {
    "sub-RESP0521": ("ses-1", "sub-RESP0521_ses-1_task-Sleep_run-030344"),
    "sub-RESP0699": ("ses-1", "sub-RESP0699_ses-1_task-sleep_run-030608"),
    "sub-RESP0724": ("ses-1", "sub-RESP0724_ses-1_task-sleep_run-022151"),
    "sub-RESP0749": ("ses-1", "sub-RESP0749_ses-1_task-Sleep_run-030241"),
    "sub-RESP0779": ("ses-1", "sub-RESP0779_ses-1_task-sleep_run-030241"),
    "sub-RESP0800": ("ses-1", "sub-RESP0800_ses-1_task-Sleep_run-030017"),
}


# ---------------------------------------------------------------- download
def _dl(sub, ses, base, ext):
    url = f"{BASE}/{sub}:{ses}:ieeg:{base}{ext}"
    dst = os.path.join(RAW, f"{base}{ext}")
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        return dst
    os.makedirs(RAW, exist_ok=True)
    tmp = dst + ".part"
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as f:
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            f.write(b)
    os.replace(tmp, dst)
    return dst


def ensure_files(sub, ses, base):
    paths = {}
    for ext in ("_channels.tsv", "_ieeg.json", "_ieeg.vhdr", "_ieeg.vmrk", "_ieeg.eeg"):
        paths[ext] = _dl(sub, ses, base, ext)
    return paths


# ---------------------------------------------------------------- channel roles
def channel_roles(channels_tsv):
    """Return dict of role -> list of row indices (0-based, = data column order), for good channels."""
    rows = list(csv.DictReader(open(channels_tsv), delimiter="\t"))
    roles = dict(ieeg=[], ecg=[], emg=[], eog=[], names=[r["name"] for r in rows])
    for i, r in enumerate(rows):
        t = (r.get("type") or "").upper()
        bad = (r.get("status") or "").lower() == "bad"
        if t in ("SEEG", "ECOG") and not bad:
            roles["ieeg"].append(i)
        elif t in ("ECG", "EKG"):
            roles["ecg"].append(i)
        elif t == "EMG" and not bad:
            roles["emg"].append(i)
        elif t == "EOG" and not bad:
            roles["eog"].append(i)
    roles["n_rows"] = len(rows)
    return roles


# ---------------------------------------------------------------- staging
def robust_z(x):
    x = np.asarray(x, float)
    m = np.nanmedian(x)
    s = np.nanmedian(np.abs(x - m)) * 1.4826 + 1e-12
    return (x - m) / s


def score_stages(ep_swa, ep_emg, ep_eog, ep_dr, ep_clean):
    """Rule-based Wake/REM/N1/N2/N3 from EMG tone, EOG movement, and iEEG slow-wave power."""
    n = len(ep_swa)
    lab = np.full(n, "", dtype=object)
    swa_z = robust_z(np.log(np.clip(ep_swa, 1e-12, None)))
    emg_z = robust_z(np.log(np.clip(ep_emg, 1e-12, None)))
    eog_z = robust_z(np.log(np.clip(ep_eog, 1e-12, None)))
    valid = np.isfinite(swa_z) & np.isfinite(emg_z) & (ep_clean >= 0.5)

    wake = valid & (emg_z > 1.0)                                   # high muscle tone
    rem = valid & ~wake & (emg_z < -0.3) & (swa_z < 0.0) & (eog_z > 0.5)   # atonia + phasic eyes + low SWA
    nrem = valid & ~wake & ~rem
    lab[wake] = "W"
    lab[rem] = "R"
    if nrem.sum() >= 40:
        # N2 vs N3 by 2-component GMM on log SWA within genuine NREM (Lecci S2/SWS axis)
        try:
            from sklearn.mixture import GaussianMixture
            v = swa_z[nrem].reshape(-1, 1)
            g = GaussianMixture(n_components=2, random_state=0, n_init=3).fit(v)
            hi = int(np.argmax(g.means_.ravel()))
            cl = g.predict(v)
            names = np.where(cl == hi, "N3", "N2")
        except Exception:
            thr = np.median(swa_z[nrem])
            names = np.where(swa_z[nrem] >= thr, "N3", "N2")
        lab[np.where(nrem)[0]] = names
    return lab, dict(n_wake=int(wake.sum()), n_rem=int(rem.sum()), n_nrem=int(nrem.sum()),
                     n_valid=int(valid.sum()))


# ---------------------------------------------------------------- main derive
def run(subject, keep_raw=False, force=False):
    fp = os.path.join(OUT, f"{subject}.npz")
    if os.path.exists(fp) and not force:
        print(f"[{subject}] cached", flush=True); return
    ses, base = SUBJECTS[subject]
    t0 = time.time()
    print(f"[{subject}] fetching {base} ...", flush=True)
    paths = ensure_files(subject, ses, base)
    roles = channel_roles(paths["_channels.tsv"])
    if not roles["ieeg"] or not roles["ecg"]:
        np.savez_compressed(fp, status="skip", reason="no iEEG or ECG")
        print(f"[{subject}] SKIP no iEEG/ECG", flush=True); return

    raw = mne.io.read_raw_brainvision(paths["_ieeg.vhdr"], preload=False, verbose="ERROR")
    sf = float(raw.info["sfreq"])
    n_samp = raw.n_times
    if len(raw.ch_names) != roles["n_rows"]:
        print(f"[{subject}] WARN channel count mismatch: raw {len(raw.ch_names)} vs tsv {roles['n_rows']}",
              flush=True)
    ie = roles["ieeg"]; ecg_i = roles["ecg"][0]
    emg_i = roles["emg"][0] if roles["emg"] else None
    eog_i = roles["eog"][0] if roles["eog"] else None
    names = roles["names"]
    ctx = [names[i] for i in ie]

    total_s = int(n_samp / sf)
    n_ep = int(total_s // EPOCH)

    # sample a 300 s window for the fast-spindle peak
    s0 = int(min(300, total_s * 0.2) * sf)
    samp = raw.get_data(picks=ie, start=s0, stop=min(s0 + int(300 * sf), n_samp)).T * 1e6  # V->uV
    fsp, fsp_real = fsp_from(samp, sf)
    print(f"[{subject}] {sf:.0f} Hz, {len(ie)} iEEG, ECG@{ecg_i}, EMG@{emg_i}, EOG@{eog_i}, "
          f"FSP {fsp:.2f} Hz, {total_s} s", flush=True)

    sig_fixed = np.full(total_s, np.nan)
    sig_fsp = np.full(total_s, np.nan)
    swa_1 = np.full(total_s, np.nan)
    ep_dr = np.full(n_ep, np.nan); ep_swa = np.full(n_ep, np.nan); ep_clean = np.zeros(n_ep)
    ep_emg = np.full(n_ep, np.nan); ep_eog = np.full(n_ep, np.nan)
    beats = []
    so_t = {c: [] for c in ctx}

    sos_fixed = band_sos(SIGMA_FIXED, sf)
    sos_fsp = band_sos((fsp - 1, fsp + 1), sf)
    sos_swa = band_sos(SWA_BAND_L, sf, 3)
    sos_so = signal.butter(3, list(SO_BAND_NAJI), btype="band", fs=sf, output="sos")
    sos_emg = signal.butter(3, [10.0, min(100.0, sf / 2 - 10)], btype="band", fs=sf, output="sos")
    sos_eog = signal.butter(3, [0.3, 6.0], btype="band", fs=sf, output="sos")

    # 50 Hz line notch (+ harmonics) as a precomputed IIR cascade -- reused per channel/chunk, far
    # cheaper than re-designing an FIR (mne.filter.notch_filter) on every call.
    notch_sos = np.vstack([signal.tf2sos(*signal.iirnotch(f0, 30.0, sf))
                           for f0 in np.arange(50.0, sf / 2 - 1, 50.0)])

    def notch50(x):
        return signal.sosfiltfilt(notch_sos, x)

    picks_all = list(ie) + [ecg_i] + ([emg_i] if emg_i is not None else []) + \
                ([eog_i] if eog_i is not None else [])
    t = 0.0
    while t < total_s:
        dur = min(CHUNK_S, total_s - t)
        a = int(t * sf); b = min(int((t + dur) * sf), n_samp)
        d = raw.get_data(picks=picks_all, start=a, stop=b) * 1e6      # (n_pick, n) in uV
        off = int(t)
        x_ie = d[:len(ie)]
        x_ecg = d[len(ie)]
        col = len(ie) + 1
        x_emg = d[col] if emg_i is not None else None
        x_eog = d[col + (1 if emg_i is not None else 0)] if eog_i is not None else None

        xa = notch50(np.nan_to_num(x_ie).mean(axis=0))

        # ECG -> beats
        try:
            cl = nk.ecg_clean(np.nan_to_num(x_ecg), sampling_rate=int(sf), method="neurokit")
            _, info = nk.ecg_peaks(cl, sampling_rate=int(sf), method="neurokit", correct_artifacts=True)
            beats.extend((np.asarray(info["ECG_R_Peaks"], int) / sf + off).tolist())
        except Exception:
            pass

        # per-channel SO troughs
        for ci, c in enumerate(ctx):
            xc = notch50(signal.detrend(np.nan_to_num(x_ie[ci])))
            if np.std(xc) < 1e-9:
                continue
            tr = detect_so_halfwaves(signal.sosfiltfilt(sos_so, xc), sf)
            if len(tr):
                so_t[c].extend((tr / sf + off).tolist())

        clean = ied_clean_mask(xa, sf)
        n_sec = int(len(xa) // int(sf))
        k = int(sf)
        for sos_b, dest in ((sos_fixed, sig_fixed), (sos_fsp, sig_fsp), (sos_swa, swa_1)):
            env2 = np.abs(signal.hilbert(signal.sosfiltfilt(sos_b, xa))) ** 2
            e2 = env2[:n_sec * k].reshape(n_sec, k)
            cm = clean[:n_sec * k].reshape(n_sec, k).astype(float)
            num, den = (e2 * cm).sum(1), cm.sum(1)
            vals = np.where(den >= 0.5 * k, num / np.maximum(den, 1e-12), np.nan)
            sl = slice(off, min(off + n_sec, total_s))
            dest[sl] = vals[:sl.stop - sl.start]

        # EMG / EOG epoch features
        emg_f = np.abs(signal.sosfiltfilt(sos_emg, notch50(np.nan_to_num(x_emg)))) if x_emg is not None else None
        eog_f = signal.sosfiltfilt(sos_eog, np.nan_to_num(x_eog)) if x_eog is not None else None

        ke = int(EPOCH * sf)
        for e in range(int(len(xa) // ke)):
            gi = int((off + e * EPOCH) // EPOCH)
            if gi >= n_ep:
                break
            seg = xa[e * ke:(e + 1) * ke]
            ep_dr[gi] = delta_ratio(seg, sf)
            fq, pp = signal.welch(seg, sf, nperseg=int(min(4 * sf, len(seg))))
            m = (fq >= SWA_BAND[0]) & (fq < SWA_BAND[1])
            ep_swa[gi] = float(np.trapezoid(pp[m], fq[m]))
            ep_clean[gi] = float(clean[e * ke:(e + 1) * ke].mean())
            if emg_f is not None:
                ep_emg[gi] = float(np.sqrt(np.mean(emg_f[e * ke:(e + 1) * ke] ** 2)))
            if eog_f is not None:
                ep_eog[gi] = float(np.var(eog_f[e * ke:(e + 1) * ke]))
        t += dur

    # heart rate grids
    beats = np.array(sorted(beats))
    hr_1 = np.full(total_s, np.nan)
    hr_4 = np.full(int(total_s * FS_RR), np.nan)
    if len(beats) > 20:
        rr = np.diff(beats); good = (rr >= 0.33) & (rr <= 1.5)
        if good.sum() > 20:
            tb, hrv = beats[1:][good], 60.0 / rr[good]
            g1 = np.arange(total_s)
            hr_1 = np.where((g1 >= tb[0]) & (g1 <= tb[-1]), np.interp(g1, tb, hrv), np.nan)
            g4 = np.arange(0, total_s, 1.0 / FS_RR)
            hr_4 = interpolate.CubicSpline(tb, hrv, extrapolate=False)(g4)

    stage_lab, stage_counts = score_stages(ep_swa, ep_emg, ep_eog, ep_dr, ep_clean)

    os.makedirs(OUT, exist_ok=True)
    payload = dict(status="ok", subject=subject, sf=sf, night_s=0.0, hours=total_s / 3600.0,
                   cortical_chans=np.array(ctx), ekg=names[ecg_i], fsp=fsp, fsp_is_real_peak=fsp_real,
                   sigma_fixed=sig_fixed, sigma_fsp=sig_fsp, swa=swa_1,
                   hr_1=hr_1, hr_4=hr_4, fs_rr=FS_RR,
                   ep_dr=ep_dr, ep_swa=ep_swa, ep_clean=ep_clean, ep_emg=ep_emg, ep_eog=ep_eog,
                   epoch_s=EPOCH, beats=beats, stage_lab=np.array(stage_lab),
                   has_eog=bool(eog_i is not None), has_emg=bool(emg_i is not None))
    for c in ctx:
        payload[f"so_t_{c}"] = np.array(sorted(so_t[c]))
    np.savez_compressed(fp, **payload)
    frac = float(np.isfinite(sig_fsp).mean())
    print(f"[{subject}] cached {time.time()-t0:.0f}s | sigma_cov {frac:.0%} | {len(beats)} beats | "
          f"stages W/R/N2/N3 = {stage_counts['n_wake']}/{stage_counts['n_rem']}/"
          f"{int((stage_lab=='N2').sum())}/{int((stage_lab=='N3').sum())} of {n_ep} ep", flush=True)
    if not keep_raw:
        try:
            os.remove(paths["_ieeg.eeg"])
            print(f"[{subject}] removed raw .eeg", flush=True)
        except OSError:
            pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default=",".join(SUBJECTS))
    ap.add_argument("--keep-raw", action="store_true")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    for s in a.subjects.split(","):
        s = s.strip()
        if s not in SUBJECTS:
            print(f"[{s}] unknown subject", flush=True); continue
        try:
            run(s, keep_raw=a.keep_raw, force=a.force)
        except Exception as e:
            import traceback
            print(f"[{s}] ERROR {type(e).__name__}: {e}", flush=True); traceback.print_exc()


if __name__ == "__main__":
    main()
