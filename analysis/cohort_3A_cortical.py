"""
Test 3A across the HUP phase-II cohort, on LATERAL NEOCORTICAL contacts.

For every subject with depth electrodes + EKG, this automates what was done by hand for HUP165:
  1. pick cortical channels  -- highest-numbered contact on each shaft (SEEG: contact 1 is deepest
     / mesial, the highest contact is the most lateral, i.e. neocortex)
  2. find a night         -- probe delta-ratio across the record, take the highest-delta stretch
  3. find a clean window  -- scan segments for NREM (delta-ratio) AND usable EKG (40-130 bpm);
                             keep the longest contiguous clean+NREM run
  4. run 3A               -- epoch-binned sigma power (FSP+/-1 Hz, IED-masked) x instantaneous HR,
                             magnitude-squared coherence, analytic significance threshold
  5. SWA control          -- same pipeline on 0.75-4 Hz, which should stay null (Lecci Fig 1)

Writes one JSON per subject as it goes (safe to interrupt / resume) plus a cohort CSV.

    .venv/bin/python analysis/cohort_3A_cortical.py [--subjects 116,130,...] [--win-min 120]
"""
import argparse, json, os, re, time, traceback
import numpy as np
from scipy import signal

from infraslow_rr_sigma_coherence import sess, get, pull_continuous, notch, ROOT
from results_3A_tutorial_style import (band_power_series, fast_spindle_peak, hr_grid, hp, msc,
                                       FS_P, INFRA, SWA, F_TARGET, EDGE_TRIM_S)

COHORT = [116, 130, 133, 138, 139, 141, 143, 150, 151, 157, 160, 165, 171, 172, 173,
          177, 178, 182, 185, 187, 191, 199, 205, 211, 212]
OUT = os.path.join(ROOT, "outputs", "cohort_3A_cortical")
EPOCH = 30.0
NREM_DR = 0.90


def delta_ratio(x, sf):
    m = np.isfinite(x)
    if m.mean() < 0.7 or np.std(x[m]) < 1e-9:
        return np.nan
    if not m.all():
        idx = np.arange(len(x)); x = x.copy(); x[~m] = np.interp(idx[~m], idx[m], x[m])
    f, p = signal.welch(x, sf, nperseg=int(min(4 * sf, len(x))))
    num = np.trapezoid(p[(f >= 0.5) & (f < 4)], f[(f >= 0.5) & (f < 4)])
    den = np.trapezoid(p[(f >= 0.5) & (f < 25)], f[(f >= 0.5) & (f < 25)])
    return float(num / den) if den > 0 else np.nan


def cortical_channels(labels, n_want=6):
    """Highest contact on each shaft = most lateral = neocortex (SEEG geometry)."""
    sh = {}
    for l in labels:
        m = re.match(r"^([A-Z]{1,3})(\d+)$", l)
        if m and not l.upper().startswith(("EKG", "ECG")):
            sh.setdefault(m.group(1), []).append(int(m.group(2)))
    # Take the highest contact on each shaft. Shaft lengths differ across subjects (HUP165 has
    # 12-contact shafts, most others 8), so an absolute contact-number cutoff would silently drop
    # whole subjects -- require only that the shaft is long enough to span from mesial to lateral.
    picks = [(pre, max(nums)) for pre, nums in sh.items() if max(nums) >= 6]
    picks.sort(key=lambda t: -t[1])               # longest shafts first = furthest lateral reach
    # Contact numbering is not always contiguous, so "prefix + max_contact" can name a channel that
    # does not exist -- validate against the real label list rather than trusting the construction.
    have = set(labels)
    return [f"{p}{n}" for p, n in picks if f"{p}{n}" in have][:n_want]


def find_night(ds, lab_idx, ch, sf, total_h, scan_h=30.0, step_min=30.0):
    step = step_min * 60
    n = int(min(scan_h, total_h) * 3600 / step)
    tr = []
    for k in range(n):
        t = k * step
        try:
            tr.append((t, delta_ratio(get(ds, [lab_idx[ch]], t, 6.0)[:, 0], sf)))
        except Exception:
            tr.append((t, np.nan))
    tr = [(t, d) for t, d in tr if np.isfinite(d)]
    if len(tr) < 6:
        return None
    win = max(1, int(3 * 3600 / step))            # 3 h window
    best_t, best_m = tr[0][0], -1
    for i in range(max(1, len(tr) - win)):
        m = np.mean([d for _, d in tr[i:i + win]])
        if m > best_m:
            best_m, best_t = m, tr[i][0]
    return best_t


def swa_power(x, sf):
    """Absolute slow-wave (0.75-4 Hz) power -- the axis that separates deep from light NREM."""
    f, p = signal.welch(np.nan_to_num(x), sf, nperseg=int(min(4 * sf, len(x))))
    m = (f >= 0.75) & (f < 4.0)
    return float(np.trapezoid(p[m], f[m]))


def scan_segments(ds, lab_idx, ch, ekg, sf, night_t, span_h=5.0, seg=300.0):
    """Per-segment NREM / EKG-usable flags plus slow-wave power, across the night."""
    import neurokit2 as nk
    out = []
    for i in range(int(span_h * 3600 / seg)):
        t = night_t + i * seg
        rec = dict(t=t, ok=False, swa=np.nan)
        try:
            d = get(ds, sorted([lab_idx[ekg], lab_idx[ch]]), t, seg)
            ei = 0 if lab_idx[ekg] < lab_idx[ch] else 1   # get_data returns ascending index order
            e, s_ = d[:, ei], d[:, 1 - ei]
            dr = delta_ratio(s_, sf)
            try:
                cl = nk.ecg_clean(np.nan_to_num(e), sampling_rate=int(sf), method="neurokit")
                _, info = nk.ecg_peaks(cl, sampling_rate=int(sf), method="neurokit", correct_artifacts=True)
                bpm = len(info["ECG_R_Peaks"]) / (seg / 60)
            except Exception:
                bpm = 0
            rec["swa"] = swa_power(s_, sf)
            rec["ok"] = bool(np.isfinite(dr) and dr >= NREM_DR and 40 <= bpm <= 130)
        except Exception:
            pass
        out.append(rec)
    return out


def longest_run(flags, seg):
    best_len = best_start = cur = cur_start = 0
    for i, g in enumerate(flags):
        if g:
            if cur == 0:
                cur_start = i
            cur += 1
            if cur > best_len:
                best_len, best_start = cur, cur_start
        else:
            cur = 0
    return best_start, best_len * seg / 60.0


def pick_stage_windows(segs, seg=300.0, need_min=45.0):
    """Split clean NREM into DEEP (high slow-wave activity, N3-like) and LIGHT (low SWA, N2-like)
    by a within-subject median split, then take the longest contiguous run of each.

    NOTE: this is NOT scored AASM staging -- no EOG/EMG exists in these recordings. It is a
    slow-wave-activity split, which is the physiological axis Lecci's S2>SWS claim rides on.
    Human SWA declines across the night, so the two runs typically fall in early vs late cycles."""
    ok = np.array([s["ok"] for s in segs])
    swa = np.array([s["swa"] for s in segs], float)
    if ok.sum() < 2 * need_min * 60 / seg:
        return None
    # Smooth SWA to the sleep-cycle timescale (~25 min) first. Raw segment-to-segment SWA is
    # noisy and oscillates across the median, which fragments the split into many short runs
    # instead of the long deep/light blocks that sleep cycles actually produce.
    w = 5
    sm = np.full(len(swa), np.nan)
    for i in range(len(swa)):
        s = swa[max(0, i - w // 2):i + w // 2 + 1]
        s = s[np.isfinite(s)]
        if len(s):
            sm[i] = s.mean()
    thr = np.nanmedian(sm[ok])
    deep = ok & (sm >= thr)
    light = ok & (sm < thr)
    out = {}
    for name, flags in (("deep", deep), ("light", light)):
        i0, mins = longest_run(flags, seg)
        if mins < need_min:
            return None
        sel = swa[i0:i0 + int(mins * 60 / seg)]
        out[name] = dict(start=segs[i0]["t"], minutes=mins, mean_swa=float(np.nanmean(sel)))
    return out


SO_BAND = (0.5, 1.25)


def detect_so(x, sf):
    """SO troughs: band-pass 0.5-1.25 Hz, take the largest-amplitude negative peaks."""
    sos = signal.butter(3, list(SO_BAND), btype="band", fs=sf, output="sos")
    f = signal.sosfiltfilt(sos, np.nan_to_num(x))
    neg, _ = signal.find_peaks(-f, distance=int(0.5 * sf))
    if len(neg) < 20:
        return np.array([], int), f
    thr = np.percentile(-f[neg], 75)              # top quartile = slow oscillations
    return neg[-f[neg] >= thr], f


def tort_mi(phase, amp, nbins=18):
    bins = np.linspace(-np.pi, np.pi, nbins + 1)
    m = np.array([amp[(phase >= bins[i]) & (phase < bins[i + 1])].mean() if
                  np.any((phase >= bins[i]) & (phase < bins[i + 1])) else np.nan
                  for i in range(nbins)])
    if np.isnan(m).any() or m.sum() <= 0:
        return np.nan
    p = m / m.sum()
    return float((np.log(nbins) + np.sum(p * np.log(p + 1e-12))) / np.log(nbins))


def test_3D_so_spindle(x_avg, sf, sigma_band, n_surr=100):
    """SO-phase -> sigma-amplitude coupling (Tort MI, surrogate-corrected z).
    Surrogates ARE valid here: shifting breaks the phase-amplitude relation for a
    broadband/event measure (unlike narrowband coherence, where a shift is invisible)."""
    x = np.nan_to_num(x_avg)
    sos_so = signal.butter(3, list(SO_BAND), btype="band", fs=sf, output="sos")
    ph = np.angle(signal.hilbert(signal.sosfiltfilt(sos_so, x)))
    sos_sp = signal.butter(4, list(sigma_band), btype="band", fs=sf, output="sos")
    am = np.abs(signal.hilbert(signal.sosfiltfilt(sos_sp, x)))
    mi = tort_mi(ph, am)
    if not np.isfinite(mi):
        return dict(mi=None, mi_z=None)
    rng = np.random.RandomState(0)
    null = [tort_mi(ph, np.roll(am, rng.randint(int(sf), len(am) - int(sf)))) for _ in range(n_surr)]
    null = np.array([v for v in null if np.isfinite(v)])
    z = float((mi - null.mean()) / (null.std() + 1e-12)) if len(null) > 10 else None
    return dict(mi=float(mi), mi_z=z, preferred_phase=float(np.angle(np.mean(np.exp(1j * ph) * am))))


def test_3B_so_heartbeat(x_avg, sf, hr_1hz, n_surr=200, half=10):
    """SO-trough-triggered heart rate: does HR systematically dip/rise around the SO trough?
    Null = random trigger times (destroys the SO-locking, keeps the HR autocorrelation)."""
    troughs, _ = detect_so(x_avg, sf)
    if len(troughs) < 30:
        return dict(n_so=int(len(troughs)), modulation=None, z=None)
    t_s = (troughs / sf).astype(int)
    t_s = t_s[(t_s >= half) & (t_s < len(hr_1hz) - half)]
    if len(t_s) < 30:
        return dict(n_so=int(len(t_s)), modulation=None, z=None)
    seg = np.stack([hr_1hz[t - half:t + half] for t in t_s])
    curve = seg.mean(axis=0) - seg.mean()
    mod = float(curve.max() - curve.min())
    rng = np.random.RandomState(0)
    null = []
    for _ in range(n_surr):
        r = rng.randint(half, len(hr_1hz) - half, size=len(t_s))
        s2 = np.stack([hr_1hz[t - half:t + half] for t in r])
        c2 = s2.mean(axis=0) - s2.mean()
        null.append(c2.max() - c2.min())
    null = np.array(null)
    return dict(n_so=int(len(t_s)), modulation=mod, z=float((mod - null.mean()) / (null.std() + 1e-12)),
                curve=[float(v) for v in curve])


def run_subject(n, win_min):
    name = f"HUP{n}_phaseII"
    fp = os.path.join(OUT, f"{name}.json")
    if os.path.exists(fp):
        print(f"[{name}] cached, skip", flush=True); return json.load(open(fp))
    rec = dict(subject=name)
    s = sess(); ds = s.open_dataset(name)
    labels = ds.get_channel_labels()
    lab_idx = {l: i for i, l in enumerate(labels)}
    d0 = ds.get_time_series_details(labels[0])
    sf = d0.sample_rate; total_h = (getattr(d0, "duration", 0) or 0) / 3.6e9
    ekg = next((l for l in labels if l.upper().startswith(("EKG", "ECG"))), None)
    ctx = cortical_channels(labels)
    if ekg is None or len(ctx) < 3:
        rec.update(status="skip", reason=f"ekg={ekg} n_cortical={len(ctx)}")
        json.dump(rec, open(fp, "w"), indent=2); print(f"[{name}] SKIP {rec['reason']}", flush=True); return rec
    print(f"[{name}] {sf:.0f} Hz, {total_h:.0f} h, ekg={ekg}, cortical={ctx}", flush=True)

    night = find_night(ds, lab_idx, ctx[0], sf, total_h)
    if night is None:
        rec.update(status="skip", reason="no night found")
        json.dump(rec, open(fp, "w"), indent=2); print(f"[{name}] SKIP no night", flush=True); return rec
    segs = scan_segments(ds, lab_idx, ctx[0], ekg, sf, night)
    wins = pick_stage_windows(segs)
    if wins is None:
        rec.update(status="skip", reason="no deep+light NREM runs of >=45 min each")
        json.dump(rec, open(fp, "w"), indent=2); print(f"[{name}] SKIP {rec['reason']}", flush=True); return rec
    print(f"[{name}] night {night/3600:.1f} h; "
          f"deep {wins['deep']['minutes']:.0f} min @ {wins['deep']['start']/3600:.2f} h; "
          f"light {wins['light']['minutes']:.0f} min @ {wins['light']['start']/3600:.2f} h", flush=True)

    idx = [lab_idx[c] for c in ctx] + [lab_idx[ekg]]
    # EQUAL window length for deep and light: magnitude-squared coherence is biased by the number
    # of Welch segments (its noise floor is ~1/K), so unequal windows would make the stage contrast
    # partly a window-length artifact rather than a stage effect.
    use_min = float(min(win_min, wins["deep"]["minutes"], wins["light"]["minutes"]))
    rec.update(status="ok", sf=sf, cortical_chans=ctx, ekg=ekg, night_h=night / 3600,
               matched_window_min=use_min)
    for stage in ("deep", "light"):
        w = wins[stage]
        data = pull_continuous(ds, idx, w["start"], use_min * 60)
        x_ctx, x_ekg = data[:, :len(ctx)], data[:, len(ctx)]
        dur = len(x_ekg) / sf
        fsp = fast_spindle_peak(x_ctx, sf)
        sig_p, ied = band_power_series(x_ctx, sf, (fsp - 1, fsp + 1), dur)
        swa_p, _ = band_power_series(x_ctx, sf, SWA, dur)
        hr, med_hr, nbeats, frac_bad = hr_grid(x_ekg, sf, len(sig_p))
        k = int(EDGE_TRIM_S * FS_P)
        C = msc(hp(hr)[k:-k], hp(sig_p)[k:-k])
        Cs = msc(hp(hr)[k:-k], hp(swa_p)[k:-k])
        # 3B and 3D on the SAME pulled data (the pull is the expensive part, these are ~free)
        x_avg = np.nan_to_num(x_ctx.astype(float)).mean(axis=1)
        d3 = test_3D_so_spindle(x_avg, sf, (fsp - 1, fsp + 1))
        b3 = test_3B_so_heartbeat(x_avg, sf, hr)

        rec[stage] = dict(start_h=w["start"] / 3600, dur_min=dur / 60, run_min=w["minutes"],
                          mean_swa=w["mean_swa"], fast_spindle_peak_hz=fsp, ied_masked_frac=ied,
                          nbeats=nbeats, median_hr_bpm=med_hr, rr_dropped_frac=frac_bad,
                          coh_at_0p02=C["at"], coh_band_max=C["bmax"], band_peak_hz=C["fpk"],
                          crit_bin=C["crit"], K=C["K"],
                          significant_at_0p02=bool(C["at"] > C["crit"]),
                          swa_coh_at_0p02=Cs["at"],
                          so_spindle_mi=d3["mi"], so_spindle_mi_z=d3["mi_z"],
                          so_hr_n=b3["n_so"], so_hr_modulation=b3["modulation"], so_hr_z=b3["z"])
        print(f"[{name}]   {stage:5s} {dur/60:5.0f} min  3A coh@0.02={C['at']:.3f} "
              f"(crit {C['crit']:.3f}, K={C['K']})  swa-ctrl={Cs['at']:.3f} | "
              f"3D MI_z={d3['mi_z'] if d3['mi_z'] is None else round(d3['mi_z'],1)} | "
              f"3B z={b3['z'] if b3['z'] is None else round(b3['z'],1)} (n_SO={b3['n_so']})", flush=True)
    rec["light_minus_deep"] = rec["light"]["coh_at_0p02"] - rec["deep"]["coh_at_0p02"]
    json.dump(rec, open(fp, "w"), indent=2)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default=",".join(map(str, COHORT)))
    ap.add_argument("--win-min", type=float, default=120.0)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    subs = [int(x) for x in a.subjects.split(",") if x.strip()]
    rows = []
    for n in subs:
        try:
            rows.append(run_subject(n, a.win_min))
        except Exception as e:
            print(f"[HUP{n}] ERROR {type(e).__name__}: {e}", flush=True)
            traceback.print_exc()
            rows.append(dict(subject=f"HUP{n}_phaseII", status="error", reason=str(e)[:200]))
    ok = [r for r in rows if r.get("status") == "ok"]
    if ok:
        import csv
        flat = []
        for r in ok:
            row = dict(subject=r["subject"], sf=r["sf"], night_h=r.get("night_h"),
                       cortical_chans="|".join(r["cortical_chans"]),
                       light_minus_deep=r["light_minus_deep"])
            for st in ("deep", "light"):
                for k, v in r[st].items():
                    row[f"{st}_{k}"] = v
            flat.append(row)
        keys = sorted({k for r in flat for k in r})
        with open(os.path.join(OUT, "cohort_3A_cortical.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(flat)

        d = np.array([r["deep"]["coh_at_0p02"] for r in ok])
        l = np.array([r["light"]["coh_at_0p02"] for r in ok])
        sw = np.array([r["light"]["swa_coh_at_0p02"] for r in ok])
        print(f"\n=== COHORT n={len(ok)} · coherence @0.02 Hz (cortical) ===")
        print(f"  LIGHT NREM (N2-like) : median {np.median(l):.3f}  mean {l.mean():.3f}  range {l.min():.3f}-{l.max():.3f}")
        print(f"  DEEP  NREM (N3-like) : median {np.median(d):.3f}  mean {d.mean():.3f}  range {d.min():.3f}-{d.max():.3f}")
        print(f"  SWA control (light)  : median {np.median(sw):.3f}")
        # Lecci's prediction: the 0.02 Hz rhythm is stronger in lighter NREM (S2) than in SWS
        try:
            from scipy import stats
            w_st, w_p = stats.wilcoxon(l, d)
            t_st, t_p = stats.ttest_rel(l, d)
            print(f"  paired light>deep    : Wilcoxon p={w_p:.4f}  |  paired t p={t_p:.4f}  "
                  f"({int((l > d).sum())}/{len(ok)} subjects light>deep)")
        except Exception as e:
            print(f"  paired test failed: {e}")
        # chance is per-window; report how many individual windows cleared their own threshold
        print(f"  windows over own crit: light {sum(r['light']['significant_at_0p02'] for r in ok)}/{len(ok)}, "
              f"deep {sum(r['deep']['significant_at_0p02'] for r in ok)}/{len(ok)}")
    print(f"\nper-subject JSON -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
