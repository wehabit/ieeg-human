"""
3D done properly: EVENT-BASED SO->spindle coupling, PER CHANNEL, split by N2-like / N3-like.

Replaces the continuous Tort-MI version in cohort_stages_3ABD.py, which did not follow the source
papers and very likely manufactured its own null:

  Helfrich 2018: "we detected SO (0.16-1.25 Hz) and sleep spindle (12-16 Hz) EVENTS ... phase during
                  the PEAK of the detected sleep spindle events ... Rayleigh z"
  Staresina 2015: "in an Event-locked analysis ..." (V-test on preferred phase)

The old version binned EVERY timepoint by SO phase into a modulation index. Most of a night is
neither spindle nor slow oscillation, so that averages the events together with hours of nothing and
dilutes the estimate toward zero regardless of the true coupling. It also averaged 6 channels into
one signal first, which partially cancels the slow oscillation before anything is measured.

This implementation follows master/staresina_style.py (which reproduces the published method):
  per channel -> detect SO troughs and spindle peaks by amplitude threshold -> take the SO phase at
  each spindle peak -> Rayleigh test on that phase distribution.

Reported per stage: number of channels with significant phase clustering, resultant vector length R,
and the preferred phase. Staging is reused from cohort_stages_3ABD (GMM on log slow-wave power).

    .venv/bin/python analysis/event_3D_by_stage.py [--subjects 165,...] [--hours 7]
"""
import argparse, json, os, traceback
import numpy as np
from scipy import signal

from infraslow_rr_sigma_coherence import sess, pull_continuous, notch, ROOT
from cohort_3A_cortical import COHORT, cortical_channels, delta_ratio, find_night
from cohort_stages_3ABD import (fsp_from, band_sos, stage_epochs, EPOCH, CHUNK_S, SWA_BAND, SO_BAND)

OUT = os.path.join(ROOT, "outputs", "event_3D_by_stage")
SP_THR = 1.5        # spindle-peak threshold, z of the spindle envelope (same as staresina_style)
SP_DIST = 0.3       # s, min spacing between spindle peaks
SO_DIST = 0.8       # s, min spacing between SO troughs
MIN_EVENTS = 20     # per channel per stage, for a Rayleigh test


def rayleigh(ph):
    """Resultant vector length R, and p for non-uniformity (as in staresina_style.py)."""
    n = len(ph)
    if n < MIN_EVENTS:
        return np.nan, np.nan, np.nan
    C, S = np.mean(np.cos(ph)), np.mean(np.sin(ph))
    R = float(np.hypot(C, S)); z = n * R ** 2
    p = float(np.exp(-z) * (1 + (2 * z - z ** 2) / (4 * n)))
    return R, p, float(np.arctan2(S, C))


def run(n, hours):
    name = f"HUP{n}_phaseII"
    fp = os.path.join(OUT, f"{name}.json")
    if os.path.exists(fp):
        print(f"[{name}] cached", flush=True); return
    s = sess(); ds = s.open_dataset(name)
    labels = ds.get_channel_labels(); lab_idx = {l: i for i, l in enumerate(labels)}
    d0 = ds.get_time_series_details(labels[0]); sf = d0.sample_rate
    total_h = (getattr(d0, "duration", 0) or 0) / 3.6e9
    ekg = next((l for l in labels if l.upper().startswith(("EKG", "ECG"))), None)
    ctx = cortical_channels(labels)
    if ekg is None or len(ctx) < 3:
        json.dump(dict(subject=name, status="skip"), open(fp, "w")); print(f"[{name}] SKIP", flush=True); return
    night = find_night(ds, lab_idx, ctx[0], sf, total_h)
    if night is None:
        json.dump(dict(subject=name, status="skip"), open(fp, "w")); print(f"[{name}] SKIP no night", flush=True); return

    idx = [lab_idx[c] for c in ctx] + [lab_idx[ekg]]
    fsp, real = fsp_from(pull_continuous(ds, idx, night, 300.0)[:, :len(ctx)], sf)
    sos_so = band_sos(SO_BAND, sf, 3)
    sos_sp = band_sos((fsp - 1, fsp + 1), sf)
    print(f"[{name}] {sf:.0f} Hz | {len(ctx)} cortical ch | FSP {fsp:.2f} Hz | streaming {hours} h", flush=True)

    total_s = int(hours * 3600); n_ep = int(total_s // EPOCH)
    ep = dict(dr=np.full(n_ep, np.nan), swa=np.full(n_ep, np.nan), clean=np.ones(n_ep))
    # per channel: list of (epoch_index, SO phase at spindle peak)
    phases = {c: {"ep": [], "ph": []} for c in ctx}
    n_so = {c: 0 for c in ctx}

    t = 0.0
    while t < total_s:
        dur = min(CHUNK_S, total_s - t)
        try:
            d = pull_continuous(ds, idx, night + t, dur)
        except Exception:
            t += dur; continue
        off = int(t)
        for ci, c in enumerate(ctx):
            x = notch(signal.detrend(np.nan_to_num(d[:, ci].astype(float))), sf)
            if np.std(x) < 1e-9:
                continue
            so = signal.sosfiltfilt(sos_so, x)
            so_ph = np.angle(signal.hilbert(so))
            sp = signal.sosfiltfilt(sos_sp, x)
            env = np.abs(signal.hilbert(sp))
            spz = (env - env.mean()) / (env.std() + 1e-12)
            # DISCRETE EVENTS -- this is the step the Tort-MI version omitted
            so_tr, _ = signal.find_peaks(-so, height=so.std(), distance=int(SO_DIST * sf))
            sp_pk, _ = signal.find_peaks(spz, height=SP_THR, distance=int(SP_DIST * sf))
            n_so[c] += int(len(so_tr))
            if len(sp_pk):
                phases[c]["ph"].extend(so_ph[sp_pk].tolist())
                phases[c]["ep"].extend(((off + sp_pk / sf) // EPOCH).astype(int).tolist())
        # staging features from the channel mean (staging only, not coupling)
        xa = notch(np.nan_to_num(d[:, :len(ctx)].astype(float)).mean(axis=1), sf)
        ke = int(EPOCH * sf)
        for e in range(int(len(xa) // ke)):
            gi = int((off + e * EPOCH) // EPOCH)
            if gi >= n_ep:
                break
            seg = xa[e * ke:(e + 1) * ke]
            ep["dr"][gi] = delta_ratio(seg, sf)
            fq, pp = signal.welch(seg, sf, nperseg=int(min(4 * sf, len(seg))))
            m = (fq >= SWA_BAND[0]) & (fq < SWA_BAND[1])
            ep["swa"][gi] = float(np.trapezoid(pp[m], fq[m]))
        t += dur

    lab, nrem, sep = stage_epochs(ep)
    rec = dict(subject=name, status="ok", sf=sf, cortical_chans=ctx, fsp=fsp, fsp_real=real,
               n_nrem=int(nrem.sum()), n_N2=int((lab == "N2").sum()), n_N3=int((lab == "N3").sum()),
               gmm_separation=sep, n_so_total=int(sum(n_so.values())))
    for stage in ("N2", "N3", "ALL"):
        keep = set(np.where(lab == stage)[0].tolist()) if stage != "ALL" else set(np.where(nrem)[0].tolist())
        per_ch, n_ev = [], 0
        for c in ctx:
            e = np.array(phases[c]["ep"], int); p = np.array(phases[c]["ph"], float)
            if len(e) == 0:
                continue
            sel = np.array([x in keep for x in e], bool)
            R, pv, mu = rayleigh(p[sel])
            n_ev += int(sel.sum())
            if np.isfinite(R):
                per_ch.append(dict(ch=c, n=int(sel.sum()), R=R, p=pv, preferred_phase_deg=np.degrees(mu)))
        sig = [d for d in per_ch if d["p"] < 0.05]
        rec[stage] = dict(n_spindle_events=n_ev, n_channels_tested=len(per_ch),
                          n_channels_significant=len(sig),
                          median_R=float(np.median([d["R"] for d in per_ch])) if per_ch else None,
                          mean_preferred_phase_deg=(float(np.degrees(np.angle(np.mean(
                              [np.exp(1j * np.radians(d["preferred_phase_deg"])) for d in sig]))))
                              if sig else None),
                          per_channel=per_ch)
        print(f"[{name}]   {stage:3s}: {len(sig)}/{len(per_ch)} channels significant | "
              f"{n_ev} spindle events | median R={rec[stage]['median_R']}", flush=True)
    os.makedirs(OUT, exist_ok=True)
    json.dump(rec, open(fp, "w"), indent=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default=",".join(map(str, COHORT)))
    ap.add_argument("--hours", type=float, default=7.0)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    for n in [int(x) for x in a.subjects.split(",") if x.strip()]:
        try:
            run(n, a.hours)
        except Exception as e:
            print(f"[HUP{n}] ERROR {type(e).__name__}: {e}", flush=True); traceback.print_exc()
    print(f"\nJSON -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
