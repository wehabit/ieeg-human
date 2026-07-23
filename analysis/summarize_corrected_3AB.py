"""Cohort summary of the CORRECTED 3A (Lecci-faithful) and 3B (stage-matched null) results.

Reads the 3A and 3B per-subject JSON dirs (default HUP; pass --a-dir/--b-dir for ds003848).

    .venv/bin/python analysis/summarize_corrected_3AB.py
    .venv/bin/python analysis/summarize_corrected_3AB.py \
        --a-dir outputs/ds003848_3A --b-dir outputs/ds003848_3B --label "ds003848 replication"
"""
import argparse, glob, json, os
import numpy as np
from scipy import stats

rng_np = np.random.RandomState(0)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_ap = argparse.ArgumentParser()
_ap.add_argument("--a-dir", default=os.path.join(ROOT, "outputs", "lecci_faithful_3A"))
_ap.add_argument("--b-dir", default=os.path.join(ROOT, "outputs", "event_3B_cached"))
_ap.add_argument("--label", default="HUP phaseII")
_args, _ = _ap.parse_known_args()
A_DIR = _args.a_dir if os.path.isabs(_args.a_dir) else os.path.join(ROOT, _args.a_dir)
B_DIR = _args.b_dir if os.path.isabs(_args.b_dir) else os.path.join(ROOT, _args.b_dir)
print(f"COHORT: {_args.label}   (3A: {A_DIR}, 3B: {B_DIR})")

LECCI_PEAK, LECCI_SD = 0.019, 0.0052        # human, 0.019 +/- 0.001 SEM over n=27
NAJI = {"N2": 12.09, "N3": 3.35}


def load_dir(d):
    return [json.load(open(f)) for f in sorted(glob.glob(os.path.join(d, "*.json")))]


def hr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


# ------------------------------------------------------------------ 3A
A = [r for r in load_dir(A_DIR) if r.get("status") == "ok"]
hr(f"3A -- LECCI-FAITHFUL, n = {len(A)} subjects")

if A:
    print("\nSTEP 1 (Lecci Fig 1G): does sigma power oscillate infraslow, and where is each")
    print("subject's own peak? Duration-weighted Morlet spectrum over all NREM bouts >= 120 s.\n")
    print(f"{'subject':18s} {'bouts':>5s} {'sec':>6s} {'peak Hz':>8s} {'prom':>6s} "
          f"{'SWA peak':>9s} {'SWA prom':>8s}")
    pk, prom, pk_swa, prom_swa = [], [], [], []
    for r in A:
        p, s = r["peak"], r["peak_swa"]
        if p.get("peak_hz"):
            pk.append(p["peak_hz"]); prom.append(p.get("prominence_over_background", np.nan))
        if s.get("peak_hz"):
            pk_swa.append(s["peak_hz"]); prom_swa.append(s.get("prominence_over_background", np.nan))
        print(f"{r['subject']:18s} {r['n_bouts']:5d} {int(r['bout_seconds']):6d} "
              f"{p.get('peak_hz') or float('nan'):8.4f} {p.get('prominence_over_background', np.nan):6.2f} "
              f"{s.get('peak_hz') or float('nan'):9.4f} {s.get('prominence_over_background', np.nan):8.2f}")
    pk, prom = np.array(pk, float), np.array(prom, float)
    prom_swa = np.array(prom_swa, float)
    print(f"\n  sigma peak: mean {np.nanmean(pk):.4f} Hz, SD {np.nanstd(pk, ddof=1):.4f}, "
          f"median {np.nanmedian(pk):.4f}  (Lecci human: {LECCI_PEAK:.3f} +/- {LECCI_SD:.4f})")
    within = np.abs(pk - LECCI_PEAK) <= LECCI_SD
    print(f"  subjects whose own peak is within 1 SD of Lecci's 0.019 Hz: {within.sum()}/{len(pk)}")
    t = stats.ttest_1samp(pk, LECCI_PEAK)
    print(f"  one-sample t vs 0.019 Hz: t = {t.statistic:+.2f}, p = {t.pvalue:.3f}")

    # THE DECISIVE COHORT TEST. The per-subject prominence control is weak (on synthetic data a
    # planted rhythm scored p~0.21, never reaching 0.05 -- see the module docstring). What DOES
    # discriminate is peak CLUSTERING: with a real shared rhythm the fitted peaks bunch (simulated
    # SD ~0.0008 Hz, 100% in 0.015-0.025 Hz); pure 1/f^2 noise scatters (SD ~0.013 Hz, ~48% in band
    # by chance). This is also Lecci's own form of evidence (0.019 +/- 0.001 across n=27).
    frac_band = np.mean((pk >= 0.015) & (pk <= 0.025))
    print(f"\n  CLUSTERING TEST (the discriminating one -- see docstring):")
    print(f"    observed peak SD = {np.nanstd(pk, ddof=1):.4f} Hz  "
          f"(planted-rhythm sim ~0.0008 Hz; 1/f noise sim ~0.013 Hz)")
    print(f"    fraction in 0.015-0.025 Hz = {frac_band:.2f}  (1/f noise ~0.48 by chance)")
    # surrogate scatter aggregated across subjects, if the per-subject control was run
    all_surr = [p for r in A if r.get("peak_null") for p in r["peak_null"].get("surrogate_peaks", [])]
    if all_surr:
        all_surr = np.array(all_surr)
        surr_sd = np.std(all_surr, ddof=1)
        surr_band = np.mean((all_surr >= 0.015) & (all_surr <= 0.025))
        print(f"    THIS COHORT's own scale-free surrogates: peak SD {surr_sd:.4f} Hz, "
              f"{surr_band:.2f} in band ({len(all_surr)} surrogate peaks)")
        # bootstrap p: is the real peak SD tighter than surrogate draws of the same size?
        bs = np.array([np.std(rng_np.choice(all_surr, len(pk), replace=False), ddof=1)
                       for _ in range(2000)]) if len(all_surr) >= len(pk) else None
        if bs is not None:
            p_tight = float((1 + np.sum(bs <= np.nanstd(pk, ddof=1))) / (1 + len(bs)))
            print(f"    bootstrap: real peaks are tighter than surrogate draws with p = {p_tight:.3f}")
    print(f"\n  LECCI'S NEGATIVE CONTROL -- sigma should be MORE prominent than SWA:")
    print(f"    sigma prominence median {np.nanmedian(prom):.2f} | SWA {np.nanmedian(prom_swa):.2f}")
    m = np.isfinite(prom) & np.isfinite(prom_swa)
    if m.sum() >= 5:
        w = stats.wilcoxon(prom[m], prom_swa[m], alternative="greater")
        print(f"    Wilcoxon sigma > SWA: p = {w.pvalue:.4f}  (n = {int(m.sum())}, "
              f"sigma higher in {int((prom[m] > prom_swa[m]).sum())})")

    print("\nSTEP 2 (Lecci Fig 6): does heart rate track it?")
    co = [r["coherence"] for r in A if r.get("coherence")]
    if co:
        K = np.array([c["K"] for c in co])
        e_l = np.array([c["sig_at_lecci"] for c in co])
        e_o = np.array([bool(c.get("sig_at_own_peak")) for c in co])
        print(f"\n  gap-aware coherence over ALL NREM: K median {int(np.median(K))} "
              f"(was 7-23 on one 55-min block), median crit {np.median([c['crit'] for c in co]):.4f}")
        print(f"  significant at 0.02 Hz            : {e_l.sum()}/{len(co)}")
        print(f"  significant at own fitted peak    : {e_o.sum()}/{len(co)}")
        f = np.array(co[0]["f"])
        C = np.array([c["cxy"] for c in co])
        crit = np.array([c["crit"] for c in co])
        exc = C > crit[:, None]
        band = (f >= 0.005) & (f <= 0.25)
        cnt = exc[:, band].sum(0)
        i02 = int(np.argmin(np.abs(f - 0.02)))
        c02 = int(exc[:, i02].sum())
        bg = cnt.mean() / len(co)
        print(f"\n  FREQUENCY-SPECIFICITY CONTROL (the test that retracted the original 3A):")
        print(f"    exceedance at 0.02 Hz      : {c02}/{len(co)} = {c02/len(co):.3f}")
        print(f"    mean across {band.sum()} band bins : {bg:.3f}")
        print(f"    rank of 0.02 Hz            : {int((cnt > c02).sum())+1} of {band.sum()}")
        if bg > 0:
            p = stats.binomtest(c02, len(co), min(max(bg, 1e-6), 0.999), alternative="greater").pvalue
            print(f"    binomial vs empirical background: p = {p:.4f}")
    xc = [r["xcorr"] for r in A if r.get("xcorr")]
    if xc:
        lag = np.array([x["lag_s"] for x in xc])[0]
        M = np.array([x["xcorr"] for x in xc])
        g = M.mean(0)
        i = int(np.argmax(np.abs(g)))
        r_pk = np.array([x["peak_r"] for x in xc])
        l_pk = np.array([x["peak_lag_s"] for x in xc])
        print(f"\n  CROSS-CORRELATION (Lecci's actual coupling statistic, never previously run):")
        print(f"    group mean correlogram peak |r| = {abs(g[i]):.4f} at lag {lag[i]:+.0f} s")
        print(f"    per-subject |peak r|: median {np.median(np.abs(r_pk)):.4f}, "
              f"max {np.abs(r_pk).max():.4f}")
        print(f"    per-subject peak lag: median {np.median(l_pk):+.1f} s, "
              f"IQR [{np.percentile(l_pk,25):+.0f}, {np.percentile(l_pk,75):+.0f}]")
        t = stats.ttest_1samp(r_pk, 0.0)
        print(f"    one-sample t on peak r vs 0: t = {t.statistic:+.2f}, p = {t.pvalue:.3f}")
        print(f"    (Lecci report a clear correlogram peak with a consistent lag; a null result is")
        print(f"     |r| ~ 0 with lags scattered uniformly across subjects)")

# ------------------------------------------------------------------ 3B
B = [r for r in load_dir(B_DIR) if r.get("status") == "ok"]
hr(f"3B -- NAJI 2019 WITH A STAGE-MATCHED NULL, n = {len(B)} subjects (was n = 1)")
if B:
    print(f"\n{'subject':18s} {'N2 %':>7s} {'N2 z':>7s} {'N3 %':>7s} {'N3 z':>7s} "
          f"{'N2 lag':>7s} {'N3 lag':>7s}")
    for r in B:
        v2, v3 = r.get("N2"), r.get("N3")
        print(f"{r['subject']:18s} "
              f"{(v2 or {}).get('pct_above_stage_mean', float('nan')):7.3f} "
              f"{(v2 or {}).get('z', float('nan')):7.2f} "
              f"{(v3 or {}).get('pct_above_stage_mean', float('nan')):7.3f} "
              f"{(v3 or {}).get('z', float('nan')):7.2f} "
              f"{(v2 or {}).get('peak_lag_s', float('nan')):7.2f} "
              f"{(v3 or {}).get('peak_lag_s', float('nan')):7.2f}")
    for st in ("N2", "N3"):
        v = np.array([r[st]["pct_above_stage_mean"] for r in B if r.get(st)], float)
        z = np.array([r[st]["z"] for r in B if r.get(st)], float)
        nb = np.array([r[st]["null_mean_pct"] for r in B if r.get(st)], float)
        if not len(v):
            continue
        t = stats.ttest_1samp(z, 0.0)
        print(f"\n  {st}: HR peak {np.mean(v):+.3f}% (SD {np.std(v, ddof=1):.3f}), "
              f"Naji {NAJI[st]:+.2f}%  -> {NAJI[st]/max(np.mean(v),1e-9):.0f}x weaker")
        print(f"      z: mean {np.mean(z):+.2f}, median {np.median(z):+.2f}, "
              f"significant (z>1.96) in {(z > 1.96).sum()}/{len(z)}")
        print(f"      one-sample t on z vs 0: t = {t.statistic:+.2f}, p = {t.pvalue:.4f}")
        print(f"      null baseline now {np.mean(nb):+.3f}% (was contaminated by the "
              f"stage-vs-night HR offset)")
    p2 = [(r["N2"]["pct_above_stage_mean"], r["N3"]["pct_above_stage_mean"])
          for r in B if r.get("N2") and r.get("N3")]
    if len(p2) >= 5:
        a2, a3 = np.array([x[0] for x in p2]), np.array([x[1] for x in p2])
        w = stats.wilcoxon(a2, a3)
        print(f"\n  Naji predicts N2 >> SWS (12.09% vs 3.35%, ~3.6x). Here: N2 {np.mean(a2):+.3f}% "
              f"vs N3 {np.mean(a3):+.3f}%, paired p = {w.pvalue:.3f} (n = {len(p2)})")

print()
