"""Validation for the 3B surrogate-null fix.

Ground truth by construction: heart rate contains NO slow-oscillation locking whatsoever. The only
structure is that the two stages have slightly different mean heart rates (N2 92 bpm, N3 94 bpm) --
exactly the situation in HUP165, where N2 sat below and N3 above the whole-night mean of ~93.2 bpm.

A correct test must return z ~ 0 for BOTH stages. The old whole-night surrogate pool instead returns
a large positive z for one stage and a large negative z for the other.

    .venv/bin/python analysis/test_3B_null.py
"""
import numpy as np

from event_3B_mednick import so_triggered, FS_RR, HALF_WIN

rng = np.random.RandomState(0)
EPOCH = 30.0
TOTAL_S = 3600
hw = int(HALF_WIN * FS_RR)
n = int(TOTAL_S * FS_RR)


def old_so_triggered(hr, trough_times_s, stage_mean_hr, n_sur=200, rng=None):
    """The ORIGINAL: surrogate triggers drawn from the WHOLE NIGHT."""
    idx = np.round(np.asarray(trough_times_s) * FS_RR).astype(int)
    idx = idx[(idx >= hw) & (idx < len(hr) - hw)]
    if len(idx) < 30:
        return None
    seg = np.stack([hr[i - hw:i + hw] for i in idx])
    curve = seg.mean(axis=0)
    lag = (np.arange(len(curve)) - hw) / FS_RR
    post = lag >= 0
    pct = 100.0 * (float(curve[post].max()) - stage_mean_hr) / stage_mean_hr
    null = []
    for _ in range(n_sur):
        r = rng.randint(hw, len(hr) - hw, size=len(idx))
        c2 = np.stack([hr[i - hw:i + hw] for i in r]).mean(axis=0)
        null.append(100.0 * (c2[post].max() - stage_mean_hr) / stage_mean_hr)
    null = np.array(null)
    return dict(pct=pct, z=float((pct - null.mean()) / (null.std() + 1e-12)),
                null_mean_pct=float(null.mean()))


# --- build a night with two stages and NO SO->HR coupling -------------------------------------
n_ep = int(TOTAL_S // EPOCH)
stage = np.where(rng.rand(n_ep) < 0.5, "N2", "N3")
hr = np.zeros(n)
for e in range(n_ep):
    s0, s1 = int(e * EPOCH * FS_RR), int((e + 1) * EPOCH * FS_RR)
    hr[s0:s1] = (92.0 if stage[e] == "N2" else 94.0)
hr += np.convolve(rng.randn(n), np.ones(40) / 40, mode="same") * 6.0   # smooth physiological noise

pools, means, troughs = {}, {}, {}
for st in ("N2", "N3"):
    m = np.zeros(n, bool)
    eps = np.where(stage == st)[0]
    for e in eps:
        m[int(e * EPOCH * FS_RR):int((e + 1) * EPOCH * FS_RR)] = True
    pools[st] = np.where(m)[0]
    means[st] = float(hr[m].mean())
    # SO troughs placed at RANDOM times inside the stage -> zero locking by construction
    cand = pools[st][(pools[st] >= hw) & (pools[st] < n - hw)]
    troughs[st] = np.sort(rng.choice(cand, size=900, replace=False)) / FS_RR

whole_night_mean = float(hr[hw:n - hw].mean())
print(f"\nwhole-night mean HR = {whole_night_mean:.3f} bpm")
print(f"stage mean HR: N2 = {means['N2']:.3f}, N3 = {means['N3']:.3f} bpm")
print("SO troughs are placed at RANDOM times within each stage -> TRUE effect is exactly zero.\n")

print(f"{'stage':6s} {'':>4s} {'observed %':>11s} {'null mean %':>12s} {'z':>8s}")
res = {}
for st in ("N2", "N3"):
    o = old_so_triggered(hr, troughs[st], means[st], rng=np.random.RandomState(1))
    nw = so_triggered(hr, troughs[st], means[st], pools[st], rng=np.random.RandomState(1))
    res[st] = (o, nw)
    print(f"{st:6s} {'OLD':>4s} {o['pct']:11.3f} {o['null_mean_pct']:12.3f} {o['z']:8.2f}")
    print(f"{st:6s} {'NEW':>4s} {nw['pct_above_stage_mean']:11.3f} {nw['null_mean_pct']:12.3f} {nw['z']:8.2f}")

fails = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not cond:
        fails.append(name)


print("\nchecks:")
check("OLD produces a spurious significant result despite zero true effect",
      max(abs(res['N2'][0]['z']), abs(res['N3'][0]['z'])) > 3,
      f"|z| up to {max(abs(res['N2'][0]['z']), abs(res['N3'][0]['z'])):.1f}")
check("OLD flips sign between stages (the HUP165 signature)",
      res['N2'][0]['z'] * res['N3'][0]['z'] < 0,
      f"N2 z={res['N2'][0]['z']:.2f}, N3 z={res['N3'][0]['z']:.2f}")
check("NEW returns z ~ 0 for N2", abs(res['N2'][1]['z']) < 3, f"z = {res['N2'][1]['z']:.2f}")
check("NEW returns z ~ 0 for N3", abs(res['N3'][1]['z']) < 3, f"z = {res['N3'][1]['z']:.2f}")
check("NEW null baseline is ~0 in both stages (no stage/night offset leaking in)",
      max(abs(res['N2'][1]['null_mean_pct']), abs(res['N3'][1]['null_mean_pct'])) < 0.30,
      f"N2 {res['N2'][1]['null_mean_pct']:+.3f}%, N3 {res['N3'][1]['null_mean_pct']:+.3f}%")
check("OLD null baselines invert to ONE common absolute HR (proves it is stage-blind)",
      abs(means['N2'] * (1 + res['N2'][0]['null_mean_pct'] / 100)
          - means['N3'] * (1 + res['N3'][0]['null_mean_pct'] / 100)) < 0.5,
      f"{means['N2']*(1+res['N2'][0]['null_mean_pct']/100):.3f} vs "
      f"{means['N3']*(1+res['N3'][0]['null_mean_pct']/100):.3f} bpm")

print("\n" + ("ALL CHECKS PASSED" if not fails else f"{len(fails)} FAILURE(S): {fails}"))
raise SystemExit(1 if fails else 0)
