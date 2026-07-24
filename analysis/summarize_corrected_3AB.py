"""Cohort summary of corrected 3A (Lecci-aligned approximation) and 3B results.

Reads the 3A and 3B per-subject JSON dirs (default HUP; pass --a-dir/--b-dir for ds003848).

    .venv/bin/python analysis/summarize_corrected_3AB.py
    .venv/bin/python analysis/summarize_corrected_3AB.py \
        --a-dir outputs/ds003848_3A --b-dir outputs/ds003848_3B --label "ds003848 replication"
"""
import argparse, glob, json, os
import numpy as np
from scipy import stats
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    cache_code_sha256,
    file_sha256,
    runtime_versions,
    source_tree_sha256,
)
from lecci_faithful_3A import (
    cache_lineage_entry,
    LECCI_XCORR_LAG_WINDOW_S,
    peak_location_null_is_adequate,
    verify_cache_lineage,
)

rng_np = np.random.RandomState(0)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_ap = argparse.ArgumentParser()
_ap.add_argument("--a-dir", default=os.path.join(ROOT, "outputs", "lecci_faithful_3A"))
_ap.add_argument("--b-dir", default=os.path.join(ROOT, "outputs", "event_3B_cached"))
_ap.add_argument("--label", default="HUP phaseII")
_ap.add_argument(
    "--expected-subjects",
    help="comma-separated exact subject IDs; inferred for the standard HUP/RESPect directories")
_ap.add_argument(
    "--min-completed", type=int, default=5,
    help="minimum estimable participants required for a cohort summary (default: 5)")
_args = _ap.parse_args()
A_DIR = _args.a_dir if os.path.isabs(_args.a_dir) else os.path.join(ROOT, _args.a_dir)
B_DIR = _args.b_dir if os.path.isabs(_args.b_dir) else os.path.join(ROOT, _args.b_dir)
HUP_IDS = {
    f"HUP{n}_phaseII" for n in
    (116, 130, 133, 138, 139, 141, 143, 150, 151, 157, 160, 165, 171, 172, 173,
     177, 178, 182, 185, 187, 191, 199, 205, 211, 212)
}
RESPECT_IDS = {
    "sub-RESP0521", "sub-RESP0699", "sub-RESP0724",
    "sub-RESP0749", "sub-RESP0779", "sub-RESP0800"
}
if _args.expected_subjects:
    EXPECTED = {value.strip() for value in _args.expected_subjects.split(",") if value.strip()}
elif "ds003848" in os.path.basename(A_DIR).lower():
    EXPECTED = RESPECT_IDS
elif os.path.basename(A_DIR) == "lecci_faithful_3A":
    EXPECTED = HUP_IDS
else:
    raise SystemExit("--expected-subjects is required for nonstandard output directories")
print(f"COHORT: {_args.label}   (3A: {A_DIR}, 3B: {B_DIR})")

# Distinct quantities in Lecci: main-text values are mean +/- SEM, whereas Fig. 1G's 0.008 Hz is
# the within-spectrum Gaussian width. It is not the across-participant SD of fitted locations.
LECCI_PEAK, LECCI_PEAK_SEM, LECCI_N = 0.019, 0.001, 27
LECCI_SPECTRAL_SD = 0.008
NAJI = {"N2": 12.09, "N3": 3.35}


def load_dir(d):
    rows = []
    for f in sorted(glob.glob(os.path.join(d, "*.json"))):
        if os.path.basename(f) == "RUN_MANIFEST.json":
            continue
        value = json.load(open(f))
        if value.get("subject"):
            value["_source_file"] = f
            rows.append(value)
    return rows


def validate_run(directory, records, label):
    """Reject partial/mixed production runs instead of globbing whatever happens to exist."""
    current = [r for r in records if r.get("analysis_version") == ANALYSIS_VERSION]
    subject_ids = [r.get("subject") for r in current]
    duplicates = sorted({
        subject for subject in subject_ids if subject_ids.count(subject) > 1
    })
    if duplicates:
        raise SystemExit(f"{label}: duplicate subject result files: {duplicates}")
    bad_filenames = [
        os.path.basename(r.get("_source_file", ""))
        for r in current
        if os.path.basename(r.get("_source_file", "")) != f"{r.get('subject')}.json"
    ]
    if bad_filenames:
        raise SystemExit(
            f"{label}: result filename must exactly match <subject>.json: {bad_filenames}")
    manifest_path = os.path.join(directory, "RUN_MANIFEST.json")
    if not os.path.exists(manifest_path):
        raise SystemExit(f"{label}: RUN_MANIFEST.json is missing; no complete production run")
    manifest = json.load(open(manifest_path))
    if manifest.get("run_state") != "complete" or not manifest.get("run_id"):
        raise SystemExit(f"{label}: run manifest is interrupted/incomplete")
    if manifest.get("analysis_version") != ANALYSIS_VERSION:
        raise SystemExit(f"{label}: manifest analysis version mismatch")
    if manifest.get("cache_schema_version") != CACHE_SCHEMA_VERSION:
        raise SystemExit(f"{label}: manifest cache schema mismatch")
    if manifest.get("runtime_versions") != runtime_versions():
        raise SystemExit(f"{label}: runtime/dependency versions differ from the completed run")
    expected_pipelines = {
        "3A": {"lecci_faithful_3A", "ds003848_lecci_3A"},
        "3B": {"event_3B_cached", "ds003848_event_3B"},
    }[label]
    pipeline = manifest.get("pipeline")
    if pipeline not in expected_pipelines:
        raise SystemExit(
            f"{label}: manifest pipeline {pipeline!r} is not a production {label} pipeline")
    config = manifest.get("config")
    if not isinstance(config, dict):
        raise SystemExit(f"{label}: manifest config is missing or malformed")
    common_required = dict(
        analysis_version=ANALYSIS_VERSION,
        cache_schema_version=CACHE_SCHEMA_VERSION,
        cache_code_sha256=cache_code_sha256(ROOT),
    )
    endpoint_required = (
        dict(band="fixed", smooth_4s=True, n_sur=200)
        if label == "3A"
        else dict(
            tachogram_domain="rr",
            n_surrogates=1000,
            null_method="shared circular shift in eligible stage-time",
            contact_qc=(
                "cache stable >=80%-coverage plus optional Destrieux frontal ROI intersection"),
            minimum_event_channels=2,
            stable_stage_minimum_s=180,
        )
    )
    wrong_config = {
        key: (config.get(key), expected)
        for key, expected in {**common_required, **endpoint_required}.items()
        if config.get(key) != expected
    }
    if wrong_config:
        raise SystemExit(f"{label}: non-production manifest configuration: {wrong_config}")
    if manifest.get("failed"):
        raise SystemExit(f"{label}: run manifest contains failed subjects")
    current_tree = source_tree_sha256(ROOT)
    if manifest.get("source_tree_sha256") != current_tree:
        raise SystemExit(
            f"{label}: executable source tree differs from the completed run; rerun analysis")
    requested_values = manifest.get("requested", [])
    requested = set(requested_values)
    if len(requested_values) != len(requested):
        raise SystemExit(f"{label}: manifest requested list contains duplicate subjects")
    if requested != EXPECTED:
        raise SystemExit(
            f"{label}: requested subjects differ from the prespecified expected set "
            f"(missing={sorted(EXPECTED-requested)}, unexpected={sorted(requested-EXPECTED)})")
    completed_values = manifest.get("completed", [])
    completed = set(completed_values)
    if len(completed_values) != len(completed):
        raise SystemExit(f"{label}: manifest completed list contains duplicate subjects")
    skipped_entries = manifest.get("skipped", [])
    skipped_values = [
        value.get("subject") if isinstance(value, dict) else value
        for value in skipped_entries
    ]
    skipped = set(skipped_values)
    if any(not isinstance(value, dict) or not value.get("subject") or not value.get("reason")
           for value in skipped_entries):
        raise SystemExit(f"{label}: every skipped subject needs an explicit reason")
    if len(skipped_values) != len(skipped):
        raise SystemExit(f"{label}: manifest skipped list contains duplicate subjects")
    if completed | skipped != requested or completed & skipped:
        raise SystemExit(f"{label}: requested subjects are not exactly completed or explicitly skipped")
    observed = {r["subject"] for r in current if r.get("status") == "ok"}
    if observed != completed:
        raise SystemExit(
            f"{label}: output subjects do not equal manifest completed set "
            f"(missing={sorted(completed-observed)}, extra={sorted(observed-completed)})")
    file_subjects = {r["subject"] for r in current}
    if file_subjects != completed | skipped:
        raise SystemExit(
            f"{label}: version-current result/partial files do not cover the manifest "
            f"(missing={sorted((completed | skipped)-file_subjects)}, "
            f"extra={sorted(file_subjects-(completed | skipped))})")
    actual_result_hashes = {
        r["subject"]: file_sha256(r["_source_file"]) for r in current
    }
    if manifest.get("result_files_sha256") != actual_result_hashes:
        raise SystemExit(
            f"{label}: result file hashes do not match the completed run manifest")
    schemas = {r.get("cache_schema_version") for r in current}
    if schemas != {CACHE_SCHEMA_VERSION}:
        raise SystemExit(f"{label}: mixed/missing cache schemas: {schemas}")
    if {r.get("run_id") for r in current} != {manifest["run_id"]}:
        raise SystemExit(f"{label}: result files are not linked to this completed run")
    if {r.get("source_tree_sha256") for r in current} != {current_tree}:
        raise SystemExit(f"{label}: result source digests do not match the completed run")
    record_cache_inputs = {
        r["subject"]: cache_lineage_entry(r) for r in current
    }
    if manifest.get("config", {}).get("cache_inputs") != record_cache_inputs:
        raise SystemExit(
            f"{label}: analysis manifest cache-input lineage does not match subject records")
    if pipeline in {"lecci_faithful_3A", "event_3B_cached"}:
        bad_hours = [
            r["subject"] for r in current
            if r.get("cache_hours") is None
            or not np.isclose(float(r["cache_hours"]), 7.0)
        ]
        if bad_hours:
            raise SystemExit(
                f"{label}: HUP production cache duration must be 7 h for every subject: "
                f"{bad_hours}")
    configured_records = [
        record for record in current if record.get("status") in {"ok", "partial"}
    ]
    if label == "3A":
        bad_record_config = [
            record["subject"] for record in configured_records
            if (
                record.get("band") != "fixed"
                or record.get("smooth_4s") is not True
                or record.get("n_surrogates_requested") != 200
            )
        ]
    else:
        bad_record_config = [
            record["subject"] for record in configured_records
            if (
                record.get("tachogram_domain") != "rr"
                or record.get("n_surrogates") != 1000
                or record.get("null_method")
                != "shared circular shift in eligible stage-time"
                or record.get("minimum_event_channels") != 2
                or record.get("stable_stage_minimum_s") != 180
            )
        ]
        for record in configured_records:
            for stage in ("N2", "N3"):
                endpoint = record.get(stage)
                if endpoint is None:
                    continue
                identifiers = endpoint.get("event_contact_ids")
                if (
                    endpoint.get("tachogram_domain") != "rr"
                    or endpoint.get("n_surrogates") != 1000
                    or endpoint.get("n_channels", 0) < 2
                    or not isinstance(identifiers, list)
                    or len(identifiers) != endpoint.get("n_channels")
                ):
                    bad_record_config.append(f"{record['subject']}:{stage}")
    if bad_record_config:
        raise SystemExit(
            f"{label}: subject record does not agree with the production configuration: "
            f"{sorted(set(bad_record_config))}")
    for record in current:
        try:
            verify_cache_lineage(record)
        except Exception as exc:
            raise SystemExit(f"{label}: cache lineage validation failed: {exc}") from exc
    return current


def hr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


# ------------------------------------------------------------------ 3A
all_A = load_dir(A_DIR)
current_A = validate_run(A_DIR, all_A, "3A")
A_spectrum = [
    r for r in current_A
    if (r.get("endpoint_availability") or {}).get("spectrum")
]
A_xcorr = [
    r for r in current_A
    if (r.get("endpoint_availability") or {}).get("cross_correlation")
]
A_coherence_fixed = [
    r for r in current_A
    if (r.get("endpoint_availability") or {}).get("coherence_fixed_0p02")
]
A_coherence_own = [
    r for r in current_A
    if (r.get("endpoint_availability") or {}).get("coherence_own_peak")
]
A_has_cohort_endpoint = max(len(A_spectrum), len(A_xcorr)) >= _args.min_completed
if not A_has_cohort_endpoint:
    reasons = [
        f"{r['subject']}: {r.get('reason', 'unspecified')}"
        for r in current_A if r.get("status") != "ok"
    ]
    print(
        f"3A: spectrum n={len(A_spectrum)}, cross-correlation n={len(A_xcorr)}; at least "
        f"{_args.min_completed} are required for one primary endpoint. "
        f"3A cohort result unavailable. Excluded/partial: {reasons}")
hr("3A -- LECCI-ALIGNED APPROXIMATION")
partial_A = [r for r in current_A if r.get("status") != "ok"]
print(f"\nEndpoint denominators: spectrum {len(A_spectrum)}/{len(EXPECTED)}; "
      f"cross-correlation {len(A_xcorr)}/{len(EXPECTED)}; "
      f"fixed-0.02-Hz coherence {len(A_coherence_fixed)}/{len(EXPECTED)}; "
      f"own-peak coherence {len(A_coherence_own)}/{len(EXPECTED)}; "
      f"{len(partial_A)} subject records explicitly skipped/partial.")
for r in partial_A:
    print(f"  {r['subject']}: {r.get('reason', 'unspecified')}")

if A_spectrum:
    print("\nSTEP 1 (Lecci Fig 1G): does sigma power oscillate infraslow, and where is each")
    print("subject's own peak? Duration-weighted Morlet spectrum over all NREM bouts >= 120 s.\n")
    print(f"{'subject':18s} {'bouts':>5s} {'sec':>6s} {'peak Hz':>8s} {'prom':>6s} "
          f"{'SWA peak':>9s} {'SWA prom':>8s}")
    pk, prom, paired_peak_values = [], [], []
    for r in A_spectrum:
        p, s = r["peak"], r["peak_swa"]
        if p.get("peak_hz"):
            pk.append(p["peak_hz"]); prom.append(p.get("prominence_over_background", np.nan))
        control = r.get("negative_control") or {}
        if (np.isfinite(control.get("sigma_window_mean", np.nan))
                and np.isfinite(control.get("swa_same_window_mean", np.nan))):
            paired_peak_values.append(
                (r["subject"], control["sigma_window_mean"],
                 control["swa_same_window_mean"]))
        print(f"{r['subject']:18s} {r['n_bouts']:5d} {int(r['bout_seconds']):6d} "
              f"{p.get('peak_hz') or float('nan'):8.4f} {p.get('prominence_over_background', np.nan):6.2f} "
              f"{s.get('peak_hz') or float('nan'):9.4f} {s.get('prominence_over_background', np.nan):8.2f}")
    pk, prom = np.array(pk, float), np.array(prom, float)
    if len(pk):
        print(f"\n  accepted sigma peaks: {len(pk)}/{len(A_spectrum)}")
        print(f"  peak locations: mean {np.nanmean(pk):.4f} Hz, "
              f"SD {np.nanstd(pk, ddof=1):.4f}, median {np.nanmedian(pk):.4f}")
        print(f"  Lecci reports group mean {LECCI_PEAK:.3f} +/- {LECCI_PEAK_SEM:.3f} Hz "
              f"(SEM, n={LECCI_N}); Fig. 1G's {LECCI_SPECTRAL_SD:.3f} Hz is a "
              "within-spectrum Gaussian width, not between-person SD.")
        if len(pk) >= 2:
            sem = float(np.std(pk, ddof=1) / np.sqrt(len(pk)))
            ci = stats.t.interval(0.95, len(pk) - 1, loc=float(np.mean(pk)), scale=sem)
            print(f"  accepted-peak subset mean difference from Lecci point estimate: "
                  f"{np.mean(pk) - LECCI_PEAK:+.4f} Hz; conditional 95% CI "
                  f"[{ci[0]:.4f}, {ci[1]:.4f}]")
            print("  no one-sample t test: Lecci's 0.019 Hz is an estimated mean (SEM 0.001), "
                  "not a fixed population constant, the subset was selected for accepted peaks, "
                  "and the cohorts/methods differ.")

    accepted_records = [r for r in A_spectrum if r["peak"].get("peak_hz")]
    surrogate_by_subject = []
    if len(accepted_records) >= _args.min_completed:
        for r in accepted_records:
            peak_null = r.get("peak_null") or {}
            values = np.asarray(peak_null.get("surrogate_peaks", []), float)
            if peak_location_null_is_adequate(peak_null):
                surrogate_by_subject.append((r["subject"], values))
        print("\n  PEAK-TIGHTNESS-ONLY CONTROL AGAINST SUBJECT-MATCHED SCALE-FREE SURROGATES:")
        print(f"    observed accepted-peak SD = {np.std(pk, ddof=1):.4f} Hz")
        if len(surrogate_by_subject) == len(accepted_records):
            null_sd = np.array([
                np.std([values[rng_np.randint(len(values))]
                        for _, values in surrogate_by_subject], ddof=1)
                for _ in range(10000)
            ])
            p_tight = float(
                (1 + np.sum(null_sd <= np.std(pk, ddof=1))) / (1 + len(null_sd)))
            print(f"    matched Monte-Carlo p for tighter real clustering = {p_tight:.4f}")
            print("    This tests clustering anywhere in the accepted search band; it cannot show")
            print("    compatibility with Lecci's 0.019-Hz location. No equivalence margin was")
            print("    prespecified, so this p value is not evidence of Lecci-location replication.")
        else:
            print("    unavailable: at least one accepted subject has fewer than 20 accepted "
                  "surrogate peak locations or <5% acceptance")
    elif accepted_records:
        print("\n  PEAK-LOCATION CLUSTERING unavailable: "
              f"{len(accepted_records)} accepted peaks; at least "
              f"{_args.min_completed} are required.")

    print("\n  LECCI NEGATIVE CONTROL -- both bands in the same sigma-defined peak window:")
    if paired_peak_values:
        sigma_values = np.asarray([value[1] for value in paired_peak_values], float)
        swa_values = np.asarray([value[2] for value in paired_peak_values], float)
        print(f"    sigma median {np.median(sigma_values):.3f} | "
              f"SWA median {np.median(swa_values):.3f} | n={len(sigma_values)}")
        print("    descriptive only: the window was selected for a high/accepted sigma peak, so an")
        print("    ordinary paired sigma>SWA p value would be selection-biased. Valid inference")
        print("    must repeat sigma-peak selection inside a joint null or use a prespecified window.")
    else:
        print("    unavailable: no participant has an accepted sigma peak and both spectra")

print("\nSTEP 2 (Lecci Fig 6): does heart rate track it?")
co = [r["coherence"] for r in A_coherence_fixed]
if co:
    K = np.array([c["K"] for c in co])
    e_l = np.array([c["sig_at_lecci"] for c in co])
    print(f"\n  gap-aware coherence over ALL NREM (descriptive; n={len(co)}): "
          f"K median {int(np.median(K))}, "
          f"median analytic crit {np.median([c['crit'] for c in co]):.4f}")
    print(f"  significant at 0.02 Hz            : {e_l.sum()}/{len(co)}")
    own = [r["coherence"] for r in A_coherence_own]
    if own:
        e_o = np.array([bool(c.get("sig_at_own_peak")) for c in own])
        print(f"  significant at own fitted peak    : {e_o.sum()}/{len(own)}")
    else:
        print("  significant at own fitted peak    : unavailable (0 accepted own peaks)")
    f = np.array(co[0]["f"])
    C = np.asarray([c["cxy"] for c in co], float)
    crit = np.array([c["crit"] for c in co])
    finite = np.isfinite(C)
    exc = finite & (C > crit[:, None])
    support = finite.sum(0)
    band = (f >= 0.005) & (f <= 0.25)
    minimum_bin_support = max(_args.min_completed, int(np.ceil(0.80 * len(co))))
    supported_band = band & (support >= minimum_bin_support)
    proportion = np.divide(
        exc.sum(0), support, out=np.full(len(f), np.nan), where=support > 0)
    i02 = int(np.argmin(np.abs(f - 0.02)))
    c02 = int(exc[:, i02].sum())
    c02_support = int(support[i02])
    bg = float(np.mean(proportion[supported_band])) if supported_band.any() else np.nan
    c02_proportion = c02 / c02_support
    rank = (
        int(np.sum(proportion[supported_band] > c02_proportion)) + 1
        if supported_band.any() else None)
    print(f"\n  FREQUENCY-SPECIFICITY CONTROL (the test that retracted the original 3A):")
    print(f"    exceedance at 0.02 Hz      : {c02}/{c02_support} = {c02_proportion:.3f}")
    print(f"    mean across {supported_band.sum()} supported band bins : {bg:.3f}")
    print(f"    rank of 0.02 Hz            : "
          f"{rank if rank is not None else 'n/a'} of {supported_band.sum()}")
    print(f"    per-bin support gate       : >={minimum_bin_support}/{len(co)} participants")
    print("    descriptive only: neighbouring frequency bins are correlated and cannot serve")
    print("    as independent binomial trials or as a pre-specified inferential null.")
else:
    print("\n  gap-aware coherence unavailable for every participant.")

xc = [r["xcorr"] for r in A_xcorr]
if xc:
    lag_arrays = [np.asarray(x["lag_s"], float) for x in xc]
    if any(not np.array_equal(values, lag_arrays[0]) for values in lag_arrays[1:]):
        raise SystemExit("3A cross-correlation lag grids differ between participants")
    lag = lag_arrays[0]
    M = np.array([x["xcorr"] for x in xc])
    g = M.mean(0)
    i = int(np.argmax(np.abs(g)))
    follows = (
        (lag >= LECCI_XCORR_LAG_WINDOW_S[0])
        & (lag <= LECCI_XCORR_LAG_WINDOW_S[1]))
    i_directional = np.where(follows)[0][int(np.argmax(g[follows]))]
    i_opposite = np.where(follows)[0][int(np.argmin(g[follows]))]
    r_pk = np.array([x["peak_r"] for x in xc])
    l_pk = np.array([x["peak_lag_s"] for x in xc])
    print(f"\n  CROSS-CORRELATION (Lecci's actual coupling statistic; n={len(xc)}):")
    print(f"    signed max-|r| sensitivity peak = {g[i]:+.4f} at lag {lag[i]:+.0f} s")
    print(f"    Lecci-direction/timing-window peak (positive r, sigma follows HR; "
          f"{LECCI_XCORR_LAG_WINDOW_S[0]:.0f} to "
          f"{LECCI_XCORR_LAG_WINDOW_S[1]:.0f} s) = "
          f"{g[i_directional]:+.4f} at lag {lag[i_directional]:+.0f} s")
    print(f"    opposite-direction extremum at nonnegative lag = "
          f"{g[i_opposite]:+.4f} at lag {lag[i_opposite]:+.0f} s")
    print(f"    per-subject signed max-|r|: median {np.median(r_pk):+.4f}; "
          f"median |r| {np.median(np.abs(r_pk)):.4f}")
    print(f"    per-subject peak lag: median {np.median(l_pk):+.1f} s, "
          f"IQR [{np.percentile(l_pk,25):+.0f}, {np.percentile(l_pk,75):+.0f}]")
    if len(xc) >= _args.min_completed:
        # A t-test on each subject's independently selected lag does not test Lecci's common group
        # correlogram peak. A sign-flip max-statistic asks whether one group-level lag survives
        # selection over all lags.
        obs_directional = float(np.max(g[follows]))
        obs_omnibus = float(np.max(np.abs(g)))
        if len(M) <= 15:
            assignments = np.arange(1 << len(M), dtype=np.uint64)
            bit_positions = np.arange(len(M), dtype=np.uint64)
            sign_matrix = (
                2.0 * ((assignments[:, None] >> bit_positions) & 1).astype(float) - 1.0)
            null_curves = sign_matrix @ M / len(M)
            null_directional = np.max(null_curves[:, follows], axis=1)
            null_omnibus = np.max(np.abs(null_curves), axis=1)
            p_directional = float(np.mean(null_directional >= obs_directional))
            p_omnibus = float(np.mean(null_omnibus >= obs_omnibus))
            inference_method = f"exact enumeration of {len(null_directional)} sign assignments"
        else:
            null_directional = np.empty(20000)
            null_omnibus = np.empty(20000)
            for j in range(len(null_directional)):
                signs = rng_np.choice((-1.0, 1.0), size=len(M))
                null_curve = (M * signs[:, None]).mean(0)
                null_directional[j] = np.max(null_curve[follows])
                null_omnibus[j] = np.max(np.abs(null_curve))
            p_directional = float(
                (1 + np.sum(null_directional >= obs_directional))
                / (1 + len(null_directional)))
            p_omnibus = float(
                (1 + np.sum(null_omnibus >= obs_omnibus))
                / (1 + len(null_omnibus)))
            inference_method = "20,000 Monte Carlo sign assignments"
        print(f"    prespecified Lecci-direction max-positive-r sign-flip test "
              f"({LECCI_XCORR_LAG_WINDOW_S[0]:.0f} to "
              f"{LECCI_XCORR_LAG_WINDOW_S[1]:.0f} s): p = {p_directional:.4f}")
        print(f"    two-sided max-|group r| sensitivity across all lags: p = {p_omnibus:.4f}")
        print(f"    sign-flip calculation: {inference_method}")
        print("    Opposite-sign coupling is not counted as support for Lecci's human direction.")
    else:
        print(f"    group inference unavailable: n={len(xc)}; "
              f"at least {_args.min_completed} participants are required.")
else:
    print("\n  cross-correlation unavailable for every participant.")

# ------------------------------------------------------------------ 3B
all_B = load_dir(B_DIR)
current_B = validate_run(B_DIR, all_B, "3B")
B_by_stage = {
    stage: [
        r for r in current_B
        if (r.get("endpoint_availability") or {}).get(stage) and r.get(stage)
    ]
    for stage in ("N2", "N3")
}
B = [
    r for r in current_B
    if any(r in B_by_stage[stage] for stage in ("N2", "N3"))
]
B_has_cohort_endpoint = (
    max(len(B_by_stage["N2"]), len(B_by_stage["N3"])) >= _args.min_completed)
if not B_has_cohort_endpoint:
    reasons = [
        f"{r['subject']}: {r.get('reason', 'unspecified')}"
        for r in current_B if r.get("status") != "ok"
    ]
    print(
        f"3B: N2 n={len(B_by_stage['N2'])}, N3 n={len(B_by_stage['N3'])}; at least "
        f"{_args.min_completed} estimable participants are required for one stage. "
        f"3B cohort result unavailable. Excluded/partial: {reasons}")
domains = {
    value.get("tachogram_domain")
    for stage in ("N2", "N3") for r in B_by_stage[stage] if (value := r.get(stage))
}
if B and domains != {"rr"}:
    raise SystemExit(f"3B: production summary requires the RR-domain estimator; found {domains}")
hr("3B -- NAJI-ALIGNED DESCRIPTIVE SO-RR ENDPOINTS")
partial_B = [r for r in current_B if r.get("status") != "ok"]
print(f"\nEndpoint denominators: N2 {len(B_by_stage['N2'])}/{len(EXPECTED)}; "
      f"N3 {len(B_by_stage['N3'])}/{len(EXPECTED)}; "
      f"{len(partial_B)} explicitly skipped/partial.")
for r in partial_B:
    print(f"  {r['subject']}: {r.get('reason', 'unspecified')}")
if B:
    print(f"\n{'subject':18s} {'N2 raw%':>8s} {'N2 local%':>9s} "
          f"{'N3 raw%':>8s} {'N3 local%':>9s} "
          f"{'N2 lag':>7s} {'N3 lag':>7s}")
    for r in B:
        v2, v3 = r.get("N2"), r.get("N3")
        print(f"{r['subject']:18s} "
              f"{(v2 or {}).get('pct_above_stage_mean', float('nan')):7.3f} "
              f"{(v2 or {}).get('event_locked_local_change_pct', float('nan')):9.3f} "
              f"{(v3 or {}).get('pct_above_stage_mean', float('nan')):7.3f} "
              f"{(v3 or {}).get('event_locked_local_change_pct', float('nan')):9.3f} "
              f"{(v2 or {}).get('peak_lag_s', float('nan')):7.2f} "
              f"{(v3 or {}).get('peak_lag_s', float('nan')):7.2f}")
    for st in ("N2", "N3"):
        stage_records = B_by_stage[st]
        if len(stage_records) < _args.min_completed:
            print(f"\n  {st}: cohort summary unavailable (n={len(stage_records)}; "
                  f"minimum {_args.min_completed})")
            continue
        v = np.array([r[st]["pct_above_stage_mean"] for r in stage_records], float)
        local = np.array([
            r[st]["event_locked_local_change_pct"] for r in stage_records
        ], float)
        nb = np.array([r[st]["null_mean_pct"] for r in stage_records], float)
        mean_v = float(np.mean(v))
        if mean_v > 0:
            ratio = mean_v / NAJI[st]
            direction = "below" if mean_v < NAJI[st] else "above"
            comparison = (
                f"observed/reference ratio {ratio:.2f}; "
                f"{abs(mean_v - NAJI[st]):.3f} percentage points {direction}")
        else:
            comparison = (
                f"ratio undefined because the cohort mean is nonpositive "
                f"(difference from Naji {mean_v - NAJI[st]:+.3f} percentage points)")
        print(f"\n  {st} (n={len(v)}): HR peak {mean_v:+.3f}% "
              f"(SD {np.std(v, ddof=1):.3f}), "
              f"Naji {NAJI[st]:+.2f}%  -> {comparison}")
        print(f"      local post-peak vs pre-event mean: median {np.median(local):+.3f}% "
              f"(descriptive)")
        excess = v - nb
        print(f"      raw maximum minus whole-stage-shift mean: median "
              f"{np.median(excess):+.3f}% (diagnostic only)")
        print("      event-locking inference disabled: whole-stage shifts do not preserve local "
              "nonstationary trends/event-density clustering.")
    # Raw post-trough maxima are count/noise biased: with fewer SOs the noisier average curve has a
    # larger expected maximum even under zero coupling.  Compare each stage only after subtracting
    # its own count- and estimator-matched surrogate mean.
    p2 = [(
              r["N2"].get(
                  "excess_over_null_pct",
                  r["N2"]["pct_above_stage_mean"] - r["N2"]["null_mean_pct"]),
              r["N3"].get(
                  "excess_over_null_pct",
                  r["N3"]["pct_above_stage_mean"] - r["N3"]["null_mean_pct"]),
          )
          for r in B if r.get("N2") and r.get("N3")]
    if len(p2) >= _args.min_completed:
        a2, a3 = np.array([x[0] for x in p2]), np.array([x[1] for x in p2])
        print("\n  N2-like vs N3-like diagnostic contrast (observed minus whole-stage-shift "
              "mean; no p value):")
        print(f"    N2 {np.mean(a2):+.3f}% vs N3 {np.mean(a3):+.3f}% (n = {len(p2)})")
        print("    Stage and event-locking inference are disabled pending validated staging and a "
              "local-trend/dependence-preserving null.")

print()
if not A_has_cohort_endpoint and not B_has_cohort_endpoint:
    raise SystemExit(
        "No 3A or 3B endpoint reached the prespecified cohort minimum; "
        "the unavailable-endpoint report above is the result.")
