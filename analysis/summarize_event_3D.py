"""Cohort summary for the EVENT-BASED 3D (SO->spindle), N2-like vs N3-like.

Same staging-quality filter as summarize_cohort_stages.py so the two are comparable:
a stage needs >= MIN_EVENTS spindle events and >= MIN_CH channels tested.
"""
import argparse, glob, json, os
import numpy as np
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    file_sha256,
    runtime_versions,
    source_tree_sha256,
)
from event_3d_estimators import (
    EVENT_FS,
    EVENT_PERCENTILE,
    IED_PAD_S,
    MIN_EVENTS as MIN_EVENTS_PER_CONTACT,
    MIN_POOLED_CONTACTS,
    MIN_POOLED_EVENTS,
    MIN_POOLED_NREM_COVERAGE,
    MIN_POOLED_VALID_S_PER_CONTACT,
    PRODUCTION_3D_INFERENCE_ENABLED,
    SO_BAND,
    SO_DUR,
    SPINDLE_BAND,
    SP_DUR,
    pooled_endpoint_passes_qc,
)
from event_3D_by_stage import production_config

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
parser = argparse.ArgumentParser()
parser.add_argument("--out", default=os.path.join(ROOT, "outputs", "event_3D_by_stage"))
parser.add_argument("--expected-subjects",
                    help="comma-separated exact IDs; defaults to the prespecified HUP cohort")
parser.add_argument("--min-completed", type=int, default=5,
                    help="minimum estimable participants required for cohort inference")
args = parser.parse_args()
OUT = args.out if os.path.isabs(args.out) else os.path.join(ROOT, args.out)
MIN_EVENTS = 200      # spindle events per stage
MIN_CH = 3            # channels tested per stage
EXPECTED = {
    f"HUP{n}_phaseII" for n in
    (116, 130, 133, 138, 139, 141, 143, 150, 151, 157, 160, 165, 171, 172, 173,
     177, 178, 182, 185, 187, 191, 199, 205, 211, 212)
}
if args.expected_subjects:
    EXPECTED = {
        value.strip() for value in args.expected_subjects.split(",") if value.strip()
    }

records = []
for f in sorted(glob.glob(os.path.join(OUT, "*.json"))):
    if os.path.basename(f) == "RUN_MANIFEST.json":
        continue
    d = json.load(open(f))
    if d.get("analysis_version") != ANALYSIS_VERSION or not d.get("subject"):
        continue
    d["_source_file"] = f
    records.append(d)

subject_ids = [d["subject"] for d in records]
duplicates = sorted({
    subject for subject in subject_ids if subject_ids.count(subject) > 1
})
if duplicates:
    raise SystemExit(f"3D duplicate subject result files: {duplicates}")
bad_filenames = [
    os.path.basename(d["_source_file"])
    for d in records
    if os.path.basename(d["_source_file"]) != f"{d['subject']}.json"
]
if bad_filenames:
    raise SystemExit(
        f"3D result filename must exactly match <subject>.json: {bad_filenames}")

manifest_path = os.path.join(OUT, "RUN_MANIFEST.json")
if not os.path.exists(manifest_path):
    raise SystemExit("3D RUN_MANIFEST.json is missing; no complete production run")
manifest = json.load(open(manifest_path))
if (manifest.get("run_state") != "complete" or not manifest.get("run_id")
        or manifest.get("pipeline") != "event_3D_by_stage"
        or manifest.get("analysis_version") != ANALYSIS_VERSION
        or manifest.get("cache_schema_version") != CACHE_SCHEMA_VERSION
        or manifest.get("failed")):
    raise SystemExit("3D run manifest is incomplete, failed, or version-mismatched")
if manifest.get("config") != production_config(7.0):
    raise SystemExit("3D run manifest does not contain the exact production configuration")
if manifest.get("runtime_versions") != runtime_versions():
    raise SystemExit("3D runtime/dependency versions differ from the completed run")
current_tree = source_tree_sha256(ROOT)
if manifest.get("source_tree_sha256") != current_tree:
    raise SystemExit("3D executable source tree differs from the completed run")
requested_values = manifest.get("requested", [])
requested = set(requested_values)
if len(requested_values) != len(requested):
    raise SystemExit("3D manifest requested list contains duplicate subjects")
if requested != EXPECTED:
    raise SystemExit(
        f"3D manifest requested set mismatch (missing={sorted(EXPECTED-requested)}, "
        f"unexpected={sorted(requested-EXPECTED)})")
completed_values = manifest.get("completed", [])
completed = set(completed_values)
if len(completed_values) != len(completed):
    raise SystemExit("3D manifest completed list contains duplicate subjects")
skipped_entries = manifest.get("skipped", [])
skipped_values = [
    value.get("subject") if isinstance(value, dict) else value
    for value in skipped_entries
]
skipped = set(skipped_values)
if any(not isinstance(value, dict) or not value.get("subject") or not value.get("reason")
       for value in skipped_entries):
    raise SystemExit("3D every skipped/partial subject needs an explicit reason")
if len(skipped_values) != len(skipped):
    raise SystemExit("3D manifest skipped list contains duplicate subjects")
if completed | skipped != requested or completed & skipped:
    raise SystemExit("3D requested subjects are not exactly completed or explicitly skipped")
record_subjects = {d["subject"] for d in records}
if record_subjects != completed | skipped:
    raise SystemExit(
        f"3D version-current files do not cover the manifest "
        f"(missing={sorted((completed | skipped)-record_subjects)}, "
        f"extra={sorted(record_subjects-(completed | skipped))})")
actual_result_hashes = {
    d["subject"]: file_sha256(d["_source_file"]) for d in records
}
if manifest.get("result_files_sha256") != actual_result_hashes:
    raise SystemExit("3D result file hashes do not match the completed run manifest")
if {d["subject"] for d in records if d.get("status") == "ok"} != completed:
    raise SystemExit("3D status=ok records do not equal the manifest completed set")
if {d.get("run_id") for d in records} != {manifest["run_id"]}:
    raise SystemExit("3D result/partial files are not linked to this completed run")
if {d.get("source_tree_sha256") for d in records} != {current_tree}:
    raise SystemExit("3D result/partial source digests do not match the completed run")

rows = []
for d in records:
    if d.get("status") != "ok":
        continue
    canonical = production_config(7.0)
    staging_qc = canonical["staging_qc"]
    staging_mask = d.get("staging_selected_contact_mask")
    staging_counts = d.get("staging_contact_count")
    staging_feature_coverage = d.get("staging_per_contact_feature_coverage")
    staging_observed_fraction = d.get("staging_candidate_observed_fraction")
    staging_clean_fraction = d.get("staging_candidate_clean_fraction")
    contact_candidates = d.get("lateral_contact_candidates")
    record_config_ok = (
        np.isclose(d.get("hours", np.nan), 7.0)
        and d.get("so_band_hz") == list(SO_BAND)
        and d.get("spindle_band_hz") == list(SPINDLE_BAND)
        and d.get("so_duration_s") == list(SO_DUR)
        and d.get("spindle_duration_s") == list(SP_DUR)
        and d.get("event_percentile") == EVENT_PERCENTILE
        and d.get("event_sampling_hz") == EVENT_FS
        and d.get("minimum_events_per_contact") == MIN_EVENTS_PER_CONTACT
        and d.get("threshold_scope") == canonical["threshold_scope"]
        and d.get("pairing") == canonical["pairing"]
        and d.get("ied_mask_padding_s") == IED_PAD_S
        and d.get("artifact_rejection")
        == f"IED/high-amplitude mask padded +/-{IED_PAD_S:.1f} s"
        and d.get("production_inference_enabled") is False
        and d.get("anatomy_selection_method")
        == canonical["anatomy_selection_method"]
        and isinstance(contact_candidates, list)
        and isinstance(staging_mask, list)
        and isinstance(staging_counts, list)
        and isinstance(staging_feature_coverage, list)
        and isinstance(staging_observed_fraction, list)
        and isinstance(staging_clean_fraction, list)
        and len(staging_mask) == len(contact_candidates)
        and len(staging_feature_coverage) == len(contact_candidates)
        and len(staging_observed_fraction) == len(contact_candidates)
        and len(staging_clean_fraction) == len(contact_candidates)
        and len(staging_counts) == int(round(float(d.get("hours", 0)) * 3600 / 30.0))
        and d.get("staging_n_selected_contacts")
        == int(np.sum(np.asarray(staging_mask, bool)))
        and d.get("staging_n_selected_contacts", 0)
        >= staging_qc["minimum_contacts"]
        and d.get("staging_required_contact_count", 0)
        >= staging_qc["minimum_contacts"]
        and d.get("staging_minimum_contact_feature_coverage")
        == staging_qc["minimum_contact_feature_coverage"]
        and all(
            not keep or coverage >= staging_qc["minimum_contact_feature_coverage"]
            for keep, coverage in zip(staging_mask, staging_feature_coverage))
        and all(
            not candidate
            or (
                observed >= staging_qc["minimum_candidate_observed_fraction"]
                and clean >= staging_qc["minimum_candidate_clean_fraction"]
            )
            for candidate, observed, clean in zip(
                staging_mask, staging_observed_fraction, staging_clean_fraction))
        and d.get("staging_swa_normalization") == staging_qc["swa_normalization"]
        and (d.get("ALL") or {}).get("inference_status", "").startswith("disabled:")
    )
    if not record_config_ok:
        raise SystemExit(
            f"3D {d['subject']} record does not agree with the production configuration")
    if not pooled_endpoint_passes_qc(d.get("ALL")):
        raise SystemExit(
            f"3D {d['subject']} is status=ok but fails the pooled-NREM quality gate")
    r = dict(sub=d["subject"].replace("_phaseII", ""), fsp=d.get("fsp"),
             n2_ep=d.get("n_N2"), n3_ep=d.get("n_N3"),
             run_id=d.get("run_id"), source_tree_sha256=d.get("source_tree_sha256"))
    for st in ("N2", "N3", "ALL"):
        s = d.get(st) or {}
        r[f"{st}_tested"] = s.get("n_channels_tested")
        r[f"{st}_R"] = s.get("median_R")
        r[f"{st}_participant_R"] = s.get("participant_R")
        r[f"{st}_vector_real"] = s.get("participant_vector_real")
        r[f"{st}_vector_imag"] = s.get("participant_vector_imag")
        r[f"{st}_ev"] = s.get("n_spindle_events")
        r[f"{st}_phase"] = s.get("participant_preferred_phase_deg")
    rows.append(r)

if len(rows) < args.min_completed:
    reasons = [
        f"{d['subject']}: {d.get('reason', 'unspecified')}"
        for d in records if d.get("status") != "ok"
    ]
    raise SystemExit(
        f"3D: only {len(rows)} completed pooled-NREM endpoints; at least "
        f"{args.min_completed} are required. Excluded/partial: {reasons}")

print(f"subjects with results: {len(rows)}\n")
print("pooled-NREM quality gate: "
      f">={MIN_POOLED_CONTACTS} contacts, >={MIN_POOLED_EVENTS} paired events, "
      f"and each included contact has >={MIN_POOLED_VALID_S_PER_CONTACT:g} s plus "
      f">={MIN_POOLED_NREM_COVERAGE:.0%} valid NREM coverage\n")
if skipped_entries:
    print("explicitly skipped/partial endpoints:")
    for value in skipped_entries:
        print(f"  {value['subject']}: {value['reason']}")
    print()
hdr = (f"{'sub':9s} | {'N2 pR':>7s}{'N2 chR':>8s}{'N2 ev':>8s} | "
       f"{'N3 pR':>7s}{'N3 chR':>8s}{'N3 ev':>8s} | {'ALL pR':>8s}{'ALL chR':>9s}")
print(hdr); print("-" * len(hdr))
for r in rows:
    f2 = lambda v, p=3: (f"{v:.{p}f}" if isinstance(v, float) else "  n/a")
    ok = (r["N2_ev"] or 0) >= MIN_EVENTS and (r["N3_ev"] or 0) >= MIN_EVENTS
    print(f"{r['sub']:9s} | {f2(r['N2_participant_R']):>7s}"
          f"{f2(r['N2_R']):>8s}{r['N2_ev'] or 0:8d} | "
          f"{f2(r['N3_participant_R']):>7s}{f2(r['N3_R']):>8s}"
          f"{r['N3_ev'] or 0:8d} | {f2(r['ALL_participant_R']):>8s}"
          f"{f2(r['ALL_R']):>9s}" + ("" if ok else "   <- stage contrast excluded"))

use = [r for r in rows if (r["N2_ev"] or 0) >= MIN_EVENTS and (r["N3_ev"] or 0) >= MIN_EVENTS
       and (r["N2_tested"] or 0) >= MIN_CH and (r["N3_tested"] or 0) >= MIN_CH]
print(f"\npassing filter (>= {MIN_EVENTS} spindle events and >= {MIN_CH} channels per stage): {len(use)}/{len(rows)}")

# Pooled-NREM vectors use one vector per participant. They are descriptive until the complete
# pairing estimator is rerun inside a calibrated time-shift/block null.
allr = [r for r in rows if r["ALL_vector_real"] is not None and r["ALL_vector_imag"] is not None]
if not allr:
    print("\nPooled-NREM descriptive endpoint unavailable: the complete run contains no estimable "
          "participant vectors. See the explicit endpoint reasons above.")
    raise SystemExit(0)
vectors = np.asarray([
    complex(r["ALL_vector_real"], r["ALL_vector_imag"]) for r in allr
])
print("\n=== Pooled-NREM participant vectors (descriptive only) ===")
if PRODUCTION_3D_INFERENCE_ENABLED:
    raise RuntimeError("unexpected: production 3D inference must remain disabled")
group_vector = complex(np.mean(vectors))
print(f"  participants passing QC                : {len(vectors)}")
print(f"  descriptive group vector length        : {abs(group_vector):.4f}")
print(f"  descriptive group preferred SO phase   : "
      f"{np.degrees(np.angle(group_vector)):+.1f} deg")
print(f"  median participant vector length       : "
      f"{np.median([r['ALL_participant_R'] for r in allr]):.3f}")
print(f"  median descriptive contact R           : "
      f"{np.median([r['ALL_R'] for r in allr]):.3f}")
print("  INFERENCE DISABLED: the finite +/-2 s SO-centered pairing window creates a")
print("  shared preferred phase even for independent spindle/SO trains. A valid null must")
print("  shift or block-resample complete spindle trains relative to SOs, then repeat event")
print("  pairing, contact aggregation, and the cohort statistic.")

# N2 vs N3
print("\n=== N2-like vs N3-like contrast ===")
print("  INFERENCE DISABLED. Stage estimates can use different contacts and event counts, and")
print("  serial dependence invalidates the iid finite-n R correction. A valid contrast requires")
print("  the same-contact intersection plus a dependence-preserving within-participant block/stage")
print("  permutation (and validated staging); the current per-stage values are descriptive only.")
