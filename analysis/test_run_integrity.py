"""Regression checks for cache lineage and fail-closed production summaries."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

import numpy as np

from event_3D_by_stage import production_config
from lecci_faithful_3A import CACHE, CacheSubjectSkipped, load, set_cache
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    atomic_json_dump,
    atomic_savez,
    cache_code_sha256,
    cache_is_current,
    file_sha256,
    source_tree_sha256,
    validated_complete_run_exists,
    write_run_manifest,
)


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


def artifact_hashes(directory, subjects, suffix):
    return {
        subject: file_sha256(os.path.join(directory, f"{subject}{suffix}"))
        for subject in subjects
    }


def build_cache(cache_dir, subjects, hours=7.0):
    """Create current-schema synthetic cache artifacts and a complete cache manifest."""
    digest = cache_code_sha256(ROOT)
    os.makedirs(cache_dir, exist_ok=True)
    for subject in subjects:
        atomic_savez(
            os.path.join(cache_dir, f"{subject}.npz"),
            subject=subject,
            status="ok",
            cache_schema_version=CACHE_SCHEMA_VERSION,
            cache_code_sha256=digest,
            hours=float(hours),
        )
    write_run_manifest(
        cache_dir,
        pipeline="cache_lc_series",
        requested=subjects,
        completed=subjects,
        skipped=[],
        failed=[],
        config={
            "hours": float(hours),
            "cache_schema_version": CACHE_SCHEMA_VERSION,
            "cache_code_sha256": digest,
        },
        run_id="synthetic-cache-run",
        run_state="complete",
        result_files_sha256=artifact_hashes(cache_dir, subjects, ".npz"),
    )


def lineage(cache_dir, subject):
    with open(os.path.join(cache_dir, "RUN_MANIFEST.json")) as handle:
        manifest = json.load(handle)
    with np.load(os.path.join(cache_dir, f"{subject}.npz"), allow_pickle=False) as cached:
        hours = float(np.asarray(cached["hours"]).item())
    return dict(
        cache_directory_relative=os.path.relpath(cache_dir, ROOT).replace(os.sep, "/"),
        cache_manifest_run_id=manifest["run_id"],
        cache_manifest_sha256=file_sha256(os.path.join(cache_dir, "RUN_MANIFEST.json")),
        cache_file_sha256=file_sha256(os.path.join(cache_dir, f"{subject}.npz")),
        cache_hours=hours,
    )


def config_3a(cache_inputs, *, smooth_4s=True):
    return dict(
        band="fixed",
        smooth_4s=smooth_4s,
        n_sur=200,
        analysis_version=ANALYSIS_VERSION,
        cache_schema_version=CACHE_SCHEMA_VERSION,
        cache_code_sha256=cache_code_sha256(ROOT),
        cache_inputs=cache_inputs,
    )


def config_3b(cache_inputs):
    return dict(
        tachogram_domain="rr",
        n_surrogates=1000,
        analysis_version=ANALYSIS_VERSION,
        cache_schema_version=CACHE_SCHEMA_VERSION,
        cache_code_sha256=cache_code_sha256(ROOT),
        null_method="shared circular shift in eligible stage-time",
        contact_qc="cache stable >=80%-coverage sigma-contact intersection",
        minimum_event_channels=2,
        cache_inputs=cache_inputs,
    )


with tempfile.TemporaryDirectory() as tmp:
    cache_path = os.path.join(tmp, "cache.npz")
    current_digest = cache_code_sha256(ROOT)
    atomic_savez(
        cache_path,
        status="ok",
        cache_schema_version=CACHE_SCHEMA_VERSION,
        cache_code_sha256=current_digest,
    )
    check("cache accepts the current schema and exact cache-building source digest",
          cache_is_current(cache_path, ROOT))
    atomic_savez(
        cache_path,
        status="ok",
        cache_schema_version=CACHE_SCHEMA_VERSION,
        cache_code_sha256="forged-or-stale",
    )
    check("cache rejects a matching schema with stale/forged builder lineage",
          not cache_is_current(cache_path, ROOT))

    strict_path = os.path.join(tmp, "strict.json")
    atomic_json_dump({"finite": 1.0}, strict_path)
    try:
        atomic_json_dump({"nonfinite": np.nan}, strict_path)
    except ValueError:
        strict_rejected = True
    else:
        strict_rejected = False
    with open(strict_path) as handle:
        preserved = json.load(handle)
    check("JSON writer rejects NaN and atomically preserves the prior valid artifact",
          strict_rejected and preserved == {"finite": 1.0})

with tempfile.TemporaryDirectory(prefix=".reuse-cache-", dir=ROOT) as tmp:
    subjects = ["HUP900_phaseII", "HUP901_phaseII"]
    build_cache(tmp, subjects)
    reuse_config = {
        "hours": 7.0,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "cache_code_sha256": cache_code_sha256(ROOT),
    }
    check("complete cache reuse verifies the prior manifest and every NPZ hash",
          validated_complete_run_exists(
              tmp, pipeline="cache_lc_series", requested=subjects,
              config=reuse_config, suffix=".npz", require_current_source_tree=False))
    try:
        validated_complete_run_exists(
            tmp, pipeline="cache_lc_series", requested=[subjects[0], subjects[0]],
            config=reuse_config, suffix=".npz", require_current_source_tree=False)
    except RuntimeError as exc:
        duplicate_request_rejected = "duplicates" in str(exc)
    else:
        duplicate_request_rejected = False
    check("exact-run reuse rejects a duplicate requested subject",
          duplicate_request_rejected)
    atomic_savez(
        os.path.join(tmp, f"{subjects[0]}.npz"),
        subject=subjects[0],
        status="ok",
        cache_schema_version=CACHE_SCHEMA_VERSION,
        cache_code_sha256=cache_code_sha256(ROOT),
        hours=7.0,
        tampered_before_reuse=True,
    )
    try:
        validated_complete_run_exists(
            tmp, pipeline="cache_lc_series", requested=subjects,
            config=reuse_config, suffix=".npz", require_current_source_tree=False)
    except RuntimeError as exc:
        tampered_reuse_rejected = "differ" in str(exc)
    else:
        tampered_reuse_rejected = False
    check("rerun cannot bless a cache modified after its prior terminal manifest",
          tampered_reuse_rejected)
    set_cache(tmp)
    try:
        try:
            load(subjects[0])
        except RuntimeError as exc:
            preanalysis_tamper_rejected = "terminal cache manifest" in str(exc)
        else:
            preanalysis_tamper_rejected = False
    finally:
        set_cache(CACHE)
    check("downstream analysis rejects cache bytes altered before analysis",
          preanalysis_tamper_rejected)

with tempfile.TemporaryDirectory(prefix=".swapped-cache-", dir=ROOT) as tmp:
    requested = "HUP130_phaseII"
    digest = cache_code_sha256(ROOT)
    atomic_savez(
        os.path.join(tmp, f"{requested}.npz"),
        subject="HUP116_phaseII",
        status="ok",
        cache_schema_version=CACHE_SCHEMA_VERSION,
        cache_code_sha256=digest,
        hours=1.0,
    )
    write_run_manifest(
        tmp,
        pipeline="cache_lc_series",
        requested=[requested],
        completed=[requested],
        skipped=[],
        failed=[],
        config={"cache_code_sha256": digest},
        run_id="swapped-cache-run",
        run_state="complete",
        result_files_sha256=artifact_hashes(tmp, [requested], ".npz"),
    )
    set_cache(tmp)
    try:
        try:
            load(requested)
        except RuntimeError as exc:
            swapped_rejected = "embeds subject" in str(exc)
        else:
            swapped_rejected = False
    finally:
        set_cache(CACHE)
    check("cache load rejects a current-schema file copied under another participant name",
          swapped_rejected)

with tempfile.TemporaryDirectory(prefix=".skipped-cache-", dir=ROOT) as tmp:
    subject = "HUP999_phaseII"
    reason = "synthetic prespecified exclusion"
    digest = cache_code_sha256(ROOT)
    atomic_savez(
        os.path.join(tmp, f"{subject}.npz"),
        subject=subject,
        status="skip",
        reason=reason,
        cache_schema_version=CACHE_SCHEMA_VERSION,
        cache_code_sha256=digest,
        hours=1.0,
    )
    write_run_manifest(
        tmp,
        pipeline="cache_lc_series",
        requested=[subject],
        completed=[],
        skipped=[dict(subject=subject, reason=reason)],
        failed=[],
        config={"cache_code_sha256": digest},
        run_id="skipped-cache-run",
        run_state="complete",
        result_files_sha256=artifact_hashes(tmp, [subject], ".npz"),
    )
    set_cache(tmp)
    try:
        try:
            load(subject)
        except CacheSubjectSkipped as exc:
            skip_propagated = str(exc) == reason
        else:
            skip_propagated = False
    finally:
        set_cache(CACHE)
    check("manifest-matched cache exclusion propagates as a structured subject skip",
          skip_propagated)

with tempfile.TemporaryDirectory(prefix=".run-integrity-", dir=ROOT) as tmp:
    a_dir = os.path.join(tmp, "a")
    b_dir = os.path.join(tmp, "b")
    d_dir = os.path.join(tmp, "d")
    os.makedirs(a_dir)
    os.makedirs(b_dir)
    os.makedirs(d_dir)

    no_manifest_3ab = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_corrected_3AB.py"),
            "--a-dir", a_dir,
            "--b-dir", b_dir,
            "--expected-subjects", "synthetic-subject",
        ],
        text=True,
        capture_output=True,
    )
    check("3A/3B summary rejects an empty directory without a complete manifest",
          no_manifest_3ab.returncode != 0 and "RUN_MANIFEST" in (
              no_manifest_3ab.stdout + no_manifest_3ab.stderr))

    no_manifest_3d = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_event_3D.py"),
            "--out", d_dir,
            "--expected-subjects", "synthetic-subject",
        ],
        text=True,
        capture_output=True,
    )
    check("3D summary rejects an empty directory without a complete manifest",
          no_manifest_3d.returncode != 0 and "RUN_MANIFEST" in (
              no_manifest_3d.stdout + no_manifest_3d.stderr))

    write_run_manifest(
        a_dir,
        pipeline="lecci_faithful_3A",
        requested=["synthetic-subject"],
        completed=[],
        skipped=[],
        failed=[],
        config={},
        run_id="synthetic-run",
        run_state="in_progress",
    )
    interrupted_3ab = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_corrected_3AB.py"),
            "--a-dir", a_dir,
            "--b-dir", b_dir,
            "--expected-subjects", "synthetic-subject",
        ],
        text=True,
        capture_output=True,
    )
    check("production summary rejects an interrupted run even when no result files exist",
          interrupted_3ab.returncode != 0 and "incomplete" in (
              interrupted_3ab.stdout + interrupted_3ab.stderr).lower())

    tree_digest = source_tree_sha256(ROOT)
    cache_dir = os.path.join(tmp, "cache")
    build_cache(cache_dir, ["synthetic-subject"])
    partial_lineage = lineage(cache_dir, "synthetic-subject")
    partial = dict(
        subject="synthetic-subject",
        status="partial",
        reason="synthetic unavailable endpoint",
        analysis_version=ANALYSIS_VERSION,
        cache_schema_version=CACHE_SCHEMA_VERSION,
        run_id="partial-run",
        source_tree_sha256=tree_digest,
        band="fixed",
        smooth_4s=True,
        n_surrogates_requested=200,
        **partial_lineage,
    )
    atomic_json_dump(partial, os.path.join(a_dir, "synthetic-subject.json"))
    write_run_manifest(
        a_dir,
        pipeline="lecci_faithful_3A",
        requested=["synthetic-subject"],
        completed=[],
        skipped=[dict(subject="synthetic-subject", reason=partial["reason"])],
        failed=[],
        config=config_3a({"synthetic-subject": partial_lineage}),
        run_id="partial-run",
        run_state="complete",
        result_files_sha256=artifact_hashes(
            a_dir, ["synthetic-subject"], ".json"),
    )
    zero_endpoint_3ab = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_corrected_3AB.py"),
            "--a-dir", a_dir,
            "--b-dir", b_dir,
            "--expected-subjects", "synthetic-subject",
        ],
        text=True,
        capture_output=True,
    )
    check("3A/3B publication summary rejects a complete manifest with zero estimable endpoints",
          zero_endpoint_3ab.returncode != 0 and "spectrum n=0" in (
              zero_endpoint_3ab.stdout + zero_endpoint_3ab.stderr).lower())

    with open(os.path.join(a_dir, "RUN_MANIFEST.json")) as handle:
        wrong_runtime_manifest = json.load(handle)
    wrong_runtime_manifest["runtime_versions"]["python"] = "0.0-incompatible"
    atomic_json_dump(
        wrong_runtime_manifest, os.path.join(a_dir, "RUN_MANIFEST.json"))
    wrong_runtime_3ab = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_corrected_3AB.py"),
            "--a-dir", a_dir,
            "--b-dir", b_dir,
            "--expected-subjects", "synthetic-subject",
        ],
        text=True,
        capture_output=True,
    )
    check("publication summary rejects results from a different runtime/dependency set",
          wrong_runtime_3ab.returncode != 0 and "runtime" in (
              wrong_runtime_3ab.stdout + wrong_runtime_3ab.stderr).lower())
    write_run_manifest(
        a_dir,
        pipeline="lecci_faithful_3A",
        requested=["synthetic-subject"],
        completed=[],
        skipped=[dict(subject="synthetic-subject", reason=partial["reason"])],
        failed=[],
        config=config_3a({"synthetic-subject": partial_lineage}),
        run_id="partial-run",
        run_state="complete",
        result_files_sha256=artifact_hashes(
            a_dir, ["synthetic-subject"], ".json"),
    )

    partial["reason"] = "tampered downstream statistic/metadata"
    atomic_json_dump(partial, os.path.join(a_dir, "synthetic-subject.json"))
    changed_result_3ab = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_corrected_3AB.py"),
            "--a-dir", a_dir,
            "--b-dir", b_dir,
            "--expected-subjects", "synthetic-subject",
        ],
        text=True,
        capture_output=True,
    )
    check("3A/3B summary rejects a result JSON modified after the terminal manifest",
          changed_result_3ab.returncode != 0 and "result file hashes" in (
              changed_result_3ab.stdout + changed_result_3ab.stderr).lower())
    partial["reason"] = "synthetic unavailable endpoint"
    atomic_json_dump(partial, os.path.join(a_dir, "synthetic-subject.json"))

    atomic_savez(
        os.path.join(cache_dir, "synthetic-subject.npz"),
        subject="synthetic-subject",
        status="ok",
        cache_schema_version=CACHE_SCHEMA_VERSION,
        cache_code_sha256=cache_code_sha256(ROOT),
        hours=7.0,
        tampered_after_analysis=True,
    )
    changed_cache_3ab = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_corrected_3AB.py"),
            "--a-dir", a_dir,
            "--b-dir", b_dir,
            "--expected-subjects", "synthetic-subject",
        ],
        text=True,
        capture_output=True,
    )
    check("publication summary rejects a cache file changed after downstream analysis",
          changed_cache_3ab.returncode != 0 and "cache file content changed" in (
              changed_cache_3ab.stdout + changed_cache_3ab.stderr).lower())

    atomic_json_dump(partial, os.path.join(d_dir, "synthetic-subject.json"))
    write_run_manifest(
        d_dir,
        pipeline="event_3D_by_stage",
        requested=["synthetic-subject"],
        completed=[],
        skipped=[dict(subject="synthetic-subject", reason=partial["reason"])],
        failed=[],
        config=production_config(7.0),
        run_id="partial-run",
        run_state="complete",
        result_files_sha256=artifact_hashes(
            d_dir, ["synthetic-subject"], ".json"),
    )
    zero_endpoint_3d = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_event_3D.py"),
            "--out", d_dir,
            "--expected-subjects", "synthetic-subject",
        ],
        text=True,
        capture_output=True,
    )
    check("3D publication summary rejects a complete manifest with zero estimable endpoints",
          zero_endpoint_3d.returncode != 0 and "only 0 completed" in (
              zero_endpoint_3d.stdout + zero_endpoint_3d.stderr).lower())

    write_run_manifest(
        d_dir,
        pipeline="event_3D_by_stage",
        requested=["synthetic-subject"],
        completed=[],
        skipped=[dict(subject="synthetic-subject", reason=partial["reason"])],
        failed=[],
        config=production_config(6.0),
        run_id="partial-run",
        run_state="complete",
        result_files_sha256=artifact_hashes(
            d_dir, ["synthetic-subject"], ".json"),
    )
    invalid_3d_config = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_event_3D.py"),
            "--out", d_dir,
            "--expected-subjects", "synthetic-subject",
        ],
        text=True,
        capture_output=True,
    )
    check("3D summary rejects a non-production duration/configuration",
          invalid_3d_config.returncode != 0 and "exact production configuration" in (
              invalid_3d_config.stdout + invalid_3d_config.stderr).lower())
    write_run_manifest(
        d_dir,
        pipeline="event_3D_by_stage",
        requested=["synthetic-subject"],
        completed=[],
        skipped=[dict(subject="synthetic-subject", reason=partial["reason"])],
        failed=[],
        config=production_config(7.0),
        run_id="partial-run",
        run_state="complete",
        result_files_sha256=artifact_hashes(
            d_dir, ["synthetic-subject"], ".json"),
    )

    partial["reason"] = "tampered 3D result"
    atomic_json_dump(partial, os.path.join(d_dir, "synthetic-subject.json"))
    changed_result_3d = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_event_3D.py"),
            "--out", d_dir,
            "--expected-subjects", "synthetic-subject",
        ],
        text=True,
        capture_output=True,
    )
    check("3D summary rejects a result JSON modified after the terminal manifest",
          changed_result_3d.returncode != 0 and "result file hashes" in (
              changed_result_3d.stdout + changed_result_3d.stderr).lower())
    partial["reason"] = "synthetic unavailable endpoint"
    atomic_json_dump(partial, os.path.join(d_dir, "synthetic-subject.json"))

    dup_a = os.path.join(tmp, "duplicate-a")
    dup_b = os.path.join(tmp, "duplicate-b")
    dup_d = os.path.join(tmp, "duplicate-d")
    os.makedirs(dup_a)
    os.makedirs(dup_b)
    os.makedirs(dup_d)
    current = dict(
        subject="synthetic-subject",
        status="ok",
        analysis_version=ANALYSIS_VERSION,
        cache_schema_version=CACHE_SCHEMA_VERSION,
        run_id="duplicate-run",
        source_tree_sha256=tree_digest,
    )
    for directory in (dup_a, dup_d):
        atomic_json_dump(current, os.path.join(directory, "synthetic-subject.json"))
        atomic_json_dump(current, os.path.join(directory, "duplicate.json"))
        write_run_manifest(
            directory,
            pipeline="synthetic",
            requested=["synthetic-subject"],
            completed=["synthetic-subject"],
            skipped=[],
            failed=[],
            config={},
            run_id="duplicate-run",
            run_state="complete",
        )
    duplicate_3ab = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_corrected_3AB.py"),
            "--a-dir", dup_a,
            "--b-dir", dup_b,
            "--expected-subjects", "synthetic-subject",
            "--min-completed", "0",
        ],
        text=True,
        capture_output=True,
    )
    check("3A/3B summary rejects duplicate files for one subject",
          duplicate_3ab.returncode != 0 and "duplicate subject" in (
              duplicate_3ab.stdout + duplicate_3ab.stderr).lower())
    duplicate_3d = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_event_3D.py"),
            "--out", dup_d,
            "--expected-subjects", "synthetic-subject",
            "--min-completed", "0",
        ],
        text=True,
        capture_output=True,
    )
    check("3D summary rejects duplicate files for one subject",
          duplicate_3d.returncode != 0 and "duplicate subject" in (
              duplicate_3d.stdout + duplicate_3d.stderr).lower())

    endpoint_a = os.path.join(tmp, "endpoint-a")
    endpoint_b = os.path.join(tmp, "endpoint-b")
    endpoint_cache = os.path.join(tmp, "endpoint-cache")
    os.makedirs(endpoint_a)
    os.makedirs(endpoint_b)
    subjects = [f"synthetic-{i}" for i in range(5)]
    build_cache(endpoint_cache, subjects)
    cache_inputs_a, cache_inputs_b = {}, {}
    skipped_a, skipped_b = [], []
    for subject in subjects:
        own_peak_available = subject in subjects[:2]
        cache_input = lineage(endpoint_cache, subject)
        reason_a = "unavailable primary endpoint(s): sigma-HR cross-correlation"
        record_a = dict(
            subject=subject,
            status="partial",
            reason=reason_a,
            analysis_version=ANALYSIS_VERSION,
            cache_schema_version=CACHE_SCHEMA_VERSION,
            run_id="endpoint-a-run",
            source_tree_sha256=tree_digest,
            band="fixed",
            smooth_4s=True,
            n_surrogates_requested=200,
            n_bouts=1,
            bout_seconds=120,
            peak=(
                {
                    "peak_hz": 0.02,
                    "prominence_over_background": 2.0,
                    "spectral_sd_hz": 0.008,
                }
                if own_peak_available else {"peak_hz": None}
            ),
            peak_swa={"peak_hz": None},
            negative_control=None,
            peak_null=None,
            coherence={
                "K": 10,
                "crit": 0.3,
                "n_valid": 1000,
                "filled_frac": 0.0,
                "at_lecci": 0.5,
                "sig_at_lecci": True,
                "at_own_peak": 0.5 if own_peak_available else None,
                "sig_at_own_peak": True if own_peak_available else None,
                "f": [0.01, 0.02, 0.03],
                "cxy": [0.2, 0.5, None if subject == subjects[-1] else 0.2],
            },
            xcorr=None,
            endpoint_availability={
                "spectrum": True,
                "cross_correlation": False,
                "coherence": True,
                "coherence_fixed_0p02": True,
                "coherence_own_peak": own_peak_available,
            },
            **cache_input,
        )
        atomic_json_dump(record_a, os.path.join(endpoint_a, f"{subject}.json"))
        cache_inputs_a[subject] = cache_input
        skipped_a.append(dict(subject=subject, reason=reason_a))

        reason_b = "no estimable N3-like SO-RR endpoint"
        record_b = dict(
            subject=subject,
            status="partial",
            reason=reason_b,
            analysis_version=ANALYSIS_VERSION,
            cache_schema_version=CACHE_SCHEMA_VERSION,
            run_id="endpoint-b-run",
            source_tree_sha256=tree_digest,
            tachogram_domain="rr",
            n_surrogates=1000,
            null_method="shared circular shift in eligible stage-time",
            minimum_event_channels=2,
            N2={
                "tachogram_domain": "rr",
                "n_surrogates": 1000,
                "n_channels": 2,
                "event_contact_ids": ["C1", "C2"],
                "pct_above_stage_mean": 1.0,
                "event_locked_local_change_pct": 0.1,
                "peak_lag_s": 1.0,
                "null_mean_pct": 0.2,
            },
            N3=None,
            endpoint_availability={"N2": True, "N3": False},
            **cache_input,
        )
        atomic_json_dump(record_b, os.path.join(endpoint_b, f"{subject}.json"))
        cache_inputs_b[subject] = cache_input
        skipped_b.append(dict(subject=subject, reason=reason_b))
    write_run_manifest(
        endpoint_a,
        pipeline="lecci_faithful_3A",
        requested=subjects,
        completed=[],
        skipped=skipped_a,
        failed=[],
        config=config_3a(cache_inputs_a),
        run_id="endpoint-a-run",
        run_state="complete",
        result_files_sha256=artifact_hashes(endpoint_a, subjects, ".json"),
    )
    write_run_manifest(
        endpoint_b,
        pipeline="event_3B_cached",
        requested=subjects,
        completed=[],
        skipped=skipped_b,
        failed=[],
        config=config_3b(cache_inputs_b),
        run_id="endpoint-b-run",
        run_state="complete",
        result_files_sha256=artifact_hashes(endpoint_b, subjects, ".json"),
    )
    endpoint_specific = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_corrected_3AB.py"),
            "--a-dir", endpoint_a,
            "--b-dir", endpoint_b,
            "--expected-subjects", ",".join(subjects),
            "--min-completed", "5",
        ],
        text=True,
        capture_output=True,
    )
    check("partial records contribute only to their available 3A/3B endpoint denominators",
          endpoint_specific.returncode == 0
          and "spectrum 5/5; cross-correlation 0/5" in endpoint_specific.stdout
          and "fixed-0.02-Hz coherence 5/5; own-peak coherence 2/5"
          in endpoint_specific.stdout
          and "significant at own fitted peak    : 2/2" in endpoint_specific.stdout
          and "endpoint denominators: n2 5/5; n3 0/5"
          in endpoint_specific.stdout.lower())

    first_a_path = os.path.join(endpoint_a, f"{subjects[0]}.json")
    with open(first_a_path) as handle:
        first_a = json.load(handle)
    first_a["n_surrogates_requested"] = 199
    atomic_json_dump(first_a, first_a_path)
    write_run_manifest(
        endpoint_a,
        pipeline="lecci_faithful_3A",
        requested=subjects,
        completed=[],
        skipped=skipped_a,
        failed=[],
        config=config_3a(cache_inputs_a),
        run_id="endpoint-a-run",
        run_state="complete",
        result_files_sha256=artifact_hashes(endpoint_a, subjects, ".json"),
    )
    record_manifest_disagreement = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_corrected_3AB.py"),
            "--a-dir", endpoint_a,
            "--b-dir", endpoint_b,
            "--expected-subjects", ",".join(subjects),
            "--min-completed", "5",
        ],
        text=True,
        capture_output=True,
    )
    check("production manifest cannot hide a non-production subject-record estimator",
          record_manifest_disagreement.returncode != 0
          and "subject record" in (
              record_manifest_disagreement.stdout
              + record_manifest_disagreement.stderr).lower())
    first_a["n_surrogates_requested"] = 200
    atomic_json_dump(first_a, first_a_path)
    write_run_manifest(
        endpoint_a,
        pipeline="lecci_faithful_3A",
        requested=subjects,
        completed=[],
        skipped=skipped_a,
        failed=[],
        config=config_3a(cache_inputs_a),
        run_id="endpoint-a-run",
        run_state="complete",
        result_files_sha256=artifact_hashes(endpoint_a, subjects, ".json"),
    )

    first_a["peak"]["peak_hz"] = 0.055
    atomic_json_dump(first_a, first_a_path)
    statistical_tamper = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_corrected_3AB.py"),
            "--a-dir", endpoint_a,
            "--b-dir", endpoint_b,
            "--expected-subjects", ",".join(subjects),
            "--min-completed", "5",
        ],
        text=True,
        capture_output=True,
    )
    check("terminal hashes protect a downstream numerical statistic, not only metadata",
          statistical_tamper.returncode != 0 and "result file hashes" in (
              statistical_tamper.stdout + statistical_tamper.stderr).lower())
    first_a["peak"]["peak_hz"] = 0.02
    atomic_json_dump(first_a, first_a_path)

    write_run_manifest(
        endpoint_a,
        pipeline="lecci_faithful_3A",
        requested=subjects,
        completed=[],
        skipped=skipped_a,
        failed=[],
        config=config_3a(cache_inputs_a, smooth_4s=False),
        run_id="endpoint-a-run",
        run_state="complete",
        result_files_sha256=artifact_hashes(endpoint_a, subjects, ".json"),
    )
    bad_production_config = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "analysis", "summarize_corrected_3AB.py"),
            "--a-dir", endpoint_a,
            "--b-dir", endpoint_b,
            "--expected-subjects", ",".join(subjects),
            "--min-completed", "5",
        ],
        text=True,
        capture_output=True,
    )
    check("publication summary rejects a homogeneous but non-production 3A configuration",
          bad_production_config.returncode != 0 and "non-production" in (
              bad_production_config.stdout + bad_production_config.stderr).lower())

print("ALL RUN-INTEGRITY CHECKS PASSED")
