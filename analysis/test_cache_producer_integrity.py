"""Focused regressions for HUP acquisition and cache provenance hardening."""
from __future__ import annotations

import json
import os
import shutil
import tempfile

import numpy as np

from cache_lc_series import (
    conservative_resampled_interval,
    detect_so_candidates,
    event_3d_so_candidates,
    retain_complete_clean_event_3d_candidates,
)
from cache_paired_scalp import (
    IEEG_CACHE_SCHEMA,
    PIPELINE as SCALP_PIPELINE,
    SCALP_CACHE_SCHEMA,
    SCALP_CHANNEL_PLAN,
    _DEPENDENCY_FILES as SCALP_DEPENDENCY_FILES,
    _manifest_config,
    _write_in_progress_manifest,
    assert_pinned_ieeg_inputs_unchanged,
    cache_dependency_sha256,
    freeze_pinned_ieeg_inputs,
    validate_reusable_sidecar_run,
)
from hup_portal import (
    PortalSampleCountMismatch,
    find_night,
    pull_continuous_exact,
    validate_hup_series_geometry,
)
from finalize_retried_cache_run import finalize
from pipeline_version import (
    _CACHE_SOURCE_FILES,
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    atomic_json_dump,
    atomic_savez,
    cache_code_sha256,
    file_sha256,
    require_integer_sample_rate,
    runtime_versions,
    write_run_manifest,
)
from run_qc_grid import _cache_manifest


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


check(
    "cache digest includes the active HUP acquisition/source helper",
    "analysis/hup_portal.py" in _CACHE_SOURCE_FILES,
)

check(
    "cache digest includes the shared 3D estimator source",
    "analysis/event_3d_estimators.py" in _CACHE_SOURCE_FILES,
)

CANONICAL_HELPER_SOURCES = {
    "analysis/signal_qc.py",
    "analysis/staging_helpers.py",
}
WITHDRAWN_WRAPPER_SOURCES = {
    "analysis/results_3A_tutorial_style.py",
    "analysis/cohort_stages_3ABD.py",
}
check(
    "primary cache digest follows canonical helper implementations, not withdrawn wrappers",
    CANONICAL_HELPER_SOURCES <= set(_CACHE_SOURCE_FILES)
    and not WITHDRAWN_WRAPPER_SOURCES & set(_CACHE_SOURCE_FILES),
)
check(
    "paired-scalp digest follows canonical helper implementations, not withdrawn wrappers",
    CANONICAL_HELPER_SOURCES <= set(SCALP_DEPENDENCY_FILES)
    and not WITHDRAWN_WRAPPER_SOURCES & set(SCALP_DEPENDENCY_FILES),
)

with tempfile.TemporaryDirectory() as digest_root:
    for relative_path in _CACHE_SOURCE_FILES:
        source = os.path.join(ROOT, relative_path)
        destination = os.path.join(digest_root, relative_path)
        os.makedirs(os.path.dirname(destination), exist_ok=True)
        shutil.copyfile(source, destination)
    digest_before_3d_change = cache_code_sha256(digest_root)
    estimator_copy = os.path.join(
        digest_root, "analysis", "event_3d_estimators.py")
    with open(estimator_copy, "ab") as handle:
        handle.write(b"\n# synthetic estimator change\n")
    digest_after_3d_change = cache_code_sha256(digest_root)
    helper_digest_changes = []
    digest_before_helper_change = digest_after_3d_change
    for relative_path in sorted(CANONICAL_HELPER_SOURCES):
        helper_copy = os.path.join(digest_root, relative_path)
        with open(helper_copy, "ab") as handle:
            handle.write(b"\n# synthetic helper change\n")
        digest_after_helper_change = cache_code_sha256(digest_root)
        helper_digest_changes.append(
            digest_after_helper_change != digest_before_helper_change)
        digest_before_helper_change = digest_after_helper_change
check(
    "changing the shared 3D estimator invalidates the cache digest",
    digest_before_3d_change != digest_after_3d_change,
)
check(
    "changing either canonical helper source invalidates the cache digest",
    all(helper_digest_changes),
)

check(
    "SO extent downsampling conservatively preserves a stage-boundary crossing",
    conservative_resampled_interval(
        29_999, 30_001, 1000.0, 20.0
    ) == (599, 601)
    and conservative_resampled_interval(
        30_000, 30_050, 1000.0, 20.0
    ) == (600, 601),
)


class _PortalBlockDataset:
    def __init__(self, sample_rate_hz=10.0, short_start_s=None):
        self.sample_rate_hz = float(sample_rate_hz)
        self.short_start_s = short_start_s

    def get_data(self, start_us, duration_us, channels):
        start_s = start_us / 1e6
        duration_s = duration_us / 1e6
        count = int(round(duration_s * self.sample_rate_hz))
        if (
            self.short_start_s is not None
            and np.isclose(start_s, self.short_start_s)
        ):
            count -= 1
        ordered = sorted(channels)
        return np.column_stack([
            np.full(count, channel, float) for channel in ordered
        ])


exact_records = []
exact_block = pull_continuous_exact(
    _PortalBlockDataset(),
    [2, 0],
    100.0,
    5.0,
    10.0,
    records=exact_records,
    max_request_s=2.0,
)
check(
    "every bounded portal subrequest is counted and requested order is restored",
    exact_block.shape == (50, 2)
    and np.all(exact_block[:, 0] == 2)
    and np.all(exact_block[:, 1] == 0)
    and [record["requested_sample_count"] for record in exact_records]
    == [20, 20, 10]
    and all(
        record["requested_sample_count"] == record["returned_sample_count"]
        and record["status"] == "ok"
        for record in exact_records
    ),
)

short_records = []
try:
    pull_continuous_exact(
        _PortalBlockDataset(short_start_s=102.0),
        [0],
        100.0,
        5.0,
        10.0,
        records=short_records,
        max_request_s=2.0,
    )
    short_rejected = False
except PortalSampleCountMismatch:
    short_rejected = True
check(
    "a one-sample-short portal subrequest fails closed and records both counts",
    short_rejected
    and short_records[-1]["requested_sample_count"] == 20
    and short_records[-1]["returned_sample_count"] == 19
    and short_records[-1]["status"] == "sample_count_mismatch",
)

try:
    pull_continuous_exact(
        _PortalBlockDataset(),
        [0],
        0.0,
        np.inf,
        10.0,
    )
    nonfinite_rejected = False
except ValueError:
    nonfinite_rejected = True
check(
    "nonfinite portal geometry is rejected before entering the request loop",
    nonfinite_rejected,
)
try:
    require_integer_sample_rate(512.5, source="synthetic")
    fractional_rate_rejected = False
except RuntimeError:
    fractional_rate_rejected = True
check(
    "fractional source rates fail before integer one-second binning",
    fractional_rate_rejected
    and require_integer_sample_rate(512.0, source="synthetic") == 512,
)


class _NightProbeDataset:
    def get_data(self, start_us, duration_us, channels):
        del channels
        start_s = start_us / 1e6
        count = int(round(duration_us / 1e6 * 100.0))
        if start_s == 0:
            count -= 1
        time = np.arange(count) / 100.0
        waveform = (
            3.0 * np.sin(2 * np.pi * time)
            + np.sin(2 * np.pi * 10.0 * time)
        )
        return waveform[:, None]


night, night_qc = find_night(
    _NightProbeDataset(),
    {"CTX": 0},
    "CTX",
    100.0,
    total_h=8.0,
    required_h=3.0,
    probe_workers=1,
    return_diagnostics=True,
)
check(
    "six-second night probes reject and record a short response without splicing it",
    night is not None
    and night_qc["n_probe_failures"] == 1
    and night_qc["probe_sample_counts"][0]["requested_sample_count"] == 600
    and night_qc["probe_sample_counts"][0]["returned_sample_count"] == 599
    and night_qc["probe_sample_counts"][0]["status"] == "failed",
)


def _geometry(sample_count=100):
    return {
        "start_time_us": 0,
        "end_time_us": 1_000_000,
        "duration_us": 1_000_000.0,
        "number_of_samples": sample_count,
        "sample_rate_hz": 100.0,
    }


shared_identity = {
    "channels": {
        "A1": _geometry(),
        "B1": _geometry(),
        "ECG1": _geometry(),
    },
}
geometry = validate_hup_series_geometry(
    shared_identity, ["A1", "B1"], "ECG1")
check(
    "HUP geometry is derived from the first selected cortical contact",
    geometry["reference_channel"] == "A1"
    and geometry["sample_rate_hz"] == 100.0,
)
mismatched_identity = {
    "channels": {
        **shared_identity["channels"],
        "ECG1": _geometry(sample_count=99),
    },
}
try:
    validate_hup_series_geometry(
        mismatched_identity, ["A1", "B1"], "ECG1")
    mismatch_rejected = False
except RuntimeError:
    mismatch_rejected = True
check("a selected ECG/cortical geometry mismatch fails closed", mismatch_rejected)


half_wave = np.asarray(
    [1.0] * 3
    + [-1.0, -1.0, -1.0, -1.0, -3.0]
    + [1.0, 1.0, 1.0, 1.0, 5.0]
    + [-1.0] * 5
    + [1.0] * 5
)
half_candidates = detect_so_candidates(half_wave, sf=10.0, edge_s=0)
check(
    "3B half-wave extrema include the final sample before each zero crossing",
    len(half_candidates) == 1
    and half_candidates[0] == (7, 3.0, 5.0, 8.0),
)
zero_upstate = np.asarray(
    [1.0] * 3 + [-1.0] * 5 + [0.0] * 5 + [-1.0] * 5)
check(
    "malformed SO candidates without a positive up-state are rejected",
    detect_so_candidates(zero_upstate, sf=10.0, edge_s=0) == [],
)

event_time = np.arange(600) / 100.0
event_candidates = event_3d_so_candidates(
    np.cos(2 * np.pi * event_time), 100.0)
check(
    "neutral 3D candidates retain trough, amplitude, and half-open cycle extent",
    event_candidates.shape[1] == 4
    and np.all(event_candidates[:, 2] <= event_candidates[:, 0])
    and np.all(event_candidates[:, 0] < event_candidates[:, 3]),
)
candidate = np.asarray([[50.0, 2.0, 20.0, 80.0]])
valid_cycle = np.ones(100, bool)
invalid_cycle = valid_cycle.copy()
invalid_cycle[25] = False
check(
    "3D candidate validity covers the complete SO cycle rather than its trough",
    len(retain_complete_clean_event_3d_candidates(
        candidate, valid_cycle, 0, 100)) == 1
    and len(retain_complete_clean_event_3d_candidates(
        candidate, invalid_cycle, 0, 100)) == 0,
)


FIXTURE_NIGHT_S = 10.0
FIXTURE_DURATION_S = 1.0
FIXTURE_HOURS = FIXTURE_DURATION_S / 3600.0
FIXTURE_SF = 100


def _acquisition_records(subrequest_purpose, core_purpose, *, channels):
    return [
        {
            "purpose": subrequest_purpose,
            "request_start_s": FIXTURE_NIGHT_S,
            "request_duration_s": FIXTURE_DURATION_S,
            "requested_sample_count": FIXTURE_SF,
            "returned_sample_count": FIXTURE_SF,
            "requested_channel_count": int(channels),
            "returned_channel_count": int(channels),
            "status": "ok",
        },
        {
            "purpose": core_purpose,
            "analysis_start_s": 0.0,
            "analysis_duration_s": FIXTURE_DURATION_S,
            "requested_sample_count": FIXTURE_SF,
            "returned_sample_count": FIXTURE_SF,
            "status": "ok",
        },
    ]


_MISSING_LEDGER = object()


def _write_grid_cache(
        directory, subject, *, pipeline, ledger=_MISSING_LEDGER):
    path = os.path.join(directory, f"{subject}.npz")
    payload = {
        "subject": subject,
        "status": "ok",
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "cache_code_sha256": cache_code_sha256(ROOT),
        "failed_chunks_json": "[]",
        "ecg_failures_json": "[]",
        "night_s": FIXTURE_NIGHT_S,
        "hours": FIXTURE_HOURS,
        "sf": float(FIXTURE_SF),
    }
    if ledger is not _MISSING_LEDGER:
        field = (
            "acquisition_sample_counts_json"
            if pipeline == "cache_lc_series"
            else "acquisition_chunk_sample_counts_json"
        )
        payload[field] = ledger
    atomic_savez(path, **payload)
    write_run_manifest(
        directory,
        pipeline=pipeline,
        requested=[subject],
        completed=[subject],
        skipped=[],
        failed=[],
        config={
            "cache_code_sha256": cache_code_sha256(ROOT),
            "cache_schema_version": CACHE_SCHEMA_VERSION,
        },
        result_files_sha256={subject: file_sha256(path)},
    )
    return path


with tempfile.TemporaryDirectory() as valid_hup_grid_dir:
    subject = "HUPGRID_phaseII"
    _write_grid_cache(
        valid_hup_grid_dir,
        subject,
        pipeline="cache_lc_series",
        ledger=json.dumps(_acquisition_records(
            "analysis_subrequest", "analysis_core", channels=2)),
    )
    _cache_manifest(valid_hup_grid_dir)
    check(
        "publication-grid HUP fixture accepts an exact complete acquisition ledger",
        True,
    )

for ledger, expected_message, label in (
    (
        _MISSING_LEDGER,
        "acquisition_sample_counts_json",
        "missing",
    ),
    (
        "not-json",
        "malformed JSON",
        "malformed",
    ),
):
    with tempfile.TemporaryDirectory() as hup_grid_dir:
        subject = "HUPGRID_phaseII"
        _write_grid_cache(
            hup_grid_dir,
            subject,
            pipeline="cache_lc_series",
            ledger=ledger,
        )
        try:
            _cache_manifest(hup_grid_dir)
        except RuntimeError as exc:
            ledger_rejected = expected_message in str(exc)
        else:
            ledger_rejected = False
        check(
            "publication grid rejects a manifest-hash-matched current HUP "
            f"cache with a {label} acquisition ledger",
            ledger_rejected,
        )

with tempfile.TemporaryDirectory() as respect_grid_root:
    respect_grid_dir = os.path.join(respect_grid_root, "ds003848")
    os.makedirs(respect_grid_dir)
    subject = "sub-grid"
    _write_grid_cache(
        respect_grid_dir,
        subject,
        pipeline="stage_ds003848",
        ledger=_MISSING_LEDGER,
    )
    try:
        _cache_manifest(respect_grid_dir)
    except RuntimeError as exc:
        local_ledger_rejected = (
            "acquisition_chunk_sample_counts_json" in str(exc))
    else:
        local_ledger_rejected = False
    check(
        "publication grid enforces the local acquisition ledger for RESPect caches",
        local_ledger_rejected,
    )


def _write_pinned_ieeg(directory, subject):
    source_identity = {
        "snapshot_id": "snapshot-1",
        "cortical_channels": ["A1"],
        "channels": {"A1": _geometry()},
    }
    path = os.path.join(directory, f"{subject}.npz")
    atomic_savez(
        path,
        subject=subject,
        status="ok",
        cache_schema_version=IEEG_CACHE_SCHEMA,
        cache_code_sha256=cache_code_sha256(ROOT),
        failed_chunks_json="[]",
        ecg_failures_json="[]",
        acquisition_sample_counts_json=json.dumps(
            _acquisition_records(
                "analysis_subrequest", "analysis_core", channels=2)),
        night_s=FIXTURE_NIGHT_S,
        hours=FIXTURE_HOURS,
        sf=float(FIXTURE_SF),
        source_identity_json=json.dumps(source_identity),
    )
    write_run_manifest(
        directory,
        pipeline="cache_lc_series",
        requested=[subject],
        completed=[subject],
        skipped=[],
        failed=[],
        config={
            "cache_code_sha256": cache_code_sha256(ROOT),
            "cache_schema_version": IEEG_CACHE_SCHEMA,
        },
        result_files_sha256={subject: file_sha256(path)},
    )
    return source_identity


def _write_sidecar(directory, cache_dir, subject, source_identity):
    dependency = cache_dependency_sha256()
    base_manifest = os.path.join(cache_dir, "RUN_MANIFEST.json")
    base_cache = os.path.join(cache_dir, f"{subject}.npz")
    with open(base_manifest) as handle:
        base_run_id = json.load(handle)["run_id"]
    roles = SCALP_CHANNEL_PLAN[subject]
    source = {
        "dataset_name": subject,
        "snapshot_id": source_identity["snapshot_id"],
        "scalp_role_to_channel": roles,
    }
    path = os.path.join(directory, f"{subject}.npz")
    atomic_savez(
        path,
        status="ok",
        cache_schema_version=SCALP_CACHE_SCHEMA,
        cache_dependency_sha256=dependency,
        subject=subject,
        source_dataset=subject,
        source_kind="iEEG.org API",
        ieeg_cache_sha256=file_sha256(base_cache),
        ieeg_cache_manifest_sha256=file_sha256(base_manifest),
        ieeg_cache_manifest_run_id=base_run_id,
        ieeg_cache_schema_version=IEEG_CACHE_SCHEMA,
        failed_chunks_json="[]",
        acquisition_sample_counts_json=json.dumps(
            _acquisition_records(
                "paired_scalp_subrequest",
                "paired_scalp_core",
                channels=len(roles),
            )),
        night_s=FIXTURE_NIGHT_S,
        hours=FIXTURE_HOURS,
        sf=float(FIXTURE_SF),
        channel_roles_json=json.dumps(roles),
        source_identity_json=json.dumps(source),
    )
    atomic_json_dump({
        "run_id": "synthetic-sidecar-run",
        "schema_version": SCALP_CACHE_SCHEMA,
        "pipeline": SCALP_PIPELINE,
        "run_state": "complete",
        "runtime_versions": runtime_versions(),
        "cache_dependency_sha256": dependency,
        "requested": [subject],
        "completed": [subject],
        "reused": [],
        "failed": [],
        "channel_plan": {subject: roles},
        "config": _manifest_config(cache_dir),
        "result_files_sha256": {subject: file_sha256(path)},
    }, os.path.join(directory, "RUN_MANIFEST.json"))
    return path


with tempfile.TemporaryDirectory() as base_dir, tempfile.TemporaryDirectory() as sidecar_dir:
    paired_subject = "HUP160_phaseII"
    pinned_identity = _write_pinned_ieeg(base_dir, paired_subject)
    sidecar_path = _write_sidecar(
        sidecar_dir, base_dir, paired_subject, pinned_identity)
    reuse = validate_reusable_sidecar_run(
        sidecar_dir, [paired_subject], cache_dir=base_dir)
    check(
        "paired sidecar reuse requires and accepts an exact terminal hash",
        reuse["subjects"][paired_subject]["sha256"] == file_sha256(sidecar_path),
    )
    with np.load(sidecar_path, allow_pickle=False) as cache:
        tampered_payload = {key: np.asarray(cache[key]) for key in cache.files}
    tampered_payload["tampered_scientific_value"] = np.asarray([1.0])
    atomic_savez(sidecar_path, **tampered_payload)
    try:
        validate_reusable_sidecar_run(
            sidecar_dir, [paired_subject], cache_dir=base_dir)
        tampered_rejected = False
    except RuntimeError:
        tampered_rejected = True
    check(
        "paired reuse cannot rehash and bless modified sidecar bytes",
        tampered_rejected,
    )

with tempfile.TemporaryDirectory() as base_dir, tempfile.TemporaryDirectory() as sidecar_dir:
    paired_subject = "HUP160_phaseII"
    _write_pinned_ieeg(base_dir, paired_subject)
    frozen = freeze_pinned_ieeg_inputs([paired_subject], base_dir)
    config = _manifest_config(
        base_dir,
        pinned_manifest_sha256=frozen["manifest_sha256"],
    )
    atomic_json_dump(
        {"run_state": "complete"},
        os.path.join(sidecar_dir, "RUN_MANIFEST.json"),
    )
    _write_in_progress_manifest(
        sidecar_dir,
        requested=[paired_subject],
        config=config,
        run_id="synthetic-in-progress",
    )
    with open(os.path.join(sidecar_dir, "RUN_MANIFEST.json")) as handle:
        invalidated = json.load(handle)
    check(
        "paired force/rebuild invalidates a prior terminal manifest before writes",
        invalidated["run_state"] == "in_progress"
        and invalidated["completed"] == []
        and invalidated["result_files_sha256"] == {},
    )

    with open(frozen["manifest_path"]) as handle:
        changed_manifest = json.load(handle)
    changed_manifest["generated_at_utc"] = "changed-during-sidecar-run"
    atomic_json_dump(changed_manifest, frozen["manifest_path"])
    try:
        assert_pinned_ieeg_inputs_unchanged(frozen)
        changed_base_rejected = False
    except RuntimeError:
        changed_base_rejected = True
    check(
        "a base HUP manifest change during sidecar work fails closed",
        changed_base_rejected,
    )

with tempfile.TemporaryDirectory() as contradictory_dir:
    paired_subject = "HUP160_phaseII"
    _write_pinned_ieeg(contradictory_dir, paired_subject)
    manifest_path = os.path.join(contradictory_dir, "RUN_MANIFEST.json")
    with open(manifest_path) as handle:
        contradictory = json.load(handle)
    contradictory["failed"] = [{
        "subject": paired_subject,
        "error": "synthetic contradictory failure",
    }]
    atomic_json_dump(contradictory, manifest_path)
    try:
        freeze_pinned_ieeg_inputs([paired_subject], contradictory_dir)
        contradictory_rejected = False
    except RuntimeError:
        contradictory_rejected = True
    check(
        "paired input validation rejects a contradictory complete base manifest",
        contradictory_rejected,
    )


def _write_current_cache(
        path, subject, *, generated_at="2999-01-01T00:00:00+00:00",
        ecg_failures=None, truncate_acquisition=False):
    acquisition_records = _acquisition_records(
        "analysis_subrequest", "analysis_core", channels=2)
    if truncate_acquisition:
        acquisition_records[0]["request_duration_s"] = 0.5
        acquisition_records[0]["requested_sample_count"] = 50
        acquisition_records[0]["returned_sample_count"] = 50
        acquisition_records[1]["analysis_duration_s"] = 0.5
        acquisition_records[1]["requested_sample_count"] = 50
        acquisition_records[1]["returned_sample_count"] = 50
    atomic_savez(
        path,
        subject=subject,
        status="ok",
        cache_schema_version=CACHE_SCHEMA_VERSION,
        cache_code_sha256=cache_code_sha256(ROOT),
        generated_at_utc=generated_at,
        failed_chunks_json="[]",
        ecg_failures_json=json.dumps(ecg_failures or []),
        acquisition_sample_counts_json=json.dumps(acquisition_records),
        night_s=FIXTURE_NIGHT_S,
        hours=FIXTURE_HOURS,
        sf=float(FIXTURE_SF),
    )


def _write_failed_batch(
        directory, *, unchanged_retry=False, fatal_retry=False,
        truncate_retry_acquisition=False):
    stable = "HUPSTABLE_phaseII"
    retried = "HUPRETRY_phaseII"
    stable_path = os.path.join(directory, f"{stable}.npz")
    retried_path = os.path.join(directory, f"{retried}.npz")
    _write_current_cache(stable_path, stable)
    _write_current_cache(
        retried_path,
        retried,
        ecg_failures=(
            [{"error": "synthetic detector failure"}] if fatal_retry else []
        ),
        truncate_acquisition=truncate_retry_acquisition,
    )
    retried_hash = file_sha256(retried_path)
    baseline_hash = retried_hash if unchanged_retry else "a" * 64
    config = {
        "cache_code_sha256": cache_code_sha256(ROOT),
        "cache_schema_version": CACHE_SCHEMA_VERSION,
    }
    manifest = write_run_manifest(
        directory,
        pipeline="cache_lc_series",
        requested=[stable, retried],
        completed=[stable],
        skipped=[],
        failed=[{
            "subject": retried,
            "error": "ConnectionError: transient",
            "pre_run_output": {
                "existed": True,
                "sha256": baseline_hash,
            },
        }],
        config=config,
        run_state="failed",
        result_files_sha256={stable: file_sha256(stable_path)},
    )
    return stable, retried, manifest, file_sha256(
        os.path.join(directory, "RUN_MANIFEST.json"))


with tempfile.TemporaryDirectory() as retry_dir:
    stable, retried, failed_manifest, failed_manifest_hash = (
        _write_failed_batch(retry_dir)
    )
    stable_hash = file_sha256(
        os.path.join(retry_dir, f"{stable}.npz"))
    finalized = finalize(retry_dir, [retried])
    archive = finalized["recovery_provenance"][
        "prior_batch_manifest_archive"]
    check(
        "retry finalization preserves the failed manifest and proves rebuilt bytes",
        finalized["run_state"] == "complete"
        and finalized["failed"] == []
        and file_sha256(os.path.join(retry_dir, archive["relative_path"]))
        == failed_manifest_hash
        and os.path.basename(archive["relative_path"]).startswith(
            failed_manifest_hash)
        and finalized["recovery_provenance"][
            "retry_proof_by_subject"][retried][
                "validated_acquisition_count_records"] == 2
        and finalized["recovery_provenance"][
            "retry_proof_by_subject"][retried][
                "validated_analysis_core_count"] == 1,
    )
    check(
        "retry finalization keeps every non-retried result byte-identical",
        file_sha256(os.path.join(retry_dir, f"{stable}.npz")) == stable_hash
        and finalized["result_files_sha256"][stable] == stable_hash,
    )

with tempfile.TemporaryDirectory() as unchanged_dir:
    _stable, unchanged_subject, _manifest, _manifest_hash = (
        _write_failed_batch(unchanged_dir, unchanged_retry=True)
    )
    try:
        finalize(unchanged_dir, [unchanged_subject])
        unchanged_rejected = False
    except RuntimeError:
        unchanged_rejected = True
    check(
        "retry finalization rejects a failed subject whose bytes were not rebuilt",
        unchanged_rejected,
    )

with tempfile.TemporaryDirectory() as fatal_dir:
    _stable, fatal_subject, _manifest, _manifest_hash = (
        _write_failed_batch(fatal_dir, fatal_retry=True)
    )
    try:
        finalize(fatal_dir, [fatal_subject])
        fatal_rejected = False
    except RuntimeError:
        fatal_rejected = True
    check(
        "retry finalization rejects a rebuilt cache with fatal ECG metadata",
        fatal_rejected,
    )

with tempfile.TemporaryDirectory() as truncated_dir:
    _stable, truncated_subject, _manifest, _manifest_hash = (
        _write_failed_batch(
            truncated_dir, truncate_retry_acquisition=True)
    )
    try:
        finalize(truncated_dir, [truncated_subject])
        truncated_rejected = False
    except RuntimeError:
        truncated_rejected = True
    check(
        "retry finalization rejects exact-looking records that cover only "
        "part of the cached interval",
        truncated_rejected,
    )

check(
    "cache producer tests execute against the current analysis/schema contract",
    ANALYSIS_VERSION.endswith("-v9")
    and CACHE_SCHEMA_VERSION == IEEG_CACHE_SCHEMA,
)
print("ALL CACHE-PRODUCER INTEGRITY CHECKS PASSED")
