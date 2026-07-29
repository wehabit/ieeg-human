"""Focused regressions for HUP acquisition and cache provenance hardening."""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile

import numpy as np

import cache_paired_scalp as paired_scalp_cache
from cache_lc_series import (
    FILTER_EDGE_S,
    SIGMA_FIXED,
    SO_BAND_NAJI,
    SWA_BAND_L,
    conservative_resampled_interval,
    detect_so_candidates,
    event_3d_so_candidates,
    retain_complete_clean_event_3d_candidates,
)
from cache_paired_scalp import (
    HISTORICAL_POWER_SUPPORT,
    IEEG_CACHE_SCHEMA,
    PIPELINE as SCALP_PIPELINE,
    SCALP_CACHE_SCHEMA,
    SCALP_CHANNEL_PLAN,
    TERMINAL_SIDECAR_MANIFEST_FIELDS,
    _DEPENDENCY_FILES as SCALP_DEPENDENCY_FILES,
    _manifest_config,
    _write_in_progress_manifest,
    assert_pinned_ieeg_inputs_unchanged,
    cache_dependency_sha256,
    freeze_pinned_ieeg_inputs,
    validate_reusable_sidecar_run,
    validate_terminal_sidecar_manifest,
)
from hup_portal import (
    PortalSampleCountMismatch,
    expected_portal_sample_count,
    expected_portal_sample_offset,
    find_night,
    portal_core_sample_geometry,
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
    git_revision,
    require_integer_sample_rate,
    runtime_versions,
    utc_now,
    write_run_manifest,
)
from run_qc_grid import _cache_manifest
from staging_helpers import CHUNK_S


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
    "paired scalp plan covers the v9 spectrum-eligible HUP138 candidate",
    SCALP_CHANNEL_PLAN.get("HUP138_phaseII")
    == {"c3": "C3", "f3": "F3", "f4": "F4", "fz": "Fz"},
)

check(
    "a first-chunk zero offset maps to sample zero without becoming a zero-duration request",
    expected_portal_sample_offset(0.0, 1024.0) == 0,
)
check(
    "both portal cache producers share first-chunk core geometry",
    portal_core_sample_geometry(0.0, 0.0, 600.0, 1024.0)
    == (0, 614_400),
)
try:
    expected_portal_sample_count(0.0, 1024.0)
except ValueError:
    zero_request_rejected = True
else:
    zero_request_rejected = False
check(
    "a zero-duration portal request remains invalid",
    zero_request_rejected,
)
try:
    expected_portal_sample_count(1e-12, 1.0)
except ValueError:
    sub_sample_request_rejected = True
else:
    sub_sample_request_rejected = False
check(
    "a positive duration that rounds to zero samples remains invalid",
    sub_sample_request_rejected,
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

SCALP_VALIDATION_SOURCE = "analysis/paired_scalp_sidecar_validation.py"
captured_scalp_dependencies = []
real_scalp_hash_files = paired_scalp_cache._hash_files


def capture_scalp_dependencies(relative_paths):
    captured_scalp_dependencies.extend(relative_paths)
    return "a" * 64


paired_scalp_cache._hash_files = capture_scalp_dependencies
try:
    captured_scalp_digest = cache_dependency_sha256()
finally:
    paired_scalp_cache._hash_files = real_scalp_hash_files

check(
    "validator-only source is outside the paired-scalp producer digest",
    captured_scalp_digest == "a" * 64
    and tuple(captured_scalp_dependencies) == SCALP_DEPENDENCY_FILES
    and SCALP_VALIDATION_SOURCE not in captured_scalp_dependencies,
)
sidecar_contract = paired_scalp_cache._validation_contract()
check(
    "byte-affecting sidecar contract values remain producer hash-bound",
    "analysis/cache_paired_scalp.py" in captured_scalp_dependencies
    and _manifest_config.__module__ == "cache_paired_scalp"
    and sidecar_contract.channel_plan is SCALP_CHANNEL_PLAN
    and sidecar_contract.scalp_cache_schema == SCALP_CACHE_SCHEMA
    and sidecar_contract.ieeg_cache_schema == IEEG_CACHE_SCHEMA
    and sidecar_contract.filter_edge_s == FILTER_EDGE_S
    and sidecar_contract.chunk_s == CHUNK_S
    and sidecar_contract.historical_power_support
    == HISTORICAL_POWER_SUPPORT
    and sidecar_contract.sigma_band_hz == tuple(SIGMA_FIXED)
    and sidecar_contract.swa_band_hz == tuple(SWA_BAND_L)
    and sidecar_contract.so_band_hz == tuple(SO_BAND_NAJI),
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
        "revision_id": "synthetic-revision",
        "data_check": "synthetic-data-check",
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
        "reference_ieeg_channel": "A1",
        "reference_ieeg_identity": source_identity["channels"]["A1"],
        "scalp_role_to_channel": roles,
        "scalp_channels": {
            channel: _geometry(FIXTURE_SF)
            for channel in roles.values()
        },
    }
    n_channels = len(roles)
    shape = (n_channels, int(FIXTURE_DURATION_S))
    denominator = np.full(shape, FIXTURE_SF, dtype=np.int64)
    numerator = np.full(shape, 2.0, dtype=float)
    power = numerator / denominator
    minimum = np.zeros(n_channels, dtype=float)
    maximum = np.ones(n_channels, dtype=float)
    tolerance = np.full(
        n_channels, 64 * np.finfo(float).eps, dtype=float)
    candidates = {
        f"so_candidate_{field}_{channel}": np.asarray([], dtype=float)
        for channel in roles.values()
        for field in ("t", "down", "up", "p2p")
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
        scalp_chans=np.asarray(list(roles.values()), dtype="<U16"),
        filter_edge_seconds=float(FILTER_EDGE_S),
        chunk_seconds=float(CHUNK_S),
        sigma_band_hz=np.asarray(SIGMA_FIXED, float),
        swa_band_hz=np.asarray(SWA_BAND_L, float),
        so_band_hz=np.asarray(SO_BAND_NAJI, float),
        sigma_fixed_by_channel=power,
        swa_by_channel=power,
        sigma_fixed_power_numerator_by_channel=numerator,
        sigma_fixed_clean_sample_count_by_channel=denominator,
        swa_power_numerator_by_channel=numerator,
        swa_clean_sample_count_by_channel=denominator,
        power_samples_per_second=FIXTURE_SF,
        power_historical_minimum_clean_fraction_per_second=(
            HISTORICAL_POWER_SUPPORT),
        scalp_signal_nonflat_mask=np.ones(n_channels, dtype=bool),
        scalp_signal_raw_minimum=minimum,
        scalp_signal_raw_maximum=maximum,
        scalp_signal_raw_dynamic_range=maximum - minimum,
        scalp_signal_numerical_flat_tolerance=tolerance,
        scalp_signal_finite_sample_count=np.full(
            n_channels, FIXTURE_SF, dtype=np.int64),
        ecg_reused_not_redetected=True,
        staging_reused_not_recomputed=True,
        **candidates,
    )
    atomic_json_dump({
        "run_id": "11111111-1111-4111-8111-111111111111",
        "schema_version": SCALP_CACHE_SCHEMA,
        "pipeline": SCALP_PIPELINE,
        "run_state": "complete",
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": IEEG_CACHE_SCHEMA,
        "generated_at_utc": utc_now(),
        "code_revision": git_revision(ROOT),
        "code_dirty": False,
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


def _rewrite_sidecar_and_rehash(sidecar_path, mutate):
    with np.load(sidecar_path, allow_pickle=False) as cache:
        payload = {key: np.asarray(cache[key]) for key in cache.files}
    mutate(payload)
    atomic_savez(sidecar_path, **payload)
    manifest_path = os.path.join(
        os.path.dirname(sidecar_path), "RUN_MANIFEST.json")
    with open(manifest_path) as handle:
        manifest = json.load(handle)
    subject = str(np.asarray(payload["subject"]).item())
    manifest["result_files_sha256"][subject] = file_sha256(sidecar_path)
    atomic_json_dump(manifest, manifest_path)


# The implementation under test is necessarily uncommitted in a developer run.
# Simulate its eventual source commit only for the current-revision fixture;
# every other revision below is hashed from its real stored Git bytes.
_real_dependency_sha256_at_revision = (
    paired_scalp_cache._dependency_sha256_at_revision)
_fixture_revision = git_revision(ROOT)
_fixture_dependency_sha256 = cache_dependency_sha256()


def _fixture_dependency_sha256_at_revision(revision):
    if revision == _fixture_revision:
        return _fixture_dependency_sha256
    return _real_dependency_sha256_at_revision(revision)


paired_scalp_cache._dependency_sha256_at_revision = (
    _fixture_dependency_sha256_at_revision)


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
    sidecar_manifest_path = os.path.join(
        sidecar_dir, "RUN_MANIFEST.json")
    with open(sidecar_manifest_path) as handle:
        valid_sidecar_manifest = json.load(handle)
    check(
        "paired sidecar producer fixture uses the exact terminal schema",
        set(valid_sidecar_manifest) == TERMINAL_SIDECAR_MANIFEST_FIELDS,
    )

    def terminal_manifest_rejected(mutate):
        candidate = json.loads(json.dumps(valid_sidecar_manifest))
        mutate(candidate)
        try:
            validate_terminal_sidecar_manifest(
                candidate,
                source="synthetic paired-scalp manifest",
                output_dir=sidecar_dir,
                cache_dir=base_dir,
                expected_requested=[paired_subject],
            )
        except RuntimeError:
            return True
        return False

    check(
        "paired sidecar terminal validation rejects every missing field",
        all(
            terminal_manifest_rejected(
                lambda candidate, field=field: candidate.pop(field))
            for field in TERMINAL_SIDECAR_MANIFEST_FIELDS
        ),
    )
    check(
        "paired sidecar terminal validation rejects extra fields",
        terminal_manifest_rejected(
            lambda candidate: candidate.update({"unexpected": True})),
    )
    check(
        "paired sidecar terminal validation rejects dirty production",
        terminal_manifest_rejected(
            lambda candidate: candidate.update({"code_dirty": True})),
    )
    stale_mutations = (
        lambda candidate: candidate.update({"analysis_version": "stale"}),
        lambda candidate: candidate.update({"cache_schema_version": "stale"}),
        lambda candidate: candidate.update({"schema_version": "stale"}),
        lambda candidate: candidate.update({
            "cache_dependency_sha256": "0" * 64}),
        lambda candidate: candidate.update({"runtime_versions": {}}),
        lambda candidate: candidate["config"].update({"chunk_seconds": -1}),
    )
    check(
        "paired sidecar terminal validation rejects stale provenance",
        all(
            terminal_manifest_rejected(mutate)
            for mutate in stale_mutations
        ),
    )
    malformed_identity_mutations = (
        lambda candidate: candidate.update({"run_id": "not-a-uuid"}),
        lambda candidate: candidate.update({
            "generated_at_utc": "not-a-timestamp"}),
        lambda candidate: candidate.update({"code_revision": "f" * 40}),
    )
    check(
        "paired sidecar terminal validation rejects fabricated run identity",
        all(
            terminal_manifest_rejected(mutate)
            for mutate in malformed_identity_mutations
        ),
    )
    revision_list = subprocess.run(
        ["git", "-C", ROOT, "rev-list", "--all"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
        text=True,
    ).stdout.splitlines()
    other_revision = None
    for revision in revision_list:
        if revision == _fixture_revision:
            continue
        recorded_digest = _real_dependency_sha256_at_revision(revision)
        if (
            recorded_digest is not None
            and recorded_digest != _fixture_dependency_sha256
        ):
            other_revision = revision
            break
    check(
        "paired sidecar terminal validation rejects a real commit whose "
        "dependency bytes do not match",
        other_revision is not None
        and terminal_manifest_rejected(
            lambda candidate: candidate.update({
                "code_revision": other_revision})),
    )

    contradictory_manifest = json.loads(json.dumps(valid_sidecar_manifest))
    contradictory_manifest["reused"] = [paired_subject]
    atomic_json_dump(contradictory_manifest, sidecar_manifest_path)
    try:
        validate_reusable_sidecar_run(
            sidecar_dir, [paired_subject], cache_dir=base_dir)
    except RuntimeError:
        contradictory_manifest_rejected = True
    else:
        contradictory_manifest_rejected = False
    check(
        "shared sidecar manifest validator rejects overlapping partitions",
        contradictory_manifest_rejected,
    )
    sidecar_path = _write_sidecar(
        sidecar_dir, base_dir, paired_subject, pinned_identity)

    def corrupt_ledger(payload):
        ledger = json.loads(str(np.asarray(
            payload["acquisition_sample_counts_json"]).item()))
        ledger[0]["returned_sample_count"] -= 1
        payload["acquisition_sample_counts_json"] = np.asarray(
            json.dumps(ledger))

    def corrupt_roles(payload):
        roles = json.loads(str(np.asarray(
            payload["channel_roles_json"]).item()))
        roles["f3"] = "F7"
        payload["channel_roles_json"] = np.asarray(json.dumps(roles))
        source = json.loads(str(np.asarray(
            payload["source_identity_json"]).item()))
        source["scalp_role_to_channel"] = roles
        payload["source_identity_json"] = np.asarray(json.dumps(source))

    def remove_scientific_payload(payload):
        payload.pop("sigma_fixed_by_channel")

    def split_band_support(payload):
        payload["swa_clean_sample_count_by_channel"] = np.full(
            np.asarray(
                payload["swa_clean_sample_count_by_channel"]).shape,
            FIXTURE_SF - 1,
            dtype=np.int64,
        )
        payload["swa_power_numerator_by_channel"] = np.full(
            np.asarray(
                payload["swa_power_numerator_by_channel"]).shape,
            1.98,
        )

    def exceed_raw_finite_support(payload):
        n_channels = len(SCALP_CHANNEL_PLAN[paired_subject])
        payload["scalp_signal_finite_sample_count"] = np.ones(
            n_channels, dtype=np.int64)
        payload["scalp_signal_nonflat_mask"] = np.zeros(
            n_channels, dtype=bool)

    def duplicate_so_time(payload):
        for field, values in (
            ("t", [0.25, 0.25]),
            ("down", [1.0, 1.0]),
            ("up", [2.0, 2.0]),
            ("p2p", [3.0, 3.0]),
        ):
            payload[f"so_candidate_{field}_C3"] = np.asarray(values)

    def remove_channel_revision(payload):
        source = json.loads(str(np.asarray(
            payload["source_identity_json"]).item()))
        source["scalp_channels"]["C3"].pop("revision_id")
        payload["source_identity_json"] = np.asarray(json.dumps(source))

    for label, mutate in (
        ("acquisition ledger", corrupt_ledger),
        ("role/source identity", corrupt_roles),
        ("required scientific payload", remove_scientific_payload),
        ("split sigma/SWA support", split_band_support),
        ("power support exceeding raw finite support",
         exceed_raw_finite_support),
        ("duplicate slow-oscillation times", duplicate_so_time),
        ("channel identity without revision", remove_channel_revision),
    ):
        sidecar_path = _write_sidecar(
            sidecar_dir, base_dir, paired_subject, pinned_identity)
        _rewrite_sidecar_and_rehash(sidecar_path, mutate)
        try:
            validate_reusable_sidecar_run(
                sidecar_dir, [paired_subject], cache_dir=base_dir)
        except RuntimeError:
            exploit_rejected = True
        else:
            exploit_rejected = False
        check(
            f"paired reuse rejects a rehashed malformed {label}",
            exploit_rejected,
        )

    sidecar_path = _write_sidecar(
        sidecar_dir, base_dir, paired_subject, pinned_identity)

    def within_geometry_tolerance(payload):
        source = json.loads(str(np.asarray(
            payload["source_identity_json"]).item()))
        source["reference_ieeg_identity"]["duration_us"] += 0.5
        for identity in source["scalp_channels"].values():
            identity["duration_us"] += 0.5
        payload["source_identity_json"] = np.asarray(json.dumps(source))

    _rewrite_sidecar_and_rehash(sidecar_path, within_geometry_tolerance)
    tolerated = validate_reusable_sidecar_run(
        sidecar_dir, [paired_subject], cache_dir=base_dir)
    check(
        "sidecar validator preserves the producer's geometry tolerances",
        tolerated["subjects"][paired_subject]["sha256"]
        == file_sha256(sidecar_path),
    )

    sidecar_path = _write_sidecar(
        sidecar_dir, base_dir, paired_subject, pinned_identity)
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
    four_role_subject = "HUP138_phaseII"
    pinned_identity = _write_pinned_ieeg(base_dir, four_role_subject)
    sidecar_path = _write_sidecar(
        sidecar_dir, base_dir, four_role_subject, pinned_identity)
    four_role_reuse = validate_reusable_sidecar_run(
        sidecar_dir, [four_role_subject], cache_dir=base_dir)
    check(
        "four-role HUP138 C3/F3/F4/Fz sidecar validates end to end",
        four_role_reuse["subjects"][four_role_subject]["sha256"]
        == file_sha256(sidecar_path),
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
        ecg_failures=None, truncate_acquisition=False,
        hours=FIXTURE_HOURS):
    acquisition_records = _acquisition_records(
        "analysis_subrequest", "analysis_core", channels=2)
    duration_s = float(hours) * 3600.0
    sample_count = int(round(duration_s * FIXTURE_SF))
    acquisition_records[0]["request_duration_s"] = duration_s
    acquisition_records[0]["requested_sample_count"] = sample_count
    acquisition_records[0]["returned_sample_count"] = sample_count
    acquisition_records[1]["analysis_duration_s"] = duration_s
    acquisition_records[1]["requested_sample_count"] = sample_count
    acquisition_records[1]["returned_sample_count"] = sample_count
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
        hours=hours,
        sf=float(FIXTURE_SF),
    )


def _write_failed_batch(
        directory, *, unchanged_retry=False, fatal_retry=False,
        truncate_retry_acquisition=False,
        truncate_stable_acquisition=False,
        manifest_hours=FIXTURE_HOURS,
        retried_hours=FIXTURE_HOURS):
    stable = "HUPSTABLE_phaseII"
    retried = "HUPRETRY_phaseII"
    stable_path = os.path.join(directory, f"{stable}.npz")
    retried_path = os.path.join(directory, f"{retried}.npz")
    _write_current_cache(
        stable_path,
        stable,
        truncate_acquisition=truncate_stable_acquisition,
    )
    _write_current_cache(
        retried_path,
        retried,
        hours=retried_hours,
        ecg_failures=(
            [{"error": "synthetic detector failure"}] if fatal_retry else []
        ),
        truncate_acquisition=truncate_retry_acquisition,
    )
    retried_hash = file_sha256(retried_path)
    baseline_hash = retried_hash if unchanged_retry else "a" * 64
    config = {
        "hours": manifest_hours,
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

with tempfile.TemporaryDirectory() as truncated_stable_dir:
    _stable, stable_retry_subject, _manifest, _manifest_hash = (
        _write_failed_batch(
            truncated_stable_dir,
            truncate_stable_acquisition=True,
        )
    )
    try:
        finalize(truncated_stable_dir, [stable_retry_subject])
        truncated_stable_rejected = False
    except RuntimeError:
        truncated_stable_rejected = True
    check(
        "retry finalization validates full acquisition for non-retried "
        "completed caches",
        truncated_stable_rejected,
    )

with tempfile.TemporaryDirectory() as wrong_retry_duration_dir:
    _stable, wrong_duration_subject, _manifest, _manifest_hash = (
        _write_failed_batch(
            wrong_retry_duration_dir,
            # The retry NPZ and its ledger consistently describe two seconds,
            # but the failed batch requested one second. This
            # internally consistent mismatch passed before the finalizer
            # explicitly bound caches back to manifest config.hours.
            retried_hours=2.0 / 3600.0,
        )
    )
    try:
        finalize(wrong_retry_duration_dir, [wrong_duration_subject])
        wrong_retry_duration_rejected = False
    except RuntimeError:
        wrong_retry_duration_rejected = True
    check(
        "retry finalization binds embedded and acquired duration exactly to "
        "prior manifest config.hours",
        wrong_retry_duration_rejected,
    )

check(
    "cache producer tests execute against the current analysis/schema contract",
    ANALYSIS_VERSION.endswith("-v9")
    and CACHE_SCHEMA_VERSION == IEEG_CACHE_SCHEMA,
)
print("ALL CACHE-PRODUCER INTEGRITY CHECKS PASSED")
