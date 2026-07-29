"""Synthetic regressions for the bounded HUP138 bilateral scalp-only runner."""
from __future__ import annotations

import copy
import json
import os
import tempfile

import numpy as np

from artifact_contracts import SCALP_F3_F4_EXPLORATORY_SCHEMA
from event_3b_estimators import rr_baseline_hr, subject_so_triggered
from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    file_sha256,
    runtime_versions,
    source_tree_sha256,
)
from qc_profiles import (
    load_qc_profile,
    profile_file_sha256,
    qc_profile_sha256,
)
from scalp_f3_f4_artifact_validation import (
    _validate_figure_bytes,
    _validate_result_contract,
    validate_exploratory_artifacts,
)
from scalp_f3_f4_exploratory import (
    PIPELINE,
    PROFILE_ID,
    RESULT_PNG,
    RESULT_JSON,
    STAGES,
    SUBJECT,
    _analyse_bilateral,
    _canonical_json,
    _stage_pool,
    _stage_seed,
    _write_artifacts,
)


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


def _scalp_npz(path, f3_times, f4_times):
    payload = {
        "channel_roles_json": json.dumps(
            {"f3": "F3", "f4": "F4"}, sort_keys=True),
        "scalp_chans": np.asarray(["F3", "F4"]),
        "scalp_signal_nonflat_mask": np.asarray([True, True]),
    }
    for channel, times in (("F3", f3_times), ("F4", f4_times)):
        times = np.asarray(times, float)
        payload[f"so_candidate_t_{channel}"] = times
        payload[f"so_candidate_up_{channel}"] = np.ones(len(times))
        payload[f"so_candidate_p2p_{channel}"] = np.full(len(times), 2.0)
    np.savez_compressed(path, **payload)


profile = load_qc_profile(PROFILE_ID)
labels = np.asarray(["NREM"] * 8)
rr = np.ones(8 * 30 * 4, float)
times = np.arange(10.0, 226.0, 5.0)
for time_s in times:
    centre = int(round(time_s * 4))
    rr[centre:centre + 8] = 0.92
base = {
    "subject": SUBJECT,
    "rr_4": rr,
    "stage_lab": labels,
    "hr_meets_profile": True,
}

with tempfile.TemporaryDirectory() as directory:
    cache_path = os.path.join(directory, "scalp.npz")
    _scalp_npz(cache_path, times, times + 0.25)
    with np.load(cache_path, allow_pickle=False) as scalp:
        first = _analyse_bilateral(scalp, base, profile)
        second = _analyse_bilateral(scalp, base, profile)

check(
    "bilateral estimator is deterministic",
    _canonical_json(first) == _canonical_json(second),
)
check(
    "N2/N3 fail closed without stable proxy-stage support",
    all(
        first["stages"][stage]["estimate"] is None
        and "no uninterrupted 180-s stage run"
        in first["stages"][stage]["support_reasons"]
        for stage in ("N2", "N3")
    ),
)
nrem = first["stages"]["NREM"]
check(
    "pooled NREM is explicit, non-Naji, and two-channel only",
    nrem["pooled_nrem_exploratory"]
    and "no Naji stage counterpart" in nrem["stage_interpretation"]
    and nrem["estimate"]["n_channels"] == 2
    and nrem["estimate"]["z"] is None
    and nrem["estimate"]["p_upper"] is None,
)
pool = _stage_pool(rr, np.arange(len(labels)))
direct = subject_so_triggered(
    rr,
    [times, times + 0.25],
    rr_baseline_hr(rr[pool]),
    pool,
    n_sur=199,
    rng=np.random.RandomState(_stage_seed("NREM")),
    domain="rr",
    minimum_channels=2,
    channel_ids=["F3", "F4"],
    minimum_events_per_channel=30,
    minimum_surrogate_pool_samples=100,
)
check(
    "runner uses the exact existing joint RR estimator",
    nrem["estimate"]["peak_lag_s"] == direct["peak_lag_s"]
    and nrem["estimate"]["pct_above_stage_mean"]
    == direct["pct_above_stage_mean"]
    and nrem["estimate"]["curve"] == direct["curve"],
)

with tempfile.TemporaryDirectory() as directory:
    cache_path = os.path.join(directory, "scalp.npz")
    _scalp_npz(cache_path, times, times[:29])
    with np.load(cache_path, allow_pickle=False) as scalp:
        insufficient = _analyse_bilateral(scalp, base, profile)
check(
    "one under-supported electrode does not downgrade to unilateral output",
    insufficient["stages"]["NREM"]["estimate"] is None
    and any(
        reason.startswith("F4 29 eligible SOs")
        for reason in insufficient["stages"]["NREM"]["support_reasons"]
    ),
)

def _synthetic_payload(root, result):
    public = os.path.join(root, "public")
    os.makedirs(public)
    pins = {}
    for name, value in (
        ("inventory", {"inventory": True}),
        ("inventory_manifest", {"inventory_manifest": True}),
        ("qc_grid", {"qc_grid": True}),
    ):
        path = os.path.join(public, f"{name}.json")
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(value, handle, sort_keys=True)
        pins[f"{name}_path_relative"] = os.path.relpath(path, root)
        pins[f"{name}_sha256"] = file_sha256(path)
    pins.update({
        "ieeg_cache_path_relative": "private/HUP138_phaseII.npz",
        "ieeg_cache_sha256": "1" * 64,
        "ieeg_cache_manifest_path_relative": "private/RUN_MANIFEST.json",
        "ieeg_cache_manifest_sha256": "2" * 64,
        "ieeg_cache_manifest_run_id": "synthetic-ieeg-run",
        "scalp_sidecar_path_relative": "private/scalp/HUP138_phaseII.npz",
        "scalp_sidecar_sha256": "3" * 64,
        "scalp_sidecar_manifest_path_relative": (
            "private/scalp/RUN_MANIFEST.json"),
        "scalp_sidecar_manifest_sha256": "4" * 64,
        "scalp_sidecar_cache_dependency_sha256": "a" * 64,
        "qc_grid_profile_sha256": "9" * 64,
        "rr_4_sha256": "5" * 64,
        "stage_lab_sha256": "6" * 64,
        "rr_4_n_samples": len(rr),
        "stage_lab_n_epochs": len(labels),
        "channel_roles_json_sha256": "7" * 64,
        "source_identity_json_sha256": "8" * 64,
        "validated_sidecar_roles": {"f3": "F3", "f4": "F4"},
        "portal_source_identity": {"synthetic": True},
        "inventory_channel_identity": {
            "f3": {"label": "F3"},
            "f4": {"label": "F4"},
        },
    })
    return {
        "schema_version": SCALP_F3_F4_EXPLORATORY_SCHEMA,
        "pipeline": PIPELINE,
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "subject": SUBJECT,
        "analysis_scope": "participant_level_scalp_only",
        "exploratory": True,
        "paired_comparison_performed": False,
        "group_inference_performed": False,
        "event_locking_inference_available": False,
        "provenance": {
            "code_revision": "a" * 40,
            "code_dirty": False,
            "source_tree_sha256": source_tree_sha256(root),
            "runtime_versions": runtime_versions(),
        },
        "input_lineage": pins,
        "method": {
            "estimator": "event_3b_estimators.subject_so_triggered",
            "source_profile_id": PROFILE_ID,
            "source_profile_sha256": qc_profile_sha256(profile),
            "source_profile_file_sha256": profile_file_sha256(),
            "tachogram_domain": "rr",
            "channels": ["F3", "F4"],
            "minimum_channels": 2,
            "n_surrogates_for_diagnostic_only": 199,
            "event_locking_z_p": None,
        },
        "locked_ieeg_context": {
            "available_3b_stages": [],
            "used_as_comparison_arm": False,
            "note": "synthetic context is not a comparison arm",
        },
        "result": result,
        "claim_limits": ["synthetic descriptive-only contract"],
    }


def _set_path(value, path, replacement):
    target = value
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = replacement


def _check_semantic_tamper(payload, label, mutate):
    tampered = copy.deepcopy(payload)
    mutate(tampered)
    try:
        _validate_result_contract(tampered)
    except RuntimeError:
        rejected = True
    else:
        rejected = False
    check(f"semantic validator rejects {label} tamper", rejected)


with tempfile.TemporaryDirectory() as root:
    absent_dir = os.path.join(root, "absent")
    absent = validate_exploratory_artifacts(
        root, absent_dir, mode="offline", allow_absent=True)
    check(
        "offline validation represents an unpublished artifact explicitly",
        absent["status"] == "absent",
    )

    output_dir = os.path.join(root, "outputs")
    second_output_dir = os.path.join(root, "outputs-repeat")
    payload = _synthetic_payload(root, first)
    _validate_result_contract(payload)
    semantic_tampers = (
        (
            "result subject",
            lambda value: _set_path(
                value, ("result", "subject"), "HUP999_phaseII"),
        ),
        (
            "result scope",
            lambda value: _set_path(
                value, ("result", "analysis_scope"), "paired"),
        ),
        (
            "role order",
            lambda value: _set_path(
                value, ("result", "roles"), ["f4", "f3"]),
        ),
        (
            "channel identity",
            lambda value: _set_path(
                value,
                ("result", "channels_by_role", "f4"),
                "F04",
            ),
        ),
        (
            "event-contact identity",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM", "estimate",
                    "event_contact_ids",
                ),
                ["F3", "F3"],
            ),
        ),
        (
            "negative event count",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM",
                    "eligible_so_count_by_role", "f3",
                ),
                -1,
            ),
        ),
        (
            "unreconciled retained-event count",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM", "estimate",
                    "n_so_total",
                ),
                10_000,
            ),
        ),
        (
            "locked event minimum",
            lambda value: _set_path(
                value,
                ("result", "minimum_events_per_channel"),
                29,
            ),
        ),
        (
            "locked surrogate minimum",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM", "estimate",
                    "n_surrogate_pool",
                ),
                99,
            ),
        ),
        (
            "nonflat gate bypass",
            lambda value: _set_path(
                value,
                ("result", "sensor_nonflat_by_role", "f4"),
                False,
            ),
        ),
        (
            "stable-stage gate bypass",
            lambda value: (
                _set_path(
                    value,
                    ("result", "stages", "NREM", "stable_epochs"),
                    5,
                ),
                _set_path(
                    value,
                    ("result", "stages", "NREM", "stable_seconds"),
                    150,
                ),
            ),
        ),
        (
            "finite-stage gate bypass",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM",
                    "finite_stage_rr_samples",
                ),
                99,
            ),
        ),
        (
            "spurious available-stage support reason",
            lambda value: value["result"]["stages"]["NREM"][
                "support_reasons"].append(
                    "cardiac coverage below locked profile"),
        ),
        (
            "fixed monotone lag grid",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM", "estimate",
                    "lag_s", 1,
                ),
                -5.0,
            ),
        ),
        (
            "finite curve grid",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM", "estimate",
                    "curve", 0,
                ),
                float("nan"),
            ),
        ),
        (
            "RR-to-HR curve reconciliation",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM", "estimate",
                    "curve", 0,
                ),
                value["result"]["stages"]["NREM"][
                    "estimate"]["curve"][0] + 1.0,
            ),
        ),
        (
            "channel peak off grid",
            lambda value: (
                _set_path(
                    value,
                    (
                        "result", "stages", "NREM", "estimate",
                        "channel_peak_lag_s", 0,
                    ),
                    0.1,
                ),
                _set_path(
                    value,
                    (
                        "result", "stages", "NREM", "estimate",
                        "channel_peak_lag_by_role_s", "f3",
                    ),
                    0.1,
                ),
            ),
        ),
        (
            "bilateral peak derivation",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM", "estimate",
                    "peak_lag_s",
                ),
                4.75,
            ),
        ),
        (
            "participant-curve peak derivation",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM", "estimate",
                    "participant_curve_peak_lag_s",
                ),
                4.75,
            ),
        ),
        (
            "stage-relative magnitude",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM", "estimate",
                    "pct_above_stage_mean",
                ),
                999.0,
            ),
        ),
        (
            "local baseline",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM", "estimate",
                    "local_pre_event_mean_hr",
                ),
                999.0,
            ),
        ),
        (
            "local-change magnitude",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM", "estimate",
                    "event_locked_local_change_pct",
                ),
                999.0,
            ),
        ),
        (
            "peak-to-peak magnitude",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM", "estimate",
                    "event_curve_peak_to_peak_pct",
                ),
                999.0,
            ),
        ),
        (
            "diagnostic z derivation",
            lambda value: _set_path(
                value,
                (
                    "result", "stages", "NREM", "estimate",
                    "stage_shift_z_diagnostic",
                ),
                999.0,
            ),
        ),
    )
    for label, mutate in semantic_tampers:
        _check_semantic_tamper(payload, label, mutate)
    first_manifest = _write_artifacts(payload, output_dir)
    second_manifest = _write_artifacts(payload, second_output_dir)
    check(
        "all result-file bytes are deterministic across repeated writes",
        first_manifest["result_files_sha256"]
        == second_manifest["result_files_sha256"],
    )
    offline = validate_exploratory_artifacts(
        root, output_dir, mode="offline", allow_absent=False)
    check(
        "offline artifact validation proves bytes and skips only private input",
        offline["status"] == "complete"
        and offline["skipped_private_lineage"]
        == ["private/scalp/HUP138_phaseII.npz"],
    )
    stale_sidecar_path = os.path.join(
        root, "private", "scalp", "HUP138_phaseII.npz")
    os.makedirs(os.path.dirname(stale_sidecar_path), exist_ok=True)
    with open(stale_sidecar_path, "wb") as handle:
        handle.write(b"stale private bytes")
    offline_with_private = validate_exploratory_artifacts(
        root, output_dir, mode="offline", allow_absent=False)
    check(
        "offline validation never consumes opportunistically present private data",
        offline_with_private["status"] == "complete"
        and offline_with_private["skipped_private_lineage"]
        == ["private/scalp/HUP138_phaseII.npz"],
    )
    _validate_figure_bytes(payload, output_dir)
    with open(os.path.join(output_dir, RESULT_PNG), "ab") as handle:
        handle.write(b"rendering drift")
    try:
        _validate_figure_bytes(payload, output_dir)
    except RuntimeError:
        rendering_drift_rejected = True
    else:
        rendering_drift_rejected = False
    check(
        "publication-runtime rendering comparison rejects figure drift",
        rendering_drift_rejected,
    )
    _write_artifacts(payload, output_dir)
    try:
        validate_exploratory_artifacts(
            root, output_dir, mode="publication", allow_absent=False)
    except (FileNotFoundError, RuntimeError):
        private_required = True
    else:
        private_required = False
    check(
        "publication validation requires the ignored private sidecar",
        private_required,
    )

    portable_runtime = copy.deepcopy(payload)
    portable_runtime["provenance"]["runtime_versions"]["python"] += (
        "+synthetic-patch")
    portable_dir = os.path.join(root, "portable-runtime")
    _write_artifacts(portable_runtime, portable_dir)
    portable_offline = validate_exploratory_artifacts(
        root, portable_dir, mode="offline", allow_absent=False)
    try:
        validate_exploratory_artifacts(
            root, portable_dir, mode="publication", allow_absent=False)
    except RuntimeError as error:
        portable_publication_rejected = (
            "runtime_versions" in str(error)
        )
    else:
        portable_publication_rejected = False
    check(
        "offline accepts a valid producing-runtime patch while publication "
        "requires that runtime",
        portable_offline["status"] == "complete"
        and portable_publication_rejected,
    )

    for label, mutate in (
        (
            "source tree",
            lambda value: value["provenance"].__setitem__(
                "source_tree_sha256", "0" * 64),
        ),
        (
            "runtime metadata",
            lambda value: value["provenance"].__setitem__(
                "runtime_versions", {"python": "stale"}),
        ),
        (
            "profile",
            lambda value: value["method"].__setitem__(
                "source_profile_sha256", "0" * 64),
        ),
    ):
        stale = copy.deepcopy(payload)
        mutate(stale)
        stale_dir = os.path.join(root, f"stale-{label.replace(' ', '-')}")
        _write_artifacts(stale, stale_dir)
        try:
            validate_exploratory_artifacts(
                root, stale_dir, mode="offline", allow_absent=False)
        except RuntimeError:
            stale_rejected = True
        else:
            stale_rejected = False
        check(
            f"offline validation rejects stale {label} provenance",
            stale_rejected,
        )

    manifest_path = os.path.join(output_dir, "RUN_MANIFEST.json")
    with open(manifest_path, encoding="utf-8") as handle:
        manifest = json.load(handle)
    manifest["unexpected"] = True
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, sort_keys=True)
    try:
        validate_exploratory_artifacts(
            root, output_dir, mode="offline", allow_absent=False)
    except RuntimeError:
        extra_rejected = True
    else:
        extra_rejected = False
    check("terminal manifest rejects extra schema fields", extra_rejected)

with tempfile.TemporaryDirectory() as root:
    output_dir = os.path.join(root, "outputs")
    dirty_payload = _synthetic_payload(root, first)
    dirty_payload["provenance"]["code_dirty"] = True
    try:
        _write_artifacts(dirty_payload, output_dir)
    except RuntimeError:
        dirty_rejected = True
    else:
        dirty_rejected = False
    check(
        "dirty source state fails before publication writes",
        dirty_rejected and not os.path.exists(output_dir),
    )

print("\nAll HUP138 bilateral scalp-only regressions passed.")
