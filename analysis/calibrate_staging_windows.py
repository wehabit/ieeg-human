"""Calibrate clean 4-s staging-window support without inspecting LC/HR outcomes.

The neutral cache stores fourteen fixed Welch-window band powers for every contact/epoch.
This script takes real artifact-window patterns and applies each pattern to a fully clean reference
epoch from the same participant/contact.  It then measures how closely the partial-window SWA and
delta ratio reproduce that reference epoch's all-window values.

This is a measurement-stability calibration, not validation of the automated artifact detector or
of N2/N3 labels.  LC/heart-rate outcomes are never inputs to the calibration.
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
from scipy import stats

from pipeline_version import (
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    atomic_json_dump,
    cache_code_sha256,
    file_sha256,
    npz_scalar_text,
    runtime_versions,
    source_tree_sha256,
    utc_now,
)
from qc_profiles import PROFILE_PATH, PROFILE_SET_SCHEMA


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_CACHE = os.path.join(ROOT, "data", "derived", "ds003848")
DEFAULT_OUTPUT = os.path.join(
    ROOT, "outputs", "qc_calibration", "staging_window_calibration.json")
SUPPORT_LEVELS = tuple(range(1, 15))
REFERENCES_PER_MASK = 3

# These are engineering measurement-error targets, not paper requirements.  They are independent
# of 0.02-Hz peaks, HR coupling, SO timing, and spindle phase.
CALIBRATION_TARGETS = {
    "minimum_records": 100,
    "minimum_spearman": 0.95,
    "maximum_median_absolute_log_swa_error": 0.10,
    "maximum_p95_absolute_log_swa_error": 0.30,
    "maximum_median_absolute_delta_ratio_error": 0.02,
    "maximum_p95_absolute_delta_ratio_error": 0.05,
}
CALIBRATION_SOURCE_FILES = (
    "analysis/calibrate_staging_windows.py",
    "analysis/pipeline_version.py",
)


def _reject_json_constant(value):
    raise ValueError(f"non-standard JSON constant {value!r}")


def _load_json_object(path, label):
    with open(path, encoding="utf-8") as handle:
        value = json.load(handle, parse_constant=_reject_json_constant)
    if not isinstance(value, dict):
        raise RuntimeError(f"{label} must contain one JSON object")
    return value


def update_staging_calibration_pin(
    artifact_path,
    profile_path=PROFILE_PATH,
    *,
    root=ROOT,
    expected_recommendation=None,
):
    """Atomically pin one completed calibration without changing QC rules.

    The recommendation is a method decision, not a value that this provenance
    helper may silently tune. By default it must equal the recommendation
    already locked in the profile set. A changed recommendation therefore
    stops for explicit scientific review rather than rewriting thresholds.
    """
    root = os.path.abspath(root)
    artifact_path = os.path.abspath(artifact_path)
    profile_path = os.path.abspath(profile_path)
    try:
        artifact_relative = os.path.relpath(
            artifact_path, root).replace(os.sep, "/")
        if artifact_relative == ".." or artifact_relative.startswith("../"):
            raise ValueError
    except ValueError as error:
        raise RuntimeError(
            "calibration artifact must live inside the repository root"
        ) from error

    artifact = _load_json_object(
        artifact_path, "staging calibration artifact")
    profile_set = _load_json_object(profile_path, "QC profile set")
    if profile_set.get("schema_version") != PROFILE_SET_SCHEMA:
        raise RuntimeError("QC profile set schema is not current")
    method_config = profile_set.get("method_config")
    prior = (
        method_config.get("staging_calibration")
        if isinstance(method_config, dict)
        else None
    )
    if not isinstance(prior, dict):
        raise RuntimeError(
            "QC profile set lacks a staging-calibration contract")
    for key in ("calibration_scope", "hup_transport_status"):
        if not isinstance(prior.get(key), str) or not prior[key]:
            raise RuntimeError(
                f"locked staging-calibration context lacks {key}")
    if prior.get("not_a_paper_requirement") is not True:
        raise RuntimeError(
            "staging calibration must remain labeled non-paper-derived")

    for key, expected in (
        ("analysis_version", ANALYSIS_VERSION),
        ("cache_schema_version", CACHE_SCHEMA_VERSION),
    ):
        if artifact.get(key) != expected:
            raise RuntimeError(f"calibration artifact differs at {key}")
    required = {
        "cache_manifest_sha256",
        "cache_pipeline",
        "cache_run_id",
        "recommended_minimum_valid_windows",
        "support_results",
    }
    missing = sorted(required - set(artifact))
    if missing:
        raise RuntimeError(
            f"calibration artifact lacks required fields: {missing}")
    recommendation = artifact["recommended_minimum_valid_windows"]
    if (
        isinstance(recommendation, bool)
        or not isinstance(recommendation, int)
        or recommendation not in SUPPORT_LEVELS
    ):
        raise RuntimeError(
            "calibration artifact has no valid exact-support recommendation")
    locked_recommendation = (
        prior.get("recommended_minimum_valid_windows")
        if expected_recommendation is None
        else expected_recommendation
    )
    if recommendation != locked_recommendation:
        raise RuntimeError(
            "calibration recommendation changed from the locked method; "
            "review and update the QC method explicitly before pinning")
    exact_row = next(
        (
            row for row in artifact["support_results"]
            if row.get("minimum_valid_windows") == recommendation
        ),
        None,
    )
    if (
        not isinstance(exact_row, dict)
        or not exact_row.get("exact_support", {}).get(
            "meets_calibration_targets")
    ):
        raise RuntimeError(
            "recommended exact-support stratum does not meet calibration targets")
    cache_manifest_sha = artifact.get("cache_manifest_sha256")
    if (
        not isinstance(cache_manifest_sha, str)
        or len(cache_manifest_sha) != 64
        or any(character not in "0123456789abcdef"
               for character in cache_manifest_sha)
    ):
        raise RuntimeError(
            "calibration artifact has invalid cache_manifest_sha256")
    if (
        not isinstance(artifact.get("cache_run_id"), str)
        or not artifact["cache_run_id"]
    ):
        raise RuntimeError("calibration artifact has invalid cache_run_id")
    if not isinstance(artifact.get("cache_pipeline"), str):
        raise RuntimeError("calibration artifact has invalid cache_pipeline")

    pin = {
        **prior,
        "analysis_version": artifact["analysis_version"],
        "artifact_relative_path": artifact_relative,
        "artifact_sha256": file_sha256(artifact_path),
        "cache_manifest_sha256": artifact["cache_manifest_sha256"],
        "cache_pipeline": artifact["cache_pipeline"],
        "cache_run_id": artifact["cache_run_id"],
        "cache_schema_version": artifact["cache_schema_version"],
        "recommended_minimum_valid_windows": recommendation,
    }
    method_config["staging_calibration"] = pin
    atomic_json_dump(profile_set, profile_path)
    return pin


def calibration_source_files_sha256():
    """Hash the code that can change this calibration independently of the full tree.

    ``source_tree_sha256`` is still retained as an execution snapshot, but it necessarily changes
    when the generated artifact is pinned in ``qc_profiles_v1.json``.  These purpose-specific
    hashes avoid that circularity and let later consumers prove that the calibration algorithm
    itself has not changed.
    """
    return {
        relative: file_sha256(os.path.join(ROOT, *relative.split("/")))
        for relative in CALIBRATION_SOURCE_FILES
    }


def _manifest_subjects(cache_dir):
    path = os.path.join(cache_dir, "RUN_MANIFEST.json")
    with open(path) as handle:
        manifest = json.load(handle)
    if manifest.get("run_state") != "complete":
        raise RuntimeError(f"{path} is not a terminal complete cache run")
    if manifest.get("cache_schema_version") != CACHE_SCHEMA_VERSION:
        raise RuntimeError(
            f"{path} has cache schema {manifest.get('cache_schema_version')!r}; "
            f"expected {CACHE_SCHEMA_VERSION!r}")
    if manifest.get("analysis_version") != ANALYSIS_VERSION:
        raise RuntimeError(
            f"{path} has analysis version {manifest.get('analysis_version')!r}; "
            f"expected {ANALYSIS_VERSION!r}")
    allowed_pipelines = {"stage_ds003848", "cache_lc_series"}
    if manifest.get("pipeline") not in allowed_pipelines:
        raise RuntimeError(
            f"{path} pipeline {manifest.get('pipeline')!r} is not a "
            f"neutral-cache producer")
    if manifest.get("runtime_versions") != runtime_versions():
        raise RuntimeError(
            f"{path} runtime versions differ from the calibration runtime")
    expected_digest = cache_code_sha256(ROOT)
    if manifest.get("config", {}).get("cache_code_sha256") != expected_digest:
        raise RuntimeError(f"{path} was produced by different cache-building source")
    if manifest.get("failed"):
        raise RuntimeError(f"{path} contains failed subjects")

    requested = list(manifest.get("requested", []))
    completed = list(manifest.get("completed", []))
    skipped_entries = list(manifest.get("skipped", []))
    if (
        len(requested) != len(set(requested))
        or len(completed) != len(set(completed))
        or any(
            not isinstance(value, dict)
            or not value.get("subject")
            or not value.get("reason")
            for value in skipped_entries
        )
    ):
        raise RuntimeError(f"{path} has malformed or duplicate subject records")
    skipped = [value["subject"] for value in skipped_entries]
    if len(skipped) != len(set(skipped)):
        raise RuntimeError(f"{path} has duplicate skipped subjects")
    if (
        set(completed) & set(skipped)
        or set(completed) | set(skipped) != set(requested)
    ):
        raise RuntimeError(
            f"{path} completed/skipped subjects do not exactly partition requested subjects")
    recorded_hashes = manifest.get("result_files_sha256")
    if (
        not isinstance(recorded_hashes, dict)
        or set(recorded_hashes) != set(requested)
    ):
        raise RuntimeError(f"{path} lacks exact hashes for every terminal cache")
    for subject in completed + skipped:
        cache_path = os.path.join(cache_dir, f"{subject}.npz")
        if (
            not os.path.isfile(cache_path)
            or file_sha256(cache_path) != recorded_hashes[subject]
        ):
            raise RuntimeError(
                f"{cache_path} bytes differ from the terminal cache manifest")
        with np.load(cache_path, allow_pickle=False) as cache:
            expected_status = "ok" if subject in completed else "skip"
            if (
                npz_scalar_text(cache, "subject") != subject
                or npz_scalar_text(cache, "status") != expected_status
                or npz_scalar_text(
                    cache, "cache_schema_version") != CACHE_SCHEMA_VERSION
                or npz_scalar_text(
                    cache, "cache_code_sha256") != expected_digest
            ):
                raise RuntimeError(
                    f"{cache_path} embedded identity/status/schema/source digest is invalid")
    return completed, manifest


def _finite_correlation(left, right):
    left = np.asarray(left, float)
    right = np.asarray(right, float)
    keep = np.isfinite(left) & np.isfinite(right)
    if keep.sum() < 3 or np.std(left[keep]) == 0 or np.std(right[keep]) == 0:
        return None
    return float(stats.spearmanr(left[keep], right[keep]).statistic)


def calibration_records(cache_dir, subjects, expected_hashes):
    records = []
    cache_hashes = {}
    for subject in subjects:
        path = os.path.join(cache_dir, f"{subject}.npz")
        if not os.path.isfile(path):
            raise RuntimeError(f"manifest-completed cache is missing: {path}")
        actual_hash = file_sha256(path)
        if expected_hashes.get(subject) != actual_hash:
            raise RuntimeError(
                f"{path} bytes differ from the terminal cache manifest")
        cache_hashes[subject] = actual_hash
        with np.load(path, allow_pickle=False) as cache:
            if (
                npz_scalar_text(cache, "subject") != subject
                or npz_scalar_text(cache, "status") != "ok"
                or npz_scalar_text(
                    cache, "cache_schema_version") != CACHE_SCHEMA_VERSION
                or npz_scalar_text(
                    cache, "cache_code_sha256") != cache_code_sha256(ROOT)
            ):
                raise RuntimeError(
                    f"{path} embedded identity/status/schema/source digest is invalid")
            required = {
                "ep_valid_welch_window_mask_by_contact",
                "ep_window_swa_power_by_contact",
                "ep_window_total_power_by_contact",
                "cortical_chans",
            }
            missing = sorted(required - set(cache.files))
            if missing:
                raise RuntimeError(f"{path} lacks neutral staging fields: {missing}")
            valid = np.asarray(
                cache["ep_valid_welch_window_mask_by_contact"], bool)
            window_swa = np.asarray(
                cache["ep_window_swa_power_by_contact"], float)
            window_total = np.asarray(
                cache["ep_window_total_power_by_contact"], float)
            contacts = [str(value) for value in cache["cortical_chans"]]
        if valid.shape != window_swa.shape or valid.shape != window_total.shape:
            raise RuntimeError(f"{path} staging-window arrays do not align")
        if valid.ndim != 3 or valid.shape[2] != 14:
            raise RuntimeError(f"{path} does not contain fourteen Welch windows per epoch")
        for contact_index, contact in enumerate(contacts):
            finite_windows = (
                valid[contact_index]
                & np.isfinite(window_swa[contact_index])
                & np.isfinite(window_total[contact_index])
                & (window_swa[contact_index] > 0)
                & (window_total[contact_index] > 0)
            )
            full_epochs = np.where(finite_windows.all(axis=1))[0]
            if not len(full_epochs):
                continue
            for pattern_epoch, pattern in enumerate(finite_windows):
                n_valid = int(pattern.sum())
                if n_valid < 1:
                    continue
                # Apply each empirical mask geometry to several fully clean epochs.  This avoids
                # making the recommendation hinge on one arbitrary pattern-to-reference pairing.
                n_references = min(REFERENCES_PER_MASK, len(full_epochs))
                offsets = np.linspace(
                    0, len(full_epochs) - 1, n_references, dtype=int)
                reference_positions = np.unique(
                    (pattern_epoch + offsets) % len(full_epochs))
                for reference_position in reference_positions:
                    reference_epoch = int(full_epochs[reference_position])
                    ref_swa_windows = window_swa[contact_index, reference_epoch]
                    ref_total_windows = window_total[contact_index, reference_epoch]
                    reference_swa = float(np.mean(ref_swa_windows))
                    reference_total = float(np.mean(ref_total_windows))
                    partial_swa = float(np.mean(ref_swa_windows[pattern]))
                    partial_total = float(np.mean(ref_total_windows[pattern]))
                    if min(
                            reference_swa, reference_total,
                            partial_swa, partial_total) <= 0:
                        continue
                    reference_dr = reference_swa / reference_total
                    partial_dr = partial_swa / partial_total
                    records.append(dict(
                        subject=subject,
                        contact=contact,
                        n_valid_windows=n_valid,
                        reference_epoch=reference_epoch,
                        pattern_epoch=int(pattern_epoch),
                        reference_swa=reference_swa,
                        partial_swa=partial_swa,
                        reference_delta_ratio=reference_dr,
                        partial_delta_ratio=partial_dr,
                        absolute_log_swa_error=abs(
                            np.log(partial_swa / reference_swa)),
                        absolute_delta_ratio_error=abs(
                            partial_dr - reference_dr),
                    ))
    return records, cache_hashes


def _metrics(selected):
    """Measurement-error summary for one explicitly defined record subset."""
    selected = list(selected)
    if selected:
        log_error = np.asarray([
            value["absolute_log_swa_error"] for value in selected], float)
        dr_error = np.asarray([
            value["absolute_delta_ratio_error"] for value in selected], float)
        reference_swa = np.asarray([
            value["reference_swa"] for value in selected], float)
        partial_swa = np.asarray([
            value["partial_swa"] for value in selected], float)
        reference_dr = np.asarray([
            value["reference_delta_ratio"] for value in selected], float)
        partial_dr = np.asarray([
            value["partial_delta_ratio"] for value in selected], float)
    else:
        log_error = dr_error = reference_swa = partial_swa = np.array([])
        reference_dr = partial_dr = np.array([])
    return dict(
        n_records=len(selected),
        n_subjects=len({value["subject"] for value in selected}),
        median_absolute_log_swa_error=(
            float(np.median(log_error)) if len(log_error) else None),
        p95_absolute_log_swa_error=(
            float(np.percentile(log_error, 95)) if len(log_error) else None),
        median_absolute_delta_ratio_error=(
            float(np.median(dr_error)) if len(dr_error) else None),
        p95_absolute_delta_ratio_error=(
            float(np.percentile(dr_error, 95)) if len(dr_error) else None),
        swa_spearman=_finite_correlation(reference_swa, partial_swa),
        delta_ratio_spearman=_finite_correlation(reference_dr, partial_dr),
    )


def _meets_targets(row):
    targets = CALIBRATION_TARGETS
    return bool(
        row["n_records"] >= targets["minimum_records"]
        and row["swa_spearman"] is not None
        and row["delta_ratio_spearman"] is not None
        and row["swa_spearman"] >= targets["minimum_spearman"]
        and row["delta_ratio_spearman"] >= targets["minimum_spearman"]
        and row["median_absolute_log_swa_error"]
            <= targets["maximum_median_absolute_log_swa_error"]
        and row["p95_absolute_log_swa_error"]
            <= targets["maximum_p95_absolute_log_swa_error"]
        and row["median_absolute_delta_ratio_error"]
            <= targets["maximum_median_absolute_delta_ratio_error"]
        and row["p95_absolute_delta_ratio_error"]
            <= targets["maximum_p95_absolute_delta_ratio_error"]
    )


def summarize(records):
    rows = []
    for minimum in SUPPORT_LEVELS:
        exact = _metrics(
            value for value in records
            if value["n_valid_windows"] == minimum)
        accepted = _metrics(
            value for value in records
            if value["n_valid_windows"] >= minimum)
        exact["meets_calibration_targets"] = _meets_targets(exact)
        accepted["meets_calibration_targets"] = _meets_targets(accepted)
        rows.append(dict(
            minimum_valid_windows=minimum,
            exact_support=exact,
            accepted_population_at_or_above=accepted,
        ))
    exact_passing = [
        row["minimum_valid_windows"]
        for row in rows
        if row["exact_support"]["meets_calibration_targets"]
    ]
    accepted_passing = [
        row["minimum_valid_windows"]
        for row in rows
        if row["accepted_population_at_or_above"]["meets_calibration_targets"]
    ]
    return (
        rows,
        min(exact_passing) if exact_passing else None,
        min(accepted_passing) if accepted_passing else None,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", default=DEFAULT_CACHE)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--update-profile-pin",
        action="store_true",
        help=(
            "after writing the artifact, atomically update only its provenance "
            "pin in the locked QC profile file; a changed recommendation fails"
        ),
    )
    parser.add_argument(
        "--profile-file",
        default=PROFILE_PATH,
        help="QC profile set to update with --update-profile-pin",
    )
    args = parser.parse_args()
    cache_dir = os.path.abspath(args.cache_dir)
    subjects, manifest = _manifest_subjects(cache_dir)
    records, cache_hashes = calibration_records(
        cache_dir, subjects, manifest["result_files_sha256"])
    rows, exact_recommended, accepted_recommended = summarize(records)
    result = dict(
        analysis_version=ANALYSIS_VERSION,
        cache_schema_version=CACHE_SCHEMA_VERSION,
        generated_at_utc=utc_now(),
        source_tree_sha256=source_tree_sha256(ROOT),
        calibration_source_files_sha256=calibration_source_files_sha256(),
        cache_directory_relative=os.path.relpath(cache_dir, ROOT).replace(os.sep, "/"),
        cache_manifest_sha256=file_sha256(
            os.path.join(cache_dir, "RUN_MANIFEST.json")),
        cache_run_id=manifest["run_id"],
        cache_pipeline=manifest["pipeline"],
        cache_runtime_versions=manifest["runtime_versions"],
        cache_files_sha256=cache_hashes,
        runtime_versions=runtime_versions(),
        subjects=subjects,
        n_calibration_records=len(records),
        references_per_empirical_mask=REFERENCES_PER_MASK,
        calibration_targets=CALIBRATION_TARGETS,
        support_results=rows,
        recommended_minimum_valid_windows=exact_recommended,
        recommendation_basis=(
            "first exact-support stratum meeting every reconstruction target; "
            "more conservative than pooling cleaner epochs above the threshold"),
        accepted_population_recommended_minimum_valid_windows=accepted_recommended,
        selection_blinding=(
            "calibrated only against all-window SWA/delta-ratio reconstruction; "
            "0.02-Hz, HR, SO timing, and spindle-phase outcomes were not inputs"),
        limitations=[
            "does not validate the automatic artifact detector against human labels",
            "does not validate HUP/RESPect proxy sleep stages against PSG scoring",
            "reuses real window-mask geometry on fully clean epochs; it does not model undetected artifacts",
        ],
    )
    output_path = os.path.abspath(args.output)
    atomic_json_dump(result, output_path)
    if args.update_profile_pin:
        pin = update_staging_calibration_pin(
            output_path,
            os.path.abspath(args.profile_file),
        )
        print(
            f"updated calibration pin in {args.profile_file}: "
            f"sha256={pin['artifact_sha256']}",
            flush=True,
        )
    print(
        f"wrote {args.output}: {len(records)} records; "
        f"exact-support minimum={exact_recommended}; "
        f"accepted-population minimum={accepted_recommended}", flush=True)


if __name__ == "__main__":
    main()
