"""Fail-closed provenance for participant teaching figures.

The visual generator consumes private neutral caches, one private full RESPect
grid, the public locked HUP snapshot, and two RESP0699 metadata files.  This
module freezes those exact bytes before plotting and verifies that they did not
change before a terminal visual manifest is published.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import uuid
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
sys.path.insert(0, str(ANALYSIS))

import run_qc_grid as grid  # noqa: E402
from compact_qc_grid_artifacts import validate_public_artifacts  # noqa: E402
from pipeline_version import (  # noqa: E402
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    atomic_json_dump,
    file_sha256,
    git_revision,
    runtime_versions,
    source_tree_sha256,
    utc_now,
)
from qc_profiles import (  # noqa: E402
    load_qc_profile,
    profile_file_sha256,
    qc_profile_sha256,
)
from participant_visual_narratives import (  # noqa: E402
    render_resp0699_check,
    render_visual_readme,
)


PROFILE_ID = "overlap11_endpoint_local"
VISUAL_MANIFEST_SCHEMA = "2026-07-participant-result-visuals-v1"
VISUAL_PIPELINE = "make_participant_result_visuals"
GENERATOR_PATH = ROOT / "visualization" / "make_participant_result_visuals.py"
PROVENANCE_PATH = Path(__file__).resolve()
NARRATIVE_PATH = (
    ROOT / "visualization" / "participant_visual_narratives.py"
)
VISUAL_OUTPUT = ROOT / "outputs" / "participant_result_visuals"
FIGURE_STEMS = (
    "3A_RESP0699_all_three_measurements",
    "3A_RESP0699_coherence_threshold_explained",
    "3B_HUP160_event_locked_heart_rate",
    "3B_RESP0699_naji_frontal_ecog_check",
    "3D_HUP172_SO_spindle_phase",
)
VISUAL_RESULT_NAMES = tuple(
    f"{stem}.{suffix}"
    for stem in FIGURE_STEMS
    for suffix in ("png", "svg")
) + (
    "figure_values.json",
    "README.md",
    "RESP0699_THRESHOLD_AND_NAJI_CHECK.md",
)
PUBLIC_ROOT = ROOT / "outputs" / "qc_grid_public"
PUBLIC_MANIFEST = PUBLIC_ROOT / "RUN_MANIFEST.json"
HUP_LOCKED = (
    PUBLIC_ROOT / "locked" / "overlap11_endpoint_local__hup.json"
)
RESP_PUBLIC_SUMMARY = (
    PUBLIC_ROOT
    / "staging_window_support_v1"
    / "respect_qc_grid_summary.json"
)
RESP_FULL_GRID = (
    ROOT
    / "outputs"
    / "qc_grid"
    / "staging_window_support_v1"
    / "respect_qc_grid.json"
)
HUP_CACHE_DIR = ROOT / "data" / "derived" / "lc_infraslow"
RESP_CACHE_DIR = ROOT / "data" / "derived" / "ds003848"
RESP_METADATA = (
    ROOT
    / "data"
    / "ds003848_raw"
    / "sub-RESP0699_ses-1_task-sleep_run-030608_ieeg.json"
)
RESP_CHANNELS = (
    ROOT
    / "data"
    / "ds003848_raw"
    / "sub-RESP0699_ses-1_task-sleep_run-030608_channels.tsv"
)
SELECTED_CACHE_SUBJECTS = {
    "HUP": ("HUP160_phaseII", "HUP172_phaseII"),
    "RESPect": ("sub-RESP0699",),
}
_TRACKED_INPUT_PATHS = {
    "public_qc_manifest": PUBLIC_MANIFEST,
    "hup_locked_qc_artifact": HUP_LOCKED,
    "respect_public_staging_summary": RESP_PUBLIC_SUMMARY,
}
_PRIVATE_INPUT_PATHS = {
    "respect_full_staging_grid": RESP_FULL_GRID,
    "respect_metadata": RESP_METADATA,
    "respect_channels": RESP_CHANNELS,
}


def _reject_json_constant(value: str):
    raise ValueError(f"non-standard JSON constant {value!r}")


def _load_json(path: Path, label: str) -> dict:
    if not path.is_file():
        raise RuntimeError(f"missing {label}: {path}")
    try:
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle, parse_constant=_reject_json_constant)
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"{label} is not readable strict JSON: {path}") from exc
    if not isinstance(payload, dict):
        raise RuntimeError(f"{label} must be a JSON object: {path}")
    return payload


def repository_relative(path: Path, root: Path = ROOT) -> str:
    """Return a portable repository-relative path, rejecting path escapes."""
    resolved_root = root.resolve()
    resolved = path.resolve()
    try:
        relative = resolved.relative_to(resolved_root)
    except ValueError as exc:
        raise RuntimeError(
            f"visual input/output escapes repository root: {path}"
        ) from exc
    return relative.as_posix()


def _require_sha256(value, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise RuntimeError(f"{label} is not a lowercase SHA-256")
    return value


def validate_grid_lineage(
    payload: dict,
    *,
    label: str,
    grid_id: str,
    cache_pipeline: str,
    cache_manifest: dict,
    cache_manifest_sha256: str,
    cache_files_sha256: dict,
    analysis_source_sha256: str,
    profile_sha256: str,
    profile_file_sha256_value: str,
    selected_subjects: tuple[str, ...],
) -> dict:
    """Validate a full/locked grid against one current terminal cache run."""
    expected = {
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "source_tree_sha256": analysis_source_sha256,
        "cache_pipeline": cache_pipeline,
        "cache_manifest_sha256": cache_manifest_sha256,
        "cache_manifest_run_id": cache_manifest.get("run_id"),
        "cache_files_sha256": cache_files_sha256,
    }
    for key, value in expected.items():
        if payload.get(key) != value:
            raise RuntimeError(f"{label} differs from current cache/source at {key}")
    if grid_id and payload.get("grid_id") != grid_id:
        raise RuntimeError(f"{label} is not the required {grid_id} grid")

    profiles = payload.get("profiles")
    if not isinstance(profiles, list):
        raise RuntimeError(f"{label} has no profile list")
    matches = [
        profile
        for profile in profiles
        if profile.get("qc_profile_id") == PROFILE_ID
    ]
    if len(matches) != 1:
        raise RuntimeError(f"{label} does not contain exactly one locked profile")
    profile = matches[0]
    if (
        profile.get("qc_profile_sha256") != profile_sha256
        or profile.get("qc_profile_file_sha256")
        != profile_file_sha256_value
    ):
        raise RuntimeError(
            f"{label} locked profile differs from current configuration"
        )

    subject_records = profile.get("subjects")
    if not isinstance(subject_records, list):
        raise RuntimeError(f"{label} locked profile lacks full subject records")
    subjects = [
        record.get("subject")
        for record in subject_records
        if isinstance(record, dict)
    ]
    if len(subjects) != len(subject_records) or len(subjects) != len(set(subjects)):
        raise RuntimeError(f"{label} has malformed or duplicate subject records")
    missing = sorted(set(selected_subjects) - set(subjects))
    if missing:
        raise RuntimeError(f"{label} lacks selected visual subjects: {missing}")
    return profile


def _validated_cache_run(
    cache_dir: Path,
    *,
    expected_pipeline: str,
    selected_subjects: tuple[str, ...],
) -> dict:
    """Validate a current terminal cache run and freeze selected NPZ bytes."""
    manifest, manifest_path_text, cache_hashes = grid._cache_manifest(
        str(cache_dir)
    )
    manifest_path = Path(manifest_path_text)
    if manifest.get("pipeline") != expected_pipeline:
        raise RuntimeError(
            f"{manifest_path} pipeline differs from {expected_pipeline}"
        )
    completed = set(manifest.get("completed", []))
    missing = sorted(set(selected_subjects) - completed)
    if missing:
        raise RuntimeError(
            f"{manifest_path} does not complete selected visual subjects: {missing}"
        )
    selected = {}
    for subject in selected_subjects:
        cache_path = cache_dir / f"{subject}.npz"
        actual = file_sha256(cache_path)
        if cache_hashes.get(subject) != actual:
            raise RuntimeError(
                f"{cache_path} differs from its terminal cache manifest"
            )
        selected[subject] = {
            "path_relative": repository_relative(cache_path),
            "sha256": actual,
        }
    return {
        "manifest": manifest,
        "manifest_path": manifest_path,
        "manifest_sha256": file_sha256(manifest_path),
        "cache_files_sha256": cache_hashes,
        "selected": selected,
    }


def _record_file(path: Path) -> dict:
    if not path.is_file():
        raise RuntimeError(f"required visual input is missing: {path}")
    return {
        "path_relative": repository_relative(path),
        "sha256": file_sha256(path),
    }


def _validate_grid_inputs(
    *,
    profile_sha: str,
    profile_file_sha: str,
    analysis_source_sha: str,
    hup_run: dict,
    resp_run: dict,
) -> dict:
    """Validate full/locked grids and return their exact file records."""
    _load_json(PUBLIC_MANIFEST, "public QC terminal manifest")
    hup_locked = _load_json(HUP_LOCKED, "locked HUP QC artifact")
    resp_summary = _load_json(
        RESP_PUBLIC_SUMMARY, "public RESPect staging-grid summary"
    )
    resp_full = _load_json(RESP_FULL_GRID, "full RESPect staging grid")

    public_full_sha = _require_sha256(
        resp_summary.get("full_grid_sha256"),
        "RESPect public summary full_grid_sha256",
    )
    if file_sha256(RESP_FULL_GRID) != public_full_sha:
        raise RuntimeError(
            "full RESPect staging grid bytes differ from current public evidence"
        )

    validate_grid_lineage(
        hup_locked,
        label="locked HUP QC artifact",
        grid_id="event_count_oat_v1",
        cache_pipeline="cache_lc_series",
        cache_manifest=hup_run["manifest"],
        cache_manifest_sha256=hup_run["manifest_sha256"],
        cache_files_sha256=hup_run["cache_files_sha256"],
        analysis_source_sha256=analysis_source_sha,
        profile_sha256=profile_sha,
        profile_file_sha256_value=profile_file_sha,
        selected_subjects=SELECTED_CACHE_SUBJECTS["HUP"],
    )
    validate_grid_lineage(
        resp_full,
        label="full RESPect staging grid",
        grid_id="staging_window_support_v1",
        cache_pipeline="stage_ds003848",
        cache_manifest=resp_run["manifest"],
        cache_manifest_sha256=resp_run["manifest_sha256"],
        cache_files_sha256=resp_run["cache_files_sha256"],
        analysis_source_sha256=analysis_source_sha,
        profile_sha256=profile_sha,
        profile_file_sha256_value=profile_file_sha,
        selected_subjects=SELECTED_CACHE_SUBJECTS["RESPect"],
    )

    if (
        resp_summary.get("cache_manifest_sha256")
        != resp_run["manifest_sha256"]
        or resp_summary.get("cache_files_sha256")
        != resp_run["cache_files_sha256"]
    ):
        raise RuntimeError(
            "public RESPect summary does not describe the current cache run"
        )
    return {
        "public_qc_manifest": _record_file(PUBLIC_MANIFEST),
        "hup_locked_qc_artifact": _record_file(HUP_LOCKED),
        "respect_public_staging_summary": _record_file(RESP_PUBLIC_SUMMARY),
        "respect_full_staging_grid": _record_file(RESP_FULL_GRID),
        "respect_metadata": _record_file(RESP_METADATA),
        "respect_channels": _record_file(RESP_CHANNELS),
    }


def _cache_run_record(run: dict) -> dict:
    return {
        "pipeline": run["manifest"]["pipeline"],
        "run_id": run["manifest"].get("run_id"),
        "manifest": _record_file(run["manifest_path"]),
        "selected_cache_files": run["selected"],
    }


def _freeze_file_records(
    sources: dict,
    inputs: dict,
    cache_runs: dict,
) -> dict:
    frozen = {
        record["path_relative"]: record["sha256"]
        for record in [*sources.values(), *inputs.values()]
    }
    for cache_run in cache_runs.values():
        records = [
            cache_run["manifest"],
            *cache_run["selected_cache_files"].values(),
        ]
        frozen.update(
            {
                record["path_relative"]: record["sha256"]
                for record in records
            }
        )
    return frozen


def freeze_visual_inputs(
    *,
    generator_path: Path = GENERATOR_PATH,
    provenance_path: Path = Path(__file__),
) -> dict:
    """Validate and freeze every byte consumed by the visual generator."""
    validate_public_artifacts(str(PUBLIC_ROOT))
    profile = load_qc_profile(PROFILE_ID)
    profile_sha = qc_profile_sha256(profile)
    profile_file_sha = profile_file_sha256()
    analysis_source_sha = source_tree_sha256(str(ROOT))
    hup_run = _validated_cache_run(
        HUP_CACHE_DIR,
        expected_pipeline="cache_lc_series",
        selected_subjects=SELECTED_CACHE_SUBJECTS["HUP"],
    )
    resp_run = _validated_cache_run(
        RESP_CACHE_DIR,
        expected_pipeline="stage_ds003848",
        selected_subjects=SELECTED_CACHE_SUBJECTS["RESPect"],
    )
    sources = {
        "generator": _record_file(generator_path),
        "provenance_helper": _record_file(provenance_path),
        "narrative_renderer": _record_file(NARRATIVE_PATH),
    }
    inputs = _validate_grid_inputs(
        profile_sha=profile_sha,
        profile_file_sha=profile_file_sha,
        analysis_source_sha=analysis_source_sha,
        hup_run=hup_run,
        resp_run=resp_run,
    )
    cache_runs = {
        "HUP": _cache_run_record(hup_run),
        "RESPect": _cache_run_record(resp_run),
    }
    manifest_provenance = {
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "analysis_source_tree_sha256": analysis_source_sha,
        "profile_id": PROFILE_ID,
        "profile_sha256": profile_sha,
        "profile_file_sha256": profile_file_sha,
        "source_files_sha256": sources,
        "input_artifacts": inputs,
        "cache_runs": cache_runs,
    }
    return {
        "manifest_provenance": manifest_provenance,
        "frozen_files_sha256": _freeze_file_records(
            sources, inputs, cache_runs
        ),
    }


def assert_visual_inputs_unchanged(snapshot: dict) -> None:
    """Fail if any frozen visual input changed while figures were generated."""
    for relative, expected in snapshot["frozen_files_sha256"].items():
        path = ROOT / Path(relative)
        if not path.is_file() or file_sha256(path) != expected:
            raise RuntimeError(
                f"visual input changed after provenance freeze: {relative}"
            )


def visual_code_is_dirty(root: Path = ROOT) -> bool | None:
    """Check repository changes except the authorized visual output directory."""
    excluded = "outputs/participant_result_visuals"
    try:
        output = subprocess.check_output(
            [
                "git",
                "-C",
                str(root),
                "status",
                "--porcelain=v1",
                "--untracked-files=all",
                "--",
                ".",
                f":(exclude){excluded}",
                f":(exclude){excluded}/**",
            ],
            text=True,
            stderr=subprocess.DEVNULL,
            timeout=10,
        )
    except Exception:
        return None
    return bool(output.strip())


def preflight_manifest_base(root: Path = ROOT) -> dict:
    """Capture clean-tree status before invalidating a prior terminal run."""
    return {
        "schema_version": VISUAL_MANIFEST_SCHEMA,
        "pipeline": VISUAL_PIPELINE,
        "run_id": str(uuid.uuid4()),
        "generated_at_utc": utc_now(),
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "code_revision": git_revision(str(root)),
        "code_dirty_at_start": visual_code_is_dirty(root),
        "runtime_versions": runtime_versions(),
    }


def visual_manifest_base(snapshot: dict, preflight: dict) -> dict:
    """Attach validated input lineage to the preflight run identity."""
    return {
        **preflight,
        **snapshot["manifest_provenance"],
    }


def write_preflight_manifest(output_dir: Path) -> tuple[Path, dict]:
    """Invalidate any prior terminal manifest before input validation starts."""
    base = preflight_manifest_base()
    manifest_path = output_dir / "RUN_MANIFEST.json"
    atomic_json_dump(
        {
            **base,
            "run_state": "in_progress",
            "validation_state": "input_validation_pending",
            "result_files_sha256": {},
        },
        str(manifest_path),
    )
    return manifest_path, base


def write_in_progress_manifest(output_dir: Path, base: dict) -> Path:
    """Record frozen input lineage before writing result files."""
    manifest_path = output_dir / "RUN_MANIFEST.json"
    atomic_json_dump(
        {
            **base,
            "run_state": "in_progress",
            "validation_state": "inputs_frozen",
            "result_files_sha256": {},
        },
        str(manifest_path),
    )
    return manifest_path


def write_failed_manifest(
    manifest_path: Path,
    base: dict,
    *,
    validation_state: str,
    error: Exception,
) -> None:
    """Publish an explicit non-success state after a failed validation."""
    atomic_json_dump(
        {
            **base,
            "run_state": "failed",
            "validation_state": validation_state,
            "validation_error": f"{type(error).__name__}: {error}",
            "result_files_sha256": {},
        },
        str(manifest_path),
    )


def write_terminal_manifest(
    output_dir: Path,
    base: dict,
    result_paths: list[Path],
) -> Path:
    """Publish one terminal manifest last, hashing the exact result file set."""
    if base.get("code_dirty_at_start") is not False:
        raise RuntimeError(
            "terminal visual publication requires a clean worktree at run start"
        )
    relative_paths = [repository_relative(path) for path in result_paths]
    if len(relative_paths) != len(set(relative_paths)):
        raise RuntimeError("visual result file list contains duplicates")
    result_hashes = {
        relative: file_sha256(ROOT / Path(relative))
        for relative in sorted(relative_paths)
    }
    manifest_path = output_dir / "RUN_MANIFEST.json"
    atomic_json_dump(
        {
            **base,
            "run_state": "complete",
            "result_files_sha256": result_hashes,
        },
        str(manifest_path),
    )
    return manifest_path


def _validate_visual_narratives() -> None:
    """Require exact prose regenerated from the current figure-values JSON."""
    values_path = VISUAL_OUTPUT / "figure_values.json"
    values = _load_json(values_path, "participant visual figure values")
    expected = {
        VISUAL_OUTPUT / "README.md": render_visual_readme(values),
        VISUAL_OUTPUT / "RESP0699_THRESHOLD_AND_NAJI_CHECK.md": (
            render_resp0699_check(values)
        ),
    }
    for path, text in expected.items():
        if not path.is_file():
            raise RuntimeError(f"required visual narrative is missing: {path}")
        if path.read_text(encoding="utf-8") != text:
            raise RuntimeError(
                f"{path} is not the exact narrative regenerated from "
                "figure_values.json"
            )


def _validate_recorded_runtime_versions(
    recorded,
    *,
    require_current: bool,
) -> dict:
    """Validate producer-runtime structure without making offline proof local.

    A checked-in release remains verifiable after a Python or dependency patch
    update. Publication is stricter because it recomputes private inputs and
    rendered figure bytes in the current producing environment.
    """
    current = runtime_versions()
    if (
        not isinstance(recorded, dict)
        or set(recorded) != set(current)
        or any(
            not isinstance(key, str)
            or not key
            or not isinstance(value, str)
            or not value
            for key, value in recorded.items()
        )
    ):
        raise RuntimeError(
            "visual manifest has invalid recorded runtime_versions")
    if require_current and recorded != current:
        raise RuntimeError(
            "publication runtime_versions differ from the producing runtime")
    return recorded


def _validate_terminal_visual_manifest_against_snapshot(
    manifest_path: Path,
    expected_result_paths: list[Path],
    *,
    snapshot: dict,
    validate_narratives: bool = True,
) -> dict:
    """Pure terminal contract used after production provenance is collected."""
    manifest = _load_json(manifest_path, "participant visual terminal manifest")
    expected_metadata = {
        "schema_version": VISUAL_MANIFEST_SCHEMA,
        "pipeline": VISUAL_PIPELINE,
        "run_state": "complete",
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "code_dirty_at_start": False,
    }
    for key, expected in expected_metadata.items():
        if manifest.get(key) != expected:
            raise RuntimeError(f"{manifest_path} differs at {key}")
    _validate_recorded_runtime_versions(
        manifest.get("runtime_versions"),
        require_current=True,
    )
    for key, expected in snapshot["manifest_provenance"].items():
        if manifest.get(key) != expected:
            raise RuntimeError(
                f"{manifest_path} has stale or incomplete input lineage at {key}"
            )
    assert_visual_inputs_unchanged(snapshot)
    expected_relatives = {
        repository_relative(path) for path in expected_result_paths
    }
    hashes = manifest.get("result_files_sha256")
    if not isinstance(hashes, dict) or set(hashes) != expected_relatives:
        raise RuntimeError(
            f"{manifest_path} does not name the exact visual result set"
        )
    for relative, expected in hashes.items():
        _require_sha256(expected, f"result_files_sha256[{relative!r}]")
        path = ROOT / Path(relative)
        if not path.is_file() or file_sha256(path) != expected:
            raise RuntimeError(
                f"visual result bytes differ or are missing: {relative}"
            )
    if validate_narratives:
        _validate_visual_narratives()
    return manifest


def validate_terminal_visual_manifest(
    manifest_path: Path,
    expected_result_paths: list[Path],
) -> dict:
    """Recompute production provenance, then validate one terminal manifest."""
    snapshot = freeze_visual_inputs()
    return _validate_terminal_visual_manifest_against_snapshot(
        manifest_path,
        expected_result_paths,
        snapshot=snapshot,
    )


def _validate_file_record(
    record,
    *,
    path: Path,
    label: str,
    require_bytes: bool,
) -> None:
    if not isinstance(record, dict) or set(record) != {
        "path_relative",
        "sha256",
    }:
        raise RuntimeError(f"{label} is not an exact file record")
    expected_relative = repository_relative(path)
    if record.get("path_relative") != expected_relative:
        raise RuntimeError(f"{label} records the wrong repository path")
    expected_sha = _require_sha256(record.get("sha256"), f"{label} sha256")
    if require_bytes and (
        not path.is_file() or file_sha256(path) != expected_sha
    ):
        raise RuntimeError(f"{label} bytes differ or are missing")


def _validate_offline_input_lineage(manifest: dict) -> list[str]:
    inputs = manifest.get("input_artifacts")
    expected_keys = set(_TRACKED_INPUT_PATHS) | set(_PRIVATE_INPUT_PATHS)
    if not isinstance(inputs, dict) or set(inputs) != expected_keys:
        raise RuntimeError("visual manifest input_artifacts set is incomplete")
    for key, path in _TRACKED_INPUT_PATHS.items():
        _validate_file_record(
            inputs[key],
            path=path,
            label=f"input_artifacts[{key!r}]",
            require_bytes=True,
        )
    skipped = []
    for key, path in _PRIVATE_INPUT_PATHS.items():
        _validate_file_record(
            inputs[key],
            path=path,
            label=f"input_artifacts[{key!r}]",
            require_bytes=False,
        )
        skipped.append(repository_relative(path))
    return skipped


def _validate_offline_cache_lineage(manifest: dict) -> list[str]:
    cache_runs = manifest.get("cache_runs")
    if not isinstance(cache_runs, dict) or set(cache_runs) != {
        "HUP",
        "RESPect",
    }:
        raise RuntimeError("visual manifest cache_runs set is incomplete")
    expected = {
        "HUP": (
            "cache_lc_series",
            HUP_CACHE_DIR,
            SELECTED_CACHE_SUBJECTS["HUP"],
        ),
        "RESPect": (
            "stage_ds003848",
            RESP_CACHE_DIR,
            SELECTED_CACHE_SUBJECTS["RESPect"],
        ),
    }
    skipped = []
    for cohort, (pipeline, cache_dir, subjects) in expected.items():
        run = cache_runs[cohort]
        if (
            not isinstance(run, dict)
            or set(run) != {
                "pipeline",
                "run_id",
                "manifest",
                "selected_cache_files",
            }
            or run.get("pipeline") != pipeline
            or not isinstance(run.get("run_id"), str)
            or not run["run_id"]
        ):
            raise RuntimeError(f"visual manifest {cohort} cache run is malformed")
        manifest_path = cache_dir / "RUN_MANIFEST.json"
        _validate_file_record(
            run["manifest"],
            path=manifest_path,
            label=f"cache_runs[{cohort!r}].manifest",
            require_bytes=False,
        )
        selected = run.get("selected_cache_files")
        if not isinstance(selected, dict) or set(selected) != set(subjects):
            raise RuntimeError(
                f"visual manifest {cohort} selected cache set is incomplete"
            )
        skipped.append(repository_relative(manifest_path))
        for subject in subjects:
            cache_path = cache_dir / f"{subject}.npz"
            _validate_file_record(
                selected[subject],
                path=cache_path,
                label=f"cache_runs[{cohort!r}][{subject!r}]",
                require_bytes=False,
            )
            skipped.append(repository_relative(cache_path))
    return skipped


def _require_tracked(paths: list[Path]) -> None:
    for path in paths:
        relative = repository_relative(path)
        result = subprocess.run(
            ["git", "-C", str(ROOT), "ls-files", "--error-unmatch", relative],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if result.returncode:
            raise RuntimeError(f"required visual artifact is not tracked: {relative}")


def _require_exact_tracked_visual_set(
    manifest_path: Path,
    result_paths: list[Path],
) -> None:
    output_relative = repository_relative(VISUAL_OUTPUT)
    output = subprocess.check_output(
        [
            "git",
            "-C",
            str(ROOT),
            "ls-files",
            "-z",
            "--",
            output_relative,
        ]
    )
    actual = {
        value.decode("utf-8")
        for value in output.split(b"\0")
        if value
    }
    expected = {
        repository_relative(path)
        for path in [manifest_path, *result_paths]
    }
    if actual != expected:
        raise RuntimeError(
            "tracked participant-visual artifact set differs: "
            f"missing={sorted(expected - actual)}, "
            f"extra={sorted(actual - expected)}"
        )


def _validate_figure_values(path: Path, manifest: dict) -> None:
    values = _load_json(path, "participant visual figure values")
    for key in (
        "analysis_version",
        "cache_schema_version",
        "profile_id",
        "profile_sha256",
        "profile_file_sha256",
    ):
        if values.get(key) != manifest.get(key):
            raise RuntimeError(f"{path} differs from its manifest at {key}")
    figures = values.get("figures")
    expected_stems = {
        "3A": FIGURE_STEMS[0],
        "3A_threshold_explainer": FIGURE_STEMS[1],
        "3B": FIGURE_STEMS[2],
        "3B_RESP0699_Naji_check": FIGURE_STEMS[3],
        "3D": FIGURE_STEMS[4],
    }
    if not isinstance(figures, dict) or set(figures) != set(expected_stems):
        raise RuntimeError(f"{path} does not contain the exact figure set")
    for key, stem in expected_stems.items():
        record = figures[key]
        if (
            not isinstance(record, dict)
            or record.get("png")
            != repository_relative(VISUAL_OUTPUT / f"{stem}.png")
            or record.get("svg")
            != repository_relative(VISUAL_OUTPUT / f"{stem}.svg")
        ):
            raise RuntimeError(f"{path} has invalid paths for figure {key}")


def validate_checked_in_visual_artifacts(
    mode: str = "offline",
) -> dict:
    """Validate tracked visual evidence offline or with full private lineage."""
    if mode not in {"offline", "publication"}:
        raise ValueError("visual artifact validation mode must be offline/publication")
    manifest_path = VISUAL_OUTPUT / "RUN_MANIFEST.json"
    result_paths = [
        VISUAL_OUTPUT / name for name in VISUAL_RESULT_NAMES
    ]
    if mode == "publication":
        return {
            "manifest": validate_terminal_visual_manifest(
                manifest_path, result_paths
            ),
            "private_inputs_skipped": [],
        }

    manifest = _load_json(manifest_path, "participant visual terminal manifest")
    expected_metadata = {
        "schema_version": VISUAL_MANIFEST_SCHEMA,
        "pipeline": VISUAL_PIPELINE,
        "run_state": "complete",
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "code_dirty_at_start": False,
        "analysis_source_tree_sha256": source_tree_sha256(str(ROOT)),
        "profile_id": PROFILE_ID,
        "profile_sha256": qc_profile_sha256(load_qc_profile(PROFILE_ID)),
        "profile_file_sha256": profile_file_sha256(),
    }
    for key, expected in expected_metadata.items():
        if manifest.get(key) != expected:
            raise RuntimeError(f"{manifest_path} differs at {key}")
    _validate_recorded_runtime_versions(
        manifest.get("runtime_versions"),
        require_current=False,
    )
    sources = manifest.get("source_files_sha256")
    expected_sources = {
        "generator": GENERATOR_PATH,
        "provenance_helper": PROVENANCE_PATH,
        "narrative_renderer": NARRATIVE_PATH,
    }
    if not isinstance(sources, dict) or set(sources) != set(expected_sources):
        raise RuntimeError("visual manifest source file set is incomplete")
    for key, path in expected_sources.items():
        _validate_file_record(
            sources[key],
            path=path,
            label=f"source_files_sha256[{key!r}]",
            require_bytes=True,
        )
    skipped = _validate_offline_input_lineage(manifest)
    skipped.extend(_validate_offline_cache_lineage(manifest))
    validate_public_artifacts(str(PUBLIC_ROOT))

    hashes = manifest.get("result_files_sha256")
    expected_relatives = {
        repository_relative(path) for path in result_paths
    }
    if not isinstance(hashes, dict) or set(hashes) != expected_relatives:
        raise RuntimeError("visual manifest does not name the exact result set")
    for path in result_paths:
        relative = repository_relative(path)
        expected_sha = _require_sha256(
            hashes[relative], f"result_files_sha256[{relative!r}]"
        )
        if not path.is_file() or file_sha256(path) != expected_sha:
            raise RuntimeError(f"visual result bytes differ: {relative}")
    _validate_figure_values(VISUAL_OUTPUT / "figure_values.json", manifest)
    _validate_visual_narratives()
    _require_tracked(
        [
            manifest_path,
            *result_paths,
            *expected_sources.values(),
            *_TRACKED_INPUT_PATHS.values(),
        ]
    )
    _require_exact_tracked_visual_set(manifest_path, result_paths)
    return {
        "manifest": manifest,
        "private_inputs_skipped": sorted(skipped),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--artifact-validation",
        choices=("offline", "publication"),
        default="offline",
    )
    args = parser.parse_args()
    result = validate_checked_in_visual_artifacts(args.artifact_validation)
    print(
        json.dumps(
            {
                "run_state": result["manifest"]["run_state"],
                "private_inputs_skipped": result["private_inputs_skipped"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
