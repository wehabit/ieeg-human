"""Focused regressions for participant-visual provenance and publication."""
from __future__ import annotations

import copy
import json
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "visualization"))

from participant_visual_provenance import (  # noqa: E402
    ANALYSIS_VERSION,
    CACHE_SCHEMA_VERSION,
    PROFILE_ID,
    VISUAL_RESULT_NAMES,
    _validate_recorded_runtime_versions,
    _validate_terminal_visual_manifest_against_snapshot,
    assert_visual_inputs_unchanged,
    file_sha256,
    preflight_manifest_base,
    repository_relative,
    runtime_versions,
    validate_grid_lineage,
    visual_manifest_base,
    visual_code_is_dirty,
    write_failed_manifest,
    write_in_progress_manifest,
    write_preflight_manifest,
    write_terminal_manifest,
)
from participant_visual_narratives import (  # noqa: E402
    render_resp0699_check,
    render_visual_readme,
)
from pipeline_version import atomic_json_dump  # noqa: E402


def check(name: str, condition: bool) -> None:
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


def raises_runtime(callable_value) -> bool:
    try:
        callable_value()
    except RuntimeError:
        return True
    return False


source_sha = "a" * 64
manifest_sha = "b" * 64
profile_sha = "c" * 64
profile_file_sha = "d" * 64
cache_hashes = {"S1": "e" * 64, "S2": "f" * 64}
cache_manifest = {"run_id": "synthetic-run"}
grid_payload = {
    "analysis_version": ANALYSIS_VERSION,
    "cache_schema_version": CACHE_SCHEMA_VERSION,
    "source_tree_sha256": source_sha,
    "grid_id": "staging_window_support_v1",
    "cache_pipeline": "stage_ds003848",
    "cache_manifest_sha256": manifest_sha,
    "cache_manifest_run_id": "synthetic-run",
    "cache_files_sha256": cache_hashes,
    "profiles": [
        {
            "qc_profile_id": PROFILE_ID,
            "qc_profile_sha256": profile_sha,
            "qc_profile_file_sha256": profile_file_sha,
            "subjects": [
                {"subject": "S1", "result_3a": {}},
                {"subject": "S2", "result_3a": {}},
            ],
        }
    ],
}


def validate_synthetic_grid(payload, *subjects):
    return validate_grid_lineage(
        payload,
        label="synthetic grid",
        grid_id="staging_window_support_v1",
        cache_pipeline="stage_ds003848",
        cache_manifest=cache_manifest,
        cache_manifest_sha256=manifest_sha,
        cache_files_sha256=cache_hashes,
        analysis_source_sha256=source_sha,
        profile_sha256=profile_sha,
        profile_file_sha256_value=profile_file_sha,
        selected_subjects=subjects,
    )


validated = validate_synthetic_grid(grid_payload, "S1")
check(
    "current grid/cache/profile lineage is accepted",
    validated["qc_profile_id"] == PROFILE_ID,
)
check(
    "historical grid analysis version is rejected",
    raises_runtime(
        lambda: validate_synthetic_grid(
            {**grid_payload, "analysis_version": "historical-v8"},
            "S1",
        )
    ),
)
check(
    "grid with a different terminal cache hash is rejected",
    raises_runtime(
        lambda: validate_grid_lineage(
            grid_payload,
            label="wrong-cache grid",
            grid_id="staging_window_support_v1",
            cache_pipeline="stage_ds003848",
            cache_manifest=cache_manifest,
            cache_manifest_sha256="0" * 64,
            cache_files_sha256=cache_hashes,
            analysis_source_sha256=source_sha,
            profile_sha256=profile_sha,
            profile_file_sha256_value=profile_file_sha,
            selected_subjects=("S1",),
        )
    ),
)
check(
    "grid missing a selected participant is rejected",
    raises_runtime(lambda: validate_synthetic_grid(grid_payload, "S3")),
)
check(
    "absolute/outside paths cannot be serialized as repository-relative",
    raises_runtime(lambda: repository_relative(Path("/tmp/outside-visual.png"))),
)
check(
    "offline result contract includes both written visual narratives",
    {
        "README.md",
        "RESP0699_THRESHOLD_AND_NAJI_CHECK.md",
    }.issubset(VISUAL_RESULT_NAMES),
)

portable_runtime = runtime_versions()
portable_runtime["python"] += "+synthetic-patch"
check(
    "offline visual verification accepts a valid producing-runtime patch",
    _validate_recorded_runtime_versions(
        portable_runtime, require_current=False
    ) == portable_runtime,
)
check(
    "publication visual verification requires the producing runtime",
    raises_runtime(
        lambda: _validate_recorded_runtime_versions(
            portable_runtime, require_current=True
        )
    ),
)
check(
    "offline visual verification rejects malformed runtime metadata",
    raises_runtime(
        lambda: _validate_recorded_runtime_versions(
            {"python": "stale"}, require_current=False
        )
    ),
)

with (
    ROOT
    / "outputs"
    / "participant_result_visuals"
    / "figure_values.json"
).open(encoding="utf-8") as handle:
    narrative_fixture = json.load(handle)
narrative_fixture.setdefault("analysis_version", ANALYSIS_VERSION)
narrative_fixture.setdefault("profile_id", PROFILE_ID)
rendered_readme = render_visual_readme(narrative_fixture)
rendered_check = render_resp0699_check(narrative_fixture)
changed_fixture = copy.deepcopy(narrative_fixture)
changed_fixture["figures"]["3A"]["peak_hz"] += 0.001
check(
    "visual narratives are deterministic functions of figure values",
    "v9" in rendered_readme
    and "Naji et al. 2019" in rendered_check
    and render_visual_readme(changed_fixture) != rendered_readme,
)

with tempfile.TemporaryDirectory() as git_directory:
    git_root = Path(git_directory)
    (git_root / "analysis").mkdir()
    visual_output = git_root / "outputs" / "participant_result_visuals"
    visual_output.mkdir(parents=True)
    source_path = git_root / "analysis" / "source.py"
    visual_manifest = visual_output / "RUN_MANIFEST.json"
    source_path.write_text("# clean source\n", encoding="utf-8")
    visual_manifest.write_text('{"run_state": "complete"}\n', encoding="utf-8")
    subprocess.run(["git", "init", "-q", str(git_root)], check=True)
    subprocess.run(
        ["git", "-C", str(git_root), "config", "user.email", "test@example.com"],
        check=True,
    )
    subprocess.run(
        ["git", "-C", str(git_root), "config", "user.name", "Test"],
        check=True,
    )
    subprocess.run(["git", "-C", str(git_root), "add", "."], check=True)
    subprocess.run(
        ["git", "-C", str(git_root), "commit", "-qm", "fixture"],
        check=True,
    )
    check(
        "clean synthetic repository passes the scoped code check",
        visual_code_is_dirty(git_root) is False,
    )
    visual_manifest.write_text('{"run_state": "failed"}\n', encoding="utf-8")
    check(
        "a failed visual manifest does not poison the next retry",
        visual_code_is_dirty(git_root) is False,
    )
    source_path.write_text("# dirty source\n", encoding="utf-8")
    check(
        "a source change still fails the scoped code check",
        visual_code_is_dirty(git_root) is True,
    )

output_parent = ROOT / "outputs"
with tempfile.TemporaryDirectory(
    prefix=".participant-visual-test-",
    dir=output_parent,
) as directory:
    output_dir = Path(directory)
    first = output_dir / "first.json"
    second = output_dir / "second.svg"
    generator = output_dir / "generator.py"
    helper = output_dir / "helper.py"
    frozen_input = output_dir / "input.bin"
    first.write_text('{"value": 1}\n', encoding="utf-8")
    second.write_text("<svg/>", encoding="utf-8")
    generator.write_text("# generator\n", encoding="utf-8")
    helper.write_text("# helper\n", encoding="utf-8")
    frozen_input.write_bytes(b"input")

    def record(path):
        return {
            "path_relative": repository_relative(path),
            "sha256": file_sha256(path),
        }

    source_records = {
        "generator": record(generator),
        "provenance_helper": record(helper),
    }
    input_records = {"synthetic_input": record(frozen_input)}
    frozen_records = {
        value["path_relative"]: value["sha256"]
        for value in [*source_records.values(), *input_records.values()]
    }
    snapshot = {
        "manifest_provenance": {
            "analysis_version": ANALYSIS_VERSION,
            "cache_schema_version": CACHE_SCHEMA_VERSION,
            "profile_id": PROFILE_ID,
            "profile_sha256": profile_sha,
            "profile_file_sha256": profile_file_sha,
            "analysis_source_tree_sha256": source_sha,
            "source_files_sha256": source_records,
            "input_artifacts": input_records,
            "cache_runs": {},
        },
        "frozen_files_sha256": frozen_records,
    }
    clean_preflight = {
        **preflight_manifest_base(),
        "code_dirty_at_start": False,
    }
    base = visual_manifest_base(snapshot, clean_preflight)

    old_terminal = {
        **base,
        "run_state": "complete",
        "result_files_sha256": {"obsolete": "0" * 64},
    }
    manifest_path = output_dir / "RUN_MANIFEST.json"
    atomic_json_dump(old_terminal, str(manifest_path))
    _, captured_preflight = write_preflight_manifest(output_dir)
    with manifest_path.open(encoding="utf-8") as handle:
        invalidated = json.load(handle)
    check(
        "preflight invalidates an old terminal manifest before input validation",
        invalidated["run_state"] == "in_progress"
        and invalidated["validation_state"] == "input_validation_pending"
        and invalidated["result_files_sha256"] == {},
    )

    dirty_base = {
        **visual_manifest_base(snapshot, captured_preflight),
        "code_dirty_at_start": True,
    }
    check(
        "dirty-at-start runs cannot publish terminal visual evidence",
        raises_runtime(
            lambda: write_terminal_manifest(
                output_dir, dirty_base, [first, second]
            )
        ),
    )

    write_in_progress_manifest(output_dir, base)
    check(
        "in-progress manifest cannot validate as terminal",
        raises_runtime(
            lambda: _validate_terminal_visual_manifest_against_snapshot(
                manifest_path,
                [first, second],
                snapshot=snapshot,
                validate_narratives=False,
            )
        ),
    )
    simulated_plot_error = RuntimeError("synthetic plot failure")
    write_failed_manifest(
        manifest_path,
        base,
        validation_state="terminal_validation_failed",
        error=simulated_plot_error,
    )
    with manifest_path.open(encoding="utf-8") as handle:
        failed = json.load(handle)
    check(
        "plot/JSON failures publish failed rather than lingering in-progress",
        failed["run_state"] == "failed"
        and failed["validation_state"] == "terminal_validation_failed"
        and "synthetic plot failure" in failed["validation_error"],
    )

    write_terminal_manifest(output_dir, base, [first, second])
    terminal = _validate_terminal_visual_manifest_against_snapshot(
        manifest_path,
        [first, second],
        snapshot=snapshot,
        validate_narratives=False,
    )
    check(
        "terminal manifest uses relative paths and exact output hashes",
        terminal["run_state"] == "complete"
        and len(terminal["result_files_sha256"]) == 2
        and all(
            not Path(relative).is_absolute()
            for relative in terminal["result_files_sha256"]
        ),
    )

    stale_lineage = dict(terminal)
    stale_lineage["analysis_source_tree_sha256"] = "9" * 64
    atomic_json_dump(stale_lineage, str(manifest_path))
    check(
        "manifest-rehashed stale source lineage is rejected",
        raises_runtime(
            lambda: _validate_terminal_visual_manifest_against_snapshot(
                manifest_path,
                [first, second],
                snapshot=snapshot,
                validate_narratives=False,
            )
        ),
    )

    incomplete = dict(terminal)
    incomplete.pop("cache_runs")
    atomic_json_dump(incomplete, str(manifest_path))
    check(
        "terminal manifest missing frozen lineage is rejected",
        raises_runtime(
            lambda: _validate_terminal_visual_manifest_against_snapshot(
                manifest_path,
                [first, second],
                snapshot=snapshot,
                validate_narratives=False,
            )
        ),
    )

    write_terminal_manifest(output_dir, base, [first, second])
    generator.write_text("# tampered generator\n", encoding="utf-8")
    check(
        "generator/helper source tampering after freeze is rejected",
        raises_runtime(
            lambda: _validate_terminal_visual_manifest_against_snapshot(
                manifest_path,
                [first, second],
                snapshot=snapshot,
                validate_narratives=False,
            )
        ),
    )
    generator.write_text("# generator\n", encoding="utf-8")

    first.write_text('{"value": 2}\n', encoding="utf-8")
    check(
        "output modification after terminal publication is rejected",
        raises_runtime(
            lambda: _validate_terminal_visual_manifest_against_snapshot(
                manifest_path,
                [first, second],
                snapshot=snapshot,
                validate_narratives=False,
            )
        ),
    )

    frozen_input.write_bytes(b"changed input")
    check(
        "frozen input modification during generation is rejected",
        raises_runtime(lambda: assert_visual_inputs_unchanged(snapshot)),
    )

print("ALL PARTICIPANT-VISUAL PROVENANCE CHECKS PASSED")
