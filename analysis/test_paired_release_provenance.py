"""Regressions for clean, revision-bound paired public releases."""
from __future__ import annotations

import copy
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from unittest import mock

import audit_hup_scalp_inventory as inventory_producer
import paired_scalp_ieeg_comparison as paired_producer
from pipeline_version import source_tree_sha256
from release_provenance import (
    require_clean_release_provenance,
    require_release_source_unchanged,
    source_tree_sha256_at_revision,
    validate_recorded_release_provenance,
)


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


def rejected(call, text):
    try:
        call()
    except RuntimeError as error:
        return text in str(error)
    return False


def git(root, *args):
    return subprocess.check_output(
        ["git", "-C", str(root), *args],
        text=True,
        stderr=subprocess.DEVNULL,
    ).strip()


with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    (root / "analysis").mkdir()
    (root / "env").mkdir()
    (root / "outputs").mkdir()
    producer_path = root / "analysis" / "producer.py"
    python_path = root / "env" / "PYTHON_VERSION"
    producer_path.write_text("VALUE = 1\n", encoding="utf-8")
    python_path.write_text("3.11\n", encoding="utf-8")
    subprocess.run(["git", "-C", str(root), "init", "-q"], check=True)
    subprocess.run(
        ["git", "-C", str(root), "config", "user.name", "Contract Test"],
        check=True,
    )
    subprocess.run(
        [
            "git", "-C", str(root), "config", "user.email",
            "contract@example.invalid",
        ],
        check=True,
    )
    subprocess.run(["git", "-C", str(root), "add", "."], check=True)
    subprocess.run(
        ["git", "-C", str(root), "commit", "-qm", "source one"],
        check=True,
    )
    revision_one = git(root, "rev-parse", "HEAD")
    provenance_one = require_clean_release_provenance(root)
    check(
        "clean release provenance is bound to the exact committed bytes",
        provenance_one["code_revision"] == revision_one
        and provenance_one["code_dirty"] is False
        and provenance_one["source_tree_sha256"]
        == source_tree_sha256_at_revision(root, revision_one)
        == source_tree_sha256(root),
    )
    check(
        "a valid recorded release is portable across runtime environments",
        validate_recorded_release_provenance(
            provenance_one,
            root,
            label="synthetic release",
        ) == provenance_one,
    )

    dirty_manifest = copy.deepcopy(provenance_one)
    dirty_manifest["code_dirty"] = True
    check(
        "a public manifest that records dirty generation fails closed",
        rejected(
            lambda: validate_recorded_release_provenance(
                dirty_manifest, root, label="dirty manifest"),
            "records a dirty source tree",
        ),
    )

    producer_path.write_text("VALUE = 2\n", encoding="utf-8")
    missing = root / "missing.json"
    with mock.patch.object(inventory_producer, "ROOT", str(root)):
        inventory_dirty_rejected = rejected(
            lambda: inventory_producer.build_inventory(
                cache_dir=str(missing),
                sidecar_dir=str(root / "missing-sidecar"),
                qc_grid_path=str(missing),
                source_pin_path=str(missing),
            ),
            "requires a clean Git worktree",
        )
    check(
        "inventory generation rejects dirty source before reading inputs",
        inventory_dirty_rejected,
    )
    paired_output = root / "paired-output"
    paired_argv = [
        "paired_scalp_ieeg_comparison.py",
        "--ieeg-cache", str(root / "missing-ieeg"),
        "--scalp-cache", str(root / "missing-scalp"),
        "--scalp-inventory", str(missing),
        "--qc-grid", str(missing),
        "--output-dir", str(paired_output),
    ]
    with (
        mock.patch.object(paired_producer, "ROOT", str(root)),
        mock.patch.object(sys, "argv", paired_argv),
    ):
        paired_dirty_rejected = rejected(
            paired_producer.main,
            "requires a clean Git worktree",
        )
    check(
        "paired generation rejects dirty source before creating outputs",
        paired_dirty_rejected and not paired_output.exists(),
    )
    producer_path.write_text("VALUE = 1\n", encoding="utf-8")
    check(
        "the end-of-run check accepts source restored to committed bytes",
        require_release_source_unchanged(root, provenance_one)
        == provenance_one,
    )

    producer_path.write_text("VALUE = 2\n", encoding="utf-8")
    subprocess.run(
        ["git", "-C", str(root), "add", "analysis/producer.py"],
        check=True,
    )
    subprocess.run(
        ["git", "-C", str(root), "commit", "-qm", "source two"],
        check=True,
    )
    revision_two = git(root, "rev-parse", "HEAD")
    provenance_two = require_clean_release_provenance(root)
    forged_revision_pin = {
        **provenance_two,
        "code_revision": revision_one,
    }
    check(
        "a real revision with different source bytes cannot be re-blessed",
        revision_one != revision_two
        and source_tree_sha256_at_revision(root, revision_one)
        != source_tree_sha256_at_revision(root, revision_two)
        and rejected(
            lambda: validate_recorded_release_provenance(
                forged_revision_pin,
                root,
                label="forged historical pin",
            ),
            "source bytes differ from its recorded Git revision",
        ),
    )

    output_path = root / "outputs" / "result.json"
    output_path.write_text("{}\n", encoding="utf-8")
    check(
        "end-of-run source checks tolerate only generated-output dirtiness",
        require_release_source_unchanged(root, provenance_two)
        == provenance_two,
    )
    producer_path.write_text("VALUE = 3\n", encoding="utf-8")
    check(
        "end-of-run checks still reject source changes after output writes",
        rejected(
            lambda: require_release_source_unchanged(
                root, provenance_two),
            "source bytes changed",
        ),
    )

print("ALL PAIRED RELEASE-PROVENANCE CHECKS PASSED")
