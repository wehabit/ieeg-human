"""Clean, revision-bound provenance for downstream public artifacts.

This module is intentionally outside the neutral-cache producer dependency
digest. Its functions govern inventory/result publication and cannot alter raw
cache values or cache metadata. Keeping that boundary explicit prevents a
downstream release-contract edit from invalidating hundreds of megabytes of
otherwise byte-current neutral caches.
"""
from __future__ import annotations

import functools
import hashlib
import os
import re
import subprocess

from pipeline_version import (
    git_is_dirty,
    git_revision,
    source_tree_sha256,
)


_SOURCE_TREE_ROOTS = (
    "analysis",
    "env",
    ".github/workflows",
)
_SOURCE_TREE_SUFFIXES = (".py", ".json", ".txt", ".yml", ".yaml")
_SOURCE_TREE_EXCLUDED_DIRECTORIES = frozenset({
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
})
_FULL_GIT_REVISION = re.compile(r"^[0-9a-f]{40}$")


def _source_tree_relative_path(relative):
    """Whether one repository-relative path belongs to the source digest."""
    relative = str(relative).replace(os.sep, "/")
    parts = relative.split("/")
    if any(value in _SOURCE_TREE_EXCLUDED_DIRECTORIES for value in parts):
        return False
    if not any(
        relative == root or relative.startswith(root + "/")
        for root in _SOURCE_TREE_ROOTS
    ):
        return False
    name = parts[-1]
    return name == "PYTHON_VERSION" or name.endswith(_SOURCE_TREE_SUFFIXES)


@functools.lru_cache(maxsize=32)
def source_tree_sha256_at_revision(root, revision):
    """Digest the source-tree bytes stored by one exact Git commit.

    ``None`` is returned for an abbreviated/non-hex revision, a missing commit,
    or a source file that cannot be read from that commit. Requiring a full
    object name prevents a manifest from relying on a mutable symbolic ref such
    as ``HEAD`` or a branch name.
    """
    if (
        not isinstance(revision, str)
        or _FULL_GIT_REVISION.fullmatch(revision) is None
    ):
        return None
    try:
        exists = subprocess.run(
            [
                "git", "-C", root, "cat-file", "-e",
                f"{revision}^{{commit}}",
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
            timeout=10,
        )
        if exists.returncode != 0:
            return None
        listed = subprocess.run(
            [
                "git", "-C", root, "ls-tree", "-r", "-z", "--name-only",
                revision, "--", *_SOURCE_TREE_ROOTS,
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    if listed.returncode != 0:
        return None
    try:
        relative_paths = sorted(
            value.decode("utf-8")
            for value in listed.stdout.split(b"\0")
            if value and _source_tree_relative_path(value.decode("utf-8"))
        )
    except UnicodeDecodeError:
        return None
    digest = hashlib.sha256()
    for relative in relative_paths:
        try:
            completed = subprocess.run(
                ["git", "-C", root, "show", f"{revision}:{relative}"],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                check=False,
                timeout=30,
            )
        except (OSError, subprocess.TimeoutExpired):
            return None
        if completed.returncode != 0:
            return None
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(completed.stdout)
        digest.update(b"\0")
    return digest.hexdigest()


def require_clean_release_provenance(root):
    """Return release provenance only for a stable, clean Git source tree."""
    revision = git_revision(root)
    dirty = git_is_dirty(root)
    if dirty is not False:
        raise RuntimeError(
            "public artifact generation requires a clean Git worktree")
    current_digest = source_tree_sha256(root)
    revision_digest = source_tree_sha256_at_revision(root, revision)
    if revision_digest is None:
        raise RuntimeError(
            "public artifact generation requires one exact recorded Git commit")
    if revision_digest != current_digest:
        raise RuntimeError(
            "current source bytes are not contained in the recorded Git "
            "revision")
    if git_revision(root) != revision or git_is_dirty(root) is not False:
        raise RuntimeError(
            "Git source state changed while release provenance was captured")
    return {
        "code_revision": revision,
        "code_dirty": False,
        "source_tree_sha256": current_digest,
    }


def require_release_source_unchanged(root, expected):
    """Fail if source bytes or the checked-out revision changed during a run."""
    required = {"code_revision", "code_dirty", "source_tree_sha256"}
    if not isinstance(expected, dict) or set(expected) != required:
        raise RuntimeError("release provenance has a nonexact schema")
    if expected["code_dirty"] is not False:
        raise RuntimeError("release provenance records a dirty source tree")
    revision = git_revision(root)
    current_digest = source_tree_sha256(root)
    if (
        revision != expected["code_revision"]
        or current_digest != expected["source_tree_sha256"]
        or source_tree_sha256_at_revision(root, revision) != current_digest
    ):
        raise RuntimeError(
            "release source bytes changed while artifacts were generated")
    return dict(expected)


def validate_recorded_release_provenance(
        payload, root, *, label, require_current_source=True):
    """Validate a public artifact's clean, immutable Git source claim."""
    if not isinstance(payload, dict):
        raise RuntimeError(f"{label} must be an object")
    revision = payload.get("code_revision")
    digest = payload.get("source_tree_sha256")
    if payload.get("code_dirty") is not False:
        raise RuntimeError(f"{label} records a dirty source tree")
    if (
        not isinstance(digest, str)
        or re.fullmatch(r"[0-9a-f]{64}", digest) is None
    ):
        raise RuntimeError(f"{label} has an invalid source_tree_sha256")
    recorded_digest = source_tree_sha256_at_revision(root, revision)
    if recorded_digest is None:
        raise RuntimeError(
            f"{label} code_revision is not one exact available Git commit")
    if recorded_digest != digest:
        raise RuntimeError(
            f"{label} source bytes differ from its recorded Git revision")
    if require_current_source and source_tree_sha256(root) != digest:
        raise RuntimeError(
            f"{label} was not produced by the current analysis source tree")
    return {
        "code_revision": revision,
        "code_dirty": False,
        "source_tree_sha256": digest,
    }
