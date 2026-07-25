"""Shared versioning, provenance, and atomic-output helpers for the LC-proxy pipeline."""
from __future__ import annotations

import datetime as _dt
import hashlib
import importlib.metadata as _metadata
import json
import os
import platform
import subprocess
import tempfile
import uuid

import numpy as np
import scipy


# Increment both values whenever an estimator or a cached derived signal changes materially.
ANALYSIS_VERSION = "2026-07-qc-sensitivity-v8"
CACHE_SCHEMA_VERSION = (
    "2026-07-neutral-per-contact-gap-aware-source-pin-v8")

# Only files that can change cache values or cache metadata belong here.  The broader
# ``source_tree_sha256`` intentionally includes downstream analyses and tests, which would make a
# harmless summary-script edit invalidate many hours of derived data.
_CACHE_SOURCE_FILES = (
    "analysis/cache_lc_series.py",
    "analysis/stage_ds003848.py",
    "analysis/ds003848_snapshot_1.0.3_files.json",
    "analysis/hup_ieeg_source_pin.json",
    "analysis/cohort_3A_cortical.py",
    "analysis/cohort_stages_3ABD.py",
    "analysis/infraslow_rr_sigma_coherence.py",
    "analysis/results_3A_tutorial_style.py",
    "analysis/spectral_gapped.py",
    "analysis/pipeline_version.py",
    "env/requirements.txt",
    "env/PYTHON_VERSION",
)


def utc_now():
    return _dt.datetime.now(_dt.timezone.utc).isoformat()


def git_revision(root):
    """Return the checked-out revision without making a clean worktree claim."""
    try:
        return subprocess.check_output(
            ["git", "-C", root, "rev-parse", "HEAD"], text=True,
            stderr=subprocess.DEVNULL, timeout=10).strip()
    except Exception:
        return "unknown"


def git_is_dirty(root):
    """Whether tracked or untracked files differ from the recorded revision."""
    try:
        return bool(subprocess.check_output(
            ["git", "-C", root, "status", "--porcelain"], text=True,
            stderr=subprocess.DEVNULL, timeout=10).strip())
    except Exception:
        return None


def source_tree_sha256(root):
    """Digest the executable analysis/config source, including untracked files.

    A Git revision alone is ambiguous when a pipeline is run from a dirty working tree. The digest
    covers the bytes that define this analysis without hashing data or generated outputs.
    """
    roots = [
        os.path.join(root, "analysis"),
        os.path.join(root, "env"),
        os.path.join(root, ".github", "workflows"),
    ]
    files = []
    for directory in roots:
        if not os.path.isdir(directory):
            continue
        for current, dirs, names in os.walk(directory):
            dirs[:] = sorted(
                value for value in dirs
                if value not in ("__pycache__", ".pytest_cache", ".mypy_cache"))
            for name in sorted(names):
                if name.endswith((".py", ".json", ".txt", ".yml", ".yaml")) or name == "PYTHON_VERSION":
                    files.append(os.path.join(current, name))
    digest = hashlib.sha256()
    for path in sorted(files):
        rel = os.path.relpath(path, root).replace(os.sep, "/")
        digest.update(rel.encode("utf-8"))
        digest.update(b"\0")
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        digest.update(b"\0")
    return digest.hexdigest()


def cache_code_sha256(root):
    """Digest every checked-in or untracked source file that can alter a cache.

    Schema labels alone are forgeable: after a code edit, an old NPZ could otherwise be accepted
    and advertised under a newly generated run manifest.  This digest ties each cache and its
    manifest to the exact cache-producing source bytes.
    """
    digest = hashlib.sha256()
    for rel in _CACHE_SOURCE_FILES:
        path = os.path.join(root, rel)
        if not os.path.isfile(path):
            raise FileNotFoundError(f"cache source is missing: {path}")
        digest.update(rel.encode("utf-8"))
        digest.update(b"\0")
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        digest.update(b"\0")
    return digest.hexdigest()


def file_sha256(path):
    """SHA-256 of one immutable input/output artifact."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def runtime_versions():
    versions = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    try:
        import neurokit2
        versions["neurokit2"] = neurokit2.__version__
    except Exception:
        versions["neurokit2"] = "unavailable"
    for key, distribution in (
        ("mne", "mne"),
        ("scikit_learn", "scikit-learn"),
        ("ieeg", "ieeg"),
    ):
        try:
            versions[key] = _metadata.version(distribution)
        except _metadata.PackageNotFoundError:
            versions[key] = "unavailable"
    return versions


def npz_scalar_text(d, key, default=""):
    """Read a scalar string from an npz payload without relying on array repr."""
    if key not in getattr(d, "files", []):
        return default
    value = np.asarray(d[key])
    return str(value.item()) if value.shape == () else str(value.ravel()[0])


def finite_float_or_none(value):
    """Return a JSON-safe finite float, otherwise an explicit unavailable value."""
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return value if np.isfinite(value) else None


def cache_is_current(path, root=None):
    if not os.path.exists(path):
        return False
    try:
        with np.load(path, allow_pickle=False) as d:
            schema_ok = npz_scalar_text(d, "cache_schema_version") == CACHE_SCHEMA_VERSION
            digest = npz_scalar_text(d, "cache_code_sha256")
            digest_ok = root is None or digest == cache_code_sha256(root)
            return npz_scalar_text(d, "status") == "ok" and schema_ok and digest_ok
    except Exception:
        return False


def atomic_json_dump(payload, path, *, indent=2):
    """Write standards-compliant JSON completely before replacing the destination.

    Python's JSON encoder otherwise writes bare ``NaN``/``Infinity`` tokens.  Those are not valid
    JSON and can silently propagate a failed numerical endpoint through Python-only readers.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".tmp-", suffix=".json", dir=os.path.dirname(path))
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(payload, f, indent=indent, sort_keys=True, allow_nan=False)
            f.write("\n")
        os.replace(tmp, path)
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def atomic_savez(path, **payload):
    """Write a compressed npz completely before replacing the destination."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".tmp-", suffix=".npz", dir=os.path.dirname(path))
    os.close(fd)
    try:
        np.savez_compressed(tmp, **payload)
        os.replace(tmp, path)
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def write_run_manifest(out_dir, *, pipeline, requested, completed, skipped, failed, config,
                       run_id=None, run_state="complete", result_files_sha256=None):
    requested = list(requested)
    completed = list(completed)
    skipped = list(skipped)
    failed = list(failed)
    named_lists = {
        "requested": requested,
        "completed": completed,
        "skipped": [
            value.get("subject") if isinstance(value, dict) else value
            for value in skipped
        ],
        "failed": [
            value.get("subject") if isinstance(value, dict) else value
            for value in failed
        ],
    }
    for label, values in named_lists.items():
        if any(value in (None, "") for value in values):
            raise ValueError(f"manifest {label} entries require a subject identifier")
        if len(values) != len(set(values)):
            raise ValueError(f"manifest {label} contains duplicate subject identifiers")
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    manifest = {
        "run_id": run_id or str(uuid.uuid4()),
        "run_state": run_state,
        "pipeline": pipeline,
        "analysis_version": ANALYSIS_VERSION,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "generated_at_utc": utc_now(),
        "code_revision": git_revision(root),
        "code_dirty": git_is_dirty(root),
        "source_tree_sha256": source_tree_sha256(root),
        "runtime_versions": runtime_versions(),
        "requested": requested,
        "completed": completed,
        "skipped": skipped,
        "failed": failed,
        "config": config,
        "result_files_sha256": dict(result_files_sha256 or {}),
    }
    atomic_json_dump(manifest, os.path.join(out_dir, "RUN_MANIFEST.json"))
    return manifest


def start_run_manifest(out_dir, *, pipeline, requested, config):
    """Invalidate any old complete manifest before the first result is touched."""
    run_id = str(uuid.uuid4())
    write_run_manifest(
        out_dir, pipeline=pipeline, requested=requested, completed=[], skipped=[], failed=[],
        config=config, run_id=run_id, run_state="in_progress")
    return run_id


def validated_complete_run_exists(
        out_dir, *, pipeline, requested, config, suffix, require_current_source_tree):
    """Return true only when an existing complete run is byte-for-byte reusable.

    Reuse must be decided *before* ``start_run_manifest`` replaces the previous terminal manifest.
    Metadata-only reuse is unsafe: a modified NPZ/JSON with unchanged embedded labels could
    otherwise be hashed and blessed by a fresh manifest.  If any requested output already exists,
    this function therefore requires an exact, complete prior run and verifies its stored hashes.
    Partial or unverifiable output sets fail closed and require an explicit ``--force`` rebuild.
    """
    requested = list(requested)
    if len(requested) != len(set(requested)):
        raise RuntimeError(
            "requested subjects contain duplicates; provide each subject exactly once")
    paths = {
        subject: os.path.join(out_dir, f"{subject}{suffix}")
        for subject in requested
    }
    existing = {subject for subject, path in paths.items() if os.path.exists(path)}
    if not existing:
        return False
    if existing != set(requested):
        raise RuntimeError(
            "only part of the requested output set already exists; rerun with --force")

    manifest_path = os.path.join(out_dir, "RUN_MANIFEST.json")
    if not os.path.exists(manifest_path):
        raise RuntimeError(
            "existing outputs have no prior terminal manifest; rerun with --force")
    try:
        with open(manifest_path) as handle:
            manifest = json.load(handle)
    except Exception as exc:
        raise RuntimeError(
            "the prior run manifest cannot be read; rerun with --force") from exc

    requested_manifest_values = manifest.get("requested", [])
    completed_manifest_values = manifest.get("completed", [])
    skipped_entries = manifest.get("skipped", [])
    if any(not isinstance(value, dict) or not value.get("subject") or not value.get("reason")
           for value in skipped_entries):
        raise RuntimeError(
            "the prior run has malformed skip records; rerun with --force")
    skipped_values = [value["subject"] for value in skipped_entries]
    skipped = set(skipped_values)
    completed = set(completed_manifest_values)
    expected = set(requested)
    if (
        manifest.get("run_state") != "complete"
        or manifest.get("pipeline") != pipeline
        or manifest.get("analysis_version") != ANALYSIS_VERSION
        or manifest.get("cache_schema_version") != CACHE_SCHEMA_VERSION
        or manifest.get("runtime_versions") != runtime_versions()
        or manifest.get("failed")
        or len(requested_manifest_values) != len(set(requested_manifest_values))
        or len(completed_manifest_values) != len(completed)
        or len(skipped_values) != len(skipped)
        or set(requested_manifest_values) != expected
        or completed | skipped != expected
        or completed & skipped
        or manifest.get("config") != config
    ):
        raise RuntimeError(
            "existing outputs do not match an exact complete current run; rerun with --force")

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if (require_current_source_tree
            and manifest.get("source_tree_sha256") != source_tree_sha256(root)):
        raise RuntimeError(
            "executable source changed since the prior run; rerun with --force")

    recorded_hashes = manifest.get("result_files_sha256")
    if not isinstance(recorded_hashes, dict) or set(recorded_hashes) != expected:
        raise RuntimeError(
            "the prior run lacks exact per-result hashes; rerun with --force")
    changed = [
        subject for subject, path in paths.items()
        if recorded_hashes.get(subject) != file_sha256(path)
    ]
    if changed:
        raise RuntimeError(
            f"existing output bytes differ from the prior manifest for {changed}; "
            "rerun with --force")
    return True
