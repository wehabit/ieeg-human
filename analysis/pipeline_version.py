"""Shared versioning, provenance, and atomic-output helpers for the LC-proxy pipeline."""
from __future__ import annotations

import datetime as _dt
import functools
import hashlib
import importlib.metadata as _metadata
import json
import os
import platform
import re
import subprocess
import tempfile
import uuid

import numpy as np
import scipy


# Increment both values whenever an estimator or a cached derived signal changes materially.
ANALYSIS_VERSION = "2026-07-qc-sensitivity-v9"
CACHE_SCHEMA_VERSION = (
    "2026-07-neutral-per-contact-gap-aware-source-pin-v9")

# Only files that can change cache values or cache metadata belong here.  The broader
# ``source_tree_sha256`` intentionally includes downstream analyses and tests, which would make a
# harmless summary-script edit invalidate many hours of derived data.
_CACHE_SOURCE_FILES = (
    "analysis/cache_lc_series.py",
    "analysis/event_3d_estimators.py",
    "analysis/stage_ds003848.py",
    "analysis/ds003848_snapshot_1.0.3_files.json",
    "analysis/hup_ieeg_source_pin.json",
    "analysis/hup_portal.py",
    "analysis/cohort_3A_cortical.py",
    "analysis/infraslow_rr_sigma_coherence.py",
    "analysis/signal_qc.py",
    "analysis/spectral_gapped.py",
    "analysis/staging_helpers.py",
    "analysis/pipeline_version.py",
    "env/requirements.txt",
    "env/PYTHON_VERSION",
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
    files = []
    for relative_root in _SOURCE_TREE_ROOTS:
        directory = os.path.join(root, *relative_root.split("/"))
        if not os.path.isdir(directory):
            continue
        for current, dirs, names in os.walk(directory):
            dirs[:] = sorted(
                value for value in dirs
                if value not in _SOURCE_TREE_EXCLUDED_DIRECTORIES)
            for name in sorted(names):
                if (
                    name.endswith(_SOURCE_TREE_SUFFIXES)
                    or name == "PYTHON_VERSION"
                ):
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
    # Detect a checkout/edit race while the digests were being computed.
    if git_revision(root) != revision or git_is_dirty(root) is not False:
        raise RuntimeError(
            "Git source state changed while release provenance was captured")
    return {
        "code_revision": revision,
        "code_dirty": False,
        "source_tree_sha256": current_digest,
    }


def require_release_source_unchanged(root, expected):
    """Fail if source bytes or the checked-out revision changed during a run.

    Generated output files may make the worktree dirty after the first result
    write, so this end-of-run check intentionally compares only revision and
    source bytes. The initial release check remains strictly clean.
    """
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
    """Validate a public artifact's clean, immutable Git source claim.

    Offline validation is runtime-portable: it reads only the repository and
    the recorded commit. ``require_current_source`` optionally also requires
    the checkout's source bytes to equal the producing bytes.
    """
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


def require_integer_sample_rate(value, *, source):
    """Return a positive integer sampling rate or fail before sample binning.

    The cache format stores one-second power bins whose current implementation
    requires an integral number of source samples per second.  Silently
    truncating a fractional rate with ``int(sf)`` would shift those bins.
    """
    try:
        sample_rate = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise RuntimeError(
            f"{source} sampling rate is not a finite positive integer") from exc
    if not np.isfinite(sample_rate) or sample_rate <= 0:
        raise RuntimeError(
            f"{source} sampling rate is not a finite positive integer")
    integer_rate = int(round(sample_rate))
    tolerance = max(
        1e-9,
        64.0 * np.finfo(float).eps * max(1.0, abs(sample_rate)),
    )
    if abs(sample_rate - integer_rate) > tolerance:
        raise RuntimeError(
            f"{source} sampling rate {sample_rate!r} is fractional; current "
            "one-second cache bins require an integer rate")
    return integer_rate


def cache_is_current(path, root):
    """Return whether an OK cache matches the current schema and producer bytes.

    ``root`` is deliberately required.  Schema equality alone is not sufficient
    provenance because scientific code can change without changing a label.
    """
    if not os.path.exists(path):
        return False
    try:
        with np.load(path, allow_pickle=False) as d:
            schema_ok = npz_scalar_text(d, "cache_schema_version") == CACHE_SCHEMA_VERSION
            digest = npz_scalar_text(d, "cache_code_sha256")
            digest_ok = digest == cache_code_sha256(root)
            return npz_scalar_text(d, "status") == "ok" and schema_ok and digest_ok
    except Exception:
        return False


def json_list_field(d, key, *, source, required=True):
    """Decode a required JSON-list field from an NPZ cache.

    Cache failure metadata is part of the scientific completion contract, not
    an optional diagnostic.  Missing, malformed, or non-list values therefore
    fail closed for current-schema publication consumers.
    """
    if key not in getattr(d, "files", []):
        if required:
            raise RuntimeError(f"{source} lacks required cache field {key!r}")
        return []
    try:
        value = json.loads(npz_scalar_text(d, key))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"{source} has malformed JSON in {key!r}") from exc
    if not isinstance(value, list):
        raise RuntimeError(f"{source} field {key!r} must contain a JSON list")
    return value


def validate_completed_cache_failures(d, *, source, require_ecg):
    """Require empty fatal-failure arrays in a completed scientific cache."""
    failed_chunks = json_list_field(
        d, "failed_chunks_json", source=source, required=True)
    if failed_chunks:
        raise RuntimeError(
            f"{source} contains {len(failed_chunks)} failed acquisition chunks")
    ecg_failures = json_list_field(
        d, "ecg_failures_json", source=source, required=require_ecg)
    if ecg_failures:
        raise RuntimeError(
            f"{source} contains {len(ecg_failures)} failed ECG-detector chunks")
    return {
        "failed_chunks": failed_chunks,
        "ecg_failures": ecg_failures,
    }


def _manifest_record_subjects(records, *, label, source, detail_field):
    if not isinstance(records, list):
        raise RuntimeError(f"{source} manifest field {label!r} must be a list")
    subjects = []
    for index, record in enumerate(records):
        if (
            not isinstance(record, dict)
            or not isinstance(record.get("subject"), str)
            or not record["subject"]
            or not isinstance(record.get(detail_field), str)
            or not record[detail_field]
        ):
            raise RuntimeError(
                f"{source} manifest {label}[{index}] must contain nonempty "
                f"subject/{detail_field} strings")
        subjects.append(record["subject"])
    if len(subjects) != len(set(subjects)):
        raise RuntimeError(
            f"{source} manifest {label} contains duplicate subject identifiers")
    return subjects


def validate_terminal_run_manifest(manifest, *, source):
    """Validate the standard cache-producer terminal partition and hash map."""
    if not isinstance(manifest, dict):
        raise RuntimeError(f"{source} manifest must be a JSON object")
    state = manifest.get("run_state")
    if state not in {"complete", "failed"}:
        raise RuntimeError(f"{source} manifest is not terminal")
    requested = manifest.get("requested")
    completed = manifest.get("completed")
    if (
        not isinstance(requested, list)
        or not all(isinstance(value, str) and value for value in requested)
        or len(requested) != len(set(requested))
        or not isinstance(completed, list)
        or not all(isinstance(value, str) and value for value in completed)
        or len(completed) != len(set(completed))
    ):
        raise RuntimeError(
            f"{source} manifest requested/completed subjects are malformed")
    skipped = _manifest_record_subjects(
        manifest.get("skipped"), label="skipped", source=source,
        detail_field="reason")
    failed = _manifest_record_subjects(
        manifest.get("failed"), label="failed", source=source,
        detail_field="error")
    partitions = (set(completed), set(skipped), set(failed))
    if (
        any(
            partitions[left] & partitions[right]
            for left, right in ((0, 1), (0, 2), (1, 2))
        )
        or set.union(*partitions) != set(requested)
        or (state == "complete" and failed)
        or (state == "failed" and not failed)
    ):
        raise RuntimeError(
            f"{source} manifest terminal subject partition is inconsistent")
    hashes = manifest.get("result_files_sha256")
    expected_hashes = set(completed) | set(skipped)
    if (
        not isinstance(hashes, dict)
        or set(hashes) != expected_hashes
        or any(
            not isinstance(value, str)
            or len(value) != 64
            or any(character not in "0123456789abcdef" for character in value)
            for value in hashes.values()
        )
    ):
        raise RuntimeError(
            f"{source} manifest result hashes do not match its terminal outputs")
    return {
        "run_state": state,
        "requested": list(requested),
        "completed": list(completed),
        "skipped": skipped,
        "failed": failed,
        "result_files_sha256": dict(hashes),
    }


def _npz_scalar_float(d, key, *, source):
    if key not in getattr(d, "files", []):
        raise RuntimeError(f"{source} lacks required cache field {key!r}")
    try:
        value = float(np.asarray(d[key]).item())
    except (TypeError, ValueError, OverflowError) as exc:
        raise RuntimeError(
            f"{source} cache field {key!r} is not a finite scalar") from exc
    if not np.isfinite(value):
        raise RuntimeError(
            f"{source} cache field {key!r} is not a finite scalar")
    return value


def validate_full_interval_acquisition(
        d, *, source, core_purpose, subrequest_purpose, core_chunk_s,
        filter_edge_s):
    """Prove portal responses match the producer's exact padded chunk geometry.

    Records must alternate one or more bounded subrequests with the
    non-overlapping analysis core they cover.  Every core must be the next
    ``core_chunk_s`` block (or the final remainder), and its subrequests must
    exactly tile the sample-aligned, ``filter_edge_s``-padded pull clipped to
    the cached interval.
    """
    records = json_list_field(
        d, "acquisition_sample_counts_json", source=source, required=True)
    hours = _npz_scalar_float(d, "hours", source=source)
    night_s = _npz_scalar_float(d, "night_s", source=source)
    sample_rate = require_integer_sample_rate(
        _npz_scalar_float(d, "sf", source=source), source=source)
    try:
        core_chunk_s = float(core_chunk_s)
        filter_edge_s = float(filter_edge_s)
    except (TypeError, ValueError, OverflowError) as exc:
        raise RuntimeError(
            f"{source} portal chunk/filter geometry is invalid") from exc
    if (
        not np.isfinite([core_chunk_s, filter_edge_s]).all()
        or core_chunk_s <= 0
        or filter_edge_s < 0
        or night_s < 0
    ):
        raise RuntimeError(
            f"{source} portal chunk/filter geometry is invalid")

    total_duration_s = hours * 3600.0
    rounded_duration_s = round(total_duration_s)
    tolerance_s = 1e-7
    if (
        hours <= 0
        or abs(total_duration_s - rounded_duration_s) > tolerance_s
    ):
        raise RuntimeError(
            f"{source} cache duration must be a positive whole number of seconds")
    total_duration_s = float(rounded_duration_s)

    def sample_index(value, *, label):
        scaled = float(value) * sample_rate
        if not np.isfinite(scaled):
            raise RuntimeError(
                f"{source} {label} is not sample-aligned")
        rounded = int(round(scaled))
        tolerance_samples = max(
            tolerance_s * sample_rate,
            64.0 * np.finfo(float).eps * max(1.0, abs(scaled)),
        )
        if abs(scaled - rounded) > tolerance_samples:
            raise RuntimeError(
                f"{source} {label} is not sample-aligned")
        return rounded

    total_samples = sample_index(
        total_duration_s, label="cache duration")
    night_sample = sample_index(night_s, label="night start")
    core_chunk_samples = sample_index(
        core_chunk_s, label="core chunk duration")
    filter_edge_samples = sample_index(
        filter_edge_s, label="filter edge duration")
    if core_chunk_samples <= 0 or filter_edge_samples < 0:
        raise RuntimeError(
            f"{source} portal chunk/filter geometry is invalid")

    core_cursor_sample = 0
    core_count = 0
    subrequest_count = 0
    expected_channel_count = None
    pending_subrequests = []
    for index, record in enumerate(records):
        if not isinstance(record, dict):
            raise RuntimeError(
                f"{source} acquisition record {index} is not an object")
        purpose = record.get("purpose")
        requested = record.get("requested_sample_count")
        returned = record.get("returned_sample_count")
        if (
            record.get("status") != "ok"
            or not isinstance(requested, int)
            or isinstance(requested, bool)
            or requested <= 0
            or not isinstance(returned, int)
            or isinstance(returned, bool)
            or returned != requested
        ):
            raise RuntimeError(
                f"{source} acquisition record {index} is not an exact "
                "successful request/core")

        if purpose == subrequest_purpose:
            try:
                request_start = float(record["request_start_s"])
                request_duration = float(record["request_duration_s"])
            except (KeyError, TypeError, ValueError, OverflowError) as exc:
                raise RuntimeError(
                    f"{source} subrequest record {index} has invalid geometry") from exc
            channel_count = record.get("requested_channel_count")
            returned_channel_count = record.get("returned_channel_count")
            request_start_sample = sample_index(
                request_start, label=f"subrequest {index} start")
            request_sample_count = sample_index(
                request_duration, label=f"subrequest {index} duration")
            if (
                not np.isfinite([request_start, request_duration]).all()
                or request_start < 0
                or request_duration <= 0
                or request_sample_count <= 0
                or requested != request_sample_count
                or not isinstance(channel_count, int)
                or isinstance(channel_count, bool)
                or channel_count <= 0
                or not isinstance(returned_channel_count, int)
                or isinstance(returned_channel_count, bool)
                or returned_channel_count != channel_count
            ):
                raise RuntimeError(
                    f"{source} subrequest record {index} has inconsistent geometry")
            if expected_channel_count is None:
                expected_channel_count = channel_count
            elif channel_count != expected_channel_count:
                raise RuntimeError(
                    f"{source} portal subrequests change channel count")
            if pending_subrequests:
                previous = pending_subrequests[-1]
                previous_stop_sample = (
                    previous["request_start_sample"]
                    + previous["request_sample_count"]
                )
                if request_start_sample != previous_stop_sample:
                    raise RuntimeError(
                        f"{source} bounded subrequests are not contiguous")
            pending_subrequests.append({
                "request_start_s": request_start,
                "request_duration_s": request_duration,
                "request_start_sample": request_start_sample,
                "request_sample_count": request_sample_count,
            })
            subrequest_count += 1
            continue

        if purpose != core_purpose:
            raise RuntimeError(
                f"{source} acquisition record {index} has unknown purpose "
                f"{purpose!r}")
        if not pending_subrequests:
            raise RuntimeError(
                f"{source} analysis core {core_count} has no proven subrequest")
        try:
            core_start = float(record["analysis_start_s"])
            core_duration = float(record["analysis_duration_s"])
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            raise RuntimeError(
                f"{source} analysis core {core_count} has invalid geometry") from exc
        core_start_sample = sample_index(
            core_start, label=f"analysis core {core_count} start")
        core_duration_sample_count = sample_index(
            core_duration, label=f"analysis core {core_count} duration")
        expected_core_sample_count = min(
            core_chunk_samples, total_samples - core_cursor_sample)
        if (
            not np.isfinite([core_start, core_duration]).all()
            or core_duration <= 0
            or core_start_sample != core_cursor_sample
            or core_duration_sample_count != expected_core_sample_count
            or requested != expected_core_sample_count
        ):
            raise RuntimeError(
                f"{source} analysis cores do not match the exact ordered "
                "producer chunk partition")

        expected_pull_start_sample = (
            night_sample
            + max(0, core_cursor_sample - filter_edge_samples)
        )
        expected_pull_stop_sample = (
            night_sample
            + min(
                total_samples,
                core_cursor_sample
                + expected_core_sample_count
                + filter_edge_samples,
            )
        )
        pull_cursor_sample = expected_pull_start_sample
        for subrequest in pending_subrequests:
            if subrequest["request_start_sample"] != pull_cursor_sample:
                raise RuntimeError(
                    f"{source} analysis core {core_count} subrequests do not "
                    "start at and exactly tile its intended padded pull")
            pull_cursor_sample += subrequest["request_sample_count"]
            if pull_cursor_sample > expected_pull_stop_sample:
                raise RuntimeError(
                    f"{source} analysis core {core_count} subrequests extend "
                    "beyond its intended padded pull")
        if pull_cursor_sample != expected_pull_stop_sample:
            raise RuntimeError(
                f"{source} analysis core {core_count} subrequests do not "
                "exactly tile its intended padded pull")

        core_cursor_sample += expected_core_sample_count
        core_count += 1
        pending_subrequests = []

    if pending_subrequests:
        raise RuntimeError(
            f"{source} has bounded subrequests without a matching analysis core")
    if core_count == 0 or core_cursor_sample != total_samples:
        raise RuntimeError(
            f"{source} acquisition cores do not cover the complete cached interval")
    return {
        "records": records,
        "core_count": core_count,
        "subrequest_count": subrequest_count,
        "duration_s": total_duration_s,
        "sample_rate_hz": sample_rate,
        "channel_count": expected_channel_count,
        "core_chunk_s": core_chunk_s,
        "filter_edge_s": filter_edge_s,
    }


def validate_local_chunk_acquisition(d, *, source, filter_edge_s):
    """Prove local padded reads exactly cover one complete cached interval."""
    records = json_list_field(
        d, "acquisition_chunk_sample_counts_json",
        source=source, required=True)
    hours = _npz_scalar_float(d, "hours", source=source)
    sample_rate = require_integer_sample_rate(
        _npz_scalar_float(d, "sf", source=source), source=source)
    try:
        filter_edge_s = float(filter_edge_s)
    except (TypeError, ValueError, OverflowError) as exc:
        raise RuntimeError(
            f"{source} local acquisition filter padding is invalid") from exc
    if not np.isfinite(filter_edge_s) or filter_edge_s < 0:
        raise RuntimeError(
            f"{source} local acquisition filter padding is invalid")

    total_duration_s = hours * 3600.0
    rounded_duration_s = round(total_duration_s)
    tolerance_s = 1e-7
    if (
        hours <= 0
        or abs(total_duration_s - rounded_duration_s) > tolerance_s
    ):
        raise RuntimeError(
            f"{source} cache duration must be a positive whole number of seconds")
    total_duration_s = float(rounded_duration_s)

    cursor_s = 0.0
    for index, record in enumerate(records):
        if not isinstance(record, dict):
            raise RuntimeError(
                f"{source} local acquisition record {index} is not an object")
        try:
            start_s = float(record["start_s"])
            duration_s = float(record["duration_s"])
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            raise RuntimeError(
                f"{source} local acquisition record {index} has invalid "
                "core geometry") from exc
        requested = record.get("requested_pull_samples")
        returned = record.get("returned_pull_samples")
        if (
            not np.isfinite([start_s, duration_s]).all()
            or duration_s <= 0
            or abs(start_s - cursor_s) > tolerance_s
            or abs(start_s * sample_rate - round(start_s * sample_rate))
                > tolerance_s
            or abs(duration_s * sample_rate - round(duration_s * sample_rate))
                > tolerance_s
        ):
            raise RuntimeError(
                f"{source} local acquisition cores are not an ordered "
                "gap-free sample-aligned partition")
        core_stop_s = start_s + duration_s
        if core_stop_s > total_duration_s + tolerance_s:
            raise RuntimeError(
                f"{source} local acquisition core extends beyond the "
                "cached interval")
        pull_start_s = max(0.0, start_s - filter_edge_s)
        pull_stop_s = min(
            total_duration_s, core_stop_s + filter_edge_s)
        expected_pull_samples = (
            int(round(pull_stop_s * sample_rate))
            - int(round(pull_start_s * sample_rate))
        )
        if (
            not isinstance(requested, int)
            or isinstance(requested, bool)
            or not isinstance(returned, int)
            or isinstance(returned, bool)
            or requested <= 0
            or returned != requested
            or requested != expected_pull_samples
        ):
            raise RuntimeError(
                f"{source} local acquisition record {index} does not prove "
                "the exact padded pull sample count")
        cursor_s = core_stop_s

    if not records or abs(cursor_s - total_duration_s) > tolerance_s:
        raise RuntimeError(
            f"{source} local acquisition cores do not cover the complete "
            "cached interval")
    return {
        "records": records,
        "chunk_count": len(records),
        "duration_s": total_duration_s,
        "sample_rate_hz": sample_rate,
        "filter_edge_s": filter_edge_s,
    }


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
                       run_id=None, run_state="complete", result_files_sha256=None,
                       extra_fields=None):
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
    if run_state not in {"in_progress", "complete", "failed"}:
        raise ValueError(
            "manifest run_state must be 'in_progress', 'complete', or 'failed'")
    requested_set = set(named_lists["requested"])
    completed_set = set(named_lists["completed"])
    skipped_set = set(named_lists["skipped"])
    failed_set = set(named_lists["failed"])
    terminal_sets = (completed_set, skipped_set, failed_set)
    if any(values - requested_set for values in terminal_sets):
        raise ValueError("manifest terminal subjects must all be requested")
    if (
        completed_set & skipped_set
        or completed_set & failed_set
        or skipped_set & failed_set
    ):
        raise ValueError(
            "manifest completed, skipped, and failed partitions must be disjoint")
    hashes = dict(result_files_sha256 or {})
    if run_state == "in_progress":
        if completed_set or skipped_set or failed_set or hashes:
            raise ValueError(
                "an in-progress manifest cannot claim terminal subjects or result hashes")
    else:
        if completed_set | skipped_set | failed_set != requested_set:
            raise ValueError(
                "a terminal manifest must account for every requested subject exactly once")
        if run_state == "complete" and failed_set:
            raise ValueError("a complete manifest cannot contain failed subjects")
        if run_state == "failed" and not failed_set:
            raise ValueError("a failed manifest must contain at least one failed subject")
        expected_hashes = completed_set | skipped_set
        if set(hashes) != expected_hashes:
            raise ValueError(
                "terminal manifest hashes must cover completed and skipped outputs exactly")
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
        "result_files_sha256": hashes,
    }
    if run_state != "in_progress":
        validate_terminal_run_manifest(
            manifest, source=f"{pipeline} terminal output")
    extra_fields = dict(extra_fields or {})
    reserved = set(manifest) & set(extra_fields)
    if reserved:
        raise ValueError(
            f"extra manifest fields cannot replace reserved keys: {sorted(reserved)}")
    manifest.update(extra_fields)
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
