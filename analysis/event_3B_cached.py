"""Reusable 3B estimators for Naji-aligned SO/heart-rate diagnostics.

Two problems with the previous 3B are fixed here:

  1. SCOPE. `event_3B_mednick.py` re-streamed each night and had only ever been run on HUP165, yet
     the summary reported "SO->heartbeat coupling is present but ~25x weaker than published" as a
     cohort-level finding. It was n=1. This runs every cached subject.
  2. THE NULL. Surrogate triggers were drawn from the whole night while the statistic was expressed
     as "% above THIS STAGE's mean HR", so the null measured the stage-vs-night mean-HR offset
     rather than SO-locking. On zero-effect synthetic data that produced z = -21.8 / +24.0 in the
     two stages (test_3B_null.py). Surrogates are now drawn from the same stage as the real
     triggers. The resulting whole-stage shift is retained only as a diagnostic: it is still invalid
     for event-locking inference when local HR trends and SO density are nonstationary.

The corrected Naji-aligned descriptive path uses per-channel 0.15-4 Hz down/up-state events, a 4 Hz RR
tachogram, one RR curve per channel over +/-5 s, and RR-minimum timing. HR conversion is used only
to report percent change. ECG detection is NeuroKit's default, not Naji's Pan-Tompkins plus visual
confirmation. This remains an adapted physiological comparison, not an LC assay.

The standalone writer is withdrawn. Current endpoint execution and artifact
materialization go through ``run_qc_grid.py``.
"""
import hashlib
import numpy as np

from event_3b_estimators import subject_so_triggered, rr_baseline_hr, FS_RR
from lecci_faithful_3A import (
    ANALYSIS_VERSION,
    CacheSubjectSkipped,
    cache_lineage,
    load,
    stages_for,
)
from staging_helpers import EPOCH
from spectral_gapped import contiguous_runs
from pipeline_version import (
    CACHE_SCHEMA_VERSION,
    finite_float_or_none,
    npz_scalar_text,
)

NAJI = {"N2": 12.09, "N3": 3.35}          # % above stage mean HR, frontal scalp, healthy sleepers
MIN_EVENT_CHANNELS = 2
MIN_STABLE_STAGE_EPOCHS = 6                # Naji: uninterrupted 3-min bins at 30 s/epoch
CONTACT_QC_METHOD = "cache stable >=80%-coverage plus optional Destrieux frontal ROI intersection"


def stage_diagnostic_rng(subject, stage):
    """Return a stable RNG isolated to one participant/stage diagnostic.

    This prevents an unavailable or newly added earlier stage from changing
    the later stage's circular-shift diagnostic merely by consuming a shared
    pseudo-random stream.
    """
    key = (
        f"event_3B_cached|{subject}|{stage}|"
        "stage_shift_diagnostic|v1"
    )
    digest = hashlib.sha256(key.encode("utf-8")).digest()
    seed = int.from_bytes(digest[:4], byteorder="little", signed=False)
    return np.random.RandomState(seed)


def stable_stage_epoch_indices(labels, stage, minimum_epochs=MIN_STABLE_STAGE_EPOCHS):
    """Return epochs inside uninterrupted stage runs that satisfy Naji's stable-bin rule."""
    labels = np.asarray(labels).astype(str)
    keep = np.zeros(len(labels), bool)
    for start, stop in contiguous_runs(labels == str(stage)):
        if stop - start >= int(minimum_epochs):
            keep[start:stop] = True
    return np.where(keep)[0]


def stage_so_times(d, channel, keep_epochs, percentile=75):
    """Threshold clean SO candidates within this channel and sleep stage."""
    keys = [f"so_candidate_t_{channel}", f"so_candidate_up_{channel}",
            f"so_candidate_p2p_{channel}"]
    if any(key not in d.files for key in keys):
        return np.array([], float)
    times = np.asarray(d[keys[0]], float)
    up = np.asarray(d[keys[1]], float)
    p2p = np.asarray(d[keys[2]], float)
    keep = np.array([int(value // EPOCH) in keep_epochs for value in times], bool)
    times, up, p2p = times[keep], up[keep], p2p[keep]
    if not len(times):
        return np.array([], float)
    selected = ((up >= np.percentile(up, percentile))
                & (p2p >= np.percentile(p2p, percentile)))
    return np.sort(times[selected])


def coverage_eligible_contacts(d, contacts):
    """Return the stable cache-qualified contact subset used for the 3B event endpoint."""
    if "sigma_selected_contact_mask" not in getattr(d, "files", []):
        raise RuntimeError("cache lacks the prespecified per-contact coverage mask")
    selected = np.asarray(d["sigma_selected_contact_mask"], bool).ravel()
    contacts = [str(value) for value in contacts]
    if len(selected) != len(contacts):
        raise RuntimeError(
            "cache contact coverage mask does not align with cortical_chans")
    if "frontal_selected_contact_mask" in getattr(d, "files", []):
        frontal = np.asarray(d["frontal_selected_contact_mask"], bool).ravel()
        if len(frontal) != len(contacts):
            raise RuntimeError("cache frontal ROI mask does not align with cortical_chans")
        selected &= frontal
    return [contact for contact, keep in zip(contacts, selected) if keep], [
        contact for contact, keep in zip(contacts, selected) if not keep
    ]


def analyse(subject):
    try:
        d = load(subject)
    except CacheSubjectSkipped as exc:
        return dict(
            subject=subject, status="skip", analysis_version=ANALYSIS_VERSION,
            cache_schema_version=CACHE_SCHEMA_VERSION, **cache_lineage(subject),
            reason=f"cache exclusion: {exc}")
    if d is None:
        return None
    lineage = cache_lineage(subject, d)
    # The cache builder already applies the prespecified PCHIP and a hard mask across >5-s beat
    # gaps.  Re-filling here could reopen a 5.1-s acquisition gap whose internal 4-Hz run happens
    # to contain exactly 5.0 s of NaNs.
    rr = np.asarray(d["rr_4"], float)
    hr = 60.0 / rr
    if np.isfinite(rr).sum() < 1000:
        return dict(subject=subject, status="skip", analysis_version=ANALYSIS_VERSION,
                    cache_schema_version=CACHE_SCHEMA_VERSION, **lineage,
                    reason="no usable RR")
    lab, nrem, sep = stages_for(d)
    ctx = [str(c) for c in d["cortical_chans"]]
    eligible_ctx, coverage_excluded_ctx = coverage_eligible_contacts(d, ctx)
    stable_epochs = {
        stage: stable_stage_epoch_indices(lab, stage)
        for stage in ("N2", "N3")
    }
    rec = dict(subject=subject, status="ok", analysis_version=ANALYSIS_VERSION,
               cache_schema_version=CACHE_SCHEMA_VERSION,
               **lineage,
               anatomy_selection_method=npz_scalar_text(
                   d, "anatomy_selection_method", "missing/unvalidated"),
               endpoint_roi=(
                   "non-pathological cortical Destrieux frontal contacts"
                   if "frontal_selected_contact_mask" in d.files
                   else "legacy coverage-qualified cortical contact set"),
               n_N2_raw=int((lab == "N2").sum()),
               n_N3_raw=int((lab == "N3").sum()),
               n_N2=int(len(stable_epochs["N2"])),
               n_N3=int(len(stable_epochs["N3"])),
               stable_stage_minimum_s=int(MIN_STABLE_STAGE_EPOCHS * EPOCH),
               n_cortical_contacts_total=len(ctx),
               n_coverage_eligible_contacts=len(eligible_ctx),
               coverage_eligible_contact_ids=eligible_ctx,
               coverage_excluded_contact_ids=coverage_excluded_ctx,
               contact_qc=(
                   "conservative intersection with the cache's stable >=80%-coverage sigma "
                   "contact set and, when available, Destrieux frontal ROI; at least two contacts "
                   "with >=30 eligible SOs per stage"),
               tachogram_domain="rr",
               n_surrogates=1000,
               null_method="shared circular shift in eligible stage-time",
               minimum_event_channels=MIN_EVENT_CHANNELS,
               gmm_separation=finite_float_or_none(sep),
               whole_night_hr_from_mean_rr=finite_float_or_none(rr_baseline_hr(rr)),
               whole_night_arithmetic_mean_hr=finite_float_or_none(np.nanmean(hr)))
    endpoint_reasons = {}
    for stage in ("N2", "N3"):
        eps = stable_epochs[stage]
        if len(eps) < MIN_STABLE_STAGE_EPOCHS:
            rec[stage] = None
            endpoint_reasons[stage] = (
                f"only {len(eps)} epochs inside uninterrupted >="
                f"{int(MIN_STABLE_STAGE_EPOCHS * EPOCH)}-s {stage} runs "
                f"(<{MIN_STABLE_STAGE_EPOCHS})")
            continue
        m = np.zeros(len(hr), bool)
        for e in eps:
            m[int(e * EPOCH * FS_RR):int((e + 1) * EPOCH * FS_RR)] = True
        m &= np.isfinite(hr)
        if m.sum() < 400:
            rec[stage] = None
            endpoint_reasons[stage] = f"only {int(m.sum())} finite 4-Hz stage samples (<400)"
            continue
        # Event curves are averaged in RR space and converted to HR only afterwards.  Their
        # denominator must use the same transform; mean(60/RR) is not 60/mean(RR).
        stage_mean = rr_baseline_hr(rr[m])
        pool = np.where(m)[0]
        keep = set(eps.tolist())
        troughs_by_channel = []
        event_contact_ids = []
        for c in eligible_ctx:
            tt = stage_so_times(d, c, keep)
            if len(tt) >= 30:
                troughs_by_channel.append(tt)
                event_contact_ids.append(c)
        if len(troughs_by_channel) < MIN_EVENT_CHANNELS:
            rec[stage] = None
            endpoint_reasons[stage] = (
                f"only {len(troughs_by_channel)} coverage-qualified contacts have >=30 eligible "
                f"SOs (<{MIN_EVENT_CHANNELS})")
            continue
        result = subject_so_triggered(
            rr, troughs_by_channel, stage_mean, pool, n_sur=1000,
            rng=stage_diagnostic_rng(subject, stage), domain="rr",
            minimum_channels=MIN_EVENT_CHANNELS, channel_ids=event_contact_ids)
        if result is None:
            rec[stage] = None
            endpoint_reasons[stage] = (
                f"fewer than {MIN_EVENT_CHANNELS} contacts retain >=30 SOs after requiring "
                "complete finite in-stage RR windows")
            continue
        result.pop("curve")
        result.pop("lag_s")
        result["stage_mean_hr"] = stage_mean
        result["excess_over_null_pct"] = float(
            result["pct_above_stage_mean"] - result["null_mean_pct"])
        result["aggregation"] = (
            "one participant-level magnitude from the average channel curve; "
            "Naji latency is the mean of channel peak times")
        result["n_coverage_eligible_contacts"] = len(eligible_ctx)
        result["minimum_event_contacts"] = MIN_EVENT_CHANNELS
        rec[stage] = result
    unavailable = [stage for stage in ("N2", "N3") if rec.get(stage) is None]
    if unavailable:
        rec["status"] = "partial"
        rec["reason"] = "unavailable 3B endpoint(s): " + ", ".join(
            f"{stage} ({endpoint_reasons.get(stage, 'unspecified')})"
            for stage in unavailable)
    rec["endpoint_unavailable_reasons"] = endpoint_reasons
    rec["endpoint_availability"] = {
        stage: rec.get(stage) is not None for stage in ("N2", "N3")
    }
    return rec


def main():
    raise SystemExit(
        "WITHDRAWN standalone 3B writer: run "
        "analysis/run_qc_grid.py to execute the profile-materialized "
        "3B endpoint")


if __name__ == "__main__":
    main()
