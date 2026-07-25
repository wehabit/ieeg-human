"""Contact aggregation based on observable overlap rather than an arbitrary coverage cutoff.

For partially observed contact-by-time matrices, a changing set of raw contact values can create
false dynamics.  A fixed 80% rule prevents one version of that failure but discards recoverable
overlap and has no paper-derived boundary.  This module instead:

1. builds the bipartite observation graph (contact nodes connected to time/epoch nodes);
2. retains one deterministically selected largest connected component;
3. estimates a contact offset and common time effect by robust median polish.

Disjoint contact groups are therefore never spliced into one apparently complete series.  Partly
observed contacts can contribute when their overlap makes relative offsets identifiable.
"""
from __future__ import annotations

import warnings

import numpy as np


class _UnionFind:
    def __init__(self, n):
        self.parent = np.arange(int(n))
        self.size = np.ones(int(n), int)

    def find(self, value):
        value = int(value)
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = int(self.parent[value])
        return value

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left == right:
            return
        if self.size[left] < self.size[right]:
            left, right = right, left
        self.parent[right] = left
        self.size[left] += self.size[right]


def largest_observation_component(
        valid, eligible_contacts=None, *, minimum_contacts=1):
    """Return the observation mask for the largest identifiable contact-time component."""
    valid = np.asarray(valid, bool)
    if valid.ndim != 2:
        raise ValueError("valid must be a contact-by-time matrix")
    n_contacts, n_times = valid.shape
    eligible = (
        np.ones(n_contacts, bool)
        if eligible_contacts is None
        else np.asarray(eligible_contacts, bool).ravel()
    )
    if len(eligible) != n_contacts:
        raise ValueError("eligible_contacts must align with contacts")
    valid = valid & eligible[:, None]
    union = _UnionFind(n_contacts + n_times)
    contact_idx, time_idx = np.where(valid)
    for contact, time in zip(contact_idx, time_idx):
        union.union(int(contact), n_contacts + int(time))

    components = {}
    for contact, time in zip(contact_idx, time_idx):
        root = union.find(int(contact))
        value = components.setdefault(
            root, {"contacts": set(), "times": set(), "n_observations": 0})
        value["contacts"].add(int(contact))
        value["times"].add(int(time))
        value["n_observations"] += 1
    summaries = []
    for value in components.values():
        contacts = sorted(value["contacts"])
        times = sorted(value["times"])
        summaries.append(dict(
            contact_indices=contacts,
            time_indices=times,
            n_contacts=len(contacts),
            n_times=len(times),
            n_observations=int(value["n_observations"]),
            first_contact=contacts[0],
            first_time=times[0],
        ))
    summaries.sort(
        key=lambda value: (
            -value["n_times"],
            -value["n_observations"],
            -value["n_contacts"],
            value["first_time"],
            value["first_contact"],
        )
    )
    candidates = [
        value for value in summaries
        if value["n_contacts"] >= int(minimum_contacts)
    ]
    selected = np.zeros_like(valid)
    if candidates:
        chosen = candidates[0]
        contact_mask = np.zeros(n_contacts, bool)
        time_mask = np.zeros(n_times, bool)
        contact_mask[chosen["contact_indices"]] = True
        time_mask[chosen["time_indices"]] = True
        selected = valid & contact_mask[:, None] & time_mask[None, :]
    else:
        chosen = dict(
            contact_indices=[], time_indices=[], n_contacts=0, n_times=0,
            n_observations=0, first_contact=None, first_time=None)
    return selected, dict(
        n_components=len(summaries),
        components=summaries,
        minimum_contacts=int(minimum_contacts),
        n_components_meeting_minimum_contacts=len(candidates),
        chosen_component=chosen,
        selection=(
            "among components meeting the locked minimum contact count: "
            "most time nodes, then observations, then contacts; "
            "deterministic earliest-time/contact tie break"),
    )


def _transform(values, kind):
    values = np.asarray(values, float)
    if kind == "log":
        valid = np.isfinite(values) & (values > 0)
        transformed = np.full(values.shape, np.nan)
        transformed[valid] = np.log(values[valid])
        return transformed, valid
    if kind == "logit":
        valid = np.isfinite(values) & (values > 0) & (values < 1)
        transformed = np.full(values.shape, np.nan)
        eps = np.finfo(float).eps
        clipped = np.clip(values[valid], eps, 1.0 - eps)
        transformed[valid] = np.log(clipped / (1.0 - clipped))
        return transformed, valid
    if kind == "identity":
        valid = np.isfinite(values)
        return values.copy(), valid
    raise ValueError("transform must be 'log', 'logit', or 'identity'")


def _inverse(values, kind):
    if kind == "log":
        return np.exp(values)
    if kind == "logit":
        positive = values >= 0
        out = np.empty_like(values)
        out[positive] = 1.0 / (1.0 + np.exp(-values[positive]))
        exp_value = np.exp(values[~positive])
        out[~positive] = exp_value / (1.0 + exp_value)
        return out
    return values


def median_polish_contact_time(values, observation_mask, *, transform="log",
                               max_iter=100, tolerance=1e-8):
    """Robust additive contact-offset/time-effect fit on one connected component."""
    transformed, intrinsic_valid = _transform(values, transform)
    observed = np.asarray(observation_mask, bool) & intrinsic_valid
    if observed.shape != transformed.shape:
        raise ValueError("observation_mask must align with values")
    n_contacts, n_times = transformed.shape
    contact_used = observed.any(axis=1)
    time_used = observed.any(axis=0)
    fitted = np.full(n_times, np.nan)
    contact_offsets = np.full(n_contacts, np.nan)
    time_effects = np.full(n_times, np.nan)
    if not observed.any():
        return fitted, dict(
            converged=False, iterations=0, grand_effect=None,
            contact_offsets=contact_offsets, time_effects=time_effects)

    grand = float(np.median(transformed[observed]))
    contact = np.zeros(n_contacts, float)
    time = np.zeros(n_times, float)
    converged = False
    iteration = 0
    for iteration in range(1, int(max_iter) + 1):
        previous = np.r_[grand, contact[contact_used], time[time_used]]
        # ``n_times`` is several thousand for a full night.  Updating one column at a time made
        # an otherwise offline sensitivity grid take minutes per profile.  NaN-masked reductions
        # implement the identical median updates in vectorized NumPy code.
        contact_residual = np.where(
            observed, transformed - grand - time[None, :], np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            contact_update = np.nanmedian(contact_residual, axis=1)
        contact[contact_used] = contact_update[contact_used]
        contact_center = float(np.median(contact[contact_used]))
        contact[contact_used] -= contact_center
        grand += contact_center
        time_residual = np.where(
            observed, transformed - grand - contact[:, None], np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            time_update = np.nanmedian(time_residual, axis=0)
        time[time_used] = time_update[time_used]
        time_center = float(np.median(time[time_used]))
        time[time_used] -= time_center
        grand += time_center
        current = np.r_[grand, contact[contact_used], time[time_used]]
        if np.max(np.abs(current - previous)) <= float(tolerance):
            converged = True
            break

    contact_offsets[contact_used] = contact[contact_used]
    time_effects[time_used] = time[time_used]
    fitted[time_used] = _inverse(grand + time[time_used], transform)
    if transform == "log" and np.isfinite(fitted).any():
        fitted /= float(np.nanmedian(fitted))
    residual = np.full(transformed.shape, np.nan)
    ci, ti = np.where(observed)
    residual[ci, ti] = (
        transformed[ci, ti] - grand - contact[ci] - time[ti])
    # A finite last iterate is not an estimate with established numerical
    # stability.  Fail closed instead of allowing a profile result to depend
    # on an arbitrary iteration cap.
    if not converged:
        fitted[:] = np.nan
    return fitted, dict(
        converged=converged,
        iterations=int(iteration),
        grand_effect=float(grand),
        contact_offsets=contact_offsets,
        time_effects=time_effects,
        residual_mad=float(
            1.4826 * np.median(np.abs(residual[observed] - np.median(residual[observed])))),
        n_observations=int(observed.sum()),
        contact_support=observed.sum(axis=1),
        time_support=observed.sum(axis=0),
        transform=transform,
    )


def overlap_connected_aggregate(values, *, eligible_contacts=None, observation_mask=None,
                                transform="log", minimum_contacts=1):
    """Aggregate a contact-by-time matrix without joining nonoverlapping contact groups."""
    transformed, intrinsic_valid = _transform(values, transform)
    del transformed
    valid = intrinsic_valid
    if observation_mask is not None:
        supplied = np.asarray(observation_mask, bool)
        if supplied.shape != valid.shape:
            raise ValueError("observation_mask must align with values")
        valid &= supplied
    component_mask, component = largest_observation_component(
        valid, eligible_contacts=eligible_contacts,
        minimum_contacts=minimum_contacts)
    aggregate, fit = median_polish_contact_time(
        values, component_mask, transform=transform)
    fit_converged = bool(fit["converged"])
    details = {
        **component,
        "fit": fit,
        "selected_contact_mask": component_mask.any(axis=1),
        "selected_time_mask": component_mask.any(axis=0),
        "contact_count": component_mask.sum(axis=0),
        "per_contact_observations": component_mask.sum(axis=1),
        "aggregate_coverage": float(np.isfinite(aggregate).mean()),
        "support_passes_fit_convergence": fit_converged,
        "aggregation": "largest overlap-connected component with robust median polish",
    }
    return aggregate, details


def overlap_connected_staging(
        dr_by_contact, swa_by_contact, clean_fraction_by_contact,
        candidate_contact_mask, *, valid_window_count_by_contact,
        minimum_valid_windows, minimum_contacts=1):
    """Joint overlap-connected staging aggregation on one observation graph."""
    dr = np.asarray(dr_by_contact, float)
    swa = np.asarray(swa_by_contact, float)
    clean = np.asarray(clean_fraction_by_contact, float)
    windows = np.asarray(valid_window_count_by_contact)
    if dr.shape != swa.shape or dr.shape != clean.shape or dr.shape != windows.shape:
        raise ValueError("staging contact-by-epoch inputs must align")
    joint = (
        np.isfinite(dr) & (dr > 0) & (dr < 1)
        & np.isfinite(swa) & (swa > 0)
        & (windows >= int(minimum_valid_windows))
    )
    component_mask, component = largest_observation_component(
        joint, eligible_contacts=candidate_contact_mask,
        minimum_contacts=minimum_contacts)
    ep_dr, dr_fit = median_polish_contact_time(
        dr, component_mask, transform="logit")
    ep_swa, swa_fit = median_polish_contact_time(
        swa, component_mask, transform="log")
    ep_clean = np.full(dr.shape[1], np.nan)
    for epoch in np.where(component_mask.any(axis=0))[0]:
        ep_clean[epoch] = float(np.median(clean[component_mask[:, epoch], epoch]))
    n_selected_contacts = int(component_mask.any(axis=1).sum())
    support_passes_minimum_contacts = bool(
        n_selected_contacts >= int(minimum_contacts))
    support_passes_fit_convergence = bool(
        dr_fit["converged"] and swa_fit["converged"])
    support_passes = bool(
        support_passes_minimum_contacts and support_passes_fit_convergence)
    if not support_passes:
        ep_dr[:] = np.nan
        ep_swa[:] = np.nan
        ep_clean[:] = np.nan
    details = {
        **component,
        "dr_fit": dr_fit,
        "swa_fit": swa_fit,
        "selected_contact_mask": component_mask.any(axis=1),
        "selected_epoch_mask": component_mask.any(axis=0),
        "contact_count": component_mask.sum(axis=0),
        "per_contact_feature_coverage": joint.mean(axis=1),
        "minimum_valid_welch_windows": int(minimum_valid_windows),
        "minimum_contacts": int(minimum_contacts),
        "n_selected_contacts": n_selected_contacts,
        "support_passes_minimum_contacts": support_passes_minimum_contacts,
        "support_passes_fit_convergence": support_passes_fit_convergence,
        "support_passes": support_passes,
        "aggregate_coverage": float(
            (np.isfinite(ep_dr) & np.isfinite(ep_swa)).mean()),
        "aggregation": "joint largest overlap-connected component with robust median polish",
    }
    return ep_dr, ep_swa, ep_clean, details
