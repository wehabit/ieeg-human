"""Synthetic checks for overlap-connected contact aggregation."""
from __future__ import annotations

import numpy as np

from overlap_aggregate import (
    largest_observation_component,
    median_polish_contact_time,
    overlap_connected_aggregate,
    overlap_connected_staging,
)


def check(name, condition):
    print(f"  {'PASS' if condition else 'FAIL'}  {name}")
    if not condition:
        raise AssertionError(name)


# Full overlap with multiplicative contact gains must recover the common time course.
t = np.linspace(0, 4 * np.pi, 300)
truth = np.exp(0.3 * np.sin(t))
gains = np.array([0.2, 1.0, 8.0])
complete = gains[:, None] * truth[None, :]
aggregate, details = overlap_connected_aggregate(complete)
truth_normalized = truth / np.median(truth)
check(
    "median-polish aggregation removes stable contact gain",
    details["n_components"] == 1
    and np.max(np.abs(aggregate - truth_normalized)) < 1e-8,
)

# Staggered but overlapping availability remains identifiable.
staggered = complete.copy()
staggered[0, 180:] = np.nan
staggered[1, :80] = np.nan
staggered[2, :160] = np.nan
staggered_aggregate, staggered_details = overlap_connected_aggregate(staggered)
keep = np.isfinite(staggered_aggregate)
check(
    "overlapping partial contacts recover one connected time series",
    staggered_details["n_components"] == 1
    and keep.sum() == len(t)
    and np.corrcoef(staggered_aggregate[keep], truth_normalized[keep])[0, 1] > 0.999,
)

# The A44 construction must not become a superficially complete aggregate.
disjoint = np.full((6, 600), np.nan)
for contact in range(6):
    disjoint[contact, contact * 100:(contact + 1) * 100] = contact + 1.0
disjoint_aggregate, disjoint_details = overlap_connected_aggregate(disjoint)
check(
    "disjoint contact blocks remain separate rather than manufacturing 100% coverage",
    disjoint_details["n_components"] == 6
    and np.isfinite(disjoint_aggregate).sum() == 100
    and disjoint_details["aggregate_coverage"] == 1 / 6,
)

# Component selection is deterministic even when two components have identical size/support.
tie = np.zeros((2, 8), bool)
tie[0, :4] = True
tie[1, 4:] = True
tie_mask, tie_details = largest_observation_component(tie)
check(
    "equal-support disconnected components use the documented earliest-time tie break",
    tie_details["n_components"] == 2
    and tie_details["chosen_component"]["time_indices"] == [0, 1, 2, 3]
    and np.array_equal(np.where(tie_mask.any(axis=0))[0], np.arange(4)),
)

# An ineligible contact cannot serve as a bridge between otherwise disconnected time blocks.
bridge = np.zeros((3, 8), bool)
bridge[0, :4] = True
bridge[1, :] = True
bridge[2, 4:] = True
without_bridge, without_bridge_details = largest_observation_component(
    bridge, eligible_contacts=np.array([True, False, True]))
check(
    "ineligible contacts cannot connect observation components",
    without_bridge_details["n_components"] == 2
    and np.array_equal(np.where(without_bridge.any(axis=0))[0], np.arange(4)),
)

# Joint staging uses the same connected observation graph for delta ratio and SWA.
dr = np.tile(0.4 + 0.1 * np.sin(t), (3, 1))
swa = complete.copy()
clean = np.ones_like(dr)
windows = np.full_like(dr, 10, dtype=int)
dr[0, 180:] = np.nan
swa[0, 180:] = np.nan
dr[1, :80] = np.nan
swa[1, :80] = np.nan
dr[2, :160] = np.nan
swa[2, :160] = np.nan
ep_dr, ep_swa, ep_clean, staging_details = overlap_connected_staging(
    dr, swa, clean, np.ones(3, bool),
    valid_window_count_by_contact=windows, minimum_valid_windows=7)
check(
    "joint overlap-connected staging returns aligned finite features",
    staging_details["n_components"] == 1
    and np.array_equal(np.isfinite(ep_dr), np.isfinite(ep_swa))
    and np.isfinite(ep_clean).all()
    and np.corrcoef(ep_swa, truth_normalized)[0, 1] > 0.999,
)

masked_windows = windows.copy()
masked_windows[:, 100:130] = 6
masked_dr, masked_swa, _, masked_details = overlap_connected_staging(
    dr, swa, clean, np.ones(3, bool),
    valid_window_count_by_contact=masked_windows, minimum_valid_windows=7)
check(
    "staging observations below the locked Welch-window support are excluded",
    not np.isfinite(masked_dr[100:130]).any()
    and not np.isfinite(masked_swa[100:130]).any()
    and not masked_details["selected_epoch_mask"][100:130].any(),
)

too_few_dr, too_few_swa, too_few_clean, too_few_details = (
    overlap_connected_staging(
        dr, swa, clean, np.array([True, False, False]),
        valid_window_count_by_contact=windows,
        minimum_valid_windows=7,
        minimum_contacts=2,
    )
)
check(
    "overlap staging fails closed below the profile's minimum contact count",
    too_few_details["n_selected_contacts"] == 0
    and too_few_details["n_components_meeting_minimum_contacts"] == 0
    and not too_few_details["support_passes_minimum_contacts"]
    and not np.isfinite(too_few_dr).any()
    and not np.isfinite(too_few_swa).any()
    and not np.isfinite(too_few_clean).any(),
)

# Without contact overlap, a low-gain half and high-gain half cannot create two global stages.
changing_dr = np.full((6, 300), 0.6)
changing_swa = np.r_[np.full((3, 300), 0.5), np.full((3, 300), 5.0)]
changing_dr[:3, 150:] = np.nan
changing_swa[:3, 150:] = np.nan
changing_dr[3:, :150] = np.nan
changing_swa[3:, :150] = np.nan
changing_windows = np.where(np.isfinite(changing_swa), 14, 0)
split_dr, split_swa, _, split_details = overlap_connected_staging(
    changing_dr, changing_swa, np.ones_like(changing_swa), np.ones(6, bool),
    valid_window_count_by_contact=changing_windows, minimum_valid_windows=7)
check(
    "disjoint staging contact groups cannot manufacture a two-state full-night series",
    split_details["n_components"] == 2
    and np.isfinite(split_dr).sum() == 150
    and np.isfinite(split_swa).sum() == 150,
)

# This deterministic sparse case needs 148 iterations at the locked tolerance.
# The production cap is 100, so the last iterate must never escape as a result.
rng = np.random.RandomState(42)
for nonconvergence_case in range(18):
    nonconvergent_values = np.exp(rng.normal(size=(8, 30)))
    nonconvergent_draws = rng.rand(8, 30)
    nonconvergent_density = rng.uniform(0.1, 0.9)
    nonconvergent_mask = nonconvergent_draws < nonconvergent_density
nonconvergent_aggregate, nonconvergent_fit = median_polish_contact_time(
    nonconvergent_values,
    nonconvergent_mask,
    max_iter=100,
    tolerance=1e-8,
)
check(
    "median polish fails closed when the locked iteration cap is reached",
    nonconvergence_case == 17
    and nonconvergent_mask.sum() == 82
    and not nonconvergent_fit["converged"]
    and nonconvergent_fit["iterations"] == 100
    and not np.isfinite(nonconvergent_aggregate).any(),
)
