"""Focused tests for paired result reporting and compatibility exports."""
from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np

import paired_reporting as reporting
import paired_scalp_ieeg_comparison as comparison


class ReportingMetricTests(unittest.TestCase):
    def test_metric_and_negative_control_ratio(self):
        result = {
            "negative_control": {
                "sigma_window_mean": 6.0,
                "swa_same_window_mean": 3.0,
            },
        }
        self.assertEqual(
            reporting._metric(
                result, ("negative_control", "sigma_window_mean")),
            6.0,
        )
        self.assertEqual(reporting._negative_control_ratio(result), 2.0)
        self.assertIsNone(
            reporting._metric({"value": np.nan}, ("value",)))

    def test_original_module_reexports_reporting_helpers(self):
        self.assertIs(comparison._pair_summary, reporting._pair_summary)
        self.assertIs(comparison._group_summary, reporting._group_summary)
        self.assertIs(
            comparison._normalized_csv_rows,
            reporting._normalized_csv_rows,
        )
        self.assertIs(comparison._write_csv, reporting._write_csv)


class NajiScopeTests(unittest.TestCase):
    def test_labels_and_staged_ieeg_never_claim_exact_naji_comparison(self):
        frozen_subjects = {
            "HUP138_phaseII": {
                "result_3b": {
                    "stages": {
                        stage: {"available_under_profile": True}
                        for stage in reporting.STAGES_3B
                    },
                },
            },
        }
        records = [{
            "subject": "HUP138_phaseII",
            "lineage": {
                "scalp_sidecar_roles": {
                    "f3": "F3",
                    "f4": "F4",
                    "fz": "Fz",
                },
            },
            "scalp": {
                "result_3b": {
                    "f3": {
                        "stages": {
                            "N2": {
                                "available_under_exploratory_profile": True,
                            },
                        },
                    },
                    "fz": {
                        "stages": {
                            "NREM": {
                                "available_under_exploratory_profile": True,
                            },
                        },
                    },
                },
            },
        }]

        status = comparison._naji_label_inventory_status(
            "HUP138_phaseII",
            frozen_subjects,
            records,
        )

        self.assertFalse(status["current_exact_naji_comparison_available"])
        self.assertEqual(
            status["locked_ieeg_3b_available_stages"],
            list(reporting.STAGES_3B),
        )
        self.assertEqual(
            status["current_exploratory_scalp_3b"][
                "unilateral_f3"
            ]["available_stages"],
            ["N2"],
        )
        self.assertTrue(
            status["current_exploratory_scalp_3b"]["f4"]["acquired"]
        )
        self.assertFalse(
            status["current_exploratory_scalp_3b"]["f4"]["analyzed"]
        )
        self.assertFalse(
            status["channel_label_inventory"][
                "labels_establish_naji_reference_montage"
            ]
        )
        self.assertGreaterEqual(len(status["exact_naji_blocking_reasons"]), 3)
        scope = comparison._naji_method_scope()
        self.assertFalse(scope["exact_comparison_available"])
        self.assertIn("not analyzed", scope["f4_status"])

    def test_outcome_neutral_sidecar_status_has_no_stale_version_name(self):
        self.assertEqual(
            comparison.POLARITY_SENSITIVITY_STATUS,
            "not_computable_from_outcome_neutral_sidecar",
        )

    def test_plan_roles_do_not_masquerade_as_terminal_acquisition(self):
        record = {
            "subject": "HUP138_phaseII",
            "lineage": {"scalp_sidecar_roles": {}},
            "scalp": {"result_3b": {}},
        }
        status = comparison._naji_label_inventory_status(
            "HUP138_phaseII",
            {"HUP138_phaseII": {"result_3b": {"stages": {}}}},
            [record],
        )
        self.assertEqual(status["terminal_sidecar_roles_acquired"], [])
        self.assertFalse(
            status["current_exploratory_scalp_3b"]["f4"]["acquired"])


class EndpointLocal3ASupportTests(unittest.TestCase):
    @staticmethod
    def _materialized(signal, swa, hr, *, contact):
        return {
            "contacts": np.asarray([contact]),
            "sigma_parietal": np.asarray(signal, float),
            "swa_parietal": np.asarray(swa, float),
            "sigma_global": np.asarray(signal, float),
            "swa_global": np.asarray(swa, float),
            "hr_1": np.asarray(hr, float),
            "hr_coverage": float(np.isfinite(hr).mean()),
            "hr_meets_profile": True,
            "stage_lab": np.asarray(["N2"] * 20),
            "parietal_power_qc": {
                "sigma": {
                    "support_passes_fit_convergence": True,
                    "selected_mask": np.asarray([True]),
                    "n_selected": 1,
                },
            },
        }

    def test_hr_gap_does_not_remove_eeg_only_spectrum_support(self):
        seconds = np.arange(600, dtype=float)
        ieeg_sigma = (
            2.0
            + 0.25 * np.sin(2 * np.pi * 0.020 * seconds)
            + 0.08 * np.sin(2 * np.pi * 0.047 * seconds)
        )
        scalp_sigma = (
            3.0
            + 0.30 * np.sin(2 * np.pi * 0.020 * seconds + 0.2)
            + 0.06 * np.sin(2 * np.pi * 0.047 * seconds)
        )
        ieeg_swa = 4.0 + 0.2 * np.sin(
            2 * np.pi * 0.011 * seconds)
        scalp_swa = 5.0 + 0.3 * np.sin(
            2 * np.pi * 0.011 * seconds + 0.1)
        hr = 60.0 + np.sin(2 * np.pi * 0.020 * seconds - 0.1)
        # A complete 120-s cardiac gap splits the otherwise continuous
        # 600-s NREM series into two valid cardiac runs. It must not split the
        # EEG support used by spectra, peak fitting, or the SWA control.
        hr[240:360] = np.nan

        ieeg = self._materialized(
            ieeg_sigma, ieeg_swa, hr, contact="A1")
        scalp = self._materialized(
            scalp_sigma, scalp_swa, hr, contact="C3")
        profile = {
            "power": {
                "minimum_contacts": 1,
                "minimum_aggregate_coverage": 0.0,
            },
            "hr": {"minimum_coverage": 0.0},
            "endpoint_3a": {
                "minimum_nrem_epochs": 4,
                "minimum_cross_correlation_windows": 1,
            },
        }
        materializations, support = (
            comparison._endpoint_local_3a_materializations(
                ieeg, scalp, profile)
        )
        ieeg_eeg = comparison.analyse_3a(
            materializations["ieeg"]["eeg"], profile)
        scalp_eeg = comparison.analyse_3a(
            materializations["scalp"]["eeg"], profile)
        ieeg_cardiac = comparison.analyse_3a(
            materializations["ieeg"]["cardiac"], profile)
        scalp_cardiac = comparison.analyse_3a(
            materializations["scalp"]["cardiac"], profile)
        assembled = comparison._assemble_endpoint_local_3a_result(
            ieeg_eeg, ieeg_cardiac)
        geometry = comparison._assert_matched_3a_geometry(
            "SYNTHETIC",
            ieeg_eeg,
            scalp_eeg,
            ieeg_cardiac,
            scalp_cardiac,
            materializations["ieeg"]["cardiac"],
            materializations["scalp"]["cardiac"],
        )

        self.assertEqual(support["eeg_common"]["n_seconds"], 600)
        self.assertEqual(support["cardiac_common"]["n_seconds"], 480)
        self.assertNotEqual(
            support["eeg_common"]["support_mask_sha256"],
            support["cardiac_common"]["support_mask_sha256"],
        )
        self.assertTrue(
            assembled["endpoint_availability"]["spectrum"])
        self.assertEqual(assembled["n_bouts"], 1)
        self.assertGreater(assembled["bout_seconds"], 590.0)
        self.assertEqual(
            assembled["endpoint_support_diagnostics"][
                "coherence_cross_correlation"]["n_bouts"],
            2,
        )
        self.assertLess(
            assembled["endpoint_support_diagnostics"][
                "coherence_cross_correlation"]["bout_seconds"],
            assembled["bout_seconds"],
        )
        self.assertEqual(
            geometry["spectrum_peak_negative_control"]["status"],
            "exact_shared_eeg_geometry",
        )
        self.assertEqual(
            geometry["coherence_cross_correlation"]["status"],
            "exact_shared_cardiac_geometry",
        )

    def test_cardiac_gate_uses_shared_not_unmasked_hr_coverage(self):
        seconds = np.arange(600, dtype=float)
        signal = 2.0 + 0.1 * np.sin(2 * np.pi * 0.02 * seconds)
        swa = 4.0 + 0.1 * np.sin(2 * np.pi * 0.01 * seconds)
        hr = 60.0 + np.sin(2 * np.pi * 0.02 * seconds)
        ieeg_signal = signal.copy()
        ieeg_signal[240:360] = np.nan
        ieeg = self._materialized(
            ieeg_signal, swa, hr, contact="A1")
        scalp = self._materialized(
            signal, swa, hr, contact="C3")
        profile = {"hr": {"minimum_coverage": 0.9}}

        materializations, support = (
            comparison._endpoint_local_3a_materializations(
                ieeg, scalp, profile)
        )

        self.assertEqual(ieeg["hr_coverage"], 1.0)
        self.assertTrue(ieeg["hr_meets_profile"])
        self.assertAlmostEqual(
            support["cardiac_common"]["fraction"], 0.8)
        self.assertFalse(
            support["cardiac_common"][
                "passes_hr_coverage_profile"])
        for arm in ("ieeg", "scalp"):
            cardiac = materializations[arm]["cardiac"]
            self.assertAlmostEqual(cardiac["hr_coverage"], 0.8)
            self.assertFalse(cardiac["hr_meets_profile"])

    def test_cardiac_qc_reason_is_not_reported_as_spectral_failure(self):
        eeg = {
            "endpoint_availability": {
                "spectrum": True,
                "fixed_0p02_coherence": False,
                "cross_correlation": False,
            },
            "coherence": None,
            "cross_correlation": None,
            "aggregate_coverage": 1.0,
            "hr_coverage": 0.2,
            "support_passes_profile": True,
            "n_bouts": 1,
            "bout_seconds": 300.0,
            "support_reasons": [
                "RR coverage 0.200 < 0.800; cardiac endpoints only",
            ],
        }
        cardiac = {
            **eeg,
            "aggregate_coverage": 0.2,
            "support_reasons": [
                "RR coverage 0.200 < 0.800; cardiac endpoints only",
            ],
        }
        assembled = comparison._assemble_endpoint_local_3a_result(
            eeg, cardiac)
        self.assertEqual(assembled["support_reasons"], [])
        self.assertEqual(
            assembled["endpoint_support_diagnostics"][
                "coherence_cross_correlation"]["support_reasons"],
            cardiac["support_reasons"],
        )


class PairSummaryTests(unittest.TestCase):
    def test_fifteen_pairs_use_exact_enumeration(self):
        summary = reporting._pair_summary([
            (float(index), float(index + 1))
            for index in range(15)
        ])
        self.assertEqual(summary["signflip_method"], "exact enumeration")
        self.assertEqual(summary["signflip_draws"], 2 ** 15)
        self.assertIsNone(summary["signflip_seed"])
        self.assertEqual(
            summary["exploratory_exact_signflip_p"],
            summary["exploratory_signflip_p"],
        )

    def test_sixteen_pairs_use_deterministic_monte_carlo(self):
        pairs = [
            (float(index), float(index + (index % 3) - 1))
            for index in range(16)
        ]
        first = reporting._pair_summary(pairs)
        second = reporting._pair_summary(pairs)
        self.assertEqual(first, second)
        self.assertEqual(
            first["signflip_method"],
            "deterministic Monte Carlo sign flips",
        )
        self.assertEqual(first["signflip_draws"], 20_000)
        self.assertIsInstance(first["signflip_seed"], int)
        self.assertIsNone(first["exploratory_exact_signflip_p"])


class CsvSchemaContractTests(unittest.TestCase):
    @staticmethod
    def _canonical_row(fields):
        row = {field: None for field in fields}
        row.update({
            "subject": "S1",
            "question": "3A",
            "comparison_role": "c3",
            "stage": "NREM first 210 min",
            "modality": "iEEG",
        })
        return row

    @staticmethod
    def _write_csv(path, fields, row=None):
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(fields)
            if row is not None:
                writer.writerow([row.get(field, "") for field in fields])

    def test_writers_use_exact_per_artifact_ordered_schemas(self):
        contracts = {
            "paired_metrics.csv": reporting.PAIRED_METRICS_CSV_FIELDS,
            "role_pair_metrics.csv": reporting.ROLE_PAIR_METRICS_CSV_FIELDS,
        }
        with tempfile.TemporaryDirectory() as directory:
            for name, fields in contracts.items():
                with self.subTest(name=name):
                    path = Path(directory) / name
                    reporting._write_csv(
                        path,
                        [self._canonical_row(fields)],
                    )
                    with path.open(newline="", encoding="utf-8") as handle:
                        header = tuple(next(csv.reader(handle)))
                    self.assertEqual(header, fields)
                    self.assertEqual(len(header), len(set(header)))
                    self.assertIn(
                        "shared_eeg_support_fraction", header)
                    self.assertIn(
                        "shared_cardiac_support_fraction", header)
                    self.assertNotIn("shared_support_fraction", header)

    def test_truncated_attacker_headers_are_rejected(self):
        attacks = (
            (
                "paired_metrics.csv",
                reporting.PAIRED_METRICS_CSV_FIELDS,
                ("subject", "question", "stage", "modality"),
            ),
            (
                "role_pair_metrics.csv",
                reporting.ROLE_PAIR_METRICS_CSV_FIELDS,
                ("subject",),
            ),
        )
        with tempfile.TemporaryDirectory() as directory:
            for name, canonical_fields, attacker_fields in attacks:
                with self.subTest(name=name):
                    path = Path(directory) / name
                    expected_row = self._canonical_row(canonical_fields)
                    self._write_csv(path, attacker_fields, expected_row)
                    with self.assertRaisesRegex(
                            RuntimeError, "exact ordered canonical schema"):
                        reporting._validate_csv_artifact(
                            path,
                            [expected_row],
                        )

    def test_missing_extra_duplicate_and_reordered_headers_are_rejected(self):
        fields = reporting.PAIRED_METRICS_CSV_FIELDS
        mutations = {
            "missing": fields[:-1],
            "extra": fields + ("attacker_column",),
            "duplicate": fields + (fields[0],),
            "reordered": (fields[1], fields[0], *fields[2:]),
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "paired_metrics.csv"
            for label, attacker_fields in mutations.items():
                with self.subTest(label=label):
                    self._write_csv(path, attacker_fields)
                    with self.assertRaisesRegex(
                            RuntimeError, "exact ordered canonical schema"):
                        reporting._validate_csv_artifact(path, [])


if __name__ == "__main__":
    unittest.main()
