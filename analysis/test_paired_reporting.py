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
