"""Independent numeric and file-interface regressions; stdlib only."""
from __future__ import annotations

import contextlib
import csv
from fractions import Fraction
import gzip
import hashlib
import io
from itertools import combinations
import json
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

from extremarank import compare_t, score
from extremarank.cli import main


def exact_t_key(values):
    """Signed squared t from the original binary64 fractions, without core moments."""
    xs = [Fraction.from_float(float(value)) for value in values]
    n, s, q = len(xs), sum(xs), sum(value * value for value in xs)
    sign = (s > 0) - (s < 0)
    if not sign:
        return 0, False, Fraction(0)
    variance_numerator = n * q - s * s
    if not variance_numerator:
        return sign, True, Fraction(0)
    return sign, False, (n - 1) * s * s / variance_numerator


def exact_t_compare(left, right):
    ls, li, lv = exact_t_key(left)
    rs, ri, rv = exact_t_key(right)
    if ls != rs:
        return (ls > rs) - (ls < rs)
    if not ls:
        return 0
    if li != ri:
        return ls * ((li > ri) - (li < ri))
    if li:
        return 0
    return ls * ((lv > rv) - (lv < rv))


class CLITests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.directory = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def table(self, gene_ids, rows, donor_ids=None, extension=".tsv"):
        path = self.directory / ("effects" + extension)
        if donor_ids is None:
            donor_ids = [f"d{i + 1}" for i in range(len(rows))]
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="," if extension == ".csv" else "\t", lineterminator="\n")
            writer.writerow(["donor_id", *gene_ids])
            for donor_id, row in zip(donor_ids, rows):
                writer.writerow([donor_id, *(format(float(x), ".17g") for x in row)])
        return path

    def run_cli(self, path, *options):
        output = self.directory / "results"
        with contextlib.redirect_stdout(io.StringIO()):
            result = main([str(path), "--output", str(output), *options])
        self.assertEqual(result, 0)
        report = json.loads((output / "audit.json").read_text(encoding="utf-8"))
        with (output / "summary.tsv").open(encoding="utf-8", newline="") as handle:
            summary = list(csv.DictReader(handle, delimiter="\t"))
        with (output / "envelopes.tsv").open(encoding="utf-8", newline="") as handle:
            envelopes = list(csv.DictReader(handle, delimiter="\t"))
        return report, summary, envelopes

    def assert_invalid(self, path, *options):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as raised:
            main([str(path), "--output", str(self.directory / "invalid"), *options])
        self.assertEqual(raised.exception.code, 2)
        self.assertFalse((self.directory / "invalid").exists())

    def test_cross_cardinality_summary_and_exhaustive_witnesses(self):
        values = [1.0, 2.0, 3.0, 4.0]
        path = self.table(["positive"], [[x] for x in values])
        report, rows, envelopes = self.run_cli(path, "--budget", "2", "--top-k", "1", "--skip-topk")
        self.assertEqual(report["status"], "NOT_RUN")
        self.assertEqual(report["topk_audit"]["status"], "NOT_RUN")
        self.assertEqual(report["input"]["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
        row = rows[0]
        self.assertAlmostEqual(float(row["baseline_t"]), math.sqrt(15))
        self.assertAlmostEqual(float(row["worst_t"]), 5 / 3)
        self.assertAlmostEqual(float(row["best_t"]), 7)
        self.assertEqual(json.loads(row["worst_deleted_donors"]), ["d2", "d3"])
        self.assertEqual(json.loads(row["best_deleted_donors"]), ["d1", "d2"])
        self.assertEqual(row["worst_retained_count"], "2")
        self.assertEqual(row["best_retained_count"], "2")
        self.assertEqual(float(row["worst_mean"]), 2.5)
        self.assertEqual(float(row["worst_sample_variance"]), 4.5)
        self.assertEqual(row["robust_positive"], "True")
        self.assertEqual(row["robust_negative"], "False")
        self.assertEqual(row["baseline_topk"], "True")
        self.assertEqual(len(envelopes), 3)
        for envelope in envelopes:
            kept_count = int(envelope["retained_count"])
            candidates = [[values[i] for i in kept] for kept in combinations(range(4), kept_count)]
            for kind in ("minimum", "maximum"):
                deleted = set(json.loads(envelope[kind + "_deleted_donors"]))
                actual = [value for i, value in enumerate(values) if f"d{i + 1}" not in deleted]
                self.assertEqual(len(actual), kept_count)
                for candidate in candidates:
                    comparison = exact_t_compare(actual, candidate)
                    self.assertLessEqual(comparison, 0) if kind == "minimum" else self.assertGreaterEqual(comparison, 0)
                xs = [Fraction.from_float(value) for value in actual]
                self.assertEqual(Fraction(envelope[kind + "_sum"]), sum(xs))
                self.assertEqual(Fraction(envelope[kind + "_sum_squares"]), sum(x*x for x in xs))
            self.assertEqual(int(envelope["total_subsets"]), math.comb(4, kept_count))

    def test_directions_and_lexical_absolute_ties(self):
        path = self.table(["z_negative", "a_positive", "m_mixed", "b_zero"],
                          [[-1, 1, 2, 0], [-2, 2, -3, 0], [-3, 3, 4, 0], [-4, 4, -1, 0]])
        for direction, expected in (("up", ["a_positive", "m_mixed", "b_zero", "z_negative"]),
                                    ("down", ["z_negative", "b_zero", "m_mixed", "a_positive"]),
                                    ("absolute", ["a_positive", "z_negative", "m_mixed", "b_zero"])):
            with self.subTest(direction=direction):
                report, rows, envelopes = self.run_cli(path, "--direction", direction, "--budget", "0", "--top-k", "2")
                self.assertEqual([row["gene_id"] for row in rows], expected)
                self.assertEqual(report["topk_audit"]["baseline_topk"], expected[:2])
                self.assertEqual(report["parameters"]["direction"], direction)
                self.assertEqual(report["status"], "CERTIFIED")
                # Scalar envelopes are signed even for an absolute-rank audit.
                negative = next(row for row in envelopes if row["gene_id"] == "z_negative")
                self.assertLess(float(negative["minimum_t"]), 0)

    def test_default_topk_reduction_and_explicit_large_topk_rejected(self):
        path = self.table(["A", "B"], [[4, 0], [4, 0], [-1, 0], [-1, 0]])
        report, _, _ = self.run_cli(path)
        self.assertEqual(report["parameters"]["top_k"], 2)
        self.assertTrue(report["parameters"]["top_k_default_reduced"])
        self.assertIsNone(report["parameters"]["top_k_requested"])
        self.assertEqual(report["status"], "CERTIFIED")
        self.assert_invalid(path, "--top-k", "20")

    def test_refutation_is_an_actual_shared_deletion(self):
        values = [[4, 0], [4, 0], [-1, 0], [-1, 0]]
        path = self.table(["A", "B"], values)
        report, _, _ = self.run_cli(path, "--top-k", "1")
        self.assertEqual(report["status"], "REFUTED")
        audit = report["topk_audit"]
        self.assertEqual(audit["baseline_topk"], ["A"])
        self.assertEqual(audit["witness_topk"], ["B"])
        deleted = set(audit["deleted_indices"])
        retained_a = [row[0] for i, row in enumerate(values) if i not in deleted]
        retained_b = [row[1] for i, row in enumerate(values) if i not in deleted]
        self.assertLess(exact_t_compare(retained_a, retained_b), 0)
        self.assertEqual(audit["deleted_donors"], [f"d{i+1}" for i in audit["deleted_indices"]])
        self.assertLessEqual(len(deleted), 2)

    def test_csv_quoted_ids_constants_and_nonfinite_json_display(self):
        donors = ["donor,1", 'donor"2', "δ3", "D04"]
        path = self.table(["constant\tpositive", "constant_negative", "all_zero"],
                          [[1, -2, 0]] * 4, donors, ".csv")
        report, rows, envelopes = self.run_cli(path, "--top-k", "1")
        by_gene = {row["gene_id"]: row for row in rows}
        self.assertEqual(report["input"]["donor_ids"], donors)
        self.assertEqual(report["input"]["delimiter"], "comma")
        self.assertEqual(report["baseline_topk_scores"][0]["t"], "inf")
        self.assertEqual(by_gene["constant\tpositive"]["worst_t"], "inf")
        self.assertEqual(by_gene["constant_negative"]["best_t"], "-inf")
        self.assertEqual(by_gene["all_zero"]["worst_t"], "0")
        self.assertEqual(by_gene["all_zero"]["robust_positive"], "False")
        self.assertEqual(by_gene["constant_negative"]["robust_negative"], "True")
        self.assertEqual(by_gene["all_zero"]["robust_negative"], "False")
        self.assertEqual(by_gene["constant\tpositive"]["worst_retained_count"], "4")
        self.assertEqual(json.loads(by_gene["constant\tpositive"]["worst_deleted_donors"]), [])
        self.assertEqual(len(envelopes), 9)
        # A standards-compliant JSON parser must never see NaN/Infinity tokens.
        text = (self.directory / "results" / "audit.json").read_text(encoding="utf-8")
        json.loads(text, parse_constant=lambda token: self.fail(f"nonstandard JSON constant: {token}"))

    def test_node_limit_does_not_change_exact_scalar_outputs(self):
        # A=B+1/8 retains identical variance and strictly greater positive t
        # for every common subset, but their separate marginal ranges overlap.
        path = self.table(["A", "B"], [[x + 0.125, x] for x in (1, 2, 3, 4)])
        report, rows, envelopes = self.run_cli(path, "--top-k", "1", "--max-nodes", "1")
        self.assertEqual(report["status"], "UNRESOLVED")
        self.assertEqual(report["topk_audit"]["nodes"], 1)
        self.assertEqual(report["topk_audit"]["deleted_donors"], None)
        self.assertEqual(len(rows), 2)
        self.assertEqual(len(envelopes), 6)
        by_gene = {row["gene_id"]: row for row in rows}
        self.assertEqual(float(by_gene["B"]["best_t"]), 7)

    def test_gzip_roundtrip_and_original_file_hash(self):
        for extension in (".csv", ".tsv"):
            with self.subTest(extension=extension):
                path = self.table(["A", "B"], [[4, 0], [4, 0], [-1, 0], [-1, 0]], extension=extension)
                reference, rows, envelopes = self.run_cli(path, "--top-k", "1")
                compressed = Path(str(path) + ".gz")
                compressed.write_bytes(gzip.compress(path.read_bytes(), mtime=0))
                report, compressed_rows, compressed_envelopes = self.run_cli(compressed, "--top-k", "1")
                self.assertEqual(compressed_rows, rows)
                self.assertEqual(compressed_envelopes, envelopes)
                self.assertEqual(report["topk_audit"], reference["topk_audit"])
                self.assertEqual(report["input"]["compression"], "gzip")
                self.assertEqual(report["input"]["sha256"], hashlib.sha256(compressed.read_bytes()).hexdigest())
                self.assertNotEqual(report["input"]["sha256"], reference["input"]["sha256"])
                self.assertEqual(report["input"]["delimiter"], "comma" if extension == ".csv" else "tab")

    def test_invalid_and_truncated_gzip(self):
        path = self.directory / "effects.csv.gz"
        for raw in (b"not gzip", gzip.compress(b"donor_id,A\nd1,1\nd2,2\n", mtime=0)[:-4]):
            with self.subTest(raw=raw):
                path.write_bytes(raw)
                self.assert_invalid(path, "--budget", "0")

    def test_invalid_input_and_parameters(self):
        cases = {
            "wrong_header": "sample\tA\nd1\t1\nd2\t2\n",
            "duplicate_gene": "donor_id\tA\tA\nd1\t1\t2\nd2\t3\t4\n",
            "empty_gene": "donor_id\t\nd1\t1\nd2\t2\n",
            "duplicate_donor": "donor_id\tA\nd1\t1\nd1\t2\n",
            "empty_donor": "donor_id\tA\n\t1\nd2\t2\n",
            "short_row": "donor_id\tA\tB\nd1\t1\nd2\t2\t3\n",
            "nan": "donor_id\tA\nd1\tnan\nd2\t2\n",
            "inf": "donor_id\tA\nd1\t1e999\nd2\t2\n",
            "missing_value": "donor_id\tA\nd1\t\nd2\t2\n",
            "single_donor": "donor_id\tA\nd1\t1\n",
        }
        for name, source in cases.items():
            with self.subTest(case=name):
                path = self.directory / (name + ".tsv")
                path.write_text(source, encoding="utf-8")
                self.assert_invalid(path, "--budget", "0")
        path = self.table(["A"], [[1], [2]])
        self.assert_invalid(path)  # Default budget two is not silently changed.
        self.assert_invalid(path, "--budget", "-1")
        self.assert_invalid(path, "--budget", "0", "--top-k", "0")
        self.assert_invalid(path, "--budget", "0", "--max-nodes", "0", "--skip-topk")

    def test_compare_t_exact_across_counts_signs_and_infinities(self):
        vectors = [[0, 0], [1, 1], [1, 1, 1, 1], [-1, -1], [-2, -2, -2],
                   [-1, 1], [1, 2], [1, 2, 3, 4], [-1, -2], [-1, -2, -3, -4],
                   [math.ldexp(1.0, -1074), 0], [math.ldexp(1.0, 900), 0, -math.ldexp(1.0, 900)],
                   [1, math.nextafter(1.0, 2.0)], [1, 1, math.nextafter(1.0, 2.0)]]
        rng = random.Random(20261004)
        vectors.extend([[rng.randrange(-12, 13) / 8 for _ in range(rng.randrange(2, 7))] for _ in range(20)])
        for left in vectors:
            for right in vectors:
                with self.subTest(left=left, right=right):
                    self.assertEqual(compare_t(score(left), score(right)), exact_t_compare(left, right))

    def test_python_module_entrypoint(self):
        path = self.table(["A"], [[1], [2]])
        output = self.directory / "module_results"
        source = Path(__file__).resolve().parents[1] / "src"
        environment = dict(os.environ)
        environment["PYTHONPATH"] = str(source) + os.pathsep + environment.get("PYTHONPATH", "")
        process = subprocess.run([sys.executable, "-m", "extremarank", str(path), "--budget", "0",
                                  "--output", str(output)], env=environment, capture_output=True, text=True)
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertIn("CERTIFIED", process.stdout)
        self.assertTrue((output / "audit.json").exists())


if __name__ == "__main__":
    unittest.main()
