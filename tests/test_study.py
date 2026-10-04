import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from extremarank.cli import main
from extremarank.external import compare_rankings
from extremarank.preparation import prepare_study


class StudyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.matrix = self.root / "matrix.tsv"
        self.meta = self.root / "meta.tsv"
        self.matrix.write_text("feature_id\ta\tb\tc\td\te\tf\ng1\t4\t6\t5\t1\t2\t3\ng2\t1\t1\t2\t3\t4\t3\n")
        self.meta.write_text("sample_id\tgroup\tdonor_id\nf\tR\tf\nc\tT\tc\na\tT\ta\ne\tR\te\nb\tT\tb\nd\tR\td\n")

    def tearDown(self):
        self.tmp.cleanup()

    def test_strict_join_reorders_and_direction_is_explicit(self):
        study = prepare_study(self.matrix, self.meta, "T", "R", design="welch")
        self.assertEqual(study.groups, ("T", "T", "T", "R", "R", "R"))
        self.assertEqual(study.values[0], (4, 1))
        self.meta.write_text(self.meta.read_text().replace("f\tR\tf", "x\tR\tx"))
        with self.assertRaisesRegex(ValueError, "match exactly"):
            prepare_study(self.matrix, self.meta, "T", "R", design="welch")

    def test_paired_contrast_incomplete_and_repeated_samples(self):
        self.meta.write_text("sample_id\tgroup\tdonor_id\na\tT\t1\nb\tT\t2\nc\tT\t3\nd\tR\t1\ne\tR\t2\nf\tR\t3\n")
        study = prepare_study(self.matrix, self.meta, "T", "R")
        self.assertEqual(study.values, ((3, -2), (4, -3), (2, -1)))
        reverse = prepare_study(self.matrix, self.meta, "R", "T")
        self.assertEqual(reverse.values, tuple(tuple(-x for x in row) for row in study.values))
        with self.assertRaisesRegex(ValueError, "independent"):
            prepare_study(self.matrix, self.meta, "T", "R", design="welch")
        self.meta.write_text(self.meta.read_text().replace("f\tR\t3", "f\tOther\t3"))
        with self.assertRaisesRegex(ValueError, "incomplete"):
            prepare_study(self.matrix, self.meta, "T", "R")
        dropped = prepare_study(self.matrix, self.meta, "T", "R", incomplete_pairs="drop")
        self.assertEqual(dropped.unit_ids, ("1", "2"))
        self.assertEqual(dropped.manifest["excluded_incomplete_donors"], ["3"])

    def test_explicit_missing_policy_and_cpm_totals(self):
        self.matrix.write_text(self.matrix.read_text().replace("g2\t1", "g2\tNA"))
        with self.assertRaisesRegex(ValueError, "missing"):
            prepare_study(self.matrix, self.meta, "T", "R", design="welch")
        dropped = prepare_study(self.matrix, self.meta, "T", "R", design="welch", missing="drop-features")
        self.assertEqual(dropped.feature_ids, ("g1",))
        self.assertEqual(dropped.manifest["excluded_missing_features"], ["g2"])
        with self.assertRaisesRegex(ValueError, "library totals"):
            prepare_study(self.matrix, self.meta, "T", "R", design="welch", missing="drop-features", transform="cpm-log2")
        self.matrix.write_text("feature_id\ta\tb\tc\td\te\tf\nhigh\t100\t100\t100\t100\t100\t100\nlow\t1\t0\t0\t0\t0\t0\n")
        prepared = prepare_study(self.matrix, self.meta, "T", "R", design="welch", transform="cpm-log2", min_cpm=10000)
        self.assertEqual(prepared.feature_ids, ("high",))
        self.assertEqual(prepared.manifest["library_totals"]["a"], 101)

    def test_study_cli_welch_report_and_prepare_only(self):
        output = self.root / "audit"
        argv = ["study", str(self.matrix), "--metadata", str(self.meta), "--target", "T", "--reference", "R",
                "--design", "welch", "--budget", "2", "--top-k", "1", "--output", str(output)]
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(argv), 0)
        report = json.loads((output / "audit.json").read_text())
        self.assertEqual(report["preparation"]["contrast"], "T - R")
        self.assertEqual(report["delete_one_summary"]["n_scenarios"], 6)
        self.assertIn(report["status"], {"CERTIFIED", "REFUTED"})
        self.assertIn("Candidate rank sensitivity", (output / "report.html").read_text())
        prep = self.root / "prepared"
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(["prepare", str(self.matrix), "--metadata", str(self.meta), "--target", "T",
                "--reference", "R", "--design", "welch", "--output", str(prep)]), 0)
        self.assertTrue((prep / "values.csv").is_file())
        self.assertFalse((prep / "audit.json").exists())
        from extremarank.cli import read_effects
        with self.assertRaisesRegex(ValueError, "donor_id"):
            read_effects(prep / "values.csv")

    def test_prepare_refuses_to_overwrite_source_input(self):
        source = self.root / "values.csv"
        source.write_text("sample_id,g1,g2\na,4,1\nb,6,1\nc,5,2\nd,1,3\ne,2,4\nf,3,3\n")
        old = source.read_bytes()
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            main(["prepare", str(source), "--metadata", str(self.meta), "--target", "T", "--reference", "R",
                  "--design", "welch", "--orientation", "samples", "--output", str(self.root)])
        self.assertEqual(old, source.read_bytes())

    def test_study_cli_paired_and_small_group_diagnostics(self):
        self.meta.write_text("sample_id\tgroup\tdonor_id\na\tT\t1\nb\tT\t2\nc\tT\t3\nd\tR\t1\ne\tR\t2\nf\tR\t3\n")
        out = self.root / "paired"
        with contextlib.redirect_stdout(io.StringIO()):
            main(["study", str(self.matrix), "--metadata", str(self.meta), "--target", "T", "--reference", "R",
                  "--budget", "1", "--top-k", "1", "--output", str(out)])
        report = json.loads((out / "audit.json").read_text())
        self.assertEqual(report["delete_one_summary"]["n_scenarios"], 3)
        self.assertTrue((out / "envelopes.tsv").exists())

    def test_sample_orientation(self):
        transposed = self.root / "transpose.tsv"
        transposed.write_text("sample_id\tg1\tg2\na\t4\t1\nb\t6\t1\nc\t5\t2\nd\t1\t3\ne\t2\t4\nf\t3\t3\n")
        self.assertEqual(prepare_study(self.matrix, self.meta, "T", "R", design="welch").values,
                         prepare_study(transposed, self.meta, "T", "R", design="welch", orientation="samples").values)


class ExternalTests(unittest.TestCase):
    def test_supplied_tables_never_yield_certificate(self):
        baseline = {"a": 5, "b": 3, "c": 1}
        report = compare_rankings(baseline, {"d1": {"a": 1, "b": 4, "c": 2}}, {"d1": ["D1"]}, 1)
        self.assertEqual(report["status"], "OBSERVED_CHANGED")
        self.assertEqual(report["scenarios"][0]["entered_features"], ["b"])
        stable = compare_rankings(baseline, {"d1": baseline}, {"d1": ["D1"]}, 1)
        self.assertEqual(stable["status"], "OBSERVED_STABLE")
        self.assertIn("unobserved", stable["scope"])
        with self.assertRaises(ValueError):
            compare_rankings(baseline, {"d1": {"a": 5}}, {"d1": ["D1"]}, 1)
        with self.assertRaises(ValueError):
            compare_rankings(baseline, {"d1": baseline}, {"d1": []}, 1)

    def test_external_cli_escapes_labels_and_records_method(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/"baseline.tsv").write_text('feature_id\tscore\n<script>\t1\nb\t0\n')
            (root/"refits.tsv").write_text('scenario_id\tfeature_id\tscore\ns1\t<script>\t0\ns1\tb\t1\n')
            # Use CSV writer for the JSON field to avoid ambiguous quoting.
            import csv
            with (root/"manifest.tsv").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["scenario_id", "deleted_units"])
                writer.writerow(["s1", json.dumps(["D1"])])
            with contextlib.redirect_stdout(io.StringIO()):
                main(["compare", str(root/"baseline.tsv"), str(root/"refits.tsv"), "--manifest", str(root/"manifest.tsv"),
                      "--top-k", "1", "--method", "custom refit t", "--output", str(root/"out")])
            report = json.loads((root/"out/comparison.json").read_text())
            self.assertEqual(report["method"], "custom refit t")
            html = (root/"out/report.html").read_text()
            self.assertNotIn("<script>", html)
            self.assertIn("&lt;script&gt;", html)


if __name__ == "__main__":
    unittest.main()
