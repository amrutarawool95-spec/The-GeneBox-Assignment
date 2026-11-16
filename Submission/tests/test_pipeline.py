"""Tests for src/problemC.py. Run from the submission/ folder:  python3 -m unittest tests.test_pipeline -v
Standard library only."""
import csv, json, subprocess, sys, tempfile, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "src" / "problemC.py"
INPUT = ROOT / "problemC_somatic"
SAVED = ROOT / "out"


def run_pipeline(out_dir):
    subprocess.run([sys.executable, str(SCRIPT), str(INPUT), str(out_dir)],
                   check=True, capture_output=True)


def calls(out_dir):
    with open(Path(out_dir) / "problemC_calls.tsv", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class TestProblemC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.out = Path(cls.tmp.name)
        run_pipeline(cls.out)
        cls.rows = calls(cls.out)
        cls.by_key = {(r["chrom"], r["pos"], r["ref"], r["alt"]): r for r in cls.rows}
        cls.summary = json.loads((cls.out / "problemC_summary.json").read_text())

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_all_48_sites_and_columns(self):
        self.assertEqual(len(self.rows), 48)
        self.assertEqual(list(self.rows[0].keys()),
                         ["chrom", "pos", "ref", "alt", "depth", "alt_count", "vaf", "normal_vaf",
                          "af_used", "expected_contam_alt", "filter", "filter_rank"])

    def test_one_filter_per_site_and_rank_matches(self):
        order = ["off_target", "low_depth", "min_alt", "germline", "panel_of_normals", "normal_artifact",
                 "contamination", "orientation_bias", "strand_bias", "mapping_quality", "base_quality",
                 "short_tandem_repeat", "read_position", "PASS"]
        for r in self.rows:
            self.assertEqual(order.index(r["filter"]) + 1, int(r["filter_rank"]))

    def test_deterministic(self):
        with tempfile.TemporaryDirectory() as d2:
            run_pipeline(d2)
            for name in ("problemC_calls.tsv", "problemC_summary.json"):
                self.assertEqual((self.out / name).read_bytes(), (Path(d2) / name).read_bytes())

    def test_matches_saved_outputs(self):
        # the calls table must match byte for byte; the JSON is compared by content
        # (so a trailing newline added by an editor or by GitHub does not matter)
        self.assertEqual((self.out / "problemC_calls.tsv").read_bytes(),
                         (SAVED / "problemC_calls.tsv").read_bytes())
        self.assertEqual(json.loads((self.out / "problemC_summary.json").read_text()),
                         json.loads((SAVED / "problemC_summary.json").read_text()))

    def test_bed_is_zero_based_half_open(self):
        # 1-based pos equal to a BED start is OUTSIDE; start+1 is inside; end is inside; end+1 is outside
        self.assertEqual(self.by_key[("chr3", "2500000", "G", "A")]["filter"], "off_target")
        self.assertNotEqual(self.by_key[("chr3", "2500001", "A", "G")]["filter"], "off_target")
        self.assertNotEqual(self.by_key[("chr1", "1400000", "G", "A")]["filter"], "off_target")
        self.assertEqual(self.by_key[("chr1", "1400001", "T", "C")]["filter"], "off_target")

    def test_exact_allele_join(self):
        # resource has G>C here, pileup has G>A: must NOT inherit the 0.47 AF
        self.assertEqual(self.by_key[("chr17", "7577538", "G", "A")]["af_used"], "0.00000100")
        # exact match keeps its AF
        self.assertEqual(self.by_key[("chr17", "7577120", "C", "T")]["af_used"], "0.47000000")

    def test_two_alts_at_one_position_are_separate(self):
        self.assertEqual(self.by_key[("chr7", "5350700", "C", "T")]["filter"], "PASS")
        self.assertEqual(self.by_key[("chr7", "5350700", "C", "G")]["filter"], "germline")

    def test_priority_order(self):
        # in the panel of normals (3 of 50) but tumour depth 42 < 50: low_depth wins (rank 2 before 5)
        self.assertEqual(self.by_key[("chr1", "1050120", "C", "T")]["filter"], "low_depth")

    def test_summary_consistency(self):
        s = self.summary
        self.assertEqual(s["n_sites_evaluated"], 48)
        self.assertEqual(sum(s["filter_counts"].values()), 48)
        self.assertEqual(s["n_pass"], s["filter_counts"]["PASS"])
        self.assertEqual(s["tmb_per_mb"], round(s["n_pass"] / 1.18644, 3))


if __name__ == "__main__":
    unittest.main()
