"""
Unit tests for ganon.reassign.

These tests exercise the module functions directly with small, self-contained
fixtures (no database is built and no external data is downloaded), in the same
spirit as the standalone unit tests under tests/ganon-build/.
"""

import io
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from types import SimpleNamespace

from ganon.config import Config
from ganon.reassign import (
    HARD_CAP_ITER,
    define_all_files,
    generate_rep,
    get_top_match,
    parse_all_file,
    parse_rep,
    reassign,
    run_em,
    write_one,
    write_rep,
)


SINGLE_REP = (
    "1\tA\t3\t2\t1\tS\tSpecies A\n"
    "1\tB\t1\t0\t0\tS\tSpecies B\n"
    "#total_classified\t3\n"
    "#total_unclassified\t1\n"
)

SINGLE_ALL = "r1\tA\t10\nr2\tA\t10\nr3\tA\t10\nr3\tB\t1\n"


class ReassignTestCase(unittest.TestCase):
    """Common fixture handling for the reassign unit tests."""

    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.dir = Path(self._tmpdir.name)

    def tearDown(self):
        self._tmpdir.cleanup()

    def write(self, name, text):
        path = self.dir / name
        path.write_text(text)
        return path

    def read(self, name):
        return (self.dir / name).read_text()

    def exists(self, name):
        return (self.dir / name).is_file()

    def em_cfg(self, **kwargs):
        """Minimal config object accepted by run_em."""
        values = {
            "quiet": True,
            "threshold": 0,
            "max_iter": 10,
            "skip_one": False,
        }
        values.update(kwargs)
        return SimpleNamespace(**values)

    def run_reassign(self, **kwargs):
        """Run ganon.reassign.reassign and capture its standard error."""
        params = {"input": str(self.dir / "sample.rep"), "quiet": True}
        params.update(kwargs)
        cfg = Config("reassign", **params)
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            result = reassign(cfg)
        return result, stderr.getvalue()

    def read_rep(self, name):
        """Return (rows, totals) from a .rep file, keyed by (hierarchy, target)."""
        rows = {}
        totals = []
        for line in self.read(name).splitlines():
            if line.startswith("#"):
                totals.append(line)
            else:
                fields = line.split("\t")
                rows[(fields[0], fields[1])] = fields
        return rows, totals

    def create_single_fixture(self):
        """Consistent .rep/.all pair with two unique reads and one multi read."""
        rep = self.write("sample.rep", SINGLE_REP)
        self.write("sample.all", SINGLE_ALL)
        return rep


class TestParseRep(ReassignTestCase):
    def test_hierarchies_are_sorted_and_unique(self):
        """.rep hierarchies are returned sorted and deduplicated"""
        rep = self.write(
            "sample.rep",
            "2\tB\t1\t1\t0\tS\tSpecies B\n"
            "1\tA\t1\t1\t0\tS\tSpecies A\n"
            "2\tC\t1\t1\t0\tS\tSpecies C\n",
        )
        hierarchies, totals = parse_rep(rep)
        self.assertEqual(hierarchies, ["1", "2"], "hierarchies not sorted/unique")
        self.assertEqual(totals, [], "unexpected totals")

    def test_total_lines_are_kept_in_order(self):
        """.rep totals lines are preserved verbatim and in order"""
        rep = self.write(
            "sample.rep",
            "#total_classified\t10\n#total_unclassified\t2\n",
        )
        hierarchies, totals = parse_rep(rep)
        self.assertEqual(hierarchies, [], "unexpected hierarchies")
        self.assertEqual(
            totals,
            [["#total_classified\t10"], ["#total_unclassified\t2"]],
            "totals lines not preserved",
        )

    def test_blank_lines_are_skipped(self):
        """Empty lines do not break parsing"""
        rep = self.write(
            "sample.rep",
            "\n1\tA\t1\t1\t0\tS\tSpecies A\n\n#total_classified\t1\n\n",
        )
        hierarchies, totals = parse_rep(rep)
        self.assertEqual(hierarchies, ["1"], "blank lines produced wrong hierarchies")
        self.assertEqual(
            totals, [["#total_classified\t1"]], "blank lines produced wrong totals"
        )


class TestDefineAllFiles(ReassignTestCase):
    def test_per_hierarchy_all_files(self):
        """One .all per hierarchy is detected"""
        self.write("sample.1.all", "r1\tA\t10\n")
        self.write("sample.2.all", "r1\tB\t10\n")
        all_files = define_all_files(str(self.dir / "sample"), ["1", "2"], True)
        self.assertEqual(
            all_files,
            {
                "1": str(self.dir / "sample.1.all"),
                "2": str(self.dir / "sample.2.all"),
            },
            "per-hierarchy .all files not detected",
        )

    def test_single_all_file_for_one_hierarchy(self):
        """A concatenated .all is used when a single hierarchy is present"""
        self.write("sample.all", "r1\tA\t10\n")
        all_files = define_all_files(str(self.dir / "sample"), ["1"], True)
        self.assertEqual(
            all_files,
            {"": str(self.dir / "sample.all")},
            "single .all file not detected",
        )

    def test_single_all_file_with_multiple_hierarchies_is_unsupported(self):
        """A concatenated .all cannot be split between several hierarchies"""
        self.write("sample.all", "r1\tA\t10\n")
        all_files = define_all_files(str(self.dir / "sample"), ["1", "2"], True)
        self.assertIsNone(all_files, "concatenated .all should be rejected")

    def test_missing_all_file(self):
        """No .all file at all is reported as an error"""
        all_files = define_all_files(str(self.dir / "sample"), ["1"], True)
        self.assertIsNone(all_files, "missing .all not detected")

    def test_partial_per_hierarchy_all_files(self):
        """A missing .all for one of the hierarchies is reported as an error"""
        self.write("sample.1.all", "r1\tA\t10\n")
        all_files = define_all_files(str(self.dir / "sample"), ["1", "2"], True)
        self.assertIsNone(all_files, "missing hierarchy .all not detected")

    def test_quiet_suppresses_messages(self):
        """--quiet suppresses the missing file message"""
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            define_all_files(str(self.dir / "sample"), ["1"], True)
        self.assertEqual(stderr.getvalue(), "", "message printed in quiet mode")

        stderr = io.StringIO()
        with redirect_stderr(stderr):
            define_all_files(str(self.dir / "sample"), ["1"], False)
        self.assertIn("Could not find matching file", stderr.getvalue())


class TestParseAllFile(ReassignTestCase):
    def test_matches_are_grouped_by_read(self):
        """.all rows are grouped per read identifier"""
        all_file = self.write(
            "sample.all",
            "r1\tA\t10\nr1\tB\t5\nr2\tA\t3\n",
        )
        read_matches, targets = parse_all_file(all_file, "1")

        self.assertEqual(list(read_matches), ["r1", "r2"], "wrong reads parsed")
        self.assertEqual(
            read_matches["r1"],
            [(targets[("1", "A")], 10), (targets[("1", "B")], 5)],
            "matches not grouped by read",
        )
        self.assertEqual(
            read_matches["r2"],
            [(targets[("1", "A")], 3)],
            "matches not grouped by read",
        )

    def test_kcount_is_converted_to_integer(self):
        """The k-mer count is stored as an integer"""
        all_file = self.write("sample.all", "r1\tA\t10\n")
        read_matches, targets = parse_all_file(all_file, "1")
        target, kcount = read_matches["r1"][0]
        self.assertEqual(kcount, 10)
        self.assertIsInstance(kcount, int, "k-count is not an integer")
        self.assertEqual(target, targets[("1", "A")])

    def test_hierarchy_is_part_of_target_identity(self):
        """The same target name in different hierarchies gets distinct ids"""
        all_1 = self.write("sample.1.all", "r1\tA\t10\n")
        all_2 = self.write("sample.2.all", "r1\tA\t10\n")
        _, targets_1 = parse_all_file(all_1, "1")
        _, targets_2 = parse_all_file(all_2, "2")

        self.assertIn(("1", "A"), targets_1)
        self.assertNotIn(("1", "A"), targets_2)
        self.assertIn(("2", "A"), targets_2)

    def test_extra_columns_are_ignored(self):
        """.all rows with extra (future) columns are still parsed"""
        all_file = self.write("sample.all", "r1\tA\t10\t42\n")
        read_matches, targets = parse_all_file(all_file, "1")
        self.assertEqual(read_matches["r1"], [(targets[("1", "A")], 10)])

    def test_blank_lines_are_skipped(self):
        """Empty lines in .all do not break parsing"""
        all_file = self.write("sample.all", "\n\nr1\tA\t10\n\n")
        read_matches, _ = parse_all_file(all_file, "1")
        self.assertEqual(list(read_matches), ["r1"])


class TestRunEm(ReassignTestCase):
    def test_unique_reads_are_kept(self):
        """Reads with a single match are counted for that target only"""
        targets = {("1", "A"): 0, ("1", "B"): 1}
        read_matches = {"r1": [(0, 10)], "r2": [(0, 10)], "r3": [(1, 10)]}
        reassigned, top_matches = run_em(self.em_cfg(), read_matches, targets)

        self.assertEqual(reassigned, {0: 2, 1: 1}, "unique reads miscounted")
        self.assertEqual(top_matches, {}, "unique reads should not be reassigned")

    def test_multi_reads_are_assigned_once(self):
        """Each multi-matching read contributes exactly one count"""
        targets = {("1", "A"): 0, ("1", "B"): 1}
        read_matches = {
            "r1": [(0, 10)],
            "r2": [(0, 10)],
            "r3": [(0, 10), (1, 10)],
            "r4": [(1, 10)],
        }
        reassigned, top_matches = run_em(self.em_cfg(), read_matches, targets)

        self.assertEqual(reassigned, {0: 3, 1: 1}, "wrong EM assignment")
        self.assertEqual(top_matches, {"r3": (0, 10)}, "wrong .one assignment")
        self.assertEqual(
            sum(reassigned.values()), len(read_matches), "reads were lost/duplicated"
        )

    def test_max_iter_is_respected(self):
        """--max-iter stops the EM loop after the requested number of iterations"""
        targets = {("1", "A"): 0, ("1", "B"): 1}
        read_matches = {
            "r1": [(0, 1)],
            "r2": [(1, 1)],
            "r3": [(0, 1), (1, 100)],
        }
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            run_em(self.em_cfg(quiet=False, max_iter=2), read_matches, targets)

        iterations = [
            line for line in stderr.getvalue().splitlines() if "Iteration" in line
        ]
        self.assertEqual(len(iterations), 2, "--max-iter not respected")

    def test_kcount_breaks_a_probability_tie(self):
        """With equal probabilities the match with the highest k-count wins"""
        targets = {("1", "A"): 0, ("1", "B"): 1}
        read_matches = {
            "r1": [(0, 10)],
            "r2": [(1, 10)],
            "r3": [(0, 1), (1, 100)],
        }
        reassigned, top_matches = run_em(self.em_cfg(), read_matches, targets)

        self.assertEqual(reassigned, {0: 1, 1: 2}, "k-count was not considered")
        self.assertEqual(top_matches, {"r3": (1, 100)})

    def test_target_prior_wins_when_kcounts_are_equal(self):
        """A target with more unique evidence attracts shared reads"""
        targets = {("1", "A"): 0, ("1", "B"): 1}
        read_matches = {f"r{i}": [(0, 10)] for i in range(5)}
        read_matches["r5"] = [(1, 10)]
        read_matches["r6"] = [(0, 10), (1, 10)]

        reassigned, top_matches = run_em(self.em_cfg(), read_matches, targets)

        self.assertEqual(reassigned, {0: 6, 1: 1}, "prior was not used")
        self.assertEqual(top_matches, {"r6": (0, 10)})

    def test_counts_are_integers(self):
        """The final assignment is a hard, integer assignment"""
        targets = {("1", "A"): 0, ("1", "B"): 1}
        read_matches = {"r1": [(0, 10), (1, 10)], "r2": [(0, 10)]}
        reassigned, _ = run_em(self.em_cfg(), read_matches, targets)

        for target, count in reassigned.items():
            self.assertIsInstance(count, int, f"count for {target} is not an integer")

    def test_all_reads_multimatching(self):
        """No unique reads: priors are uniform and the run does not crash"""
        targets = {("1", "A"): 0, ("1", "B"): 1}
        read_matches = {
            "r1": [(0, 10), (1, 10)],
            "r2": [(0, 10), (1, 10)],
        }
        reassigned, top_matches = run_em(self.em_cfg(), read_matches, targets)

        self.assertEqual(sum(reassigned.values()), 2, "reads were lost")
        self.assertEqual(len(top_matches), 2, "multi reads not reassigned")

    def test_skip_one_does_not_fill_top_matches(self):
        """--skip-one avoids storing the per-read assignment"""
        targets = {("1", "A"): 0, ("1", "B"): 1}
        read_matches = {"r1": [(0, 10), (1, 10)], "r2": [(0, 10)]}
        reassigned, top_matches = run_em(
            self.em_cfg(skip_one=True), read_matches, targets
        )

        self.assertEqual(sum(reassigned.values()), 2)
        self.assertEqual(top_matches, {}, "top matches computed with --skip-one")

    def test_max_iter_zero_runs_until_convergence(self):
        """--max-iter 0 iterates until convergence, without hitting the hard cap"""
        targets = {("1", "A"): 0, ("1", "B"): 1}
        read_matches = {
            "r1": [(0, 1)],
            "r2": [(1, 1)],
            "r3": [(0, 1), (1, 100)],
        }
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            run_em(
                self.em_cfg(quiet=False, max_iter=0, threshold=1e-9),
                read_matches,
                targets,
            )
        log = stderr.getvalue()

        iterations = [line for line in log.splitlines() if "Iteration" in line]
        self.assertGreaterEqual(len(iterations), 1, "EM did not run")
        self.assertLess(len(iterations), HARD_CAP_ITER, "EM did not converge")
        self.assertNotIn("Reach max. default iterations", log)


class TestGetTopMatch(ReassignTestCase):
    def test_highest_score_wins(self):
        """The match with the highest probability*k-count wins"""
        # A higher probability can outweigh a larger k-count
        self.assertEqual(get_top_match([(0, 10), (1, 20)], {0: 0.9, 1: 0.1}), (0, 10))
        # A much larger k-count can outweigh a lower probability
        self.assertEqual(get_top_match([(0, 1), (1, 100)], {0: 0.9, 1: 0.1}), (1, 100))

    def test_kcount_breaks_equal_score(self):
        """Equal scores are broken by the highest k-count"""
        self.assertEqual(get_top_match([(0, 10), (1, 20)], {0: 0.5, 1: 0.5}), (1, 20))

    def test_first_match_wins_on_full_tie(self):
        """A full tie keeps the first match (current behaviour)"""
        self.assertEqual(get_top_match([(0, 10), (1, 10)], {0: 0.5, 1: 0.5}), (0, 10))


class TestWriteOne(ReassignTestCase):
    def test_one_row_per_read(self):
        """The .one file has one row per read, in .all order"""
        all_file = self.write("sample.all", SINGLE_ALL)
        read_matches, targets = parse_all_file(all_file, "1")
        _, top_matches = run_em(self.em_cfg(), read_matches, targets)

        write_one(str(self.dir / "sample.one"), read_matches, targets, top_matches)

        self.assertEqual(
            self.read("sample.one").splitlines(),
            ["r1\tA\t10", "r2\tA\t10", "r3\tA\t10"],
            "wrong .one output",
        )


class TestGenerateRep(ReassignTestCase):
    def test_rep_rows_keep_direct_matches_and_set_em_reads(self):
        """.rep keeps field 3 (direct matches) and writes assigned-unique as em"""
        rep = self.write(
            "sample.rep",
            "1\tA\t3\t2\t1\tS\tSpecies A\n"
            "1\tB\t1\t0\t0\tS\tSpecies B\n"
            "1\tC\t1\t0\t0\tS\tSpecies C\n",
        )
        all_file = self.write("sample.all", SINGLE_ALL)
        read_matches, targets = parse_all_file(all_file, "1")
        reassigned, _ = run_em(self.em_cfg(), read_matches, targets)

        rows = generate_rep(rep, "1", reassigned, targets)

        # C is a LCA-only row (not present in .all) and must be dropped
        self.assertEqual(len(rows), 2, "unexpected number of .rep rows")
        self.assertEqual(
            rows[0], ["1", "A", "3", 2, 1, "S", "Species A"], "wrong A row"
        )
        self.assertEqual(
            rows[1], ["1", "B", "1", 0, 0, "S", "Species B"], "wrong B row"
        )
        # direct matches must not be replaced by unique + em
        self.assertEqual(rows[0][2], "3")

    def test_rows_are_filtered_by_hierarchy(self):
        """When .all files are split, each .rep hierarchy is handled separately"""
        rep = self.write(
            "sample.rep",
            "1\tA\t2\t1\t1\tS\tSpecies A\n"
            "1\tC\t1\t0\t0\tS\tSpecies C\n"
            "2\tB\t1\t1\t0\tS\tSpecies B\n",
        )
        all_1 = self.write("sample.1.all", "r1\tA\t10\nr3\tA\t7\nr3\tC\t5\n")
        all_2 = self.write("sample.2.all", "r2\tB\t10\n")

        read_matches_1, targets_1 = parse_all_file(all_1, "1")
        reassigned_1, _ = run_em(self.em_cfg(), read_matches_1, targets_1)
        rows_1 = generate_rep(rep, "1", reassigned_1, targets_1)

        read_matches_2, targets_2 = parse_all_file(all_2, "2")
        reassigned_2, _ = run_em(self.em_cfg(), read_matches_2, targets_2)
        rows_2 = generate_rep(rep, "2", reassigned_2, targets_2)

        self.assertEqual(
            [(row[0], row[1], row[4]) for row in rows_1],
            [("1", "A", 1), ("1", "C", 0)],
            "wrong hierarchy 1 rows",
        )
        self.assertEqual(
            [(row[0], row[1], row[4]) for row in rows_2],
            [("2", "B", 0)],
            "wrong hierarchy 2 rows",
        )

    def test_short_rep_rows_have_empty_rank_and_name(self):
        """.rep rows without rank/name columns are tolerated"""
        rep = self.write("sample.rep", "1\tA\t3\t2\n")
        all_file = self.write("sample.all", SINGLE_ALL)
        read_matches, targets = parse_all_file(all_file, "1")
        reassigned, _ = run_em(self.em_cfg(), read_matches, targets)

        rows = generate_rep(rep, "1", reassigned, targets)
        self.assertEqual(rows[0][5:], ["", ""], "rank/name should be empty")


class TestWriteRep(ReassignTestCase):
    def test_rows_and_totals_are_written(self):
        """The .rep file contains the new rows followed by the totals lines"""
        write_rep(
            str(self.dir / "sample"),
            [["1", "A", "3", 2, 1, "S", "Species A"]],
            [["#total_classified\t3"], ["#total_unclassified\t1"]],
        )

        self.assertEqual(
            self.read("sample.rep").splitlines(),
            [
                "1\tA\t3\t2\t1\tS\tSpecies A",
                "#total_classified\t3",
                "#total_unclassified\t1",
            ],
            "wrong .rep output",
        )


class TestReassign(ReassignTestCase):
    def test_single_hierarchy(self):
        """ganon reassign with one .rep/.all pair"""
        self.create_single_fixture()
        result, _ = self.run_reassign(output_prefix=str(self.dir / "sample_RE"))

        self.assertTrue(result, "reassign failed")
        self.assertTrue(self.exists("sample_RE.one"), ".one was not created")
        self.assertTrue(self.exists("sample_RE.rep"), ".rep was not created")

        self.assertEqual(
            self.read("sample_RE.one").splitlines(),
            ["r1\tA\t10", "r2\tA\t10", "r3\tA\t10"],
            "wrong .one output",
        )
        rows, totals = self.read_rep("sample_RE.rep")
        self.assertEqual(rows[("1", "A")][2:5], ["3", "2", "1"], "wrong A row")
        self.assertEqual(rows[("1", "B")][2:5], ["1", "0", "0"], "wrong B row")
        self.assertEqual(
            totals,
            ["#total_classified\t3", "#total_unclassified\t1"],
            "totals were not preserved",
        )

    def test_direct_matches_are_not_reassigned(self):
        """.rep field 3 stays the pre-EM direct match count"""
        self.create_single_fixture()
        self.run_reassign(output_prefix=str(self.dir / "sample_RE"))
        rows, _ = self.read_rep("sample_RE.rep")

        for fields in rows.values():
            direct = int(fields[2])
            assigned = int(fields[3]) + int(fields[4])
            self.assertGreaterEqual(
                direct, assigned, "direct matches are not a superset of assignments"
            )

    def test_input_files_are_not_modified_with_output_prefix(self):
        """--output-prefix leaves the original .rep/.all untouched"""
        self.create_single_fixture()
        self.run_reassign(output_prefix=str(self.dir / "sample_RE"))

        self.assertEqual(self.read("sample.rep"), SINGLE_REP, ".rep was modified")
        self.assertEqual(self.read("sample.all"), SINGLE_ALL, ".all was modified")

    def test_default_output_prefix_overwrites_input_prefix(self):
        """Without --output-prefix the input .rep/.one are (re)written"""
        self.create_single_fixture()
        result, _ = self.run_reassign()

        self.assertTrue(result, "reassign failed")
        self.assertTrue(self.exists("sample.one"), ".one was not created")
        rows, totals = self.read_rep("sample.rep")
        self.assertIn(("1", "A"), rows, ".rep was not rewritten")
        self.assertEqual(len(totals), 2, "totals were not preserved")

    def test_skip_one(self):
        """--skip-one writes .rep only"""
        self.create_single_fixture()
        result, _ = self.run_reassign(
            output_prefix=str(self.dir / "sample_RE"), skip_one=True
        )

        self.assertTrue(result, "reassign failed")
        self.assertFalse(self.exists("sample_RE.one"), ".one should be skipped")
        self.assertTrue(self.exists("sample_RE.rep"), ".rep was not created")

    def test_skip_rep(self):
        """--skip-rep writes .one only"""
        self.create_single_fixture()
        result, _ = self.run_reassign(
            output_prefix=str(self.dir / "sample_RE"), skip_rep=True
        )

        self.assertTrue(result, "reassign failed")
        self.assertTrue(self.exists("sample_RE.one"), ".one was not created")
        self.assertFalse(self.exists("sample_RE.rep"), ".rep should be skipped")
        self.assertEqual(self.read("sample.rep"), SINGLE_REP, "input .rep was modified")

    def test_skip_one_and_skip_rep(self):
        """Skipping both outputs is an error"""
        self.create_single_fixture()
        result, stderr = self.run_reassign(
            output_prefix=str(self.dir / "sample_RE"),
            skip_one=True,
            skip_rep=True,
            quiet=False,
        )

        self.assertFalse(result, "reassign should fail with nothing to write")
        self.assertIn("No valid output", stderr)
        self.assertFalse(self.exists("sample_RE.one"))
        self.assertFalse(self.exists("sample_RE.rep"))

    def test_missing_input_file(self):
        """A missing .rep returns False instead of raising"""
        rep = self.create_single_fixture()
        cfg = Config(
            "reassign",
            input=str(rep),
            output_prefix=str(self.dir / "sample_RE"),
            quiet=False,
        )
        rep.unlink()

        stderr = io.StringIO()
        with redirect_stderr(stderr):
            result = reassign(cfg)

        self.assertFalse(result, "missing input should return False")
        self.assertIn("Report file not found", stderr.getvalue())

    def test_totals_only_rep(self):
        """A .rep without target rows returns False"""
        self.write("sample.rep", "#total_classified\t0\n#total_unclassified\t5\n")
        self.write("sample.all", "r1\tA\t10\n")
        result, stderr = self.run_reassign(
            output_prefix=str(self.dir / "sample_RE"), quiet=False
        )

        self.assertFalse(result, "empty .rep should return False")
        self.assertIn("No valid report entries", stderr)

    def test_missing_all_file(self):
        """A .rep without its .all file returns False"""
        self.write("sample.rep", SINGLE_REP)
        result, _ = self.run_reassign(output_prefix=str(self.dir / "sample_RE"))

        self.assertFalse(result, "missing .all should return False")
        self.assertFalse(self.exists("sample_RE.one"))
        self.assertFalse(self.exists("sample_RE.rep"))

    def test_remove_all(self):
        """--remove-all deletes the .all file after processing"""
        self.create_single_fixture()
        result, _ = self.run_reassign(
            output_prefix=str(self.dir / "sample_RE"), remove_all=True
        )

        self.assertTrue(result, "reassign failed")
        self.assertFalse(self.exists("sample.all"), ".all was not removed")
        self.assertTrue(self.exists("sample_RE.one"), ".one was not created")

    def test_multiple_hierarchies_split_files(self):
        """Split .all files produce one .one per hierarchy and a combined .rep"""
        self.write(
            "multi.rep",
            "1\tA\t2\t1\t1\tS\tSpecies A\n"
            "1\tC\t1\t0\t0\tS\tSpecies C\n"
            "2\tB\t1\t1\t0\tS\tSpecies B\n"
            "#total_classified\t3\n"
            "#total_unclassified\t0\n",
        )
        self.write("multi.1.all", "r1\tA\t10\nr3\tA\t7\nr3\tC\t5\n")
        self.write("multi.2.all", "r2\tB\t10\n")

        params = {"input": str(self.dir / "multi.rep"), "quiet": True}
        params["output_prefix"] = str(self.dir / "multi_RE")
        cfg = Config("reassign", **params)
        with redirect_stderr(io.StringIO()):
            result = reassign(cfg)

        self.assertTrue(result, "reassign failed")
        self.assertTrue(self.exists("multi_RE.1.one"), "hierarchy 1 .one missing")
        self.assertTrue(self.exists("multi_RE.2.one"), "hierarchy 2 .one missing")
        self.assertFalse(self.exists("multi_RE.one"), "unexpected combined .one")

        self.assertEqual(
            self.read("multi_RE.1.one").splitlines(),
            ["r1\tA\t10", "r3\tA\t7"],
            "wrong hierarchy 1 .one",
        )
        self.assertEqual(
            self.read("multi_RE.2.one").splitlines(),
            ["r2\tB\t10"],
            "wrong hierarchy 2 .one",
        )

        rows, totals = self.read_rep("multi_RE.rep")
        self.assertEqual(rows[("1", "A")][2:5], ["2", "1", "1"])
        self.assertEqual(rows[("1", "C")][2:5], ["1", "0", "0"])
        self.assertEqual(rows[("2", "B")][2:5], ["1", "1", "0"])
        self.assertEqual(len(totals), 2, "totals were not preserved")

    def test_multiple_hierarchies_single_all_file(self):
        """A concatenated .all cannot be used with several hierarchies"""
        self.write(
            "sample.rep",
            "1\tA\t1\t1\t0\tS\tSpecies A\n2\tB\t1\t1\t0\tS\tSpecies B\n",
        )
        self.write("sample.all", "r1\tA\t10\nr2\tB\t10\n")
        result, stderr = self.run_reassign(
            output_prefix=str(self.dir / "sample_RE"), quiet=False
        )

        self.assertFalse(result, "concatenated .all should be rejected")
        self.assertIn("not supported", stderr)
        self.assertFalse(self.exists("sample_RE.rep"))


if __name__ == "__main__":
    unittest.main()
