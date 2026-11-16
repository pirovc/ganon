import contextlib
import io
import unittest

from multitax import CustomTx

from ganon.report import correct_genome_size


class TestCorrectGenomeSize(unittest.TestCase):
    """
    Tests for the global genome size correction in ganon.report.

    Trees are built in memory so expected values can be verified by hand.
    """

    def build_tax(self, tree):
        tax = CustomTx(
            undefined_node="",
            undefined_rank="na",
            undefined_name="na",
            root_node="1",
            root_name="root",
            root_rank="root",
            empty=True,
        )
        for node, parent, name, rank in tree:
            tax.add(node, parent, name, rank)
        tax.check_consistency()
        tax.build_lineages()
        return tax

    def run_correction(
        self, tree, target_counts, genome_sizes, ranks_abundance=None, quiet=False
    ):
        """
        Run correct_genome_size, capturing the warning messages
        """
        if ranks_abundance is None:
            ranks_abundance = ["all"]
        tax = self.build_tax(tree)
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            result = correct_genome_size(
                target_counts, genome_sizes, tax, ranks_abundance, quiet
            )
        return result, tax, stderr.getvalue()

    def assert_tree_values(self, tree, expected, places=4):
        for node, value in expected.items():
            self.assertIn(node, tree, "node " + node + " not in corrected tree")
            self.assertAlmostEqual(
                tree[node],
                value,
                places=places,
                msg="wrong corrected value for node " + node,
            )

    @staticmethod
    def simple_example():
        """
        root
        +-- G1 (genus)        100 reads, 20 Mb
        +-- G2 (genus)
            +-- S (species)   100 reads,  2 Mb
        """
        tree = [
            ("7", "1", "G1", "genus"),
            ("11", "1", "G2", "genus"),
            ("12", "11", "S", "species"),
        ]
        target_counts = {"7": 100, "12": 100}
        genome_sizes = {"1": 1000, "7": 20, "12": 2}
        return tree, target_counts, genome_sizes

    @staticmethod
    def complex_example():
        """
        root
        +-- D (domain)
            +-- P (phylum)
                +-- C (class)
                    +-- O (order)
                        +-- F1 (family)
                        |   +-- SF (subfamily)   40 reads, 100 Mb
                        |       +-- G1 (genus)
                        |       |   +-- S1 (species)   20 reads,  10 Mb
                        |       |   |   +-- T1 (strain) 120 reads,  5 Mb
                        |       |   |   +-- T2 (strain)  30 reads, 15 Mb
                        |       |   +-- S2 (species)   80 reads, 20 Mb
                        |       +-- G2 (genus)    50 reads, 25 Mb
                        |           +-- S3 (species) 10 reads,  5 Mb
                        +-- F2 (family)          30 reads, 200 Mb
                            +-- G4 (genus)
                                +-- S5 (species) 70 reads, 10 Mb
        """
        tree = [
            ("2", "1", "D", "domain"),
            ("3", "2", "P", "phylum"),
            ("4", "3", "C", "class"),
            ("5", "4", "O", "order"),
            ("6", "5", "F1", "family"),
            ("7", "6", "SF", "subfamily"),
            ("8", "7", "G1", "genus"),
            ("9", "8", "S1", "species"),
            ("10", "9", "T1", "strain"),
            ("11", "9", "T2", "strain"),
            ("12", "8", "S2", "species"),
            ("13", "7", "G2", "genus"),
            ("14", "13", "S3", "species"),
            ("16", "5", "F2", "family"),
            ("17", "16", "G4", "genus"),
            ("18", "17", "S5", "species"),
        ]
        target_counts = {
            "10": 120,  # T1
            "11": 30,  # T2
            "9": 20,  # S1 direct
            "12": 80,  # S2
            "7": 40,  # SF
            "13": 50,  # G2 direct
            "14": 10,  # S3
            "16": 30,  # F2 direct
            "18": 70,  # S5
        }
        genome_sizes = {
            "1": 1000,
            "6": 100,  # F1
            "7": 100,  # SF
            "9": 10,  # S1
            "10": 5,  # T1
            "11": 15,  # T2
            "12": 20,  # S2
            "13": 25,  # G2
            "14": 5,  # S3
            "16": 200,  # F2
            "18": 10,  # S5
        }
        return tree, target_counts, genome_sizes

    @staticmethod
    def proxy_example(strain_size=None):
        """
        root
        +-- G (genus)          50 Mb
            +-- S1 (species)    10 Mb
            |   +-- T (strain)
            +-- S2 (species)    20 Mb
        """
        tree = [
            ("7", "1", "G", "genus"),
            ("8", "7", "S1", "species"),
            ("9", "8", "T", "strain"),
            ("10", "7", "S2", "species"),
        ]
        target_counts = {"9": 100, "10": 100}
        genome_sizes = {"1": 1000, "7": 50, "8": 10, "10": 20}
        if strain_size is not None:
            genome_sizes["9"] = strain_size
        return tree, target_counts, genome_sizes

    def test_global_normalization_simple(self):
        tree, target_counts, genome_sizes = self.simple_example()
        result, _, stderr = self.run_correction(
            tree, target_counts, genome_sizes, ["all"]
        )
        # total ratio: 100/20 + 100/2 = 55
        self.assert_tree_values(
            result,
            {
                "1": 200.0,
                "7": 18.1818,  # G1: 200 * (100/20) / 55
                "11": 181.8182,  # G2: cumulative from S
                "12": 181.8182,  # S: 200 * (100/2) / 55
            },
        )
        self.assertEqual(stderr, "", "unexpected warning for valid genome sizes")

    def test_global_normalization_complex(self):
        tree, target_counts, genome_sizes = self.complex_example()
        result, _, stderr = self.run_correction(
            tree, target_counts, genome_sizes, ["all"]
        )
        self.assert_tree_values(
            result,
            {
                "1": 450.0,
                "2": 450.0,
                "3": 450.0,
                "4": 450.0,
                "5": 450.0,
                "10": 247.9908,  # T1: 120/5
                "11": 20.6659,  # T2: 30/15
                "9": 289.3226,  # S1 cumulative: own + T1 + T2
                "12": 41.3318,  # S2: 80/20
                "14": 20.6659,  # S3: 10/5
                "18": 72.3307,  # S5: 70/10
                "8": 330.6544,  # G1: S1 + S2
                "13": 41.3318,  # G2: direct + S3
                "7": 376.1194,  # SF: own + G1 + G2
                "6": 376.1194,  # F1: SF only
                "16": 73.8806,  # F2: own + S5
                "17": 72.3307,  # G4: S5 only
            },
        )
        self.assertEqual(stderr, "", "unexpected warning for valid genome sizes")

    def test_total_is_preserved(self):
        for example in (self.simple_example, self.complex_example):
            tree, target_counts, genome_sizes = example()
            result, tax, _ = self.run_correction(
                tree, target_counts, genome_sizes, ["all"]
            )
            self.assertAlmostEqual(
                result[tax.root_node],
                sum(target_counts.values()),
                places=4,
                msg="total count was not preserved",
            )

    def test_parent_size_used_when_target_size_missing(self):
        # T (strain) has no valid size (1), so S1 size (10) is used
        tree, target_counts, genome_sizes = self.proxy_example(strain_size=1)
        result, _, stderr = self.run_correction(
            tree, target_counts, genome_sizes, ["all"]
        )
        # total ratio: 100/10 (T via S1) + 100/20 (S2) = 15
        self.assert_tree_values(
            result,
            {
                "1": 200.0,
                "9": 133.3333,  # T: 200 * 10 / 15
                "10": 66.6667,  # S2: 200 * 5 / 15
            },
        )
        self.assertEqual(stderr, "", "parent proxy should not warn")

    def test_root_size_used_when_no_parent_size(self):
        tree = [
            ("7", "1", "G", "genus"),
            ("8", "7", "S", "species"),
        ]
        # No size for G nor S: root size is used
        result, tax, stderr = self.run_correction(
            tree, {"8": 100}, {"1": 1000}, ["all"]
        )
        self.assert_tree_values(result, {"1": 100.0, "8": 100.0})
        self.assertIn(
            "1 target(s) were corrected with the root node genome size",
            stderr,
            "missing root fallback warning",
        )

    def test_root_fallback_warning_counts_targets(self):
        tree = [
            ("7", "1", "G", "genus"),
            ("8", "7", "S1", "species"),
            ("9", "7", "S2", "species"),
        ]
        result, _, stderr = self.run_correction(
            tree, {"8": 100, "9": 100}, {"1": 1000}, ["all"]
        )
        # both targets use root, so they are corrected equally
        self.assert_tree_values(result, {"8": 100.0, "9": 100.0})
        self.assertIn(
            "2 target(s) were corrected with the root node genome size",
            stderr,
            "wrong number of targets in root fallback warning",
        )

    def test_ranks_abundance_limits_correction(self):
        tree, target_counts, genome_sizes = self.proxy_example(strain_size=5)

        # All ranks: T uses its own size (5), S2 uses 20
        result, _, _ = self.run_correction(tree, target_counts, genome_sizes, ["all"])
        self.assert_tree_values(result, {"9": 160.0, "10": 40.0})

        # species only: T uses S1 (10), S2 uses 20
        result, _, _ = self.run_correction(
            tree, target_counts, genome_sizes, ["species"]
        )
        self.assert_tree_values(result, {"9": 133.3333, "10": 66.6667})

        # genus only: both use G (50)
        result, _, _ = self.run_correction(tree, target_counts, genome_sizes, ["genus"])
        self.assert_tree_values(result, {"9": 100.0, "10": 100.0})

    def test_ranks_abundance_all_overrides_other_ranks(self):
        tree, target_counts, genome_sizes = self.proxy_example(strain_size=5)
        result, _, _ = self.run_correction(
            tree, target_counts, genome_sizes, ["all", "species"]
        )
        # "all" disables the rank filter, so T uses its own size (5)
        self.assert_tree_values(result, {"9": 160.0, "10": 40.0})

    def test_root_fallback_warning_is_quiet(self):
        tree = [
            ("7", "1", "G", "genus"),
            ("8", "7", "S", "species"),
        ]
        result, _, stderr = self.run_correction(
            tree, {"8": 100}, {"1": 1000}, ["all"], quiet=True
        )
        self.assert_tree_values(result, {"8": 100.0})
        self.assertEqual(stderr, "", "warning should be suppressed with quiet")

    def test_ranks_abundance_higher_rank_target_uses_root(self):
        # With --ranks-abundance species, a direct genus target (7)
        # has no allowed ancestor and falls back to the root size
        tree = [
            ("7", "1", "G", "genus"),
            ("8", "7", "S1", "species"),
            ("9", "8", "T", "strain"),
            ("11", "1", "G2", "genus"),
            ("12", "11", "S2", "species"),
        ]
        genome_sizes = {"1": 1000, "7": 50, "8": 10, "9": 5, "11": 25, "12": 20}
        result, _, stderr = self.run_correction(
            tree, {"7": 100, "12": 100}, genome_sizes, ["species"]
        )
        # total ratio: 100/1000 (G -> root) + 100/20 (S2) = 5.1
        self.assert_tree_values(
            result,
            {
                "7": 3.9216,  # G: 200 * 0.1 / 5.1
                "12": 196.0784,  # S2: 200 * 5 / 5.1
            },
        )
        self.assertIn("1 target(s) were corrected", stderr)

    def test_all_sizes_unknown_keeps_counts_unchanged(self):
        tree = [
            ("7", "1", "G", "genus"),
            ("8", "7", "S", "species"),
        ]
        result, _, _ = self.run_correction(tree, {"8": 100}, {"1": 1}, ["all"])
        # a root size of 1 does not change values
        self.assert_tree_values(result, {"1": 100.0, "7": 100.0, "8": 100.0})

    def test_empty_target_counts(self):
        tree, _, _ = self.proxy_example()
        result, _, stderr = self.run_correction(tree, {}, {"1": 1000}, ["all"])
        self.assertEqual(result, {})
        self.assertEqual(stderr, "")


if __name__ == "__main__":
    unittest.main()
