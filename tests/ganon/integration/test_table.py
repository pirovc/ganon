import unittest
import sys
sys.path.append('src')
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import setup_dir
from utils import table_sanity_check_and_parse
from utils import run_ganon
data_dir = base_dir + "data/"


class TestTable(unittest.TestCase):

    results_dir = base_dir + "results/integration/table/"
    default_params = {"input": [data_dir+"table/report_reads1.tre",
                                data_dir+"table/report_reads2.tre",
                                data_dir+"table/report_reads3.tre"],
                      "rank": "species",
                      "transpose": True,
                      "verbose": True,
                      "quiet": False}

    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)

    def test_default(self):
        """
        Test ganon table with default parameters
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_default.tsv"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")

    def test_input_directory(self):
        """
        Test run with default parameters using input directory and extension
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_input_directory.tsv"
        params["input"] = data_dir+"table/"
        params["input_extension"] = ".tre"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")

    def test_rank(self):
        """
        Test ganon table with --rank
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_rank.tsv"
        params["rank"] = "superkingdom"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should output just bacteria
        self.assertEqual(res["out_pd"].columns.values.size, 1, "ganon table rank selection failed")

    def test_no_rank(self):
        """
        Test ganon table without specific --rank
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_no_rank.tsv"
        params["rank"] = ""
        params["header"] = "lineage"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        self.assertEqual(res["out_pd"].columns.values.size, 50, "ganon table without rank failed")
        # Test some outcomes
        self.assertTrue("1|2" in res["out_pd"].columns.values, "ganon table without rank failed")
        self.assertTrue("1|2|1239|909932|1843489|31977" in res["out_pd"].columns.values, "ganon table without rank failed")

        # Sum of counts should be total of all reads
        self.assertEqual(res["out_pd"].sum().sum(), 1973568, "ganon table without rank failed")

    def test_no_rank_unclassified(self):
        """
        Test ganon table without specific --rank reporting unclassified
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_no_rank_unc.tsv"
        params["rank"] = ""
        params["header"] = "lineage"
        params["unclassified_label"] = "unclassified"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        self.assertEqual(res["out_pd"].columns.values.size, 51, "ganon table without rank failed")
        # Test some outcomes
        self.assertTrue("1|2" in res["out_pd"].columns.values, "ganon table without rank failed")
        self.assertTrue("1|2|1239|909932|1843489|31977" in res["out_pd"].columns.values, "ganon table without rank failed")

        # Sum of counts should be total of all reads
        self.assertEqual(res["out_pd"].sum().sum(), 3786439, "ganon table without rank failed")

    def test_no_rank_no_root(self):
        """
        Test ganon table without specific --rank and no root
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_no_rank_no_root.tsv"
        params["rank"] = ""
        params["header"] = "lineage"
        params["no_root"] = True

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should output just bacteria
        self.assertEqual(res["out_pd"].columns.values.size, 49, "ganon table without rank failed")
        # Test some outcomes
        self.assertTrue("2" in res["out_pd"].columns.values, "ganon table without rank failed")
        self.assertTrue("2|1239|909932|1843489|31977" in res["out_pd"].columns.values, "ganon table without rank failed")
        # Sum of counts should be total of all reads with remaining root counts not reported
        self.assertEqual(res["out_pd"].sum().sum(), 1703834, "ganon table without rank failed")

    def test_no_rank_no_root_unclassified(self):
        """
        Test ganon table without specific --rank and no root reporting to unclassified
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_no_rank_no_root_unc.tsv"
        params["rank"] = ""
        params["header"] = "lineage"
        params["no_root"] = True
        params["unclassified_label"] = "unclassified"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should output just bacteria
        self.assertEqual(res["out_pd"].columns.values.size, 50, "ganon table without rank failed")
        # Test some outcomes
        self.assertTrue("2" in res["out_pd"].columns.values, "ganon table without rank failed")
        self.assertTrue("2|1239|909932|1843489|31977" in res["out_pd"].columns.values, "ganon table without rank failed")
        self.assertTrue("unclassified" in res["out_pd"].columns.values, "ganon table without rank failed")

        # Sum of counts should be total of all reads with remaining root matches counted as unclassified
        self.assertEqual(res["out_pd"].sum().sum(), 3786439, "ganon table without rank failed")


    def test_min_count(self):
        """
        Test ganon table with --min-count
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_min_count.tsv"
        params["output_value"] = "counts"
        params["min_count"] = 15000

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should output just counts higher than min_count (or zeros)
        self.assertTrue(((res["out_pd"]==0) | (res["out_pd"]>=params["min_count"])).all(axis=None) , "ganon table min count filter failed")

    def test_min_count_perc(self):
        """
        Test ganon table with --min-count below 1
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_min_count_perc.tsv"
        params["output_value"] = "percentage"
        params["min_count"] = 0.01

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should output just value higher than min_count (or zeros)
        self.assertTrue(((res["out_pd"]==0) | (res["out_pd"] >= params["min_count"])).all(axis=None) , "ganon table min count filter failed")

    def test_max_count(self):
        """
        Test ganon table with --max-count
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_max_count.tsv"
        params["output_value"] = "counts"
        params["max_count"] = 14000

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should output just counts higher than min_count (or zeros)
        self.assertTrue(((res["out_pd"]==0) | (res["out_pd"]<=params["max_count"])).all(axis=None) , "ganon table min count filter failed")

    def test_max_count_perc(self):
        """
        Test ganon table with --max-count below 1
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_max_count_perc.tsv"
        params["output_value"] = "percentage"
        params["max_count"] = 0.02

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should output just counts higher than min_count (or zeros)
        self.assertTrue(((res["out_pd"]==0) | (res["out_pd"]<=params["max_count"])).all(axis=None) , "ganon table min count filter failed")


    def test_taxids_relative(self):
        """
        Test ganon table with --taxids not on the chosen rank
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_taxids_relative.tsv"
        params["taxids"] = "838" # genus: Prevotella

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should output species of the genus 838 (Prevotella) only
        self.assertTrue(all("Prevotella" in r for r in res["out_pd"].columns.values), "ganon table taxids filter failed")

    def test_taxids_direct(self):
        """
        Test ganon table with --taxids of the chosen rank
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_taxids.tsv"
        params["taxids"] = "1110546" # species: Veillonella tobetsuensis

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should output species with the taxid 1110546 (Veillonella tobetsuensis)
        self.assertEqual(res["out_pd"].columns.values, "Veillonella tobetsuensis", "ganon table taxids filter failed")

    def test_names(self):
        """
        Test ganon table with --names
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_names.tsv"
        params["names"] = "Veillonella tobetsuensis"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should output species with the name Veillonella tobetsuensis
        self.assertEqual(res["out_pd"].columns.values, "Veillonella tobetsuensis", "ganon table names filter failed")

    def test_names_with(self):
        """
        Test ganon table with --names-with
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_names_with.tsv"
        params["names_with"] = "Prevotella"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should output species with the name starting with Veillonella
        self.assertTrue(all("Prevotella" in r for r in res["out_pd"].columns.values), "ganon table names with filter failed")

    def test_top_sample(self):
        """
        Test ganon table with --top-sample
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_top_sample.tsv"
        params["top_sample"] = 1
        params["rank"] = "genus"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should have 3 cols (each file has a different top genus)
        self.assertEqual(res["out_pd"].shape[1], 3, "ganon table top sample filter failed")

    def test_top_all(self):
        """
        Test ganon table with --top-all
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_top_all.tsv"
        params["top_all"] = 1
        params["rank"] = "genus"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should have 1 col (top genus)
        self.assertEqual(res["out_pd"].shape[1], 1, "ganon table top sample filter failed")

    def test_min_frequency(self):
        """
        Test ganon table with --min-frequency
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_min_frequency.tsv"
        params["min_frequency"] = 3
        params["rank"] = "phylum" # Fusobacteria is left out from report_reads3.tre

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should have no zero entries (3 samples/3 min frequency)
        self.assertTrue((res["out_pd"].values>0).all(), "ganon table min frequency filter failed")

    def test_min_frequency_perc(self):
        """
        Test ganon table with --min-frequency below 1
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_min_frequency_perc.tsv"
        params["min_frequency"] = 0.9
        params["rank"] = "phylum" # Fusobacteria is left out from report_reads3.tre

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should have no zero entries (3 samples/3 min frequency)
        self.assertTrue((res["out_pd"].values>0).all(), "ganon table min frequency filter failed")

    def test_extra_cols(self):
        """
        Test ganon table with --unclassified-label and --filtered-label
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_extra_cols1.tsv"
        params["min_count"] = 0.02
        params["unclassified_label"] = "UNC"
        params["filtered_label"] = "FIL"
        params["rank"] = "genus"
        cfg = Config("table", **params)
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        res1 = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res1, "ganon table has inconsistent results")
        # last 2 cols should be the fixed unclassified and filtered
        self.assertTrue(all(c in ["UNC", "FIL"] for c in res1["out_pd"].columns.values[-2:]), "ganon table extra cols failed")

        params["output_file"] = self.results_dir + "test_extra_cols2.tsv"
        params["unclassified_label"] = "UNC"
        params["filtered_label"] = None
        cfg = Config("table", **params)
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        res2 = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res2, "ganon table has inconsistent results")
        # last col should be the fixed unclassified
        self.assertEqual(res2["out_pd"].columns.values[-1], "UNC", "ganon table extra cols failed")
        # should not have the filtered
        self.assertFalse("FIL" in res2["out_pd"].columns.values, "ganon table extra cols failed")

        params["output_file"] = self.results_dir + "test_extra_cols3.tsv"
        params["unclassified_label"] = None
        params["filtered_label"] = "FIL"
        cfg = Config("table", **params)
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        res3 = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res3, "ganon table has inconsistent results")
        # last col should be the fixed filtered
        self.assertEqual(res3["out_pd"].columns.values[-1], "FIL", "ganon table extra cols failed")
        # should not have the unclassified
        self.assertFalse("UNC" in res3["out_pd"].columns.values, "ganon table extra cols failed")

        ## SAME LABEL, report together
        params["output_file"] = self.results_dir + "test_extra_cols4.tsv"
        params["unclassified_label"] = "UNASSIGNED"
        params["filtered_label"] = "UNASSIGNED"
        cfg = Config("table", **params)
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        res4 = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res4, "ganon table has inconsistent results")
        # last col should be the fixed label
        self.assertEqual(res4["out_pd"].columns.values[-1], "UNASSIGNED", "ganon table extra cols failed")

        # when reporting together, should match values
        self.assertTrue(all(res1["out_pd"][["UNC", "FIL"]].sum(axis=1)==res4["out_pd"]["UNASSIGNED"]), "ganon table extra cols failed")


    def test_skip_zeros(self):
        """
        Test ganon table with --skip-zeros
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_extra_cols.tsv"
        params["skip_zeros"] = True
        params["min_count"] = 0.02
        params["no_root"] = True

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should have 1 line with values higher than 0.02
        self.assertEqual(res["out_pd"].shape[0], 1, "ganon table skip zeros option failed")

    def test_matches(self):
        """
        Test ganon table with report type "matches" from ganon report
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_matches.tsv"
        params["input"] = [data_dir+"table/report_matches1.tre",
                           data_dir+"table/report_matches2.tre",
                           data_dir+"table/report_matches3.tre"]
        params["unclassified_label"] = "unclassified"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should have unclassified summing to 0 (not unclassified line reported)
        self.assertEqual(res["out_pd"]["unclassified"].sum(), 0, "ganon table min frequency filter failed")

    def test_header(self):
        """
        Test ganon table with different headers
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_header.tsv"
        params["header"] = "lineage"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # check if printed lineage on all headers (but one root)
        self.assertTrue(all(["|" in c for c in res["out_pd"].columns.values]), "ganon table headers are wrong (lineage)")

        params["header"] = "taxid"
        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_file"]), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # check if printed taxid (just numeric for this specific test)
        self.assertTrue(all([c.isdigit() for c in res["out_pd"].columns.values]), "ganon table headers are wrong (taxid)")

if __name__ == '__main__':
    unittest.main()
