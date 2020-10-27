import unittest, sys
sys.path.append('src')
from ganon import ganon
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import *
data_dir = base_dir + "data/"

class TestTableOffline(unittest.TestCase):

    results_dir = base_dir + "results/integration/table/"
    default_params = {"tre_files": [data_dir+"table/report_reads1.tre", 
                                    data_dir+"table/report_reads2.tre",
                                    data_dir+"table/report_reads3.tre"],
                      "quiet": True}
    
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
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")

    def test_rank(self):
        """
        Test ganon table with --ranks
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_rank.tsv"
        params["rank"] = "superkingdom"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should output just bacteria
        self.assertEqual(res["out_pd"].columns.values.size, 1, "ganon table rank selection failed")
        
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
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should output just counts higher than min_count (or zeros)
        self.assertTrue(((res["out_pd"]==0) | (res["out_pd"]>=params["min_count"])).all(axis=None) , "ganon table min count filter failed")
    
    def test_min_percentage(self):
        """
        Test ganon table with --min-percentage
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_min_percentage.tsv"
        params["min_percentage"] = 0.01

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should output just value higher than min_percentage (or zeros)
        self.assertTrue(((res["out_pd"]==0) | (res["out_pd"]>=params["min_percentage"])).all(axis=None) , "ganon table min count filter failed")

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
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
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
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
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
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
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
        params["names"] = "Prevotella"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
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
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
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
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should have 1 col (top genus)
        self.assertEqual(res["out_pd"].shape[1], 1, "ganon table top sample filter failed")

    def test_min_occurence(self):
        """
        Test ganon table with --min-occurence
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_min_occurence.tsv"
        params["min_occurence"] = 3
        params["rank"] = "phylum" # Fusobacteria is left out from report_reads3.tre

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should have no zero entries (3 samples/3 min occurence)
        self.assertTrue((res["out_pd"].values>0).all(), "ganon table min occurence filter failed")

    def test_extra_cols(self):
        """
        Test ganon table with --add-unclassified --add-unclassified-rank --add-filtered
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_extra_cols.tsv"
        params["add_unclassified"] = True
        params["add_unclassified_rank"] = True
        params["add_filtered"] = True
        params["rank"] = "genus"

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # last 3 cols should be the fixed ones
        self.assertTrue(all(c in ["unclassified", "unclassified_"+params["rank"], "filtered"] for c in  res["out_pd"].columns.values[-3:]), "ganon table extra cols failed")

    def test_skip_zeros(self):
        """
        Test ganon table with --skip-zeros
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_extra_cols.tsv"
        params["skip_zeros"] = True
        params["min_percentage"] = 0.02

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
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
        params["tre_files"] = [data_dir+"table/report_matches1.tre", 
                               data_dir+"table/report_matches2.tre",
                               data_dir+"table/report_matches3.tre"]
        params["add_unclassified"] = True

        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        # should have unclassified summing to 0 (not unclassified line reported)
        self.assertEqual(res["out_pd"]["unclassified"].sum(), 0, "ganon table min occurence filter failed")

if __name__ == '__main__':
    unittest.main()