import unittest, sys
sys.path.append('src')
from ganon import ganon
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import *
data_dir = base_dir + "data/"

class TestReportOffline(unittest.TestCase):

    results_dir = base_dir + "results/integration/report/"
    default_params = {"db_prefix": data_dir+"bacteria_assembly",
                      "rep_files": data_dir+"report/results.rep",
                      "output_format": "tsv",
                      "quiet": True}

    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)
       
    def test_default(self):
        """
        Test run with default parameters
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_default"
        
        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")

    def test_input_directory(self):
        """
        Test run with default parameters using input directory and extension
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_input_directory"
        del params["rep_files"]
        params["input_directory"] = data_dir+"report/"
        params["input_extension"] = ".rep"

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
    
    def test_min_percentage(self):
        """
        Test run with min_percentage
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_min_percentage"
        params["min_percentage"] = 0.2
        
        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if none is higher than min_percentage
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["cumulative_perc"] >= params["min_percentage"]).all(), "ganon report failed filtering with --min-percentage")

    def test_min_count(self):
        """
        Test run with min_count
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_min_count"
        params["min_count"] = 60
        
        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if none is higher than min_percentage
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["cumulative"] >= params["min_count"]).all(), "ganon report failed filtering with --min-count")


    def test_min_count_and_percentages(self):
        """
        Test run with min_percentage and min_count
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_min_count_and_percentages"
        params["min_percentage"] = 0.2
        params["min_count"] = 50
        
        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if none is higher than min_percentage
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["cumulative"] >= params["min_count"]).all(), "ganon report failed filtering with --min-count")
        # check if none is higher than min_percentage
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["cumulative_perc"] >= params["min_percentage"]).all(), "ganon report failed filtering with --min-percentage")

    def test_report_type(self):
        """
        Test run with report_type matches
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_report_type"
        params["report_type"] = "matches"

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # should not output unclassified
        self.assertFalse((res["tre_pd"]['rank'] == "unclassified").any(), "ganon report has wrong output for --report_type matches")

    def test_ranks(self):
        """
        Test run with limited ranks
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_ranks"
        params["ranks"] = ["phylum","species"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if only selected ranks were reported
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["rank"].isin(params["ranks"])).all(),"ganon report did not report the correct ranks")

    def test_ranks_all(self):
        """
        Test run with all ranks
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_ranks_all"
        params["ranks"] = "all"

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if reported any "no rank" rank
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["rank"]=="no rank").any(),"ganon report did not report the correct ranks")

    def test_skip_hierachy(self):
        """
        Test run skipping hierachies
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_skip_hierachy"
        params["skip_hierarchy"] = ["1_assembly"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # should not have any assembly reported
        self.assertFalse((res["tre_pd"][~res["idx_base"]]["rank"].isin(["assembly"])).any(),"ganon report did not skip the hierarchy")

    def test_names(self):
        """
        Test run filtering for specific names
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_names"
        params["names"] = ["Bacteria"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # should have only reported the species asked
        self.assertEqual(res["tre_pd"][~res["idx_base"]]["name"].values[0],params["names"][0], "ganon report did not filter by name")

    def test_names_with(self):
        """
        Test run filtering for names with
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_names_with"
        params["names_with"] = ["bacter"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # should have only matches with pattern
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["name"].str.contains(params["names_with"][0])).all(), "ganon report did not filter by names with")

    def test_taxids(self):
        """
        Test run filtering for names with
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_taxids"
        params["taxids"] = ["1224"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # should have only matches with pattern
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["lineage"].str.contains(params["taxids"][0])).all(), "ganon report did not filter by taxids")

    def test_taxdump_file(self):
        """
        Test run using taxdump instead of db_prefix
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_taxdump_file"
        params["db_prefix"] = ""
        params["taxdump_file"] = [data_dir+"mini_nodes.dmp", data_dir+"mini_names.dmp"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")

    def test_na(self):
        """
        Test run reporting missing taxa
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_na"
        params["db_prefix"] = ""
        params["ranks"] = "all"
        params["taxdump_file"] = [data_dir+"mini_nodes.dmp", data_dir+"mini_names.dmp"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if reported any "na" rank
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["rank"]=="na").any(),"ganon report did not report the correct ranks")


    def test_na_ranks(self):
        """
        Test run reporting missing taxa
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_na_ranks"
        params["db_prefix"] = ""
        params["ranks"] = ["genus","species","na"]
        params["taxdump_file"] = [data_dir+"mini_nodes.dmp", data_dir+"mini_names.dmp"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if only selected ranks were reported
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["rank"].isin(params["ranks"])).all(),"ganon report did not report the correct ranks")

if __name__ == '__main__':
    unittest.main()