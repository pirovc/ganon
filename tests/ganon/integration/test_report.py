from utils import setup_dir
from utils import run_ganon
from utils import report_sanity_check_and_parse
from utils import build_sanity_check_and_parse
from utils import classify_sanity_check_and_parse
from utils import list_files_folder
from ganon.config import Config
from math import ceil
import unittest
import sys
import shutil
sys.path.append('src')

base_dir = "tests/ganon/"
sys.path.append(base_dir)
data_dir = base_dir + "data/"


class TestReport(unittest.TestCase):

    results_dir = base_dir + "results/integration/report/"
    default_params = {"input": results_dir + "base_classify.rep",
                      "db_prefix": [results_dir + "base_build2", results_dir + "base_build"],
                      "verbose": True,
                      "quiet": False}

    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)
        # Build base database
        build_params = {"db_prefix": self.results_dir + "base_build",
                        "input": data_dir + "build-custom/files/",
                        "taxonomy": "ncbi",
                        "taxonomy_files": data_dir + "build-custom/taxdump.tar.gz",
                        "ncbi_file_info":  data_dir + "build-custom/assembly_summary.txt",
                        "genome_size_file": data_dir + "build-custom/species_genome_size.txt.gz",
                        "level": "species",
                        "filter_type": "ibf",
                        "threads": 1,
                        "keep_files": True,
                        "write_info_file": True,
                        "verbose": True,
                        "quiet": False}
        build_cfg = Config("build-custom", **build_params)
        self.assertTrue(run_ganon(
            build_cfg, build_params["db_prefix"]), "ganon build-custom run failed")
        self.assertIsNotNone(build_sanity_check_and_parse(
            vars(build_cfg)), "ganon build-custom sanity check failed")

        # Build second base database
        build_params = {"db_prefix": self.results_dir + "base_build2",
                        "input": data_dir + "build-custom/files/more/",
                        "taxonomy": "ncbi",
                        "taxonomy_files": data_dir + "build-custom/taxdump.tar.gz",
                        "ncbi_file_info":  data_dir + "build-custom/assembly_summary.txt",
                        "genome_size_file": data_dir + "build-custom/species_genome_size.txt.gz",
                        "level": "assembly",
                        "threads": 1,
                        "filter_type": "ibf",
                        "keep_files": True,
                        "write_info_file": True,
                        "verbose": True,
                        "quiet": False}
        build_cfg = Config("build-custom", **build_params)
        self.assertTrue(run_ganon(
            build_cfg, build_params["db_prefix"]), "ganon build-custom run failed")
        self.assertIsNotNone(build_sanity_check_and_parse(
            vars(build_cfg)), "ganon build-custom sanity check failed")

        classify_params = {"db_prefix": [self.results_dir + "base_build2",
                                         self.results_dir + "base_build"],
                           "hierarchy_labels": ["A", "B"],
                           "rel_cutoff": "0",  # no cutoff and filter to allow spurious matches everywhere
                           "rel_filter": "1",
                           "output_prefix": self.results_dir + "base_classify",
                           "paired_reads": [data_dir+"classify/sim.1.fq.gz",
                                            data_dir+"classify/sim.2.fq.gz"],
                           "verbose": True,
                           "multiple_matches": "lca",
                           "quiet": False}
        # Build config from params
        classify_cfg = Config("classify", **classify_params)
        self.assertTrue(run_ganon(
            classify_cfg, classify_params["output_prefix"]), "ganon classify exited with an error")
        self.assertIsNotNone(classify_sanity_check_and_parse(
            vars(classify_cfg)), "ganon classify sanity check failed")

        classify_params = {"db_prefix": [self.results_dir + "base_build2",
                                         self.results_dir + "base_build"],
                           "hierarchy_labels": ["C", "D"],
                           "output_prefix": self.results_dir + "base_classify2",
                           "paired_reads": [data_dir+"classify/sim.1.fq.gz",
                                            data_dir+"classify/sim.2.fq.gz"],
                           "verbose": True,
                           "quiet": False}
        # Build config from params
        classify_cfg = Config("classify", **classify_params)
        self.assertTrue(run_ganon(
            classify_cfg, classify_params["output_prefix"]), "ganon classify exited with an error")
        self.assertIsNotNone(classify_sanity_check_and_parse(
            vars(classify_cfg)), "ganon classify sanity check failed")

    def test_default(self):
        """
        Test run with default parameters
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_default"

        # Build config from params
        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")

    def test_min_count(self):
        """
        Test run with min_count
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_min_count"
        params["min_count"] = 15

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if none is higher than min_count
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["cumulative"] >= params["min_count"]).all(),
                        "ganon report failed filtering with --min-count")

    def test_min_count_perc(self):
        """
        Test run with min_count using percentages
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_min_count_perc"
        params["min_count"] = 0.2

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if none is higher than min_count
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["cumulative_perc"] >= (params["min_count"]*100)).all(),
                        "ganon report failed filtering with --min-count")

    def test_max_count(self):
        """
        Test run with --max-count
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_max_count"
        params["max_count"] = 30

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if none is higher than min_count
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["cumulative"] <= params["max_count"]).all(),
                        "ganon report failed filtering with --max-count")

    def test_max_count_perc(self):
        """
        Test run with --max-count below 1
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_max_count_perc"
        params["max_count"] = 0.2

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if none is higher than min_count
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["cumulative_perc"] <= (params["max_count"]*100)).all(),
                        "ganon report failed filtering with --max-count")

    def test_report_type_abundance(self):
        """
        Test run with report_type abundance
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + \
            "test_report_type_abundance"
        params["report_type"] = "abundance"

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")

        # Re-distribution and genome correction, shared sum bigger then 0
        self.assertTrue(res["tre_pd"][res["tre_pd"]["rank"] == "assembly"]["shared"].sum(
        ) > 0, "ganon report has wrong output for --report_type abundance")

    def test_report_type_reads(self):
        """
        Test run with report_type reads
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_report_type_reads"
        params["report_type"] = "reads"

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")

        # No re-distribution, shared sum to 0
        self.assertEqual(res["tre_pd"][res["tre_pd"]["rank"] == "assembly"]["shared"].sum(
        ), 0, "ganon report has wrong output for --report_type reads")

    def test_report_type_matches(self):
        """
        Test run with report_type matches
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_report_type_matches"
        params["report_type"] = "matches"

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # should not output unclassified
        self.assertFalse((res["tre_pd"]['rank'] == "unclassified").any(),
                         "ganon report has wrong output for --report_type matches")

    def test_report_type_corr(self):
        """
        Test run with report_type abundance
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_report_type_corr"
        params["report_type"] = "corr"

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")

        # No re-distribution, shared sum to 0
        self.assertEqual(res["tre_pd"][res["tre_pd"]["rank"] == "assembly"]["shared"].sum(
        ), 0, "ganon report has wrong output for --report_type corr")

    def test_report_type_dist(self):
        """
        Test run with report_type abundance
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_report_type_dist"
        params["report_type"] = "dist"

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")

        # Re-distribution, shared sum bigger then 0
        self.assertTrue(res["tre_pd"][res["tre_pd"]["rank"] == "assembly"]["shared"].sum(
        ) > 0, "ganon report has wrong output for --report_type dist")

    def test_ranks(self):
        """
        Test run with limited ranks
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_ranks"
        params["ranks"] = ["phylum", "species"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if only selected ranks were reported
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["rank"].isin(params["ranks"])).all(),
                        "ganon report did not report the correct ranks")

    def test_ranks_all(self):
        """
        Test run with all ranks
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_ranks_all_abundance"
        params["report_type"] = "abundance"
        params["ranks"] = "all"
        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        res = report_sanity_check_and_parse(vars(cfg))
    
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if reported any "no rank" rank
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["rank"] == "no rank").any(),
                        "ganon report did not report the correct ranks")

        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_ranks_all_reads"
        params["report_type"] = "reads"
        params["ranks"] = "all"
        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if reported any "no rank" rank
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["rank"] == "no rank").any(),
                        "ganon report did not report the correct ranks")

        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_ranks_all_matches"
        params["report_type"] = "matches"
        params["ranks"] = "all"
        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if reported any "no rank" rank
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["rank"] == "no rank").any(),
                        "ganon report did not report the correct ranks")

        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_ranks_all_corr"
        params["report_type"] = "corr"
        params["ranks"] = "all"
        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if reported any "no rank" rank
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["rank"] == "no rank").any(),
                        "ganon report did not report the correct ranks")

        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_ranks_all_dist"
        params["report_type"] = "dist"
        params["ranks"] = "all"
        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if reported any "no rank" rank
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["rank"] == "no rank").any(),
                        "ganon report did not report the correct ranks")

    def test_skip_hierachy(self):
        """
        Test run skipping hierachies
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_skip_hierachy"
        params["skip_hierarchy"] = ["A"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(
            vars(cfg), sum_full_percentage=False)
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # should not have any assembly reported (base_build)
        self.assertFalse((res["tre_pd"][~res["idx_base"]]["rank"].isin(["assembly"])).any(),
                         "ganon report did not skip the hierarchy")

    def test_keep_hierachy(self):
        """
        Test run keeping hierachies
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_keep_hierachy"
        params["keep_hierarchy"] = ["B"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(
            vars(cfg), sum_full_percentage=False)
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # should not have any assembly reported
        self.assertFalse((res["tre_pd"][~res["idx_base"]]["rank"].isin(["assembly"])).any(),
                         "ganon report did not skip the hierarchy")

    def test_split_hierachy(self):
        """
        Test run splitting hierachies
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_split_hierachy"
        params["split_hierarchy"] = True

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(
            vars(cfg), sum_full_percentage=False)
        self.assertIsNotNone(res, "ganon report has inconsistent results")

        # sum all root value
        total_root_split = 0
        for file, r in res.items():
            total_root_split += r["tre_pd"][r["tre_pd"]
                                            ['rank'] == "root"]["cumulative_perc"].values[0]
        # sum one time unclassified
        total_root_split += r["tre_pd"][r["tre_pd"]['rank']
                                        == "unclassified"]["cumulative_perc"].values[0]
        # values reported on root of splitted reports should equal 100
        self.assertEqual(int(total_root_split), 100,
                         "ganon report with wrong root counts")

    def test_multiple_rep_files(self):
        """
        Test run with multiple rep files as input
        """
        params = self.default_params.copy()
        params["input"] = [self.results_dir + "base_classify.rep",
                           self.results_dir + "base_classify2.rep"]
        params["output_prefix"] = self.results_dir + "test_multiple_rep_files"

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # should have two outputs
        self.assertEqual(len(res), len(
            params["input"]), "ganon report did not generate multiple report files")

    def test_multiple_rep_files_folder(self):
        """
        Test run with multiple rep files as input
        """
        setup_dir(self.results_dir + "multi_files")
        shutil.copy(self.results_dir + "base_classify.rep",
                    self.results_dir + "multi_files/")
        shutil.copy(self.results_dir + "base_classify2.rep",
                    self.results_dir + "multi_files/")

        params = self.default_params.copy()
        params["input"] = [self.results_dir + "multi_files"]
        params["input_extension"] = "rep"
        params["output_prefix"] = self.results_dir + \
            "test_multiple_rep_files_folder"

        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # should have two outputs
        self.assertEqual(len(res), len(list_files_folder(params["input"][0], ext=params["input_extension"])),
                         "ganon report did not generate multiple report files")

    def test_multiple_rep_files_split_hierachy(self):
        """
        Test run with multiple rep files as input
        """
        setup_dir(self.results_dir + "multi_files")
        shutil.copy(self.results_dir + "base_classify.rep",
                    self.results_dir + "multi_files/")
        shutil.copy(self.results_dir + "base_classify2.rep",
                    self.results_dir + "multi_files/")

        params = self.default_params.copy()
        params["input"] = [self.results_dir + "base_classify.rep",
                           self.results_dir + "base_classify2.rep"]
        params["output_prefix"] = self.results_dir + \
            "test_multiple_rep_files_split_hierachy_"
        params["split_hierarchy"] = True

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(
            vars(cfg), sum_full_percentage=False)
        self.assertIsNotNone(res, "ganon report has inconsistent results")

        # should have 2+2 outputs (4 hiearchies)
        self.assertEqual(
            len(res), 4, "ganon report did not generate multiple report files")

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
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # should have only reported the species asked
        self.assertEqual(res["tre_pd"][~res["idx_base"]]["name"].values[0], params["names"][0],
                         "ganon report did not filter by name")

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
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # should have only matches with pattern
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["name"].str.contains(params["names_with"][0])).all(),
                        "ganon report did not filter by names with")

    def test_taxids(self):
        """
        Test run filtering for taxids
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_taxids"
        params["taxids"] = ["1224"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # should have only matches with pattern
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["lineage"].str.contains(params["taxids"][0])).all(),
                        "ganon report did not filter by taxids")

    def test_top_percentile(self):
        """
        Test --top-percentile
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + \
            "test_top_percentile_default"
        params["top_percentile"] = 0

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")

        # count output taxa by rank
        total_by_rank = res["tre_pd"][~res["idx_base"]].groupby([
                                                                'rank']).size()

        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_top_percentile_50"
        params["top_percentile"] = 0.5
        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if only half of taxa is reported
        self.assertTrue((total_by_rank/2).apply(ceil).equals(
            res["tre_pd"][~res["idx_base"]].groupby(['rank']).size()),
            "percentile reporting wrong number of taxa")

        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_top_percentile_25"
        params["top_percentile"] = 0.25
        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if only a quarter of taxa is reported
        self.assertTrue((total_by_rank/4).apply(ceil).equals(
            res["tre_pd"][~res["idx_base"]].groupby(['rank']).size()),
            "percentile reporting wrong number of taxa")

    def test_taxdump_file(self):
        """
        Test run using taxdump instead of db_prefix
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_taxdump_file"
        params["db_prefix"] = ""
        params["taxonomy_files"] = [data_dir + "build-custom/taxdump.tar.gz"]
        params["genome_size_files"] = [data_dir +
                                       "build-custom/species_genome_size.txt.gz"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
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
        params["taxonomy_files"] = [data_dir + "build-custom/taxdump.tar.gz"]
        params["genome_size_files"] = [data_dir +
                                       "build-custom/species_genome_size.txt.gz"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if reported any "na" rank
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["rank"] == "na").any(),
                        "ganon report did not report the correct ranks")

    def test_na_ranks(self):
        """
        Test run reporting missing taxa
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_na_ranks"
        params["db_prefix"] = ""
        params["ranks"] = ["genus", "species", "na"]
        params["taxonomy_files"] = [data_dir + "build-custom/taxdump.tar.gz"]
        params["genome_size_files"] = [data_dir +
                                       "build-custom/species_genome_size.txt.gz"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if reported any "na" rank
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["rank"] == "na").any(),
                        "ganon report did not report the correct ranks")

    def test_no_orphan(self):
        """
        Test run reporting missing taxa and not reporting them
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_no_orphan"
        params["db_prefix"] = ""
        params["ranks"] = "all"
        params["no_orphan"] = True
        params["taxonomy_files"] = [data_dir + "build-custom/taxdump.tar.gz"]
        params["genome_size_files"] = [data_dir +
                                       "build-custom/species_genome_size.txt.gz"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # Should not report any "na"
        self.assertFalse((res["tre_pd"][~res["idx_base"]]["rank"] == "na").any(),
                         "ganon report did not report the correct ranks")

    def test_only_orphan(self):
        """
        Test run reporting missing taxa and reporting only orphan
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_only_orphan"
        params["db_prefix"] = ""
        params["ranks"] = "na"
        params["no_orphan"] = False
        params["taxonomy_files"] = [data_dir + "build-custom/taxdump.tar.gz"]
        params["genome_size_files"] = [data_dir +
                                       "build-custom/species_genome_size.txt.gz"]

        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # Should not report any "na"
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["rank"] == "na").all(),
                        "ganon report did not report the correct ranks")

    def test_output_format_tsv(self):
        """
        Test run with --output-format tsv
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_output_format_tsv"
        params["output_format"] = "tsv"

        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")

    def test_output_format_csv(self):
        """
        Test run with --output-format csv
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_output_format_csv"
        params["output_format"] = "csv"

        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")

    def test_output_format_text(self):
        """
        Test run with --output-format text
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_output_format_text"
        params["output_format"] = "text"

        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")

    def test_output_format_bioboxes(self):
        """
        Test run with --output-format bioboxes
        """

        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + \
            "test_output_format_bioboxes_base"
        params["output_format"] = "tsv"

        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")

        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + \
            "test_output_format_bioboxes"
        params["output_format"] = "bioboxes"
        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        taxid = []

        with open(params["output_prefix"] + ".tre", "r") as file:
            for line in file:
                if line[0] == "@":
                    continue
                taxid.append(line.rstrip().split("\t")[0])

        # Check if output has same taxids
        self.assertTrue(res["tre_pd"][~res["idx_base"]]["target"].isin(taxid).all(),
                        "ganon report inconsistend on bioboxes output format")


if __name__ == '__main__':
    unittest.main()
