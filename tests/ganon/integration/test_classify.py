import unittest
import sys
sys.path.append('src')
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import setup_dir
from utils import build_sanity_check_and_parse
from utils import classify_sanity_check_and_parse, report_sanity_check_and_parse
from utils import run_ganon
from utils import check_files
data_dir = base_dir + "data/"

from parameterized import parameterized_class

@parameterized_class([
   { "filter_type": "ibf"},
   { "filter_type": "hibf" },
])
class TestClassify(unittest.TestCase):

    default_params = {"paired_reads": [data_dir+"classify/sim.1.fq.gz",
                                       data_dir+"classify/sim.2.fq.gz"],
                      "verbose": True,
                      "quiet": False}

    @classmethod
    def setUpClass(self):
        self.results_dir = base_dir + "results/integration/classify_" + self.filter_type + "/"
        self.default_params["db_prefix"] = [self.results_dir + "base_build"]
        setup_dir(self.results_dir)
        # Build base database
        build_params = {"db_prefix": self.results_dir + "base_build",
                        "input": data_dir + "build-custom/files/",
                        "taxonomy": "ncbi",
                        "taxonomy_files": data_dir + "build-custom/taxdump.tar.gz",
                        "ncbi_file_info":  data_dir + "build-custom/assembly_summary.txt",
                        "genome_size_file": data_dir + "build-custom/species_genome_size.txt.gz",
                        "level": "assembly",
                        "threads": 1,
                        "filter_type": self.filter_type,
                        "keep_files": True,
                        "write_info_file": True,
                        "verbose": True,
                        "quiet": False}
        build_cfg = Config("build-custom", **build_params)
        self.assertTrue(run_ganon(build_cfg, build_params["db_prefix"]), "ganon build-custom run failed")
        self.assertIsNotNone(build_sanity_check_and_parse(vars(build_cfg)), "ganon build-custom sanity check failed")

        # Build second base database
        build_params = {"db_prefix": self.results_dir + "base_build2",
                        "input": data_dir + "build-custom/files/more/",
                        "taxonomy": "ncbi",
                        "taxonomy_files": data_dir + "build-custom/taxdump.tar.gz",
                        "ncbi_file_info":  data_dir + "build-custom/assembly_summary.txt",
                        "genome_size_file": data_dir + "build-custom/species_genome_size.txt.gz",
                        "level": "assembly",
                        "threads": 1,
                        "filter_type": self.filter_type,
                        "keep_files": True,
                        "write_info_file": True,
                        "verbose": True,
                        "quiet": False}
        build_cfg = Config("build-custom", **build_params)
        self.assertTrue(run_ganon(build_cfg, build_params["db_prefix"]), "ganon build-custom run failed")
        self.assertIsNotNone(build_sanity_check_and_parse(vars(build_cfg)), "ganon build-custom sanity check failed")

    def test_single_db(self):
        """
        Test ganon classify with one database
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "single_db"

        # Build config from params
        cfg = Config("classify", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_prefix"]), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")

    def test_multi_db(self):
        """
        Test ganon classify with two database
        """
        params = self.default_params.copy()
        params["db_prefix"] = [self.results_dir + "base_build", self.results_dir + "base_build2"]
        params["output_prefix"] = self.results_dir + "multi_db"

        # Build config from params
        cfg = Config("classify", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_prefix"]), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")

    def test_hierarchy_db(self):
        """
        Test ganon classify with two database
        """
        params = self.default_params.copy()
        params["db_prefix"] = [self.results_dir + "base_build2", self.results_dir + "base_build"]
        params["output_prefix"] = self.results_dir + "hierarchy_db"
        params["hierarchy_labels"] = ["A", "B"]

        # Build config from params
        cfg = Config("classify", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_prefix"]), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")

    def test_output_files(self):
        """
        Test ganon classify with one database
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "output_files"
        params["output_one"] = True
        params["output_all"] = True
        params["output_unclassified"] = True

        # Build config from params
        cfg = Config("classify", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_prefix"]), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")

    def test_multiple_matches_em(self):
        """
        Test ganon classify with --multiple-matches em
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "multiple_matches_em"
        params["multiple_matches"] = "em"
        params["output_all"] = True
        params["output_one"] = True
        params["rel_cutoff"] = 0.001
        params["rel_filter"] = 1

        # Build config from params
        cfg = Config("classify", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_prefix"]), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        
        # There are only single matches on output
        self.assertEqual(len(res["one_pd"].readid), len(res["all_pd"].readid.unique()))

    def test_multiple_matches_lca(self):
        """
        Test ganon classify with --multiple-matches lca
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "multiple_matches_lca"
        params["multiple_matches"] = "lca"
        params["output_all"] = True
        params["output_one"] = True
        params["rel_cutoff"] = 0.001
        params["rel_filter"] = 1

        # Build config from params
        cfg = Config("classify", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_prefix"]), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        
        # There are only single matches on output
        self.assertEqual(len(res["one_pd"].readid), len(res["all_pd"].readid.unique()))

    def test_multiple_matches_skip(self):
        """
        Test ganon classify with --multiple-matches skip
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "multiple_matches_skip"
        params["multiple_matches"] = "skip"
        params["output_all"] = True
        params["rel_cutoff"] = 0.001
        params["rel_filter"] = 1

        # Build config from params
        cfg = Config("classify", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_prefix"]), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        
        # There is no .one output
        self.assertFalse(check_files(params["output_prefix"], "one"))

    def test_report_params(self):
        """
        Test ganon classify with report parameters
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "report_params"
        params["min_count"] = 0.1
        params["ranks"] = ["genus", "species"]
        params["report_type"] = "matches"

        # Build config from params
        cfg = Config("classify", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_prefix"]), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")

        # General sanity check of reports
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")
        # check if only selected ranks were reported
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["rank"].isin(params["ranks"])).all(),
                        "ganon report did not report the correct ranks")
        # check if none is higher than min_count
        self.assertTrue((res["tre_pd"][~res["idx_base"]]["cumulative"] >= params["min_count"]).all(),
                        "ganon report failed filtering with --min-count")
        # should not output unclassified
        self.assertFalse((res["tre_pd"]['rank'] == "unclassified").any(),
                         "ganon report has wrong output for --report_type matches")

    def test_skip_report(self):
        """
        Test ganon classify with --skip-report
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "skip_report"
        params["skip_report"] = True

        # Build config from params
        cfg = Config("classify", **params)
        # Run
        self.assertTrue(run_ganon(cfg, params["output_prefix"]), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")
        
        # There is no .tre output
        self.assertFalse(check_files(params["output_prefix"], "tre"))

if __name__ == '__main__':
    unittest.main()
