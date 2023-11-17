import unittest
import sys
import os
sys.path.append('src')
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import setup_dir
from utils import build_sanity_check_and_parse
from utils import classify_sanity_check_and_parse
from utils import reassign_sanity_check_and_parse
from utils import parse_all_one, parse_rep
from utils import run_ganon
from utils import check_files
data_dir = base_dir + "data/"


class TestReassign(unittest.TestCase):

    results_dir = base_dir + "results/integration/reassign/"
    default_params_classify = {"db_prefix": [results_dir + "base_build", results_dir + "base_build2"],
                               "single_reads": [data_dir+"reassign/sim.fq.gz"],
                               "rel_cutoff": 0.01,
                               "rel_filter": 1,
                               "output_all": True,
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
                        "level": "assembly",
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

    def test_single(self):
        """
        Test ganon reassign with one database/output
        """
        params_classify = self.default_params_classify.copy()
        params_classify["output_prefix"] = self.results_dir + "single"

        # Build config from params
        cfg = Config("classify", **params_classify)
        self.assertTrue(run_ganon(
            cfg, params_classify["output_prefix"]), "ganon classify exited with an error")
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # There are multiple matches on output
        self.assertTrue(len(res["all_pd"].readid) > len(
            res["all_pd"].readid.unique()), "ganon classify has only unique matches")
        total_reads_classified = res["rep_pd"]["unique"].fillna(0).astype(
            int).sum() + res["rep_pd"]["lca"].fillna(0).astype(int).sum()

        # Reassign
        params = {"input_prefix": params_classify["output_prefix"],
                  "output_prefix": params_classify["output_prefix"] + "_reassigned"}
        cfg = Config("reassign", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon reassign exited with an error")
        res = reassign_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon reassign has inconsistent results")
        # There are only single matches on output
        self.assertEqual(len(res["one_pd"].readid), len(
            res["all_pd"].readid.unique()), "ganon reassign has multiple matches")
        # Check if all reads got properly reported
        self.assertEqual(total_reads_classified, res["rep_pd"]["unique"].fillna(0).astype(int).sum(
        ) + res["rep_pd"]["lca"].fillna(0).astype(int).sum(), "ganon reassign reported wrong number of reads")

    def test_hierarchy_output_single(self):
        """
        Test ganon classify with two database
        """
        params_classify = self.default_params_classify.copy()
        params_classify["output_prefix"] = self.results_dir + \
            "hierarchy_output_single"
        params_classify["hierarchy_labels"] = ["A", "B"]
        params_classify["output_single"] = True
        params_classify["rel_cutoff"] = ["0.1", "0.01"]

        # Build config from params
        cfg = Config("classify", **params_classify)
        self.assertTrue(run_ganon(
            cfg, params_classify["output_prefix"]), "ganon classify exited with an error")
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # There are multiple matches on output
        self.assertTrue(len(res["all_pd"].readid) > len(
            res["all_pd"].readid.unique()), "ganon classify has only unique matches")
        total_reads_classified = res["rep_pd"]["unique"].fillna(0).astype(
            int).sum() + res["rep_pd"]["lca"].fillna(0).astype(int).sum()

        # Reassign
        params = {"input_prefix": params_classify["output_prefix"],
                  "output_prefix": params_classify["output_prefix"] + "_reassigned"}
        cfg = Config("reassign", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon reassign exited with an error")
        res = reassign_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon reassign has inconsistent results")
        # There are only single matches on output
        self.assertEqual(len(res["one_pd"].readid), len(
            res["all_pd"].readid.unique()), "ganon reassign has multiple matches")
        # Check if all reads got properly reported
        self.assertEqual(total_reads_classified, res["rep_pd"]["unique"].fillna(0).astype(int).sum(
        ) + res["rep_pd"]["lca"].fillna(0).astype(int).sum(), "ganon reassign reported wrong number of reads")

    def test_hierarchy_output_split(self):
        """
        Test ganon classify with two database
        """
        params_classify = self.default_params_classify.copy()
        params_classify["output_prefix"] = self.results_dir + \
            "hierarchy_output_split"
        params_classify["hierarchy_labels"] = ["A", "B"]
        params_classify["rel_cutoff"] = ["0.1", "0.01"]

        # Build config from params
        cfg = Config("classify", **params_classify)
        self.assertTrue(run_ganon(
            cfg, params_classify["output_prefix"]), "ganon classify exited with an error")
        rep = parse_rep(params_classify["output_prefix"] + ".rep")
        all_A = parse_all_one(params_classify["output_prefix"] + ".A.all")
        all_B = parse_all_one(params_classify["output_prefix"] + ".B.all")
        self.assertIsNotNone(all_A, "ganon classify has inconsistent results")
        self.assertIsNotNone(all_B, "ganon classify has inconsistent results")
        total_reads_classified = rep["unique"].fillna(0).astype(
            int).sum() + rep["lca"].fillna(0).astype(int).sum()

        # There are multiple matches on output
        self.assertTrue(len(all_A.readid) > len(
            all_A.readid.unique()), "ganon classify has only unique matches")
        self.assertTrue(len(all_B.readid) > len(
            all_B.readid.unique()), "ganon classify has only unique matches")

        # Reassign
        params = {"input_prefix": params_classify["output_prefix"],
                  "output_prefix": params_classify["output_prefix"] + "_reassigned"}
        cfg = Config("reassign", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon reassign exited with an error")
        rep = parse_rep(params["output_prefix"] + ".rep")
        one_A = parse_all_one(params["output_prefix"] + ".A.one")
        one_B = parse_all_one(params["output_prefix"] + ".B.one")

        # There are only single matches on output
        self.assertEqual(len(one_A.readid), len(
            all_A.readid.unique()), "ganon reassign has multiple matches")
        # Check if all reads got properly reported
        self.assertEqual(total_reads_classified, rep["unique"].fillna(0).astype(int).sum(
        ) + rep["lca"].fillna(0).astype(int).sum(), "ganon reassign reported wrong number of reads")

    def test_no_rep(self):
        """
        Test ganon reassign with one database/output and no .rep file (just .all)
        """
        params_classify = self.default_params_classify.copy()
        params_classify["output_prefix"] = self.results_dir + "no_rep"

        # Build config from params
        cfg = Config("classify", **params_classify)
        self.assertTrue(run_ganon(
            cfg, params_classify["output_prefix"]), "ganon classify exited with an error")
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # There are multiple matches on output
        self.assertTrue(len(res["all_pd"].readid) > len(
            res["all_pd"].readid.unique()), "ganon classify has only unique matches")
        total_reads_classified = res["rep_pd"]["unique"].fillna(0).astype(
            int).sum() + res["rep_pd"]["lca"].fillna(0).astype(int).sum()

        # remove .rep before reassign
        os.remove(params_classify["output_prefix"]+".rep")

        # Reassign
        params = {"input_prefix": params_classify["output_prefix"],
                  "output_prefix": params_classify["output_prefix"] + "_reassigned"}
        cfg = Config("reassign", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon reassign exited with an error")
        res = reassign_sanity_check_and_parse(vars(cfg))
        
        self.assertIsNotNone(res["one_pd"], "ganon reassign has inconsistent results")

        # There are only single matches on output
        self.assertEqual(len(res["one_pd"].readid), len(
            res["all_pd"].readid.unique()), "ganon reassign has multiple matches")

        # Check if all reads got properly reported
        self.assertEqual(total_reads_classified, len(res["one_pd"].readid), "ganon reassign reported wrong number of reads")

    def test_remove_all(self):
        """
        Test ganon reassign with --remove-all
        """
        params_classify = self.default_params_classify.copy()
        params_classify["output_prefix"] = self.results_dir + "remove_all"

        # Build config from params
        cfg = Config("classify", **params_classify)
        self.assertTrue(run_ganon(
            cfg, params_classify["output_prefix"]), "ganon classify exited with an error")
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # There are multiple matches on output
        self.assertTrue(len(res["all_pd"].readid) > len(
            res["all_pd"].readid.unique()), "ganon classify has only unique matches")
        total_reads_classified = res["rep_pd"]["unique"].fillna(0).astype(
            int).sum() + res["rep_pd"]["lca"].fillna(0).astype(int).sum()

        # Reassign
        params = {"input_prefix": params_classify["output_prefix"],
                  "output_prefix": params_classify["output_prefix"] + "_reassigned",
                  "remove_all": True}

        cfg = Config("reassign", **params)
        self.assertTrue(run_ganon(cfg, params["output_prefix"]), "ganon reassign exited with an error")

        res = reassign_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon reassign has inconsistent results")
        
        # .all file removed 
        self.assertFalse(check_files(params["input_prefix"], "all"))

        # Check if all reads got properly reported
        self.assertEqual(total_reads_classified, res["rep_pd"]["unique"].fillna(0).astype(int).sum() + res["rep_pd"]["lca"].fillna(0).astype(int).sum(), "ganon reassign reported wrong number of reads")

    def test_skip_one(self):
        """
        Test ganon reassign with --skip-one
        """
        params_classify = self.default_params_classify.copy()
        params_classify["output_prefix"] = self.results_dir + "skip_one"

        # Build config from params
        cfg = Config("classify", **params_classify)
        self.assertTrue(run_ganon(
            cfg, params_classify["output_prefix"]), "ganon classify exited with an error")
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # There are multiple matches on output
        self.assertTrue(len(res["all_pd"].readid) > len(
            res["all_pd"].readid.unique()), "ganon classify has only unique matches")
        total_reads_classified = res["rep_pd"]["unique"].fillna(0).astype(
            int).sum() + res["rep_pd"]["lca"].fillna(0).astype(int).sum()

        # Reassign
        params = {"input_prefix": params_classify["output_prefix"],
                  "output_prefix": params_classify["output_prefix"] + "_reassigned",
                  "skip_one": True}

        cfg = Config("reassign", **params)
        self.assertTrue(run_ganon(cfg, params["output_prefix"]), "ganon reassign exited with an error")

        res = reassign_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon reassign has inconsistent results")
        
        # .one file not created
        self.assertFalse(check_files(params["output_prefix"], "one"))

        # Check if all reads got properly reported
        self.assertEqual(total_reads_classified, res["rep_pd"]["unique"].fillna(0).astype(int).sum() + res["rep_pd"]["lca"].fillna(0).astype(int).sum(), "ganon reassign reported wrong number of reads")


if __name__ == '__main__':
    unittest.main()
