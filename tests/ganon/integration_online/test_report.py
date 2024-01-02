import os
from utils import run_ganon
from utils import download_bulk_files
from utils import classify_sanity_check_and_parse
from utils import build_sanity_check_and_parse
from utils import report_sanity_check_and_parse
from utils import setup_dir
from ganon.config import Config
import unittest
import sys
sys.path.append('src')

base_dir = "tests/ganon/"
sys.path.append(base_dir)
data_dir = base_dir + "data/"


class TestReport(unittest.TestCase):

    results_dir = base_dir + "results/integration_online/report/"
    download_dir = base_dir + "downloads/"

    default_params = {"ncbi_url": "file://" + os.path.abspath(download_dir),
                      "gtdb_url": "file://" + os.path.abspath(download_dir),
                      "verbose": True,
                      "quiet": False}

    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)

        # Download full files once (and simulate download by setting ncbi_url and gtdb_url on params)
        download_bulk_files(self.download_dir)

    def test_ncbi(self):
        """
        Test run with --taxonomy ncbi, downloading .tax
        """

        # Build base database
        build_params = {"db_prefix": self.results_dir + "base_build_ncbi",
                        "input": data_dir + "build-custom/files/",
                        "taxonomy": "ncbi",
                        "taxonomy_files": data_dir + "build-custom/taxdump.tar.gz",
                        "ncbi_file_info":  data_dir + "build-custom/assembly_summary.txt",
                        "genome_size_file": data_dir + "build-custom/species_genome_size.txt.gz",
                        "level": "species",
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

        classify_params = {"db_prefix": [self.results_dir + "base_build_ncbi"],
                           "rel_cutoff": "0",  # no cutoff and filter to allow spurious matches everywhere
                           "rel_filter": "1",
                           "output_prefix": self.results_dir + "base_classify_ncbi",
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

        # build at species level and remove tax (will be re-created on report)
        os.remove(self.results_dir + "base_build_ncbi.tax")

        params = self.default_params.copy()
        params["input"] = self.results_dir + "base_classify_ncbi.rep"
        params["output_prefix"] = self.results_dir + "test_ncbi"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        # Build config from params
        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")


        # No genome size
        params = self.default_params.copy()
        params["input"] = self.results_dir + "base_classify_ncbi.rep"
        params["output_prefix"] = self.results_dir + "test_ncbi_skip_genome_size"
        params["taxonomy"] = "ncbi"
        params["taxonomy_files"] = data_dir + "build-custom/taxdump.tar.gz"
        params["skip_genome_size"] = True
        # Build config from params
        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")

    def test_gtdb(self):
        """
        Test run with --taxonomy gtdb, downloading .tax
        """

        # Build base database
        build_params = {"db_prefix": self.results_dir + "base_build_gtdb",
                        "input": data_dir + "build-custom/files/",
                        "taxonomy": "gtdb",
                        "taxonomy_files": [data_dir + "build-custom/ar53_taxonomy.tsv.gz",
                                           data_dir + "build-custom/bac120_taxonomy.tsv.gz"],
                        "genome_size_files": [data_dir + "build-custom/ar53_metadata.tsv.gz",
                                              data_dir + "build-custom/bac120_metadata.tsv.gz"],
                        "level": "species",
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

        classify_params = {"db_prefix": [self.results_dir + "base_build_gtdb"],
                           "rel_cutoff": "0",  # no cutoff and filter to allow spurious matches everywhere
                           "rel_filter": "1",
                           "output_prefix": self.results_dir + "base_classify_gtdb",
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

        # build at species level and remove tax (will be re-created on report)
        os.remove(self.results_dir + "base_build_gtdb.tax")

        params = self.default_params.copy()
        params["input"] = self.results_dir + "base_classify_gtdb.rep"
        params["output_prefix"] = self.results_dir + "test_gtdb"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [data_dir + "build-custom/ar53_taxonomy.tsv.gz",
                                    data_dir + "build-custom/bac120_taxonomy.tsv.gz"]
        # Build config from params
        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")


        # No genome size
        params = self.default_params.copy()
        params["input"] = self.results_dir + "base_classify_gtdb.rep"
        params["output_prefix"] = self.results_dir + "test_gtdb_skip_genome_size"
        params["taxonomy"] = "gtdb"
        params["taxonomy_files"] = [data_dir + "build-custom/ar53_taxonomy.tsv.gz",
                                    data_dir + "build-custom/bac120_taxonomy.tsv.gz"]
        params["skip_genome_size"] = True
        # Build config from params
        cfg = Config("report", **params)
        self.assertTrue(
            run_ganon(cfg, params["output_prefix"]), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")

if __name__ == '__main__':
    unittest.main()
