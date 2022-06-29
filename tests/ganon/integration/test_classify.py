import unittest
import sys
sys.path.append('src')
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import setup_dir
from utils import build_sanity_check_and_parse
from utils import classify_sanity_check_and_parse
from utils import run_ganon
data_dir = base_dir + "data/"


class TestClassify(unittest.TestCase):

    results_dir = base_dir + "results/integration/classify/"
    default_params = {"db_prefix": [results_dir + "base_build"],
                      "paired_reads": [data_dir+"classify/sim.1.fq.gz",
                                       data_dir+"classify/sim.2.fq.gz"],
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
                        "threads": 1,
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

if __name__ == '__main__':
    unittest.main()
