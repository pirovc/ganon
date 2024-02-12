import unittest
import sys
import os
import pickle 
from time import sleep

sys.path.append('src')
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import run_ganon
from utils import setup_dir
from utils import build_sanity_check_and_parse
data_dir = base_dir + "data/"

from parameterized import parameterized_class

@parameterized_class([
   { "filter_type": "ibf"},
   { "filter_type": "hibf" },
])
class TestUpdate(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.results_dir = base_dir + "results/integration/update_" + self.filter_type + "/"    
        # 3 small genomes from genbank (also on GTDB R207 for bac arc)
        # archaea GCA_002254805.1
        # bacteria GCA_000147015.1
        # viral GCA_004132065.1
        # Export local_dir for genome_updater, uses local folder as repository
        os.environ["local_dir"] = os.path.abspath(data_dir + "build/")
        setup_dir(self.results_dir)

    def test_output_db_prefix(self):
        """
        ganon build with --organism-group archaea bacteria and update "adding" viral with --output-db-prefix
        """

        build_params = {"db_prefix": self.results_dir + "test_output_db_prefix",
                        "organism_group": ["archaea", "bacteria"],
                        "source": ["genbank"],
                        "taxonomy": "skip",
                        "level": "assembly",
                        "threads": 1,
                        "filter_type": self.filter_type,
                        "write_info_file": True,
                        "keep_files": True,
                        "verbose": True,
                        "quiet": False}
        cfg = Config("build", **build_params)    
        # Run ganon build
        self.assertTrue(run_ganon(cfg, build_params["db_prefix"]), "ganon build run failed")
        # Load config from written file (to get all arguments generated on build for build custom)
        cfg = pickle.load(open(build_params["db_prefix"] + "_files/config.pkl", "rb"))
        res = build_sanity_check_and_parse(cfg)
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")
        
        # check if all 2 assemblies were used
        self.assertEqual(res["info"].shape[0], 2, "Wrong number of files")

        # Change genome_update history and add viral to simulate an "update"
        os.system("sed -i 's/archaea,bacteria/archaea,bacteria,viral/' " + build_params["db_prefix"] + "_files/history.tsv")
        # Sleep for a second to not overlap prefixes on genome_updater
        sleep(1)

        update_params = {"db_prefix": self.results_dir + "test_output_db_prefix",
                         "output_db_prefix": self.results_dir + "test_output_db_prefix2",
                         "threads": 1,
                         "write_info_file": True,
                         "keep_files": True,
                         "verbose": True,
                         "quiet": False}
        cfg = Config("update", **update_params)    
        # Run ganon build
        self.assertTrue(run_ganon(cfg, update_params["output_db_prefix"]), "ganon update run failed")
        # Load config from written file (to get all arguments generated on build for build custom)
        cfg = pickle.load(open(update_params["output_db_prefix"] + "_files/config.pkl", "rb"))

        res = build_sanity_check_and_parse(cfg)
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")        
        # check if all 3 assemblies were used
        self.assertEqual(res["info"].shape[0], 3, "Wrong number of files")

    def test_db_prefix(self):
        """
        ganon build with --organism-group archaea bacteria and update "adding" viral with same --db-prefix
        """

        build_params = {"db_prefix": self.results_dir + "test_db_prefix",
                        "organism_group": ["archaea", "bacteria"],
                        "source": ["genbank"],
                        "taxonomy": "skip",
                        "level": "assembly",
                        "threads": 1,
                        "filter_type": self.filter_type,
                        "write_info_file": True,
                        "keep_files": True,
                        "verbose": True,
                        "quiet": False}
        cfg = Config("build", **build_params)    
        # Run ganon build
        self.assertTrue(run_ganon(cfg, build_params["db_prefix"]), "ganon build run failed")
        # Load config from written file (to get all arguments generated on build for build custom)
        cfg = pickle.load(open(build_params["db_prefix"] + "_files/config.pkl", "rb"))
        res = build_sanity_check_and_parse(cfg)
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")
        
        # check if all 2 assemblies were used
        self.assertEqual(res["info"].shape[0], 2, "Wrong number of files")

        # Change genome_update history and add viral to simulate an "update"
        os.system("sed -i 's/archaea,bacteria/archaea,bacteria,viral/' " + build_params["db_prefix"] + "_files/history.tsv")
        # Sleep for a second to not overlap prefixes on genome_updater
        sleep(1)

        update_params = {"db_prefix": self.results_dir + "test_db_prefix",
                        "threads": 1,
                        "write_info_file": True,
                        "keep_files": True,
                        "verbose": True,
                        "quiet": False}
        cfg = Config("update", **update_params)    
        # Run ganon build
        self.assertTrue(run_ganon(cfg, update_params["db_prefix"]), "ganon update run failed")
        # Load config from written file (to get all arguments generated on build for build custom)
        cfg = pickle.load(open(update_params["db_prefix"] + "_files/config.pkl", "rb"))

        res = build_sanity_check_and_parse(cfg)
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")        
        # check if all 3 assemblies were used
        self.assertEqual(res["info"].shape[0], 3, "Wrong number of files")

if __name__ == '__main__':
    unittest.main()
