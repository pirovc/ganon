import unittest
import sys
import os
import pickle 

sys.path.append('src')
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import run_ganon
from utils import setup_dir
from utils import list_files_folder
from utils import list_sequences
from utils import build_sanity_check_and_parse
data_dir = base_dir + "data/"


class TestBuild(unittest.TestCase):

    results_dir = base_dir + "results/integration/build/"

    default_params = {"organism_group": ["archaea", "bacteria", "viral"],
                      "source": ["genbank"],
                      "taxonomy": "skip",
                      "threads": 1,
                      "write_info_file": True,
                      "keep_files": True,
                      "verbose": True,
                      "quiet": False}

    @classmethod
    def setUpClass(self):    
        # 3 small genomes from genbank (also on GTDB R207 for bac arc)
        # archaea GCA_002254805.1
        # bacteria GCA_000147015.1
        # viral GCA_004132065.1
        # Export local_dir for genome_updater, uses local folder as repository
        os.environ["local_dir"] = os.path.abspath(data_dir + "build/")
        setup_dir(self.results_dir)

    def test_og_all(self):
        """
        ganon build with --organism-group "archaea" "bacteria" "viral" 
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_og_all"

        # Run ganon build
        self.assertTrue(run_ganon(Config("build", **params), params["db_prefix"]), "ganon build run failed")
        # Load config from written file (to get all arguments generated on build)
        cfg = pickle.load(open(params["db_prefix"] + "_files/config.pkl", "rb"))
        res = build_sanity_check_and_parse(cfg)
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")

if __name__ == '__main__':
    unittest.main()
