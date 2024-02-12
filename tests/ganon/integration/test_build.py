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
from utils import build_sanity_check_and_parse
data_dir = base_dir + "data/"

from parameterized import parameterized_class

@parameterized_class([
   { "filter_type": "ibf"},
   { "filter_type": "hibf" },
])
class TestBuild(unittest.TestCase):

    default_params = {"organism_group": ["archaea", "bacteria", "viral"],
                      "source": ["genbank"],
                      "taxonomy": "skip",
                      "level": "assembly",
                      "threads": 1,
                      "filter_type": "ibf",
                      "write_info_file": True,
                      "keep_files": True,
                      "verbose": True,
                      "quiet": False}

    @classmethod
    def setUpClass(self):    
        # 3 small genomes from genbank (also on GTDB R207 for bac arc)
        # archaea GCA_002254805.1 txid 2012515
        # bacteria GCA_000147015.1 txid 871271
        # viral GCA_004132065.1 txid 2161879
        # Export local_dir for genome_updater, uses local folder as repository
        os.environ["local_dir"] = os.path.abspath(data_dir + "build/")

        # parametrization for filter type
        self.default_params["filter_type"] = self.filter_type
        self.results_dir = base_dir + "results/integration/build_" + self.filter_type + "/"
        setup_dir(self.results_dir)

    def test_og_arc_bac_vir(self):
        """
        ganon build with --organism-group "archaea" "bacteria" "viral" 
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_og_arc_bac_vir"

        cfg = Config("build", **params)
    
        # Run ganon build
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build run failed")
        # Load config from written file (to get all arguments generated on build for build custom)
        cfg = pickle.load(open(params["db_prefix"] + "_files/config.pkl", "rb"))

        res = build_sanity_check_and_parse(cfg)
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")
        
        # check if all 3 assemblies were used
        self.assertEqual(res["info"].shape[0], 3, "Wrong number of files")

    def test_taxid(self):
        """
        ganon build with --taxid 131567 (cellular organisms)
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_taxid"
        params["organism_group"] = []
        params["taxid"] = "131567"

        cfg = Config("build", **params)
    
        # Run ganon build
        self.assertTrue(run_ganon(cfg, params["db_prefix"]), "ganon build run failed")
        # Load config from written file (to get all arguments generated on build for build custom)
        cfg = pickle.load(open(params["db_prefix"] + "_files/config.pkl", "rb"))

        res = build_sanity_check_and_parse(cfg)
        self.assertIsNotNone(res, "ganon build-custom sanity check failed")
        
        # Only 2 assemblies should be part of cellular organisms (bac, arc)
        self.assertEqual(res["info"].shape[0], 2, "Wrong number of files")


if __name__ == '__main__':
    unittest.main()
