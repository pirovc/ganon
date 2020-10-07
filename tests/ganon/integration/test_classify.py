import unittest, sys
sys.path.append('src')
from ganon import ganon
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import *
data_dir = base_dir + "data/"

class TestClassifyOffline(unittest.TestCase):

    results_dir = base_dir + "results/integration/classify/"
    default_params = {"db_prefix": data_dir+"bacteria_default",
                      "single_reads": data_dir+"bac.sim.1.fq",
                      "output_all": True,
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
        cfg = Config("classify", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")

    def test_multiple_databases(self):
        """
        Test run with default parameters
        """
        params = self.default_params.copy()
        params["output_prefix"] = self.results_dir + "test_multiple_databases"
        params["db_prefix"] = [data_dir+"bacteria_assembly", data_dir+"bacteria_default"]
        params["min_kmers"] = ["0.25", "0.1"]
        params["hierarchy_labels"] = ["1_assembly", "2_default"]
        params["output_single"] = True
        
        # Build config from params
        cfg = Config("classify", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")


def classify_classify_sanity_check_and_parse(params):
    # Provide sanity checks for outputs (not specific to a test) and return loaded data

    if not check_files(params["output_prefix"], ["lca","all","rep","tre"]):
        return None

    res = {}
    # Sequence information from database to be updated
    res["tre_pd"] =  parse_tre(params["output_prefix"]+".tre")
    return res

if __name__ == '__main__':
    unittest.main()