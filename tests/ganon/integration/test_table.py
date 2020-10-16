import unittest, sys
sys.path.append('src')
from ganon import ganon
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import *
data_dir = base_dir + "data/"

class TestTableOffline(unittest.TestCase):

    results_dir = base_dir + "results/integration/table/"
    default_params = {"tre_file": data_dir+"table/results.tre",
                      "quiet": True}
    
    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)
       
    def test_default(self):
        """
        Test run with default parameters
        """
        params = self.default_params.copy()
        params["output_file"] = self.results_dir + "test_default.tsv"
        
        # Build config from params
        cfg = Config("table", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon table exited with an error")
        # General sanity check of results
        res = table_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon table has inconsistent results")


if __name__ == '__main__':
    unittest.main()