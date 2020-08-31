import unittest, sys
sys.path.append('src')
from ganon import ganon
from ganon.config import Config
sys.path.append('tests/ganon/integration/')
from utils import *

base_dir = "tests/ganon/integration/"
data_dir = base_dir + "data/"


class TestReportOffline(unittest.TestCase):

    results_dir = base_dir + "results/report/"
    default_params = {"db_prefix": data_dir+"bacteria_assembly",
                      "rep_file": data_dir+"report/results.rep",
                      "quiet": True}
    
    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)
       
    def test_default(self):
        """
        Test run with default parameters
        """
        params = self.default_params.copy()
        params["output_report"] = self.results_dir + "test_default.tre"
        
        # Build config from params
        cfg = Config("report", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon report exited with an error")
        # General sanity check of results
        res = report_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon report has inconsistent results")


if __name__ == '__main__':
    unittest.main()