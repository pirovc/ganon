import unittest, sys
sys.path.append('src')
from ganon import ganon
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import *
data_dir = base_dir + "data/"

class TestBuildOnline(unittest.TestCase):
    
    results_dir = base_dir + "results/integration_online/build/"
    default_params = {"input_files": [data_dir+"build/bacteria_NC_010333.1.fasta.gz",
                                      data_dir+"build/bacteria_NC_017164.1.fasta.gz", 
                                      data_dir+"build/bacteria_NC_017163.1.fasta.gz", 
                                      data_dir+"build/bacteria_NC_017543.1.fasta.gz"],
                      "write_seq_info_file": True,
                      "rank": "species",
                      "quiet": True}

    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)

    def test_default(self):
        """
        With default parameters online
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_default"
        
        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
       
    def test_assembly(self):
        """
        Test rank as assembly online
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_assembly"
        params["rank"] = "assembly"

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results") 



if __name__ == '__main__':
    unittest.main()
