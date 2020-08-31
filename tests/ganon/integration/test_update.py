import unittest, sys
sys.path.append('src')
from ganon import ganon
from ganon.config import Config

sys.path.append('tests/ganon/integration/')
from utils import *

base_dir = "tests/ganon/integration/"
data_dir = base_dir + "data/"

class TestUpdateOffline(unittest.TestCase):

    results_dir = base_dir + "results/update/"
    default_params = {"taxdump_file": [data_dir+"mini_nodes.dmp", 
                                       data_dir+"mini_names.dmp"],
                      "db_prefix": data_dir+"bacteria_default",
                      "input_files": [data_dir+"update/virus_NC_003676.1.fasta.gz",
                                      data_dir+"update/virus_NC_011646.1.fasta.gz", 
                                      data_dir+"update/virus_NC_032412.1.fasta.gz", 
                                      data_dir+"update/virus_NC_035470.1.fasta.gz"],
                      "seq_info_file": data_dir+"update/virus_seqinfo.txt",
                      "write_seq_info_file": True,
                      "quiet": True}
    
    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)
       
    def test_default(self):
        """
        Test run with default parameters
        """
        params = self.default_params.copy()
        params["output_db_prefix"] = self.results_dir + "test_default"
        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")

        # Classify simulated virus against updated index
        params_classify = {}
        params_classify["db_prefix"] = params["output_db_prefix"] #params["output_db_prefix"]
        params_classify["single_reads"] = data_dir+"vir.sim.1.fq"
        params_classify["max_error"] = 0
        params_classify["output_all"] = True
        params_classify["quiet"] = True
        params_classify["output_prefix"] = self.results_dir + "test_default"
        # Build config from params
        cfg_classify = Config("classify", **params_classify)
        # Run
        self.assertTrue(ganon.main(cfg=cfg_classify), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg_classify))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # Specific tes - should return Viruses matches on the updated index
        self.assertEqual(res["tre_pd"][res["tre_pd"]["rank"]=="superkingdom"]["name"].values[0], "Viruses", "classification on updated index failed")

    def test_assembly(self):
        """
        Test rank as assembly
        """
        params = self.default_params.copy()
        params["db_prefix"] = data_dir+"bacteria_assembly"
        params["output_db_prefix"] = self.results_dir + "test_assembly"

        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # Specific test - count assemblies on tax (3 bac + 4 vir)
        self.assertEqual(sum(res["tax_pd"]["rank"]=="assembly"), 7, "error updating assemblies")

class TestUpdateOnline(unittest.TestCase):
    
    results_dir = base_dir + "results/update/online/"
    default_params = {"input_files": [data_dir+"update/virus_NC_003676.1.fasta.gz",
                                      data_dir+"update/virus_NC_011646.1.fasta.gz", 
                                      data_dir+"update/virus_NC_032412.1.fasta.gz", 
                                      data_dir+"update/virus_NC_035470.1.fasta.gz"],
                      "db_prefix": data_dir+"bacteria_default",
                      "write_seq_info_file": True,
                      "quiet": True}
    
    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)

    def test_default(self):
        """
        With default parameters online
        """
        params = self.default_params.copy()
        params["output_db_prefix"] = self.results_dir + "test_default"
        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        
        # Classify simulated virus against updated index
        params_classify = {}
        params_classify["db_prefix"] = params["output_db_prefix"] #params["output_db_prefix"]
        params_classify["single_reads"] = data_dir+"vir.sim.1.fq"
        params_classify["max_error"] = 0
        params_classify["output_all"] = True
        params_classify["quiet"] = True
        params_classify["output_prefix"] = self.results_dir + "test_default"
        # Build config from params
        cfg_classify = Config("classify", **params_classify)
        # Run
        self.assertTrue(ganon.main(cfg=cfg_classify), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg_classify))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # Specific tes - should return Viruses matches on the updated index
        self.assertEqual(res["tre_pd"][res["tre_pd"]["rank"]=="superkingdom"]["name"].values[0], "Viruses", "classification on updated index failed")


    def test_assembly(self):
        """
        Test rank as assembly online
        """
        params = self.default_params.copy()
        params["db_prefix"] = data_dir+"bacteria_assembly"
        params["output_db_prefix"] = self.results_dir + "test_assembly"

        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # Specific test - count assemblies on tax (3 bac + 4 vir)
        self.assertEqual(sum(res["tax_pd"]["rank"]=="assembly"), 7, "error updating assemblies")

if __name__ == '__main__':
    unittest.main()