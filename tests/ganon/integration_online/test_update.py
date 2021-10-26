import unittest, sys, shutil
sys.path.append('src')
from ganon import ganon
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import *
data_dir = base_dir + "data/"

class TestUpdateOnline(unittest.TestCase):
    
    results_dir = base_dir + "results/integration_online/update/"
    default_params = {"input_files": [data_dir+"update/virus_NC_003676.1.fasta.gz",
                                      data_dir+"update/virus_NC_011646.1.fasta.gz", 
                                      data_dir+"update/virus_NC_032412.1.fasta.gz", 
                                      data_dir+"update/virus_NC_035470.1.fasta.gz"],
                      "taxdump_file": [data_dir+"mini_nodes.dmp", 
                                       data_dir+"mini_names.dmp"],
                      "db_prefix": data_dir+"bacteria_default",
                      "write_seq_info_file": True,
                      "quiet": True}
    
    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)

    def test_default(self):
        """
        ganon update with default parameters (online: eutils, taxdump)
        """
        params = self.default_params.copy()
        params["output_db_prefix"] = self.results_dir + "test_default"
        params["taxdump_file"] = []
        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        
        # Classify simulated virus against updated index
        params_classify = {"db_prefix": params["output_db_prefix"],
            "single_reads": [data_dir+"vir.sim.1.fq", data_dir+"bac.sim.1.fq"],
            "max_error": 0,
            "output_all": True,
            "quiet": True,
            "output_lca": True,
            "output_prefix": self.results_dir + "test_default"}

        # Build config from params
        cfg_classify = Config("classify", **params_classify)
        # Run
        self.assertTrue(ganon.main(cfg=cfg_classify), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg_classify))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # Specific tes - should return Viruses and Bacteria matches on the updated index
        self.assertTrue(res["tre_pd"][res["tre_pd"]["rank"]=="superkingdom"]["name"].isin(["Bacteria","Viruses"]).all(), "classification on updated index failed")

    def test_specialization_assembly(self):
        """
        ganon update --specialization assembly (online: eutils)
        """
        params = self.default_params.copy()
        params["db_prefix"] = data_dir+"bacteria_assembly"
        params["output_db_prefix"] = self.results_dir + "test_specialization_assembly"
        params["specialization"] = "assembly"
        
        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # Specific test - count assemblies on tax (3 bac + 4 vir)
        self.assertEqual(sum(res["tax_pd"]["rank"]=="assembly"), 7, "error updating with assemblies")

    def test_specialization_file_on_custom(self):
        """
        ganon update --specialization file with previous generated index --specialization custom (online: eutils)
        """
        params = self.default_params.copy()
        params["db_prefix"] = data_dir+"bacteria_custom"
        params["output_db_prefix"] = self.results_dir + "test_specialization_file_on_custom"
        params["specialization"] = "file"
        
        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # Specific test - count assemblies on tax (3 bac + 4 vir)
        self.assertEqual(sum(res["tax_pd"]["rank"]=="custom"), 3, "error updating")
        self.assertEqual(sum(res["tax_pd"]["rank"]=="file"), 4, "error updating")

    def test_specialization_sequence_on_custom(self):
        """
        ganon update --specialization sequence with previous generated index --specialization custom (online: eutils)
        """
        params = self.default_params.copy()
        params["db_prefix"] = data_dir+"bacteria_custom"
        params["output_db_prefix"] = self.results_dir + "test_specialization_sequence_on_custom"
        params["specialization"] = "sequence"
        
        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # Specific test - count assemblies on tax (3 bac + 4 vir)
        self.assertEqual(sum(res["tax_pd"]["rank"]=="custom"), 3, "error updating")
        self.assertEqual(sum(res["tax_pd"]["rank"]=="sequence"), 4, "error updating")



    def test_specialization_assembly_on_custom(self):
        """
        ganon update --specialization assembly with previous generated index --specialization custom (online: eutils)
        """
        params = self.default_params.copy()
        params["db_prefix"] = data_dir+"bacteria_custom"
        params["output_db_prefix"] = self.results_dir + "test_specialization_assembly_on_custom"
        params["specialization"] = "assembly"
        
        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # Specific test - count assemblies on tax (3 bac + 4 vir)
        self.assertEqual(sum(res["tax_pd"]["rank"]=="custom"), 3, "error updating")
        self.assertEqual(sum(res["tax_pd"]["rank"]=="assembly"), 4, "error updating")


if __name__ == '__main__':
    unittest.main()
