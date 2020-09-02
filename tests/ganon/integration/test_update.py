import unittest, sys, shutil
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
        # Specific - check if number of bins increased
        self.assertTrue(res["map_pd"].binid.max()>41, "no bins were added")

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

    def test_duplicated(self):
        """
        Test duplicated entries on update
        """
        params = self.default_params.copy()

        params["output_db_prefix"] = self.results_dir + "test_duplicated"
        params["input_files"].append([data_dir+"build/bacteria_NC_010333.1.fasta.gz",
                                      data_dir+"build/bacteria_NC_017164.1.fasta.gz", 
                                      data_dir+"build/bacteria_NC_017163.1.fasta.gz", 
                                      data_dir+"build/bacteria_NC_017543.1.fasta.gz"])
        params["seq_info_file"] = data_dir+"update/bacteria_virus_seqinfo.txt"

        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # Specific test - check if there are no duplicates in the datastructure
        self.assertFalse(res["tax_pd"]["taxid"].duplicated().any(), "duplicated entries on .tax after update")
        # Check if new map has any target from before in new bins
        map_before = parse_map(params["db_prefix"]+".map")
        new_targets = res["map_pd"][res["map_pd"]["binid"]>map_before.binid.max()]
        self.assertFalse(map_before["target"].isin(new_targets["target"]).any(), "duplicated entries on .map after update")


    def test_duplicated_only(self):
        """
        Test only duplicated entries on update
        """
        params = self.default_params.copy()

        params["output_db_prefix"] = self.results_dir + "test_duplicated_only"
        params["input_files"] = [data_dir+"build/bacteria_NC_010333.1.fasta.gz",
                                      data_dir+"build/bacteria_NC_017164.1.fasta.gz", 
                                      data_dir+"build/bacteria_NC_017163.1.fasta.gz", 
                                      data_dir+"build/bacteria_NC_017543.1.fasta.gz"]
        params["seq_info_file"] = data_dir+"update/bacteria_seqinfo.txt"
        # Build config from params
        cfg = Config("update", **params)
        # Should not run - nothing to update
        self.assertFalse(ganon.main(cfg=cfg), "ganon update exited with an error")

    def test_reuse_bins(self):
        """
        Test update without creating new bins
        """
        # Build database with 2 sequences at superkingdom level with large bins
        params_build = {"db_prefix": self.results_dir + "test_reuse_bins_part1",
                        "taxdump_file": [data_dir+"mini_nodes.dmp", 
                                         data_dir+"mini_names.dmp"],
                        "input_files": [data_dir+"update/virus_NC_003676.1.fasta.gz",
                                        data_dir+"update/virus_NC_011646.1.fasta.gz"],
                        "seq_info_file": data_dir+"update/virus_part1_seqinfo.txt",
                        "write_seq_info_file": True,
                        "rank": "superkingdom",
                        "bin_length": 200000,
                        "quiet": True}
        # Build config from params
        cfg_build = Config("build", **params_build)
        # Run
        self.assertTrue(ganon.main(cfg=cfg_build), "ganon update exited with an error")
        # General sanity check of results
        res_build = build_sanity_check_and_parse(vars(cfg_build))
        # Copy seqinfo to be parsed later
        shutil.copy(params_build["seq_info_file"], params_build["db_prefix"]+".seqinfo.txt")

        # Update with part2 - add virus to same bin and bacteria to new bins
        params = self.default_params.copy()
        params["db_prefix"] = params_build["db_prefix"]
        params["output_db_prefix"] = self.results_dir + "test_reuse_bins_part2"
        params["input_files"] = [data_dir+"update/virus_NC_032412.1.fasta.gz", 
                                      data_dir+"update/virus_NC_035470.1.fasta.gz",
                                      data_dir+"build/bacteria_NC_010333.1.fasta.gz",
                                      data_dir+"build/bacteria_NC_017164.1.fasta.gz", 
                                      data_dir+"build/bacteria_NC_017163.1.fasta.gz", 
                                      data_dir+"build/bacteria_NC_017543.1.fasta.gz"]
        params["seq_info_file"] = data_dir+"update/bacteria_virus_part2_seqinfo.txt"
        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")

        # Check if new virus (10239) were added to the same bin 0
        self.assertEqual(res["bins_pd"][res["bins_pd"]["taxid"]=="10239"].binid.drop_duplicates().max(), '0', "virus added to new bins")

        # Check if entries are unique (virus 10239, bacteria 2)
        self.assertEqual(res["bins_pd"].taxid.drop_duplicates().size, 2, "multiple taxids")


        # Classify to first part
        params_classify = {"db_prefix": params_build["db_prefix"],
                            "single_reads": [data_dir+"vir.sim.1.fq", data_dir+"bac.sim.1.fq"],
                            "max_error": 0,
                            "output_all": True,
                            "quiet": True,
                            "output_prefix": self.results_dir + "test_reuse_bins_classify_part1"}
        # Build config from params
        cfg_classify = Config("classify", **params_classify)
        # Run
        self.assertTrue(ganon.main(cfg=cfg_classify), "ganon classify exited with an error")
        # General sanity check of results
        res_classify1 = classify_sanity_check_and_parse(vars(cfg_classify))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # Specific test - should contain only one superkingdom (virus 10239, bacteria 2)
        self.assertEqual(res_classify1["tre_pd"][res_classify1["tre_pd"]["rank"]=="superkingdom"]["target"].shape[0], 1, "more than one superkingdom as target")
        self.assertEqual(res_classify1["tre_pd"][res_classify1["tre_pd"]["rank"]=="superkingdom"]["target"].values[0], '10239', "wrong target")

        # Classify to second part
        params_classify["db_prefix"] = params["output_db_prefix"]
        # Build config from params
        cfg_classify = Config("classify", **params_classify)
        # Run
        self.assertTrue(ganon.main(cfg=cfg_classify), "ganon classify exited with an error")
        # General sanity check of results
        res_classify2 = classify_sanity_check_and_parse(vars(cfg_classify))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # Specific
        # Classification to the second updated index has to have more matches than the first
        self.assertTrue(res_classify2["all_pd"].shape[0]>res_classify1["all_pd"].shape[0], "updated index did not improve matches")
        # should contain only two superkingdoms (virus 10239, bacteria 2)
        self.assertEqual(res_classify2["tre_pd"][res_classify2["tre_pd"]["rank"]=="superkingdom"]["target"].shape[0], 2, "more than two superkingdom as target")
        self.assertTrue(res_classify2["tre_pd"][res_classify2["tre_pd"]["rank"]=="superkingdom"]["target"].isin(['10239','2']).all(), "wrong target")


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