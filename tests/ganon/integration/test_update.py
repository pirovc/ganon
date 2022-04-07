import unittest, sys, shutil
sys.path.append('src')
from ganon import ganon
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import *
data_dir = base_dir + "data/"


class TestUpdateOffline(unittest.TestCase):

    results_dir = base_dir + "results/integration/update/"
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
        self.assertTrue(res["map_pd"].binid.max() > 41, "no bins were added")

        # Classify simulated virus against updated index
        params_classify = {"db_prefix": params["output_db_prefix"],
                           "single_reads": [data_dir + "vir.sim.1.fq", data_dir + "bac.sim.1.fq"],
                           "abs_cutoff": 0,
                           "output_lca": True,
                           "output_all": True,
                           "quiet": True,
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

    def test_minimizers(self):
        """
        ganon update with minimizers
        """
        params_build = {"taxdump_file": [data_dir + "mini_nodes.dmp",
                                         data_dir + "mini_names.dmp"],
                        "input_files": [data_dir + "build/bacteria_NC_010333.1.fasta.gz",
                                        data_dir + "build/bacteria_NC_017164.1.fasta.gz",
                                        data_dir + "build/bacteria_NC_017163.1.fasta.gz",
                                        data_dir + "build/bacteria_NC_017543.1.fasta.gz"],
                        "seq_info_file": data_dir + "build/bacteria_seqinfo.txt",
                        "write_seq_info_file": True,
                        "rank": "species",
                        "window_size": 27,
                        "quiet": True}
        params_build["db_prefix"] = self.results_dir + "test_minimizers_build"

        # Build config from params
        cfg_build = Config("build", **params_build)
        # Run
        self.assertTrue(ganon.main(cfg=cfg_build), "ganon build exited with an error")
        # General sanity check of results
        res_build = build_sanity_check_and_parse(vars(cfg_build))
        self.assertIsNotNone(res_build, "ganon build has inconsistent results")
        shutil.copy(params_build["seq_info_file"], params_build["db_prefix"]+".seqinfo.txt")

        params = self.default_params.copy()
        params["db_prefix"] = params_build["db_prefix"]
        params["output_db_prefix"] = self.results_dir + "test_minimizers_update"
        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # Specific - check if number of bins increased
        self.assertTrue(res["map_pd"].binid.max() > 41, "no bins were added")

        # Classify simulated virus against updated index
        params_classify = {"db_prefix": params["output_db_prefix"],
                           "single_reads": [data_dir + "vir.sim.1.fq", data_dir + "bac.sim.1.fq"],
                           "rel_cutoff": 0,
                           "rel_filter": 1,
                           "output_lca": True,
                           "output_all": True,
                           "quiet": True,
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


    def test_specialization_custom(self):
        """
        Test update with custom specialization
        """
        params = self.default_params.copy()
        params["db_prefix"] = data_dir+"bacteria_custom"
        params["output_db_prefix"] = self.results_dir + "test_specialization_custom"

        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # Specific test - count assemblies on tax (3 bac + 4 vir)
        self.assertEqual(sum(res["tax_pd"]["rank"]=="custom"), 7, "error updating assemblies")

    def test_specialization_custom_on_assembly(self):
        """
        ganon update --specialization custom with previous generated index with --specialization assembly
        """
        params = self.default_params.copy()
        params["db_prefix"] = data_dir+"bacteria_assembly"
        params["output_db_prefix"] = self.results_dir + "test_specialization_custom_on_custom"
        params["specialization"] = "custom"

        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # Specific test - count assemblies on tax (3 bac + 4 vir)
        self.assertEqual(sum(res["tax_pd"]["rank"]=="custom"), 4, "error updating assemblies")
        self.assertEqual(sum(res["tax_pd"]["rank"]=="assembly"), 3, "error updating assemblies")

    def test_specialization_on_default(self):
        """
        ganon update --specialization custom on previous generated index without specialiazazion
        """
        params = self.default_params.copy()
        params["db_prefix"] = data_dir+"bacteria_default"
        params["output_db_prefix"] = self.results_dir + "test_specialization_on_default"
        params["specialization"] = "custom"

        # Build config from params
        cfg = Config("update", **params)
        # Should not run
        self.assertFalse(ganon.main(cfg=cfg), "ganon update exited with an error")


    def test_repeated(self):
        """
        ganon update with some repeated sequences already in the index
        """
        params = self.default_params.copy()

        params["output_db_prefix"] = self.results_dir + "test_repeated"
        params["input_files"].extend([data_dir+"build/bacteria_NC_010333.1.fasta.gz",
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

    def test_repeated_only(self):
        """
        ganon update with only repeated sequences already in the index
        """
        params = self.default_params.copy()

        params["output_db_prefix"] = self.results_dir + "test_repeated_only"
        params["input_files"] = [data_dir+"build/bacteria_NC_010333.1.fasta.gz",
                                      data_dir+"build/bacteria_NC_017164.1.fasta.gz", 
                                      data_dir+"build/bacteria_NC_017163.1.fasta.gz", 
                                      data_dir+"build/bacteria_NC_017543.1.fasta.gz"]
        params["seq_info_file"] = data_dir+"update/bacteria_seqinfo.txt"
        # Build config from params
        cfg = Config("update", **params)
        # Should not run - nothing to update
        self.assertFalse(ganon.main(cfg=cfg), "ganon update exited with an error")

    def test_duplicated_input_files(self):
        """
        ganon update with duplicated input files. ganon-build will process all input files, but bins should be correct
        """
        params = self.default_params.copy()
        params["output_db_prefix"] = self.results_dir + "test_duplicated_input_files"
        params["input_files"] = params["input_files"] * 4
        
        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # Unique entries on bins (not duplicated)
        self.assertTrue(res["bins_pd"][["seqid","seqstart","seqend"]].equals(res["bins_pd"][["seqid","seqstart","seqend"]].drop_duplicates()), "Duplicated entries of repeated sequences on bins")

    def test_add_existing_bins(self):
        """
        Test update without creating new bins
        """
        # Build database with 2 sequences at superkingdom level with large bins
        params_build = {"db_prefix": self.results_dir + "test_add_existing_bins_part1",
                        "taxdump_file": [data_dir+"mini_nodes.dmp", 
                                         data_dir+"mini_names.dmp"],
                        "input_files": [data_dir+"update/virus_NC_003676.1.fasta.gz",
                                        data_dir+"update/virus_NC_011646.1.fasta.gz"],
                        "seq_info_file": data_dir+"update/virus_part1_seqinfo.txt",
                        "write_seq_info_file": True,
                        "rank": "superkingdom",
                        "bin_length": 200000,
                        "window_size": 0,
                        "quiet": True}
        # Build config from params
        cfg_build = Config("build", **params_build)
        # Run
        self.assertTrue(ganon.main(cfg=cfg_build), "ganon build exited with an error")
        # General sanity check of results
        res_build = build_sanity_check_and_parse(vars(cfg_build))
        self.assertIsNotNone(res_build, "ganon build has inconsistent results")
        # Copy seqinfo to be parsed later
        shutil.copy(params_build["seq_info_file"], params_build["db_prefix"]+".seqinfo.txt")

        # Update with part2 - add virus to same bin and bacteria to new bins
        params = self.default_params.copy()
        params["db_prefix"] = params_build["db_prefix"]
        params["output_db_prefix"] = self.results_dir + "test_add_existing_bins_part2"
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
        self.assertEqual(res["bins_pd"][res["bins_pd"]["taxid"]=="10239"].binid.drop_duplicates().max(), 0, "virus added to new bins")

        # Check if entries are unique (virus 10239, bacteria 2)
        self.assertEqual(res["bins_pd"].taxid.drop_duplicates().size, 2, "multiple taxids")

        # Classify to first part
        params_classify = {"db_prefix": params_build["db_prefix"],
                           "single_reads": [data_dir+"vir.sim.1.fq", data_dir+"bac.sim.1.fq"],
                           "abs_cutoff": 0,
                           "offset": 1,
                           "output_lca": True,
                           "output_all": True,
                           "quiet": True,
                           "output_prefix": self.results_dir + "test_add_existing_bins_classify_part1"}
        # Build config from params
        cfg_classify = Config("classify", **params_classify)
        # Run
        self.assertTrue(ganon.main(cfg=cfg_classify), "ganon classify exited with an error")
        # General sanity check of results
        res_classify1 = classify_sanity_check_and_parse(vars(cfg_classify))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # Specific test - should contain only one superkingdom (virus 10239, bacteria 2)
        self.assertEqual(res_classify1["tre_pd"][res_classify1["tre_pd"]["rank"] == "superkingdom"]["target"].shape[0], 1, "not one superkingdom as target")
        self.assertEqual(res_classify1["tre_pd"][res_classify1["tre_pd"]["rank"] == "superkingdom"]["target"].values[0], '10239', "wrong target")

        # Classify to second part
        params_classify["db_prefix"] = params["output_db_prefix"]
        params_classify["output_prefix"] = self.results_dir + "test_add_existing_bins_classify_part2"

        # Build config from params
        cfg_classify = Config("classify", **params_classify)

        # Run
        self.assertTrue(ganon.main(cfg=cfg_classify), "ganon classify exited with an error")
        # General sanity check of results
        res_classify2 = classify_sanity_check_and_parse(vars(cfg_classify))

        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # Specific
        # Classification to the second updated index has to have more matches than the first
        self.assertTrue(res_classify2["all_pd"].shape[0] > res_classify1["all_pd"].shape[0], "updated index did not improve matches")

        # should contain only two superkingdoms (virus 10239, bacteria 2)
        self.assertEqual(res_classify2["tre_pd"][res_classify2["tre_pd"]["rank"] == "superkingdom"]["target"].shape[0], 2, "not two superkingdom as target")
        self.assertTrue(res_classify2["tre_pd"][res_classify2["tre_pd"]["rank"] == "superkingdom"]["target"].isin(['10239','2']).all(), "wrong target")

    def test_update_complete_add(self):
        """
        Test run update complete adding sequences only
        """
        params = self.default_params.copy()
        params["output_db_prefix"] = self.results_dir + "test_update_complete_add"
        params["update_complete"] = True
        params["seq_info_file"] = data_dir + "update/bacteria_virus_seqinfo.txt"
        params["input_files"].extend([data_dir + "build/bacteria_NC_010333.1.fasta.gz",
                                      data_dir + "build/bacteria_NC_017164.1.fasta.gz",
                                      data_dir + "build/bacteria_NC_017163.1.fasta.gz",
                                      data_dir + "build/bacteria_NC_017543.1.fasta.gz"])

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
        params_classify = {"db_prefix": params["output_db_prefix"],
                           "single_reads": [data_dir + "vir.sim.1.fq", data_dir + "bac.sim.1.fq"],
                           "abs_cutoff": 0,
                           "output_lca": True,
                           "output_all": True,
                           "quiet": True,
                           "output_prefix": self.results_dir + "test_update_complete_add"}
        # Build config from params
        cfg_classify = Config("classify", **params_classify)
        # Run
        self.assertTrue(ganon.main(cfg=cfg_classify), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg_classify))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # Specific tes - should return Viruses and Bacteria matches on the updated index
        self.assertTrue(res["tre_pd"][res["tre_pd"]["rank"]=="superkingdom"]["name"].isin(["Bacteria","Viruses"]).all(), "classification on updated index failed")
    
    def test_update_complete_remove(self):
        """
        Test run update complete removing sequences only
        """
        params = self.default_params.copy()
        params["output_db_prefix"] = self.results_dir + "test_update_complete_remove"
        params["update_complete"] = True
        params["seq_info_file"] = data_dir+"update/bacteria_half_seqinfo.txt"
        params["input_files"] = [data_dir+"build/bacteria_NC_010333.1.fasta.gz",
                              data_dir+"build/bacteria_NC_017164.1.fasta.gz"]
       
        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # Specific - keep only two entries
        self.assertEqual(res["bins_pd"]["seqid"].drop_duplicates().shape[0],2, "sequences not removed from bins")
        self.assertEqual(res["map_pd"]["target"].drop_duplicates().shape[0],2, "sequences not removed from .map")

        # Classify against reduced updated index
        params_classify = {"db_prefix": params["output_db_prefix"],
                    "single_reads": data_dir+"bac.sim.1.fq",
                    "abs_cutoff": 0,
                    "output_lca": True,
                    "output_all": True,
                    "quiet": True,
                    "output_prefix": self.results_dir + "test_default"}

        # Build config from params
        cfg_classify = Config("classify", **params_classify)
        # Run
        self.assertTrue(ganon.main(cfg=cfg_classify), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg_classify))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # Specific - only matches on remaining sequences (taxids 366602 and 470)
        self.assertTrue(res["all_pd"]["target"].isin(["366602","470"]).all(), "ganon classify has inconsistent results")
        
    def test_update_complete_add_remove(self):
        """
        Test run update complete adding and removing sequences (reusing same bins for new sequences)
        """
        params = self.default_params.copy()
        params["output_db_prefix"] = self.results_dir + "test_update_complete_add_remove"
        params["update_complete"] = True
        params["seq_info_file"] = data_dir+"update/bacteria_half_virus_seqinfo.txt"
        params["input_files"].extend([data_dir+"build/bacteria_NC_010333.1.fasta.gz",
                                      data_dir+"build/bacteria_NC_017164.1.fasta.gz"])
        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # Specific tes - should have 6 taxid targets (2 bacteria, 4 viruses)
        self.assertEqual(res["bins_pd"]["taxid"].drop_duplicates().shape[0], 6, "update failed to add new sequences")
        # Should re-use bins and reach max 18 (41 before)
        self.assertEqual(res["bins_pd"]["binid"].max(), 18, "bins were not re-used")
        self.assertEqual(res["map_pd"]["binid"].max(), 18, "bins were not re-used")

        # Classify against updated index
        params_classify = {"db_prefix": params["output_db_prefix"],
                    "single_reads": [data_dir+"vir.sim.1.fq", data_dir+"bac.sim.1.fq"],
                    "abs_cutoff": 0,
                    "output_lca": True,
                    "output_all": True,
                    "quiet": True,
                    "output_prefix": self.results_dir + "test_default"}

        # Build config from params
        cfg_classify = Config("classify", **params_classify)
        # Run
        self.assertTrue(ganon.main(cfg=cfg_classify), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg_classify))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # Specific - no matches on the removed entries (taxid 1052684)
        self.assertFalse(res["all_pd"]["target"].isin(["1052684"]).any(), "ganon classify has inconsistent results")
        # should return Viruses and Bacteria matches on the updated index
        self.assertTrue(res["tre_pd"][res["tre_pd"]["rank"]=="superkingdom"]["name"].isin(["Bacteria","Viruses"]).all(), "classification on updated index failed")

    def test_update_multiple(self):
        """
        Test multiple update runs: 1) only remove 2) only add (reusing bins) 3) remove and add
        """

        #Remove only (2 bacteria entries)
        params = self.default_params.copy()
        params["output_db_prefix"] = self.results_dir + "test_update_multiple_1"
        params["update_complete"] = True
        params["seq_info_file"] = data_dir + "update/bacteria_half_seqinfo.txt"
        params["input_files"] = [data_dir + "build/bacteria_NC_010333.1.fasta.gz",
                                 data_dir + "build/bacteria_NC_017164.1.fasta.gz"]

        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # Specific - keep only two entries
        self.assertEqual(res["bins_pd"]["seqid"].drop_duplicates().shape[0],2, "sequences not removed from bins")
        self.assertEqual(res["map_pd"]["target"].drop_duplicates().shape[0],2, "sequences not removed from .map")

        # Add only (2 viruses entries, not update complete)
        params["db_prefix"] = params["output_db_prefix"]
        params["update_complete"] = False
        params["output_db_prefix"] = self.results_dir + "test_update_multiple_2"
        params["seq_info_file"] = data_dir+"update/virus_part1_seqinfo.txt"
        params["input_files"] = [data_dir+"update/virus_NC_003676.1.fasta.gz",
                                data_dir+"update/virus_NC_011646.1.fasta.gz"]
        # Copy seqinfo to be parsed later
        shutil.copy(params["seq_info_file"], params["db_prefix"]+".seqinfo.txt")
        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")

        # Add viruses and remove all bacteria
        params["db_prefix"] = params["output_db_prefix"]
        params["update_complete"] = True
        params["output_db_prefix"] = self.results_dir + "test_update_multiple_3"
        params["seq_info_file"] = data_dir+"update/virus_seqinfo.txt"
        params["input_files"] = [data_dir+"update/virus_NC_003676.1.fasta.gz",
                                  data_dir+"update/virus_NC_011646.1.fasta.gz", 
                                  data_dir+"update/virus_NC_032412.1.fasta.gz", 
                                  data_dir+"update/virus_NC_035470.1.fasta.gz"]
        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")

        # At the end, should have only viruses on the index
        params_classify = {"db_prefix": params["output_db_prefix"],
                    "single_reads": [data_dir+"vir.sim.1.fq", data_dir+"bac.sim.1.fq"],
                    "abs_cutoff": 0,
                    "output_lca": True,
                    "output_all": True,
                    "quiet": True,
                    "output_prefix": self.results_dir + "test_update_multiple_classify"}

        # Build config from params
        cfg_classify = Config("classify", **params_classify)
        # Run
        self.assertTrue(ganon.main(cfg=cfg_classify), "ganon classify exited with an error")
        # General sanity check of results
        res = classify_sanity_check_and_parse(vars(cfg_classify))
        self.assertIsNotNone(res, "ganon classify has inconsistent results")
        # should not contain any bacteria
        self.assertFalse(res["tre_pd"][res["tre_pd"]["rank"]=="superkingdom"]["name"].isin(["Bacteria"]).any(), "index was not properly updated, bacteria sequences remain")
      
    def test_duplicated_seqinfo(self):
        """
        ganon update 
        """
        params = self.default_params.copy()
        params["output_db_prefix"] = self.results_dir + "test_duplicated_seqinfo"
        params["seq_info_file"] = data_dir+"update/virus_seqinfo_duplicated.txt"

        # Build config from params
        cfg = Config("update", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon update exited with an error")
        # General sanity check of results
        res = update_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")
        # ganon should remove the duplicates and just have unique entries on bins
        self.assertTrue(res["bins_pd"][["seqid","seqstart","seqend"]].equals(res["bins_pd"][["seqid","seqstart","seqend"]].drop_duplicates()), "Duplicated entries on bins")


if __name__ == '__main__':
    unittest.main()