import unittest, sys
sys.path.append('src')
from ganon import ganon
from ganon.config import Config

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import *
data_dir = base_dir + "data/"

class TestBuildOffline(unittest.TestCase):

    results_dir = base_dir + "results/integration/build/"
    default_params = {"taxdump_file": [data_dir+"mini_nodes.dmp", 
                                       data_dir+"mini_names.dmp"],
                      "input_files": [data_dir+"build/bacteria_NC_010333.1.fasta.gz",
                                      data_dir+"build/bacteria_NC_017164.1.fasta.gz", 
                                      data_dir+"build/bacteria_NC_017163.1.fasta.gz", 
                                      data_dir+"build/bacteria_NC_017543.1.fasta.gz"],
                      "seq_info_file": data_dir+"build/bacteria_seqinfo.txt",
                      "write_seq_info_file": True,
                      "rank": "species",
                      "quiet": True}
    
    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)
       
    def test_default(self):
        """
        ganon build with default parameters
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

    def test_rank_genus(self):
        """
        ganon build --rank genus
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_rank_genus"
        params["rank"] = "genus"

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
        # Check if did not report species from .tax
        self.assertFalse((res["tax_pd"]["rank"]=="species").any(), "Species reported with rank genus")

    def test_rank_leaves(self):
        """
        ganon build --rank leaves
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_rank_genus"
        params["rank"] = "leaves"
        
        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
        # Check if taxids provides are exactly the same used in .map
        self.assertTrue(res["seq_info"]["taxid"].isin(res["map_pd"]["target"].drop_duplicates()).all(), "Did not use leaves on .map")


    def test_specialization_custom(self):
        """
        ganon build --specialization custom (with --seq-info-file)
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_specialization_custom"
        params["specialization"] = "custom"
                
        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
        # Specific test
        # Check if all assembly ids are on bins and map and tax
        self.assertTrue(res["seq_info"]["specialization"].isin(res["bins_pd"]["specialization"]).all(), "Missing assembly ids on bins")
        self.assertTrue(res["seq_info"]["specialization"].isin(res["map_pd"]["target"].drop_duplicates()).all(), "Missing assembly ids on .map")
        self.assertTrue(res["seq_info"]["specialization"].isin(res["tax_pd"]["taxid"].drop_duplicates()).all(), "Missing assembly ids on .tax")

    def test_specialization_custom(self):
        """
        ganon build --specialization custom (with --seq-info-file)
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_specialization_custom"
        params["specialization"] = "custom"
                
        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
        # Specific test
        # Check if all assembly ids are on bins and map and tax
        self.assertTrue(res["seq_info"]["specialization"].isin(res["bins_pd"]["specialization"]).all(), "Missing assembly ids on bins")
        self.assertTrue(res["seq_info"]["specialization"].isin(res["map_pd"]["target"].drop_duplicates()).all(), "Missing assembly ids on .map")
        self.assertTrue(res["seq_info"]["specialization"].isin(res["tax_pd"]["taxid"].drop_duplicates()).all(), "Missing assembly ids on .tax")

    def test_bin_length(self):
        """
        Test changing bin length
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_bin_length"
        params["bin_length"] = 10000
        
        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
        # Specific test
        # Check max size of fragments on bins
        self.assertTrue(res["bins_pd"]["length"].max()<=params["bin_length"], "Bin length greater than max.")

    def test_fragment_length(self):
        """
        Test changing fragment length
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_fragment_length"
        params["fragment_length"] = 5000
        
        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
        # Specific test
        # Check max size of fragments on bins
        self.assertTrue(max(res["bins_pd"]["length"])<=params["fragment_length"], "Bin length greater than max.")

    def test_overlap_length(self):
        """
        Test changing overlap length
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_overlap_length"
        params["bin_length"] = 10000
        params["overlap_length"] = 999
        
        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
        # Specific test
        # Check max size of fragments on bins
        self.assertTrue(max(res["bins_pd"]["length"])<=params["bin_length"]+params["overlap_length"], "Fragment bigger than max. set")

    def test_bin_fragment_overlap_length(self):
        """
        Test changing bin, fragment and overlap length
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_bin_fragment_overlap_length"
        params["bin_length"] = 5692
        params["fragment_length"] = 667
        params["overlap_length"] = 349
        
        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
        # Specific test
        # Check max size of fragments on bins
        self.assertTrue(max(res["bins_pd"]["length"])<=params["fragment_length"]+params["overlap_length"], "Fragment greater than max.")
        # Check max size of bins
        self.assertTrue(max(res["bins_pd"].groupby("binid").sum()["length"])<=params["bin_length"], "Bin length greater than max.")

    def test_duplicated_input_files(self):
        """
        ganon build with duplicated input files. ganon-build will process all input files, but bins should be correct
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_duplicated_input_files"
        params["input_files"] = params["input_files"] * 4

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
        # Specific test
        # Unique entries on bins (not duplicated)
        self.assertTrue(res["bins_pd"][["seqid","seqstart","seqend"]].equals(res["bins_pd"][["seqid","seqstart","seqend"]].drop_duplicates()), "Duplicated entries of repeated sequences on bins")

    def test_missing_entries(self):
        """
        Test missing entries
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_missing_entries"
        params["input_files"] = params["input_files"][0]

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        # Sanity check passes because ganon build (using --seq-info-file)
        # does not iterate through sequences, just works with metadata
        # ganon-build will simply leave bins empty
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")

    def test_duplicated_seqinfo(self):
        """
        ganon build with duplicated --seq-info-file entries
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_duplicated_seqinfo"
        params["seq_info_file"] = data_dir+"build/bacteria_seqinfo_duplicated.txt"

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
         # ganon should remove the duplicates and just have unique entries on bins
        self.assertTrue(res["bins_pd"][["seqid","seqstart","seqend"]].equals(res["bins_pd"][["seqid","seqstart","seqend"]].drop_duplicates()), "Duplicated entries on bins")

    def test_missing_invalid_seqinfo(self):
        """
        ganon build --seq-info-file with missing and invalid information
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_missing_invalid_seqinfo"
        params["seq_info_file"] = data_dir+"build/bacteria_seqinfo_missing_invalid.txt"

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")

        res = build_sanity_check_and_parse(vars(cfg))
        # General sanity check of results
        # Sanity check passes because ganon build (using --seq-info-file)
        # does not iterate through sequences, just works with metadata
        # ganon-build will simply add the unique provided entries
        self.assertIsNotNone(res, "ganon build has inconsistent results")

    def test_input_directory(self):
        """
        Test duplicated entries on the seqinfo file
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_input_directory"
        del params["input_files"]
        params["input_directory"] = data_dir+"build/"
        params["input_extension"] = ".fasta.gz"
       
        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")

    def test_invalid_rank(self):
        """
        ganon build --rank xyz (invalid)
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_invalid_rank"
        params["rank"] = "xyz"

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
        
if __name__ == '__main__':
    unittest.main()