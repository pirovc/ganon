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
                      "taxdump_file": [data_dir+"mini_nodes.dmp", 
                                      data_dir+"mini_names.dmp"],
                      "write_seq_info_file": True,
                      "rank": "species",
                      "quiet": True}

    @classmethod
    def setUpClass(self):
        setup_dir(self.results_dir)

    def test_default(self):
        """
        ganon build with default parameters (online: eutils, taxdump)
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_default"
        params["taxdump_file"] = []

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
       
    def test_specialization_assembly(self):
        """
        ganon build --specialization assembly (online: eutils)
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_specialization_assembly"
        params["specialization"] = "assembly"

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results") 
        # Specific test - count assemblies on tax (3 bac)
        self.assertEqual(sum(res["tax_pd"]["rank"]=="assembly"), 3, "error retrieving assembly accessions")
        # Check if all targets starts with "GCF_"
        self.assertTrue((res["map_pd"]["target"].map(lambda x: x.startswith("GCF_"))).all(), "failed to retrieve assembly accession")

    def test_specialization_sequence(self):
        """
        ganon build --specialization sequence (online: eutils)
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_specialization_sequence"
        params["specialization"] = "sequence"

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results") 
        # Specific test - count sequences
        self.assertEqual(sum(res["tax_pd"]["rank"]=="sequence"), 4, "failed to use sequence accession as specialization")
        # Check if all targets starts with "NC_"
        self.assertTrue((res["map_pd"]["target"].map(lambda x: x.startswith("NC_"))).all(), "failed to use sequence accession as specialization")

    def test_specialization_file(self):
        """
        ganon build --specialization file (online: eutils)
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_specialization_file"
        params["specialization"] = "file"

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results") 
        # Specific test - count files
        self.assertEqual(sum(res["tax_pd"]["rank"]=="file"), 4, "failed to use file name as specialization")
        # Check if all targets ends with ".fasta.gz"
        self.assertTrue((res["map_pd"]["target"].map(lambda x: x.endswith(".fasta.gz"))).all(), "failed to use file name as specialization")

    def test_specialization_file_single(self):
        """
        ganon build --specialization file (with one file only online: eutils)
        """
        params = self.default_params.copy()

        merge_gz(params["input_files"], self.results_dir + "merged_input_files.fasta.gz")
        params["input_files"] = self.results_dir + "merged_input_files.fasta.gz"
        params["db_prefix"] = self.results_dir + "test_specialization_file_single"
        params["specialization"] = "file"

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results") 
        # Specific test - count files
        self.assertEqual(sum(res["tax_pd"]["rank"]=="file"), 4, "failed to use file name as specialization")
        # Check if all targets starts with "NC_" - fails to use file specialization and replaces it with sequence accession
        self.assertTrue((res["map_pd"]["target"].map(lambda x: x.startswith("NC_"))).all(), "failed to use sequence accession as specialization")

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

    def test_invalid_ncbi_header(self):
        """
        ganon build --specialization sequence with one invalid entry (online: eutils)
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_invalid_ncbi_header"
        params["input_files"].append(data_dir+"build/invalid_ncbi_header.fasta.gz")
        params["specialization"] = "sequence"

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = build_sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results") 
        # Specific test - count sequences on tax (4 files)
        self.assertEqual(sum(res["tax_pd"]["rank"]=="sequence"), 4, "Wrong number of sequences")


if __name__ == '__main__':
    unittest.main()
