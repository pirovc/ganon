import unittest, shlex, pickle, sys, os, shutil
from pathlib import Path
import pandas as pd
sys.path.append('src')
from ganon import ganon
from ganon.bins import Bins
from ganon.config import Config
from ganon.gnn import Gnn
from ganon.tax import Tax

base_dir = "tests/ganon/integration/"
data_dir = base_dir + "data/"

class TestBuild(unittest.TestCase):
        
    results_dir = base_dir + "results/test_build/"
    default_params = {"taxdump_file": [data_dir+"mini_nodes.dmp", data_dir+"mini_names.dmp"],
                      "input_files": [data_dir+"bacteria_NC_010333.1.fasta.gz", data_dir+"bacteria_NC_017164.1.fasta.gz", data_dir+"bacteria_NC_017163.1.fasta.gz", data_dir+"bacteria_NC_017543.1.fasta.gz"],
                      "seq_info_file": data_dir+"bacteria_seqinfo.txt",
                      "write_seq_info_file": True,
                      "rank": "species",
                      "quiet": True}
    
    @classmethod
    def setUpClass(self):
        setup_build(self.results_dir)
       
    def test_default(self):
        """
        Test run with default parameters
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_default"
        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")

    def test_assembly(self):
        """
        Test rank as assembly
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_assembly"
        params["rank"] = "assembly"
        
        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = sanity_check_and_parse(vars(cfg))
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
        res = sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
        # Specific test
        # Check max size of fragments on bins
        self.assertTrue(max(res["bins_pd"]["length"])<=params["bin_length"], "Bin length greater than max.")

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
        res = sanity_check_and_parse(vars(cfg))
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
        res = sanity_check_and_parse(vars(cfg))
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
        res = sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
        # Specific test
        # Check max size of fragments on bins
        self.assertTrue(max(res["bins_pd"]["length"])<=params["fragment_length"]+params["overlap_length"], "Fragment greater than max.")
        # Check max size of bins
        self.assertTrue(max(res["bins_pd"].groupby("binid").sum()["length"])<=params["bin_length"], "Bin length greater than max.")

    def test_duplicated_entries(self):
        """
        Test duplicated entries
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_duplicated_entries"
        params["input_files"] = params["input_files"] * 4

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        res = sanity_check_and_parse(vars(cfg))
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
        res = sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")

    def test_duplicated_entries_seqinfo(self):
        """
        Test duplicated entries on the seqinfo file
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "test_duplicated_entries_seqinfo"
        params["seq_info_file"] = data_dir+"bacteria_seqinfo_duplicated.txt"
        #params["verbose"] = True

        # Build config from params
        cfg = Config("build", **params)
        # Run
        self.assertTrue(ganon.main(cfg=cfg), "ganon build exited with an error")
        # General sanity check of results
        # Sanity check passes because ganon build (using --seq-info-file)
        # does not iterate through sequences, just works with metadata
        # ganon-build will simply add the unique provided entries
        res = sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")

class TestBuildOnline(unittest.TestCase):
    
    results_dir = base_dir + "results/test_build/online/"
    default_params = {"input_files": [data_dir+"bacteria_NC_010333.1.fasta.gz", data_dir+"bacteria_NC_017164.1.fasta.gz", data_dir+"bacteria_NC_017163.1.fasta.gz", data_dir+"bacteria_NC_017543.1.fasta.gz"],
                      "write_seq_info_file": True,
                      "rank": "species",
                      "quiet": True}

    @classmethod
    def setUpClass(self):
        setup_build(self.results_dir)

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
        res = sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
        
def setup_build(results_dir):
    shutil.rmtree(results_dir, ignore_errors=True)
    os.makedirs(results_dir)

def sanity_check_and_parse(params):
    # Provide sanity checks for outputs (not specific to a test) and return loaded data

    # Check if files were created
    for ext in ["ibf", "map", "tax", "gnn"]:
        f=params["db_prefix"]+"."+ext
        if not Path(f).is_file():
            print("File (" + f +") was not created")
            return None
        elif Path(f).stat().st_size==0:
            print("File (" + f +") is empty")
            return None

    res = {}
    # Parse in and out files
    if "seq_info_file" in params and params["seq_info_file"]:
        res["seq_info"] = parse_seq_info(params["seq_info_file"])
    else:
        res["seq_info"] = parse_seq_info(params["db_prefix"]+".seqinfo.txt")
         
    res["gnn"] = Gnn(file=params["db_prefix"]+".gnn")
    res["bins"] = Bins(taxsbp_ret=res["gnn"].bins)
    res["tax_pd"] = parse_tax(params["db_prefix"]+".tax")
    res["map_pd"] = parse_map(params["db_prefix"]+".map")
    res["bins_pd"] = res["bins"].bins

    # Check number of bins
    if res["map_pd"].binid.unique().size != res["gnn"].number_of_bins:
        print("Number of bins do not match between .gnn and .map")
        return None

    # Check if all input accession made it to the bins
    if not res["seq_info"]["seqid"].isin(res["bins_pd"]["seqid"]).all():
        print("Missing sequence accessions on bins")
        return None

    # Check if all taxids/assembly on .map appear on .tax
    if res["tax_pd"]["taxid"].isin(res["map_pd"]["target"].drop_duplicates()).all():
        print("Inconsistent entries between taxonomy (.tax) and bin map (.map)")
        return None

    return res

def parse_seq_info(seq_info_file):
    colums=['seqid', 'length', 'taxid', 'specialization']
    types={'seqid': 'str', 'length': 'uint64', 'taxid': 'str', 'specialization': 'str'}
    return pd.read_table(seq_info_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def parse_map(map_file):
    colums=['target', 'binid']
    types={'target': 'str', 'binid': 'uint64'}
    return pd.read_table(map_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def parse_tax(tax_file):
    colums=['taxid', 'parent', 'rank', 'name']
    types={'taxid': 'str', 'parent': 'str', 'rank': 'str', 'name': 'str'}
    return pd.read_table(tax_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

if __name__ == '__main__':
    unittest.main()