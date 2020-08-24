import unittest, shlex, pickle, sys, os, shutil
from pathlib import Path
import pandas as pd
sys.path.append('src')
from ganon import ganon
from ganon.build_update import build
from ganon.bins import Bins
from ganon.config import Config
from ganon.gnn import Gnn
from ganon.tax import Tax

base_dir = "tests/ganon/integration/"
data_dir = base_dir + "data/"

class TestOffline(unittest.TestCase):
        
    results_dir = base_dir + "results/test_build/TestOffline/"
    default_params = {"taxdump_file": [data_dir+"mini_nodes.dmp", data_dir+"mini_names.dmp"],
                      "input_files": [data_dir+"bacteria_NC_010333.1.fasta.gz", data_dir+"bacteria_NC_017164.1.fasta.gz", data_dir+"bacteria_NC_017163.1.fasta.gz", data_dir+"bacteria_NC_017543.1.fasta.gz"],
                      "seq_info_file": data_dir+"bacteria_seqinfo.txt",
                      "write_seq_info_file": True,
                      "rank": "species"}
    
    @classmethod
    def setUpClass(self):
        setup_build(self.results_dir)
       
    def test_default(self):
        """
        With default parameters
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "TestOffline_test_default"
        
        # Build config from params
        cfg = Config("build", **params)
        # Validate
        self.assertTrue(cfg.validate(), "Invalid configuration")
        # Run
        self.assertTrue(build(cfg), "ganon build exited with an error")
        # General sanity check of results
        res = sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")

    def test_default_main(self):
        """
        With default parameters going through ganon.main()
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "TestOffline_test_default_main"
        # Run
        self.assertTrue(ganon.main("build", **params), "ganon build exited with an error")
        # General sanity check of results
        res = sanity_check_and_parse(params)
        self.assertIsNotNone(res, "ganon build has inconsistent results")

    def test_assembly(self):
        """
        Default but rank=assembly
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "TestOffline_test_assembly"
        params["rank"] = "assembly"
        
        # Build config from params
        cfg = Config("build", **params)
        # Validate
        self.assertTrue(cfg.validate(), "Invalid configuration")
        # Run
        self.assertTrue(build(cfg), "ganon build exited with an error")
        # General sanity check of results
        res = sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon build has inconsistent results")
        # Specific test
        print(res)

class TestOnline(unittest.TestCase):
    
    results_dir = base_dir + "results/test_build/TestOnline/"
    default_params = {"input_files": [data_dir+"bacteria_NC_010333.1.fasta.gz", data_dir+"bacteria_NC_017164.1.fasta.gz", data_dir+"bacteria_NC_017163.1.fasta.gz", data_dir+"bacteria_NC_017543.1.fasta.gz"],
                      "write_seq_info_file": True,
                      "rank": "species"}

    @classmethod
    def setUpClass(self):
        setup_build(self.results_dir)

    def test_default(self):
        """
        With default parameters
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "TestOnline_test_default"
        
        # Build config from params
        cfg = Config("build", **params)
        # Validate
        self.assertTrue(cfg.validate(), "Invalid configuration")
        # Run
        self.assertTrue(build(cfg), "ganon build exited with an error")
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
    res["tax"] = Tax(tax_files=[params["db_prefix"]+".tax"])
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

    # Check if all taxids (chosen rank) are present in the .tax
    if not res["tax_pd"]['taxid'].apply(lambda x: res["tax"].get_rank(x, params["rank"])).isin(res["tax_pd"]["taxid"]).all():
        print("Missing taxonomic entries")
        return None
    
    # Check if all taxids (chosen rank) are present in the .tax
    if not res["tax_pd"]['taxid'].apply(lambda x: res["tax"].get_rank(x, params["rank"])).isin(res["tax_pd"]["taxid"]).all():
        print("Missing taxonomic entries")
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