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
        Test if build work with default paramenters
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
        self.assertTrue(sanity_check(vars(cfg)), "ganon build has inconsistent results")

    def test_default_main(self):
        """
        Test if build work with default paramenters
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "TestOffline_test_default_main"
        # Run
        self.assertTrue(ganon.main("build", **params), "ganon build exited with an error")
        # General sanity check of results
        self.assertTrue(sanity_check(params), "ganon build has inconsistent results")



    def test_default_assembly(self):
        """
        Test if build work with default paramenters
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "TestOffline_test_default_assembly"
        params["rank"] = "assembly"
        # Run
        self.assertTrue(ganon.main("build", **params), "ganon build exited with an error")
        # General sanity check of results
        self.assertTrue(sanity_check(params), "ganon build has inconsistent results")
        # 

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
        Test if build work with default paramenters without aux files
        """
        params = self.default_params.copy()
        params["db_prefix"] = self.results_dir + "TestOnline_test_default"
        # Run
        self.assertTrue(ganon.main("build", **params), "ganon build exited with an error")
        # General sanity check of results
        self.assertTrue(sanity_check(params), "ganon build has inconsistent results")

def setup_build(results_dir):
    shutil.rmtree(results_dir, ignore_errors=True)
    os.makedirs(results_dir)

def sanity_check(params):
    # Check if files were created
    for ext in ["ibf", "map", "tax", "gnn"]:
        f=params["db_prefix"]+"."+ext
        if not Path(f).is_file():
            print("File (" + f +") was not created")
            return False
        elif Path(f).stat().st_size==0:
            print("File (" + f +") is empty")
            return False

    # Parse in and out files
    if "seq_info_file" in params and params["seq_info_file"]:
        in_seq_info = parse_seq_info(params["seq_info_file"])
    else:
        in_seq_info = parse_seq_info(params["db_prefix"]+".seqinfo.txt")
         
    out_map = parse_map(params["db_prefix"]+".map")
    out_gnn = Gnn(file=params["db_prefix"]+".gnn")
    out_tax = parse_tax(params["db_prefix"]+".tax")
    out_bins = Bins(taxsbp_ret=out_gnn.bins).bins
    tax = Tax(tax_files=[params["db_prefix"]+".tax"])

    # Check number of bins
    if out_map.binid.unique().size != out_gnn.number_of_bins:
        print("Number of bins do not match between .gnn and .map")
        return False

    # Check if all input accession made it to the bins
    if not in_seq_info["seqid"].isin(out_bins["seqid"]).all():
        print("Missing sequence accessions on bins")
        return False

    # Check if all taxids (chosen rank) are present in the .tax
    if not out_tax['taxid'].apply(lambda x: tax.get_rank(x, params["rank"])).isin(out_tax["taxid"]).all():
        print("Missing taxonomic entries")
        return False

    return True

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