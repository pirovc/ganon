import unittest, sys
sys.path.append('src')
from ganon import ganon
from ganon.bins import Bins
from ganon.config import Config
from ganon.gnn import Gnn
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
        res = sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")

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
        res = sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")

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
        res = sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")

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
        res = sanity_check_and_parse(vars(cfg))
        self.assertIsNotNone(res, "ganon update has inconsistent results")

def sanity_check_and_parse(params):
    # Provide sanity checks for outputs (not specific to a test) and return loaded data

    if not check_files(params["output_db_prefix"], ["ibf", "map", "tax", "gnn"]):
        return None

    res = {}
    # Sequence information from database to be updated
    res["seq_info"] =  parse_seq_info(params["db_prefix"]+".seqinfo.txt")
    # Parse in and out files
    if "seq_info_file" in params and params["seq_info_file"]:
        res["seq_info"] = res["seq_info"].append(parse_seq_info(params["seq_info_file"]), ignore_index=True)
    else:
        res["seq_info"] = res["seq_info"].append(parse_seq_info(params["output_db_prefix"]+".seqinfo.txt"), ignore_index=True)

    res["gnn"] = Gnn(file=params["output_db_prefix"]+".gnn")
    res["bins"] = Bins(taxsbp_ret=res["gnn"].bins)
    res["tax_pd"] = parse_tax(params["output_db_prefix"]+".tax")
    res["map_pd"] = parse_map(params["output_db_prefix"]+".map")
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

if __name__ == '__main__':
    unittest.main()