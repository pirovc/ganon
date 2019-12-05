import unittest
import shlex, pickle
from pathlib import Path
from src.ganon import ganon

class TestBuild(unittest.TestCase):
    def test_build(self):
        """
        Test if build on sample data is working
        """
        path_data = "tests/ganon/integration/data/"
        prefix = "test_build"
        ret = ganon.main(shlex.split("ganon build --db-prefix "+prefix+" --taxdump-file "+path_data+"mini_nodes.dmp "+path_data+"mini_names.dmp --seq-info-file "+path_data+"bacteria_acc_len_taxid.txt --input-files "+path_data+"bacteria_NC_010333.1.fasta.gz "+path_data+"bacteria_NC_017164.1.fasta.gz "+path_data+"bacteria_NC_017163.1.fasta.gz "+path_data+"bacteria_NC_017543.1.fasta.gz"))
        
        # check if ran okay
        self.assertFalse(ret, "ganon build finish with an error")
       
        # check if files were created
        for ext in ["ibf", "map", "tax", "gnn"]:
            self.assertTrue(Path(prefix+"."+ext).is_file() , "File (" + ext +") was not created") # TODO check file contents

if __name__ == '__main__':
    unittest.main()