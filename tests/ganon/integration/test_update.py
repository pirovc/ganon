import unittest
import shlex, pickle
from pathlib import Path
from src.ganon import ganon

class TestUpdate(unittest.TestCase):
    def test_build(self):
        """
        Test if update on sample data is working
        """
        path_data = "tests/ganon/integration/data/"
        prefix = "test_update"
        ret = ganon.main(shlex.split("ganon update --db-prefix "+path_data+"sample_bacteria --output-db-prefix "+prefix+" --taxdump-file "+path_data+"mini_nodes.dmp "+path_data+"mini_names.dmp --seq-info-file "+path_data+"virus_acc_len_taxid.txt --input-files "+path_data+"virus_NC_003676.1.fasta.gz "+path_data+"virus_NC_011646.1.fasta.gz "+path_data+"virus_NC_032412.1.fasta.gz "+path_data+"virus_NC_035470.1.fasta.gz"))
        
        # check if ran okay
        self.assertFalse(ret, "ganon update finish with an error")
       
        # check if files were created
        for ext in ["filter", "map", "bins", "nodes"]:
            self.assertTrue(Path(prefix+"."+ext).is_file() , "File (" + ext +") was not created") # TODO check file contents

if __name__ == '__main__':
    unittest.main()