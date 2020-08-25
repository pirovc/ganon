import unittest, shlex, pickle, sys
from pathlib import Path
sys.path.append('src')
from ganon import ganon

class TestUpdate(unittest.TestCase):
    def test_update(self):
        """
        Test if update on sample data is working
        """
        path_data = "tests/ganon/integration/data/"
        prefix = "test_update"

        ret = ganon.main("update", 
                            db_prefix=path_data+"sample_bacteria", 
                            output_db_prefix=prefix,
                            taxdump_file=[path_data+"mini_nodes.dmp", path_data+"mini_names.dmp"], 
                            seq_info_file=path_data+"virus_acc_len_taxid.txt", 
                            input_files=[path_data+"virus_NC_003676.1.fasta.gz", path_data+"virus_NC_011646.1.fasta.gz", path_data+"virus_NC_032412.1.fasta.gz", path_data+"virus_NC_035470.1.fasta.gz"],
                            quiet=True)

        # check if ran okay
        self.assertTrue(ret, "ganon update finish with an error")
       
        # check if files were created
        for ext in ["ibf", "map", "tax", "gnn"]:
            self.assertTrue(Path(prefix+"."+ext).is_file() , "File (" + ext +") was not created") # TODO check file contents

if __name__ == '__main__':
    unittest.main()