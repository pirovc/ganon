import unittest
import shlex, pickle
from pathlib import Path
from src.ganon import ganon

class TestBuild(unittest.TestCase):
    def test_build(self):
        """
        Test if build on sample data is working
        """
        prefix = "test_output"
        ret = ganon.main(shlex.split("ganon build --db-prefix " + prefix + " --taxdump-file tests/ganon/integration/data/mini_nodes.dmp tests/ganon/integration/data/mini_names.dmp --input-files tests/ganon-build/data/sequences/bacteria_NC_010333.1.fasta.gz tests/ganon-build/data/sequences/bacteria_NC_017164.1.fasta.gz tests/ganon-build/data/sequences/bacteria_NC_017163.1.fasta.gz tests/ganon-build/data/sequences/bacteria_NC_017543.1.fasta.gz"))
        
        # check if ran okay
        self.assertFalse(ret, "ganon build finish with an error")
       
        # check if files were created
        for ext in ["filter", "map", "bins", "nodes"]:
            self.assertTrue(Path(prefix+"."+ext).is_file() , "File (" + ext +") was not created")

        # TODO open pickle .nodes and compare variables, compare number of bins, etc...

if __name__ == '__main__':
    unittest.main()