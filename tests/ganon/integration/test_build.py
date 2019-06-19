import unittest

import shlex, filecmp, pickle

from src.ganon import ganon

class TestBuild(unittest.TestCase):
    def test_build(self):
        """
        Test if build on sample data is working
        """
        prefix = "test_output"
        ganon.main(shlex.split("ganon build --db-prefix " + prefix + " --taxdump-file tests/ganon/integration/data/mini_nodes.dmp tests/ganon/integration/data/mini_names.dmp --input-files tests/ganon-build/data/sequences/bacteria_NC_010333.1.fasta.gz tests/ganon-build/data/sequences/bacteria_NC_017164.1.fasta.gz tests/ganon-build/data/sequences/bacteria_NC_017163.1.fasta.gz tests/ganon-build/data/sequences/bacteria_NC_017543.1.fasta.gz"))
        
        # bins can be different given python version, compare metadata (size filter, number of bins, etc)
        for ext in ["filter", "map"]:
            self.assertTrue(filecmp.cmp(prefix+"."+ext,"tests/ganon/integration/data/sample_bacteria."+ext, shallow=False), "File (" + ext +") is not equal")

        # TODO open pickle .nodes and compare variables

if __name__ == '__main__':
    unittest.main()