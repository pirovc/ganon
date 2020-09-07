import unittest, sys
sys.path.append('src')
from ganon.tax import Tax

base_dir = "tests/ganon/"
sys.path.append(base_dir)
from utils import *
data_dir = base_dir + "data/"

class TestTax(unittest.TestCase):
    def test_default(self):
        """
        Test if read nodes is working properly
        """
        tax = Tax(ncbi_nodes=data_dir + "mini_nodes.dmp", ncbi_names=data_dir + "mini_names.dmp")
        self.assertEqual(len(tax.nodes), 36)
        self.assertEqual(tax.nodes["2"][0], "131567")
        self.assertEqual(tax.nodes["131567"][0], "1")
        self.assertEqual(tax.nodes["1"][0], "0")
        self.assertEqual(tax.nodes["2"][2], "Bacteria")
    
if __name__ == '__main__':
    unittest.main()
