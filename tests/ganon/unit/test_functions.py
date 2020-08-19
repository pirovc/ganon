import unittest

from src.ganon.tax import Tax

class TestTax(unittest.TestCase):
    def testTax(self):
        """
        Test if read nodes is working properly
        """
        tax = Tax(ncbi_nodes="tests/ganon/unit/data/mini_nodes.dmp", ncbi_names="tests/ganon/unit/data/mini_names.dmp")
        self.assertEqual(len(tax.nodes), 34)
        self.assertEqual(tax.nodes["2"][0], "131567")
        self.assertEqual(tax.nodes["131567"][0], "1")
        self.assertEqual(tax.nodes["1"][0], "0")
        self.assertEqual(tax.nodes["2"][2], "Bacteria")
    

if __name__ == '__main__':
    unittest.main()