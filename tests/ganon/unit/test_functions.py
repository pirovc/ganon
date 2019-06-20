import unittest

from src.ganon.ganon import read_nodes, read_names

class TestReadNodes(unittest.TestCase):
    def test_read_nodes(self):
        """
        Test if read nodes is working properly
        """
        nodes_file = ("tests/ganon/unit/data/mini_nodes.dmp")
        nodes, ranks = read_nodes(nodes_file)
        self.assertEqual(len(nodes), 34)
        self.assertEqual(len(ranks), 34)
        self.assertEqual(nodes["1"], "0")
    
    def test_read_names(self):
        """
        Test if read names is working properly
        """
        names_file = ("tests/ganon/unit/data/mini_names.dmp")
        names = read_names(names_file)
        self.assertEqual(len(names), 34)
        self.assertEqual(names["2"], "Bacteria")

if __name__ == '__main__':
    unittest.main()