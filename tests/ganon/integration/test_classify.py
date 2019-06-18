# import unittest

# from src.ganon.ganon import read_nodes

# class TestReadNodes(unittest.TestCase):
#     def test_read_nodes(self):
#         """
#         Test if read nodes is working properly
#         """
#         nodes_file = ("tests/ganon/unit/data/mini_nodes.dmp")
#         nodes, ranks = read_nodes(nodes_file)
#         self.assertEqual(len(nodes), 101)
#         self.assertEqual(len(ranks), 101)
#         self.assertEqual(nodes["1"], "0")

# if __name__ == '__main__':
#     unittest.main()