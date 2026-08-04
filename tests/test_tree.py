import unittest

from dnctree.tree import Tree


class TestTree(unittest.TestCase):
    def setUp(self):
        self.t = Tree()
        self.t.add_edge('a', '1')
        self.t.add_edge('b', '1')
        self.t.add_edge('1', '2')
        self.t.add_edge('d', '3')
        self.t.add_edge('3', '4')
        self.t.add_edge('2', '3')
        self.t.add_edge('4', 'e')
        self.t.add_edge('f', '4')
        self.t.add_edge('c', '2')
        self.t.set_start_node('1')
        
        self.t1 = Tree()
        self.t1.add_edge('a', '1')
        self.t1.add_edge('b', '1')
        self.t1.add_edge('1', '4')

        self.t2 = Tree()
        self.t2.add_edge('c', '2')
        self.t2.add_edge('d', '2')
        self.t2.add_edge('2', '4')

        self.t3 = Tree()
        self.t3.add_edge('e', '3')
        self.t3.add_edge('f', '3')
        self.t3.add_edge('3', '4')


    def test_building(self):
        self.assertEqual(str(self.t), '(a,b,((d,(e,f)),c));')

        self.t.set_start_node('3')
        self.assertEqual(str(self.t), '(d,(e,f),((a,b),c));')


    def test_merge(self):
        self.t1.merge(self.t2, '4')
        self.t1.merge(self.t3, '4')

        self.t1.set_start_node('4')
        self.assertEqual(str(self.t1), '((a,b),(c,d),(e,f));')
