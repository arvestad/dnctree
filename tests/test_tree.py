import unittest

from dnctree.tree import Tree


class TestTree(unittest.TestCase):
    def test_building(self):
        t = Tree()
        t.add_edge('a', '1')
        t.add_edge('b', '1')
        t.add_edge('1', '2')
        t.add_edge('d', '3')
        t.add_edge('3', '4')
        t.add_edge('2', '3')
        t.add_edge('4', 'e')
        t.add_edge('f', '4')
        t.add_edge('c', '2')
        t.set_start_node('1')
        self.assertEqual(str(t), '(a,b,((d,(e,f)),c));')

        t.set_start_node('3')
        self.assertEqual(str(t), '(d,(e,f),((a,b),c));')

    def test_merge(self):
        t1 = Tree()
        t1.add_edge('a', '1')
        t1.add_edge('b', '1')
        t1.add_edge('1', '4')

        t2 = Tree()
        t2.add_edge('c', '2')
        t2.add_edge('d', '2')
        t2.add_edge('2', '4')

        t3 = Tree()
        t3.add_edge('e', '3')
        t3.add_edge('f', '3')
        t3.add_edge('3', '4')

        t1.merge(t2, '4')
        t1.merge(t3, '4')

        t1.set_start_node('4')
        self.assertEqual(str(t1), '((a,b),(c,d),(e,f));')
