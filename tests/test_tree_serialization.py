import unittest

from dnctree.tree import Tree


class TestTreeSerialization(unittest.TestCase):
    def test_single_vertex(self):
        t = Tree()
        t.add_vertex('a')
        self.assertEqual(str(t), 'a;')

    def test_two_leaves_single_edge(self):
        # Regression test: with no internal vertex to root the traversal at,
        # whichever leaf ends up as the default start node used to have its
        # own label silently dropped from the Newick output.
        t = Tree()
        t.add_edge('a', 'b')
        s = str(t)
        self.assertIn('a', s)
        self.assertIn('b', s)

    def test_star_rooted_at_internal_vertex(self):
        t = Tree()
        t.add_edge('a', 'c')
        t.add_edge('b', 'c')
        t.add_edge('d', 'c')
        t.set_start_node('c')
        self.assertEqual(str(t), '(a,b,d);')

    def test_larger_tree_rooted_at_a_leaf(self):
        # Regression test: same underlying bug as test_two_leaves_single_edge,
        # but for a bigger tree where the default (or explicitly set) start
        # node happens to be a leaf rather than an internal vertex.
        t = Tree()
        t.add_edge('a', 'x')
        t.add_edge('b', 'x')
        t.add_edge('x', 'y')
        t.add_edge('c', 'y')
        t.add_edge('d', 'y')
        t.set_start_node('a')  # 'a' is a leaf
        s = str(t)
        for taxon in ('a', 'b', 'c', 'd'):
            self.assertIn(taxon, s)


if __name__ == '__main__':
    unittest.main()
