import unittest

import ete3

from dnctree.dnctree_multi import dnc_tree_multi, _reachable_leaves, _gini_impurity, _assign_non_core_taxa_tree
from dnctree.msa import MSA
from dnctree.partialdistancematrix import PartialDistanceMatrix
from dnctree.distances import kimura_distance
from dnctree.tree import Tree


class TestDncTreeMulti(unittest.TestCase):
    def setUp(self):
        msa = MSA.from_seq_list([
            'PEPTIDELARS', 'PEPTIDELARS', 'PEPTIDELARA', 'SEPTICELARA', 'SEPTICELARS',
        ])
        self.pdm = PartialDistanceMatrix(msa, kimura_distance)
        self.taxa = ['seq0', 'seq1', 'seq2', 'seq3', 'seq4']

    def test_dnc_tree_multi_quartet_end_to_end(self):
        # Too few taxa here for exact-reconstruction expectations given
        # random core sampling, so just check every taxon is present and the
        # tree is parseable.
        t = dnc_tree_multi(self.pdm, self.taxa, base_case_size=3, core_size=3,
                            partition_method='quartet')
        s = str(t)
        for taxon in self.taxa:
            self.assertIn(taxon, s)
        result_tree = ete3.Tree(s)
        self.assertEqual(set(result_tree.get_leaf_names()), set(self.taxa))

    def test_dnc_tree_multi_tree_end_to_end(self):
        t = dnc_tree_multi(self.pdm, self.taxa, base_case_size=3, core_size=3,
                            partition_method='tree')
        s = str(t)
        for taxon in self.taxa:
            self.assertIn(taxon, s)
        result_tree = ete3.Tree(s)
        self.assertEqual(set(result_tree.get_leaf_names()), set(self.taxa))

    def test_dnc_tree_multi_base_case_nj_vs_bionj(self):
        # len(taxa) <= base_case_size, so this only exercises _run_base_case.
        for base_case_method in ('nj', 'bionj'):
            t = dnc_tree_multi(self.pdm, self.taxa, base_case_size=5, core_size=5,
                                base_case_method=base_case_method)
            s = str(t)
            for taxon in self.taxa:
                self.assertIn(taxon, s)
            result_tree = ete3.Tree(s)
            self.assertEqual(set(result_tree.get_leaf_names()), set(self.taxa))


class TestReachableLeaves(unittest.TestCase):
    def test_partitions_by_directed_edge(self):
        # p - a - (a0, a1)
        #  \
        #   b - (b0, b1)
        t = Tree()
        t.add_edge('p', 'a')
        t.add_edge('a', 'a0')
        t.add_edge('a', 'a1')
        t.add_edge('p', 'b')
        t.add_edge('b', 'b0')
        t.add_edge('b', 'b1')

        self.assertEqual(_reachable_leaves(t, ('p', 'a')), {'a0', 'a1'})
        self.assertEqual(_reachable_leaves(t, ('p', 'b')), {'b0', 'b1'})

    def test_single_leaf(self):
        t = Tree()
        t.add_edge('p', 'a')
        self.assertEqual(_reachable_leaves(t, ('p', 'a')), {'a'})


class TestGiniImpurity(unittest.TestCase):
    def test_pure_group_is_zero(self):
        self.assertEqual(_gini_impurity([4, 0, 0]), 0.0)

    def test_even_split_is_high(self):
        self.assertAlmostEqual(_gini_impurity([2, 2, 2]), 2 / 3)

    def test_more_mixed_is_less_pure(self):
        self.assertLess(_gini_impurity([4, 0, 0]), _gini_impurity([3, 1, 0]))
        self.assertLess(_gini_impurity([3, 1, 0]), _gini_impurity([2, 1, 1]))

    def test_empty_group_is_zero(self):
        self.assertEqual(_gini_impurity([0, 0, 0]), 0.0)


class TestAssignNonCoreTaxaTree(unittest.TestCase):
    def test_picks_the_purer_sides_dominant_clade(self):
        # taxon_to_clade already known for core taxa (a0/a1 -> 0, b0/b1 -> 1,
        # c0/c1 -> 2). v attaches so that removing v and its neighbor p
        # yields one pure side {a0, a1} (clade 0) and one mixed side
        # {b0, b1, c0, c1} (clades 1 and 2) - v should be assigned clade 0.
        class FakeDm:
            def get(self, a, b):
                return 1.0  # unused by this function directly

        import dnctree.dnctree_multi as multi

        original_nj = multi.dnc_neighborjoining
        t_w = Tree()
        t_w.add_edge('p', 'v')
        t_w.add_edge('p', 'a')
        t_w.add_edge('a', 'a0')
        t_w.add_edge('a', 'a1')
        t_w.add_edge('p', 'q')
        t_w.add_edge('q', 'b0')
        t_w.add_edge('q', 'b1')
        t_w.add_edge('q', 'r')
        t_w.add_edge('r', 'c0')
        t_w.add_edge('r', 'c1')
        multi.dnc_neighborjoining = lambda dm, taxa: t_w
        try:
            taxon_to_clade = {'a0': 0, 'a1': 0, 'b0': 1, 'b1': 1, 'c0': 2, 'c1': 2}
            result = _assign_non_core_taxa_tree(
                FakeDm(), ['a0', 'a1', 'b0', 'b1', 'c0', 'c1'], ['v'], taxon_to_clade
            )
        finally:
            multi.dnc_neighborjoining = original_nj

        self.assertEqual(result, {0: ['v'], 1: [], 2: []})
        self.assertEqual(taxon_to_clade['v'], 0)


if __name__ == '__main__':
    unittest.main()
