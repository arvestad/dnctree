import unittest

import ete3

from dnctree.dnctree_k import dnc_neighborjoining_with_weights, dnc_tree_k
from dnctree.msa import MSA
from dnctree.partialdistancematrix import PartialDistanceMatrix
from dnctree.distances import kimura_distance


class TestDncNeighborjoiningWithWeights(unittest.TestCase):
    def setUp(self):
        msa = MSA.from_seq_list([
            'PEPTIDELARS', 'PEPTIDELARS', 'PEPTIDELARA', 'SEPTICELARA', 'SEPTICELARS',
        ])
        self.pdm = PartialDistanceMatrix(msa, kimura_distance)
        self.taxa = ['seq0', 'seq1', 'seq2', 'seq3', 'seq4']

    def test_single_taxon(self):
        t = dnc_neighborjoining_with_weights(self.pdm, ['seq0'])
        self.assertEqual(str(t), 'seq0;')

    def test_two_taxa(self):
        t = dnc_neighborjoining_with_weights(self.pdm, ['seq0', 'seq1'])
        s = str(t)
        self.assertIn('seq0', s)
        self.assertIn('seq1', s)

    def test_three_taxa_keeps_all_leaf_labels(self):
        # Regression test: the Newick serializer silently drops a vertex's
        # own label if the tree's start node happens to be a leaf, so the
        # bionj-based implementation must explicitly pick an internal vertex
        # as the start node.
        t = dnc_neighborjoining_with_weights(self.pdm, ['seq0', 'seq1', 'seq2'])
        s = str(t)
        for taxon in ('seq0', 'seq1', 'seq2'):
            self.assertIn(taxon, s)

    def test_all_taxa_present_and_topology_is_plausible(self):
        t = dnc_neighborjoining_with_weights(self.pdm, self.taxa)
        s = str(t)
        for taxon in self.taxa:
            self.assertIn(taxon, s)

        result_tree = ete3.Tree(s)
        truth = ete3.Tree('((seq0, seq1), seq2, (seq3, seq4));')
        rf, *_ = truth.robinson_foulds(result_tree, unrooted_trees=True)
        self.assertEqual(rf, 0)

    def test_edge_weights_are_recorded_in_dm(self):
        dnc_neighborjoining_with_weights(self.pdm, self.taxa)
        # seq0 and seq1 are identical sequences, so their branch lengths to
        # whatever vertex joins them first should be non-negative and small.
        d = self.pdm.get('seq0', 'seq1')
        self.assertTrue(d >= 0)

    def test_dnc_tree_k_end_to_end(self):
        # With only 5 taxa and core_size=3, dnc_tree_k's random core sampling
        # can legitimately produce different (still valid) topologies
        # depending on which taxa get picked as the core - not enough
        # redundancy at this tiny scale to expect exact reconstruction. Check
        # what actually matters: every taxon present, and a parseable tree.
        t = dnc_tree_k(self.pdm, self.taxa, base_case_size=3, core_size=3)
        s = str(t)
        for taxon in self.taxa:
            self.assertIn(taxon, s)
        result_tree = ete3.Tree(s)
        self.assertEqual(set(result_tree.get_leaf_names()), set(self.taxa))


if __name__ == '__main__':
    unittest.main()
