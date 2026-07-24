import unittest
from dnctree.neighbor_joining import dnc_neighborjoining
from dnctree.msa import MSA
from dnctree.partialdistancematrix import PartialDistanceMatrix
from dnctree.distances import kimura_distance

import ete3


class Test_dnctree(unittest.TestCase):
    def setUp(self):
        msa = MSA.from_seq_list(['PEPTIDELARS', 'PEPTIDELARS', 'PEPTIDELARA', 'SEPTICELARA', 'SEPTICELARS'])
        self.pdm = PartialDistanceMatrix(msa, kimura_distance)

    def test_dnc_nj(self):
        s = dnc_neighborjoining(self.pdm, ['seq0', 'seq1', 'seq2'])
        # Leaf order in the output isn't guaranteed, only that all three
        # taxa are present (there's only one possible unrooted topology
        # for 3 leaves anyway).
        for taxon in ('seq0', 'seq1', 'seq2'):
            self.assertIn(taxon, str(s))

        s_simple = dnc_neighborjoining(self.pdm, ['seq0', 'seq1', 'seq2', 'seq3', 'seq4'])
        s_tree = ete3.Tree(str(s_simple))

        s_weighted = dnc_neighborjoining(self.pdm, ['seq0', 'seq1', 'seq2', 'seq3', 'seq4'])
        s_weighted_tree = ete3.Tree(str(s_weighted))

        truth = ete3.Tree('((seq0, seq1), seq2, (seq3, seq4));')
        rf_simple, *x = truth.robinson_foulds(s_tree, unrooted_trees=True)
        rf_weighted, *x = truth.robinson_foulds(s_weighted_tree, unrooted_trees=True)
        self.assertEqual(rf_simple, 0)
        self.assertEqual(rf_weighted, 0)
        

        # s = dnc_neighborjoining(self.pdm, self.pdm.taxa)
        # print(s)
