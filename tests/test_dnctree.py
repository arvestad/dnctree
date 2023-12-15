import unittest
from dnctree import dnc_neighborjoining, nj_selection_function
from dnctree.msa import MSA
from dnctree.partialdistancematrix import PartialDistanceMatrix
from dnctree.distances import kimura_distance

import ete3


class Test_dnctree(unittest.TestCase):
    def setUp(self):
        msa = MSA.from_seq_list(['PEPTIDELARS', 'PEPTIDELARS', 'PEPTIDELARA', 'SEPTICELARA', 'SEPTICELARS'])
        self.pdm = PartialDistanceMatrix(msa, kimura_distance)

    def test_nj_selection(self):
        pair = nj_selection_function(self.pdm, ['seq0', 'seq1', 'seq2', 'seq3', 'seq4'])
        self.assertTrue('seq3' in pair)
        self.assertTrue('seq4' in pair)

        pair = nj_selection_function(self.pdm, ['seq0', 'seq1', 'seq2'])
        self.assertTrue('seq0' in pair)
        self.assertTrue('seq1' in pair)


    def test_dnc_nj(self):
        s = dnc_neighborjoining(self.pdm, ['seq0', 'seq1', 'seq2'])
        self.assertEqual(str(s), '(seq0,seq1,seq2);') # I cannot guarantee the order though...

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
