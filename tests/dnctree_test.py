import unittest
from dnctree import dnc_neighborjoining, nj_selection_function
from dnctree.msa import MSA
from dnctree.partialdistancematrix import PartialDistanceMatrix
from dnctree.distances import kimura_distance


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

        # s = dnc_neighborjoining(self.pdm, ['seq0', 'seq1', 'seq2', 'seq3'])
        # print(s)

        # s = dnc_neighborjoining(self.pdm, self.pdm.taxa)
        # print(s)
