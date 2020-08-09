import unittest
from dnctree.msa import MSA
from dnctree.partialdistancematrix import PartialDistanceMatrix
from dnctree import choose_distance_function, divide_n_conquer_tree
from dnctree.distances import kimura_distance
import dnctree.exceptions as dncexc

class Test_partialdistancematrix(unittest.TestCase):
    def setUp(self):
        msa = MSA.from_seq_list(['PEPTIDELARS', 'PEPTIDELARS', 'PEPTIDELARA'])
        self.pdm = PartialDistanceMatrix(msa, kimura_distance)

    def test_unique_ids(self):
        acc1 = self.pdm.create_unique_vertex_id()
        acc2 = self.pdm.create_unique_vertex_id()
        self.assertNotEqual(acc1, acc2)

    def test_kimura_distances(self):
        d1 = self.pdm.get('seq0', 'seq1')
        d2 = self.pdm.get('seq0', 'seq2')
        self.assertTrue(abs(d1) < 0.00001) # Close to zero at least
        self.assertTrue(d1 < d2)

    def test_taxa(self):
        self.assertEqual(len(self.pdm.taxa), 3)

    def test_set(self):
        with self.assertRaises(KeyError):
            self.pdm.set('WhateverAccession1', 'WhateverAccession2', 0)
        self.pdm.set('seq0', 'seq1', 1)
        self.assertEqual(1, self.pdm.get('seq0', 'seq1'))

    # def test_bad_data(self):
    #     msa = MSA.from_seq_list(['LARS----', '----SEAL', 'CEREALLE', 'SEALEAGE'])
    #     with self.assertRaises(dncexc.NoSharedCharactersError):
    #         divide_n_conquer_tree(msa, verbose=[])
