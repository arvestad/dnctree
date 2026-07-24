import unittest

import alv.alignment
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import dnctree.exceptions as dncexc
from dnctree.api import choose_distance_function, divide_n_conquer_tree
from dnctree.msa import MSA
from dnctree.partialdistancematrix import PartialDistanceMatrix


def dna_msa(seqs):
    biopy_msa = MultipleSeqAlignment([])
    for i, s in enumerate(seqs):
        biopy_msa.append(SeqRecord(Seq(s), id=f'seq{i}'))
    return MSA(alv.alignment.DnaAlignment(biopy_msa))


class TestDnaDistance(unittest.TestCase):
    def setUp(self):
        self.msa = dna_msa([
            'ACGTACGTACGTACGTACGT',
            'ACGTACGAACGTTCGTACGT',
            'TTTTACGTACGTACGTACGA',
        ])

    def test_choose_distance_function_default_model(self):
        fn, model = choose_distance_function('dna')
        self.assertEqual(model, 'k2p')
        self.assertTrue(getattr(fn, 'operates_on_sequences', False))

    def test_choose_distance_function_unknown_model(self):
        with self.assertRaises(dncexc.UnknownModel):
            choose_distance_function('dna', 'bogus-model')

    def test_distance_is_computed_for_each_model(self):
        for model in ('hamming', 'jc', 'k2p', 'tn93'):
            fn, _ = choose_distance_function('dna', model)
            pdm = PartialDistanceMatrix(self.msa, fn)
            d = pdm.get('seq0', 'seq1')
            self.assertTrue(d >= 0)

    def test_rna_uses_same_models(self):
        fn, model = choose_distance_function('rna', 'k2p')
        self.assertEqual(model, 'k2p')

    def test_no_informative_sites_falls_back_to_maximal_distance(self):
        degenerate = dna_msa(['AAAAAAAAAA', 'GGGGGGGGGG'])
        fn, _ = choose_distance_function('dna', 'k2p')
        pdm = PartialDistanceMatrix(degenerate, fn)
        d = pdm.get('seq0', 'seq1')
        self.assertEqual(d, pdm._maximal_distance)

    def test_divide_n_conquer_tree_on_dna(self):
        t, aux_info = divide_n_conquer_tree(self.msa, model_name='k2p', base_case_size=10)
        self.assertEqual(aux_info['model-name'], 'k2p')
        self.assertEqual(aux_info['description'], 'DNA alignment')
        newick = str(t)
        for i in range(3):
            self.assertIn(f'seq{i}', newick)


if __name__ == '__main__':
    unittest.main()
