from Bio.Alphabet import ProteinAlphabet, AlphabetEncoder, Gapped
from Bio.Align import MultipleSeqAlignment
import alv.alignment
import itertools as it
import numpy as np

class MSA:
    '''
    Storing multiple sequence alignments.
    Using the "alv" MSA class because we get easy reading.
    '''

    def __init__(self, alv_msa):
        '''
        alv_msa is the MSA retrieved from AlignIO.
        '''
        self._msa = alv_msa.al  # Actual BioPython alignment object
        self.type = alv_msa.type # 'aa', or 'dna'
        self._taxa = list(alv_msa.accessions())
        self._id2int = dict()

        # Create a lookup table, mapping id to an integer
        idx = 0
        for record in self._msa:
            self._id2int[record.id] = idx
            idx += 1

        # Create translation table for sequences so that it is easy to create the N matrices
        if alv_msa.type == 'aa':
            alphabet='ARNDCQEGHILKMFPSTWYV'
        else:
            alphabet='ACGT'
        self._pair_table = dict()
        for i in range(len(alphabet)):
            for j in range(len(alphabet)):
                self._pair_table[alphabet[i], alphabet[j]] = (i, j)


    @classmethod
    def from_seq_list(cls, seqs):
        '''
        Instantiate an MSA from the list of sequences.
        This method is primarily for testing.
        '''
        biopy_msa = MultipleSeqAlignment([], Gapped(ProteinAlphabet(), "-"))
        for i, seq in enumerate(seqs):
            acc = f'seq{i}'
            biopy_msa.add_sequence(acc, seq)
        return cls(alv.alignment.aaAlignment(biopy_msa))         # Instantiate


    def taxa(self):
        return self._taxa

    def __getitem__(self, key):
        '''
        Return the SeqRecord in the MSA corresponding to id=key.
        '''
        try:
            idx = self._id2int[key]
            return self._msa[idx]
        except KeyError:
            raise KeyError(f'"{key}" is not a known accession/id in the given MSA')


    def count_pairs(self, acc1, acc2):
        '''
        Count the amino acid pairs in two aligned sequences.
        Return NumPy array. Notice that it is a 20x20 matrix for a protein alignment,
        but 4x4 for DNA.
        '''
        s1 = self[acc1]
        s2 = self[acc2]
        return self._count_pairs(s1, s2)


    def _count_pairs(self, s1, s2):
        '''
        Helper for count_pairs(), but with actual sequences as parameters.
        Main reason is to enable easier testing.
        '''
        if self.type == 'aa':
            N = np.zeros((20,20))
        elif self.type == 'dna':
            N = np.zeros((4,4))
        else:
            raise Exception('dnctree currently only handles protein or DNA sequences.')

        for c1, c2 in zip(s1, s2):
            present = self._pair_table.get((c1, c2), None)
            if present:
                row, col = present
                N[row,col] += 1
        return N
