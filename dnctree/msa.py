from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import fastphylo
import numpy as np

import dnctree.exceptions as dnc


class MSA:
    '''
    Storing multiple sequence alignments, backed by a BioPython alignment
    internally regardless of how it was read in.
    '''

    def __init__(self, biopy_alignment, seq_type):
        '''
        biopy_alignment: a BioPython MultipleSeqAlignment.
        seq_type:        'aa' or 'dna'.
        '''
        self._msa = biopy_alignment
        self.msa_width = biopy_alignment.get_alignment_length()
        self.type = seq_type # 'aa', or 'dna'
        self._taxa = [record.id for record in biopy_alignment]
        self._id2int = dict()

        # Create a lookup table, mapping id to an integer
        idx = 0
        for record in self._msa:
            self._id2int[record.id] = idx
            idx += 1

        # Create translation table for sequences so that it is easy to create the N matrices
        if seq_type == 'aa':
            alphabet='ARNDCQEGHILKMFPSTWYV'
        else:
            alphabet='ACGT'
        self._pair_table = dict()
        for i in range(len(alphabet)):
            for j in range(len(alphabet)):
                self._pair_table[alphabet[i], alphabet[j]] = (i, j)


    @classmethod
    def from_alv(cls, alv_msa):
        '''
        Build an MSA from an alv-wrapped alignment. Used for formats
        fastphylo's readers don't support (Nexus, Clustal).
        '''
        return cls(alv_msa.al, alv_msa.type)

    @classmethod
    def from_fastphylo(cls, fp_alignment):
        '''
        Build an MSA from a fastphylo Alignment (read via fastphylo's fast
        FASTA/Phylip/Stockholm reader).
        '''
        biopy_msa = MultipleSeqAlignment([
            SeqRecord(Seq(seq.data), id=seq.accession) for seq in fp_alignment
        ])
        seq_type = 'aa' if fp_alignment[0].seq_type == fastphylo.SeqType.Protein else 'dna'
        return cls(biopy_msa, seq_type)

    @classmethod
    def from_seq_list(cls, seqs):
        '''
        Instantiate an MSA from the list of sequences.
        This method is primarily for testing.
        '''
        biopy_msa = MultipleSeqAlignment([])
        for i, s in enumerate(seqs):
            acc = f'seq{i}'
            new_record = SeqRecord(Seq(s), id=acc)
            biopy_msa.append(new_record)
        return cls(biopy_msa, 'aa')

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

        n_chars = 0
        n_diffs = 0
        for c1, c2 in zip(s1, s2):
            present = self._pair_table.get((c1, c2), None)
            if present:
                row, col = present
                N[row,col] += 1
                n_chars += 1
                if row != col:
                    n_diffs += 1
        if n_chars == 0:
            raise dnc.NoSharedCharactersError()
        if n_diffs == n_chars:
            raise dnc.AllCharactersDifferentError()
        return N

    @staticmethod
    def can_retrieve_distances() -> bool:
        """Checks if distances are readily available and do not have to be
        computed manually.
        """
        return False

    def distance(self, t1, t2) -> float:
        """
        We cannot calculate distances using the MSA-class.
        Do not call this method!
        """
        raise Exception("Bug! Could not retrieve distance. (MSA.can_retrieve_distances() == False)")
