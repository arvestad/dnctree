#from Bio.Alphabet import ProteinAlphabet, AlphabetEncoder, Gapped
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import alv.alignment
import itertools as it
import numpy as np
import sys

import dnctree.exceptions as dnc
from dnctree.pahmm_adaptor import pahmm_available
if pahmm_available():
    from pahmm import Sequences as PahmmSequences, BandingEstimator


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
        self.msa_width = alv_msa.al.get_alignment_length()
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
        biopy_msa = MultipleSeqAlignment([])
        for i, seq in enumerate(seqs):
            acc = f'seq{i}'
            biopy_msa.add_sequence(acc, seq)
        return cls(alv.alignment.AminoAcidAlignment(biopy_msa))         # Instantiate

    @classmethod
    def from_biopython(cls, msa):
        pass

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

if not pahmm_available():
    class MSApaHMM:             # Dummy class
        pass
else:
    class MSApaHMM:
        """
        Storing unaligned sequence data using the PAHMM module.
        The clas name is due to the related MSA class!

        :param sequence_type can be 'aa' (Amino-acids) or 'dna' (Nucleotides)
        """
        def __init__(self, sequences: PahmmSequences, sequence_type='aa'):
            self.type = sequence_type
            self._sequences = sequences
            self._taxa = list(map(self._sequences.get_seq_name, range(len(self._sequences))))


        @classmethod
        def from_file(cls, filename, model, verbosity=[]):
            '''
            Instantiate a holder for sequence data in the PAHMM module from the given filename.
            '''
            be = BandingEstimator()
            if verbosity:
                print('Using paHMM for distance estimation. Note: accepts only unaligned sequences.', file=sys.stderr)
                print(f"Reading '{filename}'...", file=sys.stderr)
            be.set_file_input(filename)

            if verbosity:
                print(f"Using {model} as model. Estimating parameters and rough pairwise distances.", file=sys.stderr)

            sequence_type = 'unknown'
            if model in ['HKY85', 'GTR']:
                sequence_type == 'dna'
            elif model in ['WAG', 'LG', 'JTT']:
                sequence_type == 'aa'
            else:
                raise dnc.UnknownModel()

            #        seq_data = cls(be.execute_wag_model())
            seq_data = cls(be.apply_model(model), sequence_type=sequence_type)
            if verbosity:
                print("paHMM initialization done.", file=sys.stderr)
            return seq_data

        def taxa(self):
            return self._taxa

        def sequences(self):
            """
            Return the sequences data-structure.
            """
            return self._sequences

        def __getitem__(self, key):
            """
            Return a fake SeqRecord. We don't want to construct a SeqRecord from a real sequence
            because that would involve converting the sequence to a Python-string. That's an O(n) operation
            and that's not acceptable. Distance calculations are left entirely to paHMM anyway.
            """
            return SeqRecord(Seq(""))

        @staticmethod
        def can_retrieve_distances() -> bool:
            """
            Checks if distances are readily available and do not have to be
            computed 'manually'.
            """
            return True
        
        def distance(self, t1: bytes, t2: bytes) -> float:
            """
            Retrieve the distance between two sequences.
            """
            d = self.sequences().get_distance_from_names(t1, t2)
            return d
