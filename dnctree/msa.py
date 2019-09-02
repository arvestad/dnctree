
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
        self._taxa = list(alv_msa.accessions())
        self._id2int = dict()

        # Create a lookup table, mapping id to an integer
        idx = 0
        for record in self._msa:
            self._id2int[record.id] = idx
            idx += 1

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
