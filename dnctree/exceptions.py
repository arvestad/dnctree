'''
Some exception classes specific for dnctree.
'''

class NoSharedCharactersError(Exception):
    '''
    Used to signal that two sequences actually have no homologous
    positions. There are plenty of examples of this in Pfam alignments.
    '''
    pass


class AllCharactersDifferentError(Exception):
    '''
    Signals that two aligned sequences have replacements/substitutions
    in all positions. There are plenty of such examples in Pfam.
    '''
    pass

class UnknownModel(Exception):
    '''
    Used when a user-requested model is not implemented or simply misspelled.
    '''
    pass
