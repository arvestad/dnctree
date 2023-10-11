import itertools
import random
import sys
import textwrap

from modelmatcher import RateMatrix
from dnctree.msa import MSA, MSApaHMM
from dnctree.tree import Tree
from dnctree.partialdistancematrix import PartialDistanceMatrix
from dnctree.distances import ml_distance_estimate


dna_models = ['HKY', 'JC']


def divide_n_conquer_tree(seq_data, model_name=None, max_n_attempts=100, max_clade_size=0.5,
                          base_case_size=100, first_triple=None, verbose=[]):
    '''
    Input: An object with sequence data (for the input alignment or unaligned if using PaHMM), 
           and a function of three args (two taxa names and a NumPy 20x20 matrix with amino acid counts
           for estimating evolutionary distance.
    Returns: a tree as a Newick-formatted string.

    Parameters:
    max_n_attempts: The number of times we try to find a suitable taxa triple to recurse on.
    max_clade_size: The fraction of the taxa that are required before deciding we are happy
                    with a triple. We require that the largest subproblem (clade) has at most this fraction of taxa.
    base_case_size: For clades with at most this many taxa, we just run NJ.
    verbose:        Whether to output a lot of unnecessary facts.
    '''
    assert base_case_size >= 3 # This is truly a minimum
    assert max_clade_size >= 0.5  # Pointless to try smaller values
    assert max_n_attempts >= 1    # Truly a minimum

    if isinstance(seq_data, MSA):
        distance_fcn, chosen_model_name = choose_distance_function(seq_data.type, model_name)
        data_description = f'{seq_data.type.upper()} alignment'
    elif isinstance(seq_data, MSApaHMM):
        distance_fcn = 'not_needed_strangely_enough' # Clearly a case of bad program structure. Should refactor.
        chosen_model_name = model_name
        data_description = f'Unaligned {seq_data.type.upper()} sequences'


    dm = PartialDistanceMatrix(seq_data, distance_fcn, verbose=verbose)                 # We will gradually add distances to this "matrix"

    t = dnc_tree(dm, seq_data.taxa(), max_n_attempts, max_clade_size, base_case_size, first_triple, verbose)
    actual_work, fraction_work, n = dm.estimate_computational_savings()
    comment = dm._computational_savings()
    if 'info' in verbose or 'verbose' in verbose:
        print(f'[{comment}]')
    aux_info = {'distances-computed': actual_work,
                'fraction-computed-distances': round(fraction_work,3),
                'n-taxa': n,
                'comment': comment,
                'model-name': chosen_model_name,
                'description': data_description
    }
    if isinstance(seq_data, MSA):
        aux_info['msa-width'] = seq_data.msa_width

    return t, aux_info


def testing_divide_n_conquer(model, max_n_attempts=100, max_clade_size=0.5, base_case_size=100, first_triple=None, verbose=[]):
    '''
    The purpose with this function is to be able to test the algorithm and
    see how sensitive it is to errors in the distance estimation.

    Input and parameters:

    model_tree:   A tree defining the distances.
    epsilon:      Decides the error added to distances.

    Other parameters are the same as for "divide_n_conquer".
    '''

    t = dnc_tree(model, model.taxa, max_n_attempts, max_clade_size, base_case_size, first_triple, verbose)
    return t


def choose_distance_function(seqtype, model_name=None):
    '''
    Pick a model depending on sequence type and model_name.
    If the two args are not in harmony (e.g., 'dna' and 'WAG'), throw an error.

    seqtype:    'aa' or 'dna'
    model_name: 'WAG' or any other model implemented in the module "modelmatcher".
                'kimura' is also an acceptable string.
    '''
    if seqtype == 'aa':
        chosen_model = 'WAG'    # Standard choice
        if model_name:           # Replace the standard choice
            if model_name in dna_models:
                raise Exception('You have requested a DNA model for amino acid sequences.')
            else:
                chosen_model = model_name
        model = RateMatrix.instantiate(chosen_model)
        distance_fcn = lambda N: ml_distance_estimate(model, N) # XXX Correct call?
        return distance_fcn, chosen_model
    elif seqtype == 'dna':
        raise dnc.UnknownModel('DNA models not yet implemented, sorry.')
    else:
        raise dnc.UnknownModel(f'Not a valid sequence type: {seqtype}.')


def dnc_tree(dm, taxa, max_n_attempts, max_clade_size, base_case_size, first_triple, verbose=[]):
    '''
    The core algorithm

    dm:       A PartialDistanceMatrix instance, to draw distances from. Note that
              distances are computed on-demand and "cached" with function registered
              in this object.
    taxa:     A list of identifiers (part of dm) to infer the tree for
    verbose:  boolean, whether to print extra info to stderr

    returns:  A Tree object.
    '''
    if len(taxa) <= base_case_size:
        if 'verbose' in verbose:
            print(f'Full NJ on {len(taxa)} taxa:', file=sys.stderr, flush=True)
            print(_taxa_string_helper(taxa, 2), file=sys.stderr, flush=True)
        t = dnc_neighborjoining(dm, taxa, verbose)
        return t
    else:
        if first_triple is not None:
            t = first_triple
            c = clade_sort_taxa(dm, taxa, t[0], t[1], t[2])
        else:
            t, c = sample_three_taxa(dm, taxa, max_n_attempts, max_clade_size, verbose)
        # t[i] are 3 taxa ids
        # c[i] are three lists of taxa representing what we think should be become a clade (c[i] together with t[i]).

        # Add an inner vertex v that connects c1, c2, and c3.
        v = dm.create_central_vertex(t, c)
        if verbose:
            taxa_string = _taxa_string_helper([t[0], t[1],t[2]], 0)
            print(f'Recursing on subproblems induced by {taxa_string}', file=sys.stderr)

        # Recurse on cx+tx+v, for x in {1,2,3}.
        subtrees = [None, None, None]
        for i in range(3):
            c[i].append(v)      # Ensure that the central vertex participates as a representative of the other clades
            subtrees[i] = dnc_tree(dm, c[i], max_n_attempts, max_clade_size, base_case_size, None, verbose)
            if verbose:
                dm.print_progress()
#            print(f't{i}:', subtrees[i], file=sys.stderr)
        t_connected = connect_the_trees(subtrees, v)
        return t_connected


def sample_three_taxa(dm, taxa, max_n_attempts, max_clade_size, verbose=[]):
    '''
    Sample three taxa and see if they split the set of taxa not too unevenly. Repeat until the largest
    clade (out of three) is small enough or we have used our max number of attempts.

    Parameters:
      dm              A PartialDistanceMatrix
      taxa            A list of taxa ids
      max_n_attempts  How many samples we will attempt at most
      max_clade_size  Try not to accept a taxa sample that defines one of the clades as larger than this,
                      as a fraction of the taxa. (It is a float.)
    '''
    N = len(taxa) * max_clade_size # The max clade size is given as a fraction. N is the actual number of taxa to allow for a subproblem.
    best_so_far = None
    smallest_large_clade_so_far = len(taxa) # Initialize with worst possible subproblem size
    for i in range(max_n_attempts):
        v1, v2, v3 = random.sample(taxa, 3) # Three taxa strings
        if 'verbose' in verbose:
            taxa_string = _taxa_string_helper([v1, v2, v3], 0)
            print(f'Attempt: {taxa_string}.', end='  ', file=sys.stderr, flush=True)
        c1, c2, c3 = clade_sort_taxa(dm, taxa, v1, v2, v3) # Three groups of taxa: "clade wannabees"
        largest = max(len(c1), len(c2), len(c3)) # Out of the three found clades, which one is largest?
        if verbose:
            print(f'Clade sizes: {len(c1)}, {len(c2)}, and {len(c3)}', file=sys.stderr, flush=True)

        if largest <= N:
            # Sample is good, so let's use it
            return (v1, v2, v3), (c1, c2, c3)
        elif largest < smallest_large_clade_so_far:
            # Not an ideal sample, but better than what we have seen before, so save it.
            smallest_large_clade_so_far = largest
            best_so_far = (v1, v2, v3), (c1, c2, c3)
    assert best_so_far != None
    return best_so_far


def clade_sort_taxa(dm, taxa, l1, l2, l3):
    '''
    Go through chosen taxa in dm (not all, because might be a subproblem) and
    decide which clade they belong to: l1, l2, or l3 (as represented by those
    leaves/taxa).  Return three lists of taxa identifiers corresponding to the
    clades (l1, l2, or l3) the belong to.
    '''
    clade_assignment = {l1: [l1], l2: [l2], l3: [l3]} # Keep track of clade assignment.
    for x in taxa:
        if x not in [l1, l2, l3]:
            clade = dm.select_clade(x, l1, l2, l3)
            clade_assignment[clade].append(x)

    # for cl in [l1, l2, l3]:
    #     print(f'Clade {cl}: {clade_assignment[cl]}', file=sys.stderr)
    return clade_assignment[l1], clade_assignment[l2], clade_assignment[l3]



def connect_the_trees(subtrees, connecting_node):
    '''
    Merge the subtrees. They are sharing a unique vertex, so it is guaranteed to connect the subtrees.
    '''
    first_tree = subtrees[0]
    for t in subtrees[1:]:
        first_tree.merge(t, connecting_node)
    first_tree.set_start_node(connecting_node)
    return first_tree


def nj_selection_function(dm, current_leaves):

    def sum_distances(a):
        sum_for_a = 0
        for b in current_leaves:
            if a != b:
                sum_for_a += dm.get(a, b)
        return sum_for_a

    n_minus_2 = len(current_leaves) - 2
    Q_ij_min = n_minus_2 * dm.get(current_leaves[0], current_leaves[1]) # Arbitrary but safe starting value that won't be chosen.
    best_so_far = current_leaves[0], current_leaves[1]
    for x, y in itertools.combinations(current_leaves, 2):
        Q_ij = n_minus_2 * dm.get(x, y) - sum_distances(x) - sum_distances(y)
        if Q_ij < Q_ij_min:
            Q_ij_min = Q_ij
            best_so_far = x, y

    return best_so_far


def dnc_neighborjoining(dm, taxa, verbose=[]):
    '''
    An implementation of NJ meant to compute base cases in dnctree.

    dm:       a PartialDistanceMatrix instance
    taxa:     a list of identifiers (part of dm) to infer the tree for
    verbose:  boolean, whether to print extra info to stderr
    '''
    t = Tree()
    if len(taxa) == 1:
        t.add_vertex(taxa[0])
        return t
    elif len(taxa) == 2:
        t.add_edge(taxa[0], taxa[1])
        return t
    else:
        current_leaves = taxa[:] # Slicing out a copy
        while len(current_leaves) > 3:
            x, y = nj_selection_function(dm, current_leaves)
            # if verbose:
            #     print(f'Joining {x} and {y}', file=sys.stderr)
            current_leaves.remove(x)
            current_leaves.remove(y)
            new_v = dm.create_representative_vertex(x, y, current_leaves)
            current_leaves.append(new_v)
            t.add_edge(x, new_v)
            t.add_edge(y, new_v)

        c = dm.create_unique_vertex_id()
        for leaf in current_leaves:
            t.add_edge(leaf, c)
        t.set_start_node(c)
        return t


def _taxa_string_helper(l, indent):
    s = textwrap.fill(', '.join(map(lambda t: t.decode() if isinstance(t, bytes) else t, l)), width=100 - indent)
    return textwrap.indent(s, ' ' * indent)

