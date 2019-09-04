import itertools
import random
import sys
import textwrap

from dnctree.tree import Tree

#from dnctree.dnctree import hubba
from dnctree.partialdistancematrix import PartialDistanceMatrix

def divide_n_conquer_tree(msa, max_n_attempts=100, max_clade_size=0.5, base_case_size=100, first_triple=None, verbose=False):
    '''
    Input: a Bio.Align object
    Returns: a tree

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

    dm = PartialDistanceMatrix(msa, verbose=verbose)                 # We will gradually add distances to this "matrix"

    t = dnc_tree(dm, msa.taxa(), max_n_attempts, max_clade_size, base_case_size, first_triple, verbose)
    if verbose:
        dm._print_computational_savings()
    return t


def dnc_tree(dm, taxa, max_n_attempts, max_clade_size, base_case_size, first_triple, verbose=False):
    if len(taxa) <= base_case_size:
        if verbose:
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
            print(f'Recursing on subproblems induced by {t[0]}, {t[1]}, and {t[2]}', file=sys.stderr)
#            print(f'   Central vertex: {v}', file=sys.stderr)

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


def sample_three_taxa(dm, taxa, max_n_attempts, max_clade_size, verbose):
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
    N = len(taxa) * max_clade_size
    best_so_far = None
    smallest_large_clade_so_far = len(taxa) # Initialize with worst possible subproblem size
    for i in range(max_n_attempts):
        t1, t2, t3 = random.sample(taxa, 3) # Three taxa strings
        if verbose:
            print(f'Attempt: {t1}, {t2}, {t3}.', end='  ', file=sys.stderr, flush=True)
        c1, c2, c3 = clade_sort_taxa(dm, taxa, t1, t2, t3) # Three groups of taxa: "clade wannabees"
        largest = max(len(c1), len(c2), len(c3)) # Out of the three found clades, which one is largest?
        if verbose:
            print(f'Clade sizes: {len(c1)}, {len(c2)}, and {len(c3)}', file=sys.stderr, flush=True)

        if largest <= N:
            # Sample is good, so let's use it
            return (t1, t2, t3), (c1, c2, c3)
        elif largest < smallest_large_clade_so_far:
            # Not an ideal sample, but better than what we have seen before, so save it.
            smallest_large_clade_so_far = largest
            best_so_far = (t1, t2, t3), (c1, c2, c3)
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
            clade = quartet_test(dm, l1, l2, l3, x)
            clade_assignment[clade].append(x)
    return clade_assignment[l1], clade_assignment[l2], clade_assignment[l3]


def quartet_test(dm, l1, l2, l3, x):
    '''
    Return the leaf (l1, l2, or l3) of which leaf a quartet test suggests x belongs to.
    '''
    def q_test(a, b, c, d):
        '''
        Check whether d(a,b) + d(c, d) < d(a, c) + d(b, d)
        '''
        # print(f'q_test: {a} with {b}, {c}, {d}', file=sys.stderr)
        # print(f'   d({a}, {b}) + d({c},{d}) = {dm.get(a, b) + dm.get(c, d)} \t d({a}, {c}) + d({b},{d}) = {dm.get(a, c) + dm.get(b,        d)}', file=sys.stderr)
        return dm.get(a, b) + dm.get(c, d) <= dm.get(a, c) + dm.get(b, d)

    def q_diff(a, b, c, d):
        '''
        Check whether d(a,b) + d(c, d) < d(a, c) + d(b, d)
        '''
        return dm.get(a, b) + dm.get(c, d) - dm.get(a, c) + dm.get(b, d)


    options = [(l1, q_diff(x, l1, l2, l3)),
               (l2, q_diff(x, l2, l1, l3)),
               (l3, q_diff(x, l3, l1, l2)),]
    choice = min(options, key = lambda pair: pair[1])
#    print(f'Choice: {x} belongs to {choice}. \t{options}', file=sys.stderr)
    return choice[0]

    # if q_test(x, l1, l2, l3):
    #     return l1
    # elif q_test(x, l2, l1, l3):
    #     return l2
    # elif q_test(x, l3, l1, l2):
    #     return l3
    # else:
    #     # We end up here if the distance matrix is not additive.
    #     # Let's just make an educated guess instead.
    #     options = [(l1, q_diff(l1, x, l2, l3)),
    #                (l2, q_diff(l2, x, l1, l3)),
    #                (l3, q_diff(l3, x, l1, l2)),]
    #     choice = min(options, key = lambda pair: pair[1])
    #     return choice[0]


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


def dnc_neighborjoining(dm, taxa, verbose):
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
        current_leaves = taxa
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
        else:
            c = dm.create_unique_vertex_id()
            for leaf in current_leaves:
                t.add_edge(leaf, c)
            t.set_start_node(c)
            return t

def _taxa_string_helper(l, indent):
    s = textwrap.fill(', '.join(l), width=100 - indent)
    return textwrap.indent(s, ' ' * indent)

def test_nj():
    dm = PartialDistanceMatrix()

