import logging
import random

from dnctree.neighbor_joining import dnc_neighborjoining, _taxa_string_helper


def dnc_tree_simple(
    dm,
    taxa,
    max_n_attempts,
    max_clade_size,
    base_case_size,
    first_triple=None
):
    """
    The core algorithm

    dm:       A PartialDistanceMatrix instance, to draw distances from. Note that
              distances are computed on-demand and "cached" with function registered
              in this object.
    taxa:     A list of identifiers (part of dm) to infer the tree for


    returns:  A Tree object.
    """
    if len(taxa) <= base_case_size:
        logging.debug(f"Full NJ on {len(taxa)} taxa:")
        logging.debug(_taxa_string_helper(taxa, 2))
        t = dnc_neighborjoining(dm, taxa)
        return t
    else:
        if first_triple is not None:
            t = first_triple
            c = clade_sort_taxa(dm, taxa, t[0], t[1], t[2])
        else:
            t, c = sample_three_taxa(dm, taxa, max_n_attempts, max_clade_size)
        # t[i] are 3 taxa ids
        # c[i] are three lists of taxa representing what we think should be become a clade (c[i] together with t[i]).

        # Add an inner vertex v that connects c1, c2, and c3.
        v = dm.create_central_vertex(t, c)
        taxa_string = _taxa_string_helper([t[0], t[1], t[2]], 0)
        logging.debug(f"Recursing on subproblems induced by {taxa_string}")

        # Recurse on cx+tx+v, for x in {1,2,3}.
        subtrees = [None, None, None]
        for i in range(3):
            c[i].append(v)  # Ensure that the central vertex participates as a representative of the other clades
            subtrees[i] = dnc_tree_simple(
                dm, c[i], max_n_attempts, max_clade_size, base_case_size, None
            )

            dm.log_progress()
        t_connected = connect_the_trees(subtrees, v)
        return t_connected


def sample_three_taxa(dm, taxa, max_n_attempts, max_clade_size):
    """
    Sample three taxa and see if they split the set of taxa not too unevenly. Repeat until the largest
    clade (out of three) is small enough or we have used our max number of attempts.

    Parameters:
      dm              A PartialDistanceMatrix
      taxa            A list of taxa ids
      max_n_attempts  How many samples we will attempt at most
      max_clade_size  Try not to accept a taxa sample that defines one of the clades as larger than this,
                      as a fraction of the taxa. (It is a float.)
    """
    N = (
        len(taxa) * max_clade_size
    )  # The max clade size is given as a fraction. N is the actual number of taxa to allow for a subproblem.
    best_so_far = None
    smallest_large_clade_so_far = len(
        taxa
    )  # Initialize with worst possible subproblem size
    for i in range(max_n_attempts):
        v1, v2, v3 = random.sample(taxa, 3)  # Three taxa strings
        taxa_string = _taxa_string_helper([v1, v2, v3], 0)
        logging.debug(f"Attempt: {taxa_string}.")
        c1, c2, c3 = clade_sort_taxa(
            dm, taxa, v1, v2, v3
        )  # Three groups of taxa: "clade wannabees"
        largest = max(
            len(c1), len(c2), len(c3)
        )  # Out of the three found clades, which one is largest?
        logging.debug(f"Clade sizes: {len(c1)}, {len(c2)}, and {len(c3)}")

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
    """
    Go through chosen taxa in dm (not all, because might be a subproblem) and
    decide which clade they belong to: l1, l2, or l3 (as represented by those
    leaves/taxa).  Return three lists of taxa identifiers corresponding to the
    clades (l1, l2, or l3) the belong to.
    """
    clade_assignment = {l1: [l1], l2: [l2], l3: [l3]}  # Keep track of clade assignment.
    for x in taxa:
        if x not in [l1, l2, l3]:
            clade = dm.select_clade(x, l1, l2, l3)
            clade_assignment[clade].append(x)

    return clade_assignment[l1], clade_assignment[l2], clade_assignment[l3]


def connect_the_trees(subtrees, connecting_node):
    """
    Merge the subtrees. They are sharing a unique vertex, so it is guaranteed to connect the subtrees.
    """
    first_tree = subtrees[0]
    for t in subtrees[1:]:
        first_tree.merge(t, connecting_node)
    first_tree.set_start_node(connecting_node)
    return first_tree
