import random

import fastphylo

from dnctree.fastphylo_bridge import TaxaDistanceView, tree_from_fastphylo
from dnctree.tree import Tree


################################################################################
## Neighbor joining where edge weights are computed
################################################################################

# Delegates to fastphylo's BioNJ for the actual reconstruction and branch
# lengths; distances are still computed on-demand through dm, same as
# elsewhere in dnctree.


def dnc_neighborjoining_with_weights(dm, taxa):
    """
    An implementation of NJ meant to compute base cases in dnctree.

    dm:       a PartialDistanceMatrix instance
    taxa:     a list of identifiers (part of dm) to infer the tree for
    """
    t = Tree()
    if len(taxa) == 1:
        t.add_vertex(taxa[0])
        return t
    elif len(taxa) == 2:
        t.add_edge(taxa[0], taxa[1])
        ######### Possibly add distance in dm for the single edge
        return t
    else:
        fp_tree = fastphylo.bionj(TaxaDistanceView(dm, taxa))
        return tree_from_fastphylo(fp_tree, dm)


################################################################################
## Find the center vertex of a tree
################################################################################


def get_longest_path_length_to_leaf(tree: Tree, lengths, edge):
    """
    edge=(x, y) is the directed edge from vertex x to y.
    Returns the length of the longest, simple path from x to a leaf such that
    the first edge in the path is (x, y).

    The function uses memoization to avoid recomputing the same path lengths
    multiple times. It makes recursive calls to itself. The path lenghts are
    stored in the dictionary lengths.
    """
    if edge in lengths:
        return lengths[edge]
    x, y = edge
    #        n[0]
    #       /
    # x -- y
    #       \
    #        n[1]
    n = [i for i in tree.adjacencies[y] if i != x]
    if not n:  # y is a leaf.
        lengths[edge] = 1  # The edge (x, y) itself.
    else:
        lengths[edge] = 1 + max(
            get_longest_path_length_to_leaf(tree, lengths, (y, n[0])),
            get_longest_path_length_to_leaf(tree, lengths, (y, n[1])),
        )
    return lengths[edge]


def get_center_vertex(tree):
    """
    Returns the center vertex of the tree.
    """
    lenghts = {}  # Memoization dictionary.
    center_candidate = None
    shortest_longest_length_so_far = len(tree.adjacencies)  # Will always be beaten.
    for v, neighbors in tree.adjacencies.items():
        longest_length = max(
            get_longest_path_length_to_leaf(tree, lenghts, (v, n)) for n in neighbors
        )
        if longest_length < shortest_longest_length_so_far:
            shortest_longest_length_so_far = longest_length
            center_candidate = v
    return center_candidate


################################################################################
## Clades and path weights from center vertex to leaves
################################################################################


def get_weights_to_leaves(tree: Tree, dm, edge):
    """
    For a given edge (x, y), the function returns two lists in a tuple:
    - a list of leaves that can be reached when going from x via y.
    - the sum of the edge weights on the simple path from x to each leaf,
    respectively.

    Example:
        Suppose that the leaves "f" and "g" can be reached from x via y (along a
        simple path). Suppose the sum of the edge weights in the simple path
        from x to f is 5, and the sum of the edge weights in the simple path
        from x to g is 12. Then the returned tuple will be
        (["f", "g"], [5, 12]).

    """
    x, y = edge
    xy_weight = dm.get(x, y)
    #        n[0]
    #       /
    # x -- y
    #       \
    #        n[1]
    n = [i for i in tree.adjacencies[y] if i != x]
    if not n:  # y is a leaf.
        return [y], [xy_weight]
    else:
        leaves0, weights0 = get_weights_to_leaves(tree, dm, (y, n[0]))
        leaves1, weights1 = get_weights_to_leaves(tree, dm, (y, n[1]))
        return (leaves0 + leaves1, [w + xy_weight for w in weights0 + weights1])


def get_clades_and_dists_to_center(tree: Tree, dm, center_vertex):
    """
    Splits the tree into three clades (the leaves), and returns three dictionaries:
    - clade_to_taxa: maps clade number 0, 1, and 2 to a list of taxa in that clade.
    - taxon_to_clade: maps each taxon to the clade it is in.
    - dist_to_center: maps each taxon to the distance from it to the center vertex.
    """
    clade_to_taxa = {}
    taxon_to_clade = {}
    dist_to_center = {}
    for clade_number, neighbor in enumerate(tree.adjacencies[center_vertex]):
        leaves, weights = get_weights_to_leaves(tree, dm, (center_vertex, neighbor))
        clade_to_taxa[clade_number] = leaves
        taxon_to_clade.update({taxon: clade_number for taxon in leaves})
        dist_to_center.update({taxon: weight for taxon, weight in zip(leaves, weights)})
    return clade_to_taxa, taxon_to_clade, dist_to_center


################################################################################
## The new DNC tree k algorithm
################################################################################


def dnc_tree_k(dm, taxa, base_case_size=5, core_size=5):
    """
    Construct a tree with the DNC tree k algorithm that runs NJ on a sample of
    core_size taxa, referred to as the core.
    If the number of taxa is less than or equal to base_case_size, then
    the tree is constructed using NJ on all taxa.

    Returns the tree.
    """
    if len(taxa) <= base_case_size:
        return dnc_neighborjoining_with_weights(dm, taxa)
    else:
        # Sample core taxa and build the T_{core} using NJ.
        core_taxa = random.sample(taxa, core_size) if len(taxa) > core_size else taxa
        non_core_taxa = set(taxa) - set(core_taxa)
        core_tree = dnc_neighborjoining_with_weights(dm, core_taxa)

        # Determine the center vertex and partition the core taxa into three clades.
        center_vertex = get_center_vertex(core_tree)
        (
            clade_to_core_taxa,  # Maps 0, 1, and 2 to list of taxa, respectively.
            taxon_to_clade,  # Maps each taxon to the clade it is in, e.g., "T1" -> 0.
            dist_to_center,  # Maps each taxon to the distance from it to the center vertex.
        ) = get_clades_and_dists_to_center(core_tree, dm, center_vertex)

        # Update the partial distance matrix with distances from core taxa to the center,
        # inferred from the edge weights of T_{core}.
        for w in core_taxa:
            dm.set(w, center_vertex, dist_to_center[w])

        # Assign non-core taxa to clades.

        # sum_distances[x] = \sum_{y \in V_{core}} d(x, y)
        # for any taxa x.
        sum_distances = {x: sum(dm.get(x, y) for y in core_taxa) for x in taxa}

        k = len(core_taxa)
        clade_to_non_core_taxa = {0: [], 1: [], 2: []}  # Will be filled in below.
        for v in non_core_taxa:
            smallest_q_so_far = None
            for w in core_taxa:
                q = (k - 2) * dm.get(v, w) - sum_distances[v] - sum_distances[w]
                if smallest_q_so_far is None or q < smallest_q_so_far:
                    smallest_q_so_far = q
                    taxon_to_clade[v] = taxon_to_clade[w]  # w is v's new mate.
            clade_to_non_core_taxa[taxon_to_clade[v]].append(v)

        # Estimate distance between each non-core taxon and the center vertex.

        # C[i] = \sum_[w \in core_taxa_in_clade_i] d(w, center_vertex)
        # for i \in {0, 1, 2}.
        C = {i: sum(dist_to_center[w] for w in clade_to_core_taxa[i]) for i in range(3)}

        for v in non_core_taxa:
            # Indices a and b are the two clades that v is not in.
            a, b = {0: (1, 2), 1: (0, 2), 2: (0, 1)}[taxon_to_clade[v]]
            n_taxa_in_a_and_b = len(clade_to_core_taxa[a]) + len(clade_to_core_taxa[b])
            # Estimate distance from v to center.
            sum_dist_estimates = (
                sum(dm.get(v, w) for w in clade_to_core_taxa[a])
                + sum(dm.get(v, w) for w in clade_to_core_taxa[b])
                - C[a]
                - C[b]
            )
            # The distance is the average of the estimates.
            dm.set(v, center_vertex, sum_dist_estimates / n_taxa_in_a_and_b)

        # Finally set the distance between the center vertex and iteself.
        dm.set(center_vertex, center_vertex, 0)

        # Recurse and get three subtrees, each sharing the center vertex.

        subtrees = [
            dnc_tree_k(
                dm,
                clade_to_core_taxa[i] + clade_to_non_core_taxa[i] + [center_vertex],
                base_case_size=base_case_size,
                core_size=core_size,
            )
            for i in (0, 1, 2)
        ]

        # Merge the subtrees on the center vertex.
        subtrees[0].merge(subtrees[1], center_vertex)
        subtrees[0].merge(subtrees[2], center_vertex)
        subtrees[0].set_start_node(center_vertex)
        return subtrees[0]
