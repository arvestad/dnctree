import itertools
import random

from dnctree.tree import Tree


################################################################################
## Neighbor joining where edge weights are computed
################################################################################

# Different from the version in `__init__.py` as this one calculates the edge
# weights, adding them to the distance matrix.


def nj_selection_function_with_weights(dm, current_leaves):
    def sum_distances(a):
        sum_for_a = 0
        for b in current_leaves:
            if a != b:
                sum_for_a += dm.get(a, b)
        return sum_for_a

    row_sum = {}
    for x in current_leaves:
        row_sum[x] = sum_distances(x)

    n_minus_2 = len(current_leaves) - 2
    Q_ij_min = n_minus_2 * dm.get(
        current_leaves[0], current_leaves[1]
    )  # Arbitrary but safe starting value that won't be chosen.
    best_so_far = None
    for x, y in itertools.combinations(current_leaves, 2):
        Q_ij = n_minus_2 * dm.get(x, y) - row_sum[x] - row_sum[y]
        if best_so_far is None or Q_ij < Q_ij_min:
            Q_ij_min = Q_ij
            best_so_far = (x, y, row_sum[x], row_sum[y])

    return best_so_far


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
        current_leaves = taxa[:]  # Slicing out a copy
        while len(current_leaves) > 3:
            x, y, sum_x_to_others, sum_y_to_others = nj_selection_function_with_weights(
                dm, current_leaves
            )
            other_taxa = current_leaves[:]
            other_taxa.remove(x)
            other_taxa.remove(y)
            new_v = dm.create_representative_vertex(x, y, other_taxa)

            d_x_newv = 0.5 * dm.get(x, y) + (1 / (2 * (len(current_leaves) - 2))) * (
                sum_x_to_others - sum_y_to_others
            )
            d_y_newv = dm.get(x, y) - d_x_newv
            dm.set(x, new_v, d_x_newv)
            dm.set(y, new_v, d_y_newv)

            other_taxa.append(new_v)
            t.add_edge(x, new_v)
            t.add_edge(y, new_v)
            current_leaves = other_taxa

        c = dm.create_unique_vertex_id()
        for leaf in current_leaves:
            t.add_edge(leaf, c)
        t.set_start_node(c)
        x, y, z = current_leaves
        d_xc = 0.5 * dm.get(x, y) + (1 / 2) * (dm.get(x, z) - dm.get(y, z))
        d_yc = dm.get(x, y) - d_xc
        d_zc = dm.get(x, z) - d_xc
        dm._dm[c] = dict()
        dm.set(x, c, d_xc)
        dm.set(y, c, d_yc)
        dm.set(z, c, d_zc)
        return t


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

        # In the following we will assign non-core taxa to the clades. We will do this
        # by assigning each non-core taxon v to the clade that contains the core
        # taxon w that minimizes the following quantity:
        #
        # Q(v, w) = (n - 2) d(v, w) - \sum_{x=1}^n d(v, x) - \sum_{x=1}^n d(w, x),
        #
        # where n is the number of taxa (i.e., core taxa plus the non-core taxon v).
        # Let C be the set of core taxa. Hence |C| = k. Let A = C \cup {v}, that is the
        # core taxa plus the non-core taxon v. The above quantity Q can be rewritten as
        # follows:
        #
        # Q(v, w) = (k + 1 - 2) d(v, w) - \sum_{x \in A} d(v, x) - \sum_{x \in A} d(w, x).
        #
        # We can rewrite this as follows:
        #
        # Q(v, w) = (k - 1) d(v, w) - \sum_{x \in C} d(v, x) - d(v,v) - \sum_{x \in C} d(w, x) - d(w, v),
        #
        # which is the same as
        #
        # Q(v, w) = (k - 2) d(v, w) - \sum_{x \in C} d(v, x) - \sum_{x \in C} d(w, x).
        #
        # This formula is used below.

        # sum_distances[v] = \sum_{x \in C} d(v, x)
        # for any taxa v.
        sum_distances = {v: sum(dm.get(v, x) for x in core_taxa) for v in taxa}

        k = len(core_taxa)
        clade_to_non_core_taxa = {0: [], 1: [], 2: []}  # Will be filled in below.
        for v in non_core_taxa:
            smallest_q_so_far = None
            for w in core_taxa:
                # q = Q(v, w) = (k - 2) d(v, w) - \sum_{x \in C} d(v, x) - \sum_{x \in C} d(w, x).
                q = (k - 2) * dm.get(v, w) - sum_distances[v] - sum_distances[w]
                if smallest_q_so_far is None or q < smallest_q_so_far:
                    smallest_q_so_far = q
                    taxon_to_clade[v] = taxon_to_clade[w]  # w is v's new mate.
            clade_to_non_core_taxa[taxon_to_clade[v]].append(v)

        # Estimate distance between each non-core taxon and the center vertex.
            
        # The logic is as follows. Suppose that v is a non-core taxon in clade 0.
        # Let w be a core taxon in clade 1. An estimate of the distance from v to
        # the center vertex, parameterized by w, is
        # d(v, w) - d(w, center_vertex).
        # The final estimate of the distance from v to the center vertex is the
        # average of the estimates over all core taxa w in clades 1 and 2.
        # The distance d(w, center_vertex) is obtained from the dictionary
        # dist_to_center.

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
