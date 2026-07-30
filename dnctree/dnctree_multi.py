import random
from collections import deque

from dnctree.neighbor_joining import dnc_neighborjoining
from dnctree.dnctree_k import (
    dnc_neighborjoining_with_weights,
    get_centroid_vertex,
    get_clades_and_dists_to_center,
)
from dnctree.simple_algorithm import compute_centroid_distances


################################################################################
## dnc_tree_multi: a variant of dnc_tree_k with a choice of non-core-taxa
## partition method (quartet testing vs tree testing), and a run-time choice
## of NJ vs BioNJ for the terminal base case.
##
## Unlike dnc_tree_k, the three subproblems are never joined on a synthetic
## vertex with estimated distances. Instead, T_core's centroid is used only
## to pick, per clade, the core taxon closest to it (a purely topological/
## branch-length-based SELECTION, not a distance estimate) - those three
## representatives then play the role of the three anchors in
## dnc_tree_simple_centroid, and the same distance-free merge mechanism is
## reused: never estimate a distance to a synthetic vertex, only ever
## compare real, directly estimated pairwise distances. See
## dnctree.simple_algorithm for the mechanism and rationale.
################################################################################


def _run_base_case(dm, taxa, base_case_method):
    if base_case_method == 'nj':
        return dnc_neighborjoining(dm, taxa)
    return dnc_neighborjoining_with_weights(dm, taxa)  # 'bionj', default


def _reachable_leaves(tree, edge):
    """
    edge=(x, y): the directed edge from x to y, topology only (no dm, no
    weights). Returns the set of leaves reachable from x via y.
    """
    x, y = edge
    n = [i for i in tree.adjacencies[y] if i != x]
    if not n:  # y is a leaf.
        return {y}
    leaves = set()
    for ni in n:
        leaves |= _reachable_leaves(tree, (y, ni))
    return leaves


def _clade_counts(leaves, taxon_to_clade):
    counts = [0, 0, 0]
    for leaf in leaves:
        counts[taxon_to_clade[leaf]] += 1
    return counts


def _gini_impurity(counts):
    """
    1 - sum(p_i^2) over the clade proportions in counts: 0 for a pure
    (single-clade) group, approaching 1 the more evenly the group is split
    across clades. Standard small-sample group-purity measure (the default
    CART splitting criterion), used here instead of e.g. entropy mainly
    because it needs no logarithms - the two rank the same way in practice
    for this comparison.
    """
    total = sum(counts)
    if total == 0:
        return 0.0
    return 1.0 - sum((c / total) ** 2 for c in counts)


def _select_clade_representatives(clade_to_core_taxa, dist_to_center):
    """
    x[i] = the core taxon in clade i estimated (via T_core's own branch
    lengths, purely as a selection heuristic) to be closest to T_core's
    centroid. Used both to assign non-core taxa in quartet mode, and as
    each clade's real stand-in for the shared merge point - no distance is
    ever estimated for a synthetic vertex from this.
    """
    return {
        i: min(clade_to_core_taxa[i], key=lambda w: dist_to_center[w])
        for i in range(3)
    }


def _assign_non_core_taxa_quartet(dm, non_core_taxa, taxon_to_clade, x):
    """
    Assign each non-core taxon v to a clade via a quartet test against one
    representative taxon per clade (x[0], x[1], x[2] - see
    _select_clade_representatives).
    """
    clade_to_non_core_taxa = {0: [], 1: [], 2: []}
    for v in non_core_taxa:
        chosen = dm.select_clade(v, x[0], x[1], x[2])
        clade_idx = taxon_to_clade[chosen]
        taxon_to_clade[v] = clade_idx
        clade_to_non_core_taxa[clade_idx].append(v)
    return clade_to_non_core_taxa


def _assign_non_core_taxa_tree(dm, core_taxa, non_core_taxa, taxon_to_clade):
    """
    Assign each non-core taxon v to a clade by building a fresh NJ tree on
    W = core_taxa + [v]. Removing v and its unique neighbor p splits the
    remaining core taxa into exactly two groups (p always has degree 3 in a
    binary NJ tree: v, plus two other directions). If v truly belongs to
    clade C_i and T_W correctly reflects that, one of those two groups
    should be dominated by (ideally, contain only) C_i, while the other
    contains a mix of all three clades. So: compute each side's clade-
    membership counts, pick whichever side is LESS diverse (lowest Gini
    impurity - see _gini_impurity), and assign v to that side's dominant
    clade.

    Much more expensive than quartet testing (a separate NJ tree per
    non-core taxon), and only needs topology, so it always uses plain NJ.
    """
    clade_to_non_core_taxa = {0: [], 1: [], 2: []}
    for v in non_core_taxa:
        W = core_taxa + [v]
        t_w = dnc_neighborjoining(dm, W)

        p = t_w.adjacencies[v][0]
        other_dirs = [n for n in t_w.adjacencies[p] if n != v]
        assert len(other_dirs) == 2, f"expected p to have degree 3, got {t_w.adjacencies[p]!r}"
        side_counts = [
            _clade_counts(_reachable_leaves(t_w, (p, n)), taxon_to_clade)
            for n in other_dirs
        ]

        purer = min(side_counts, key=_gini_impurity)
        clade_idx = max(range(3), key=lambda i: purer[i])

        taxon_to_clade[v] = clade_idx
        clade_to_non_core_taxa[clade_idx].append(v)
    return clade_to_non_core_taxa


def dnc_tree_multi(dm, taxa, base_case_size=5, core_size=5,
                    base_case_method='bionj', partition_method='quartet'):
    """
    A variant of dnc_tree_k with the same overall divide-and-conquer
    structure (sample core taxa, find T_core's centroid, partition core taxa
    into three clades) but differing in two ways:

    1. A choice of how non-core taxa get assigned to a clade:
       - 'quartet' (default): a 4-point quartet test against one
         representative taxon per clade (see _assign_non_core_taxa_quartet).
         Cheap - reuses PartialDistanceMatrix.select_clade and the
         already-computed dist_to_center, no new tree construction.
       - 'tree': builds a fresh NJ tree per non-core taxon and inherits the
         dominant clade of whichever side of its attachment point is less
         diverse across clades (see _assign_non_core_taxa_tree). Much more
         expensive - one NJ tree per non-core taxon.

    2. The three subproblems are joined the way dnc_tree_simple_centroid
       joins its three subproblems: never on a synthetic vertex with
       estimated distances. T_core's centroid is only ever used to *select*
       (via T_core's own branch lengths) which core taxon in each clade is
       closest to it; those three real representatives then play the role
       of the simple algorithm's three anchors, and are merged via the same
       distance-free mechanism (see dnctree.simple_algorithm). No distance
       is ever estimated for a synthetic vertex - only real, directly
       estimated pairwise distances are ever compared.

    base_case_method ('nj' or 'bionj', default 'bionj') controls only the
    terminal base-case tree construction; the core tree (and, in 'tree'
    partition mode, every per-taxon probe tree) always uses the method
    documented at their respective call sites below.

    dm:       A PartialDistanceMatrix instance.
    taxa:     A list of identifiers (part of dm) to infer the tree for.

    Returns the tree.
    """
    if len(taxa) <= base_case_size:
        return _run_base_case(dm, taxa, base_case_method)
    else:
        # Sample core taxa and build the T_{core} using NJ.
        core_taxa = random.sample(taxa, core_size) if len(taxa) > core_size else taxa
        non_core_taxa = list(set(taxa) - set(core_taxa))
        # Always BioNJ: get_clades_and_dists_to_center needs real edge
        # weights, which plain NJ doesn't currently provide (fastphylo's NJ
        # returns a -1.0 "no branch length" sentinel that isn't recorded
        # into dm - see fastphylo_bridge.tree_from_fastphylo). Revisit this
        # once fastphylo's NJ can provide branch lengths.
        core_tree = dnc_neighborjoining_with_weights(dm, core_taxa)

        # Partition the core taxa into three clades via T_core's centroid.
        # dist_to_center is only ever used below to SELECT clade
        # representatives (T_core's own branch lengths as a heuristic) -
        # never to estimate a distance for a synthetic vertex.
        center_vertex = get_centroid_vertex(core_tree)
        clade_to_core_taxa, taxon_to_clade, dist_to_center = get_clades_and_dists_to_center(
            core_tree, dm, center_vertex
        )

        # x[i] = clade i's real representative: the core taxon closest to
        # T_core's centroid. Used by quartet assignment below, and always
        # needed for the merge step further down.
        x = _select_clade_representatives(clade_to_core_taxa, dist_to_center)

        # Assign non-core taxa to clades.
        if partition_method == 'tree':
            clade_to_non_core_taxa = _assign_non_core_taxa_tree(
                dm, core_taxa, non_core_taxa, taxon_to_clade
            )
        else:
            clade_to_non_core_taxa = _assign_non_core_taxa_quartet(
                dm, non_core_taxa, taxon_to_clade, x
            )

        # Merge the three subproblems the way dnc_tree_simple_centroid
        # does: x[0], x[1], x[2] play the role of the simple algorithm's
        # three anchors. For clade i, the closer (by real pairwise
        # distance) of the OTHER two representatives stands in for where
        # the shared centroid attaches during that subproblem's own
        # recursion; once solved, it's a leaf (real taxa are never
        # anything else), and gets replaced by a fresh, purely topological
        # vertex (never assigned any distance) before merging.
        d_c = compute_centroid_distances(dm, x[0], x[1], x[2])
        c_id = dm.create_unique_vertex_id()

        subtrees = [None, None, None]
        for i in range(3):
            other_two = [x[j] for j in range(3) if j != i]
            representative = min(other_two, key=lambda w: d_c[w])
            subtree = dnc_tree_multi(
                dm,
                clade_to_core_taxa[i] + clade_to_non_core_taxa[i] + [representative],
                base_case_size=base_case_size,
                core_size=core_size,
                base_case_method=base_case_method,
                partition_method=partition_method,
            )
            subtree.replace_leaf(representative, c_id)
            subtrees[i] = subtree
            dm.log_progress()

        subtrees[0].merge(subtrees[1], c_id)
        subtrees[0].merge(subtrees[2], c_id)
        subtrees[0].set_start_node(c_id)
        return subtrees[0]
