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
# Find the centroid vertex of a tree
#
# A centroid minimizes the size (vertex count) of the largest branch hanging
# off of it.  This# guarantees that there are at most half the leaves in any
# subtree off of the centroid.


def get_subtree_size(tree: Tree, sizes, edge):
    """
    edge=(x, y) is the directed edge from vertex x to y.
    Returns the number of vertices in the branch hanging off of x in the
    direction of y (i.e. y's subtree, once x is removed) - including y
    itself.

    Memoized: sizes is a dict shared across the recursion.
    """
    if edge in sizes:
        return sizes[edge]
    x, y = edge
    #        n[0]
    #       /
    # x -- y
    #       \
    #        n[1]
    n = [i for i in tree.adjacencies[y] if i != x]
    if not n:  # y is a leaf.
        sizes[edge] = 1  # Just y itself.
    else:
        sizes[edge] = 1 + sum(get_subtree_size(tree, sizes, (y, ni)) for ni in n)
    return sizes[edge]


def get_centroid_vertex(tree):
    """
    Returns the centroid vertex of the tree: the vertex whose removal
    leaves branches that are as balanced (in vertex count) as possible,
    i.e. the largest branch is as small as possible.
    """
    sizes = {}  # Memoization dictionary.
    centroid_candidate = None
    smallest_largest_branch_so_far = len(tree.adjacencies)  # Will always be beaten.
    for v, neighbors in tree.adjacencies.items():
        largest_branch = max(
            get_subtree_size(tree, sizes, (v, n)) for n in neighbors
        )
        if largest_branch < smallest_largest_branch_so_far:
            smallest_largest_branch_so_far = largest_branch
            centroid_candidate = v
    return centroid_candidate


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
