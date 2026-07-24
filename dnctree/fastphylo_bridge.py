"""
Glue code between dnctree's PartialDistanceMatrix/Tree and fastphylo's
distance/tree types, used to delegate NJ-family base-case reconstruction to
fastphylo's C++ implementations.
"""
from dnctree.tree import Tree


class TaxaDistanceView:
    """
    Adapts a PartialDistanceMatrix, restricted to a specific ordered list of
    taxa, to fastphylo's DistanceProtocol (__len__, __getitem__((i, j)),
    names()). Distances are still computed on-demand through dm.get(), same
    as anywhere else in dnctree - this only changes who builds the tree from
    them, not when they're computed.
    """
    def __init__(self, dm, taxa):
        self._dm = dm
        self._taxa = list(taxa)

    def __len__(self):
        return len(self._taxa)

    def __getitem__(self, key):
        i, j = key
        if i == j:
            return 0.0
        return self._dm.get(self._taxa[i], self._taxa[j])

    def names(self):
        return self._taxa


def tree_from_fastphylo(fp_tree, dm):
    """
    Convert a fastphylo.Tree into a dnctree Tree. Leaf vertices keep their
    real taxon names (as already recorded in fp_tree.leaves); internal
    vertices get a fresh dnctree-style id via dm.create_unique_vertex_id().
    Edge weights are recorded into dm so that downstream code (e.g.
    get_weights_to_leaves) can look them up via dm.get(u, v), the same way
    it already does for trees built by dnctree's own NJ. The -1.0 sentinel
    fastphylo uses for "no branch length computed" (plain NJ/FNJ) is not
    recorded, since it isn't a real distance.

    Tree.__str__ starts its Newick traversal from an arbitrary vertex unless
    told otherwise, and silently drops that vertex's own label if it happens
    to be a leaf (its label is only ever emitted when some *other* vertex
    visits it). So the start node is explicitly set to one of the internal
    vertices fastphylo created, never a leaf.
    """
    t = Tree()
    id_to_name = dict(fp_tree.leaves)
    internal_name = None
    for u, v, weight in fp_tree.edges:
        for vid in (u, v):
            if vid not in id_to_name:
                name = dm.create_unique_vertex_id()
                id_to_name[vid] = name
                dm._dm[name] = dict()
                internal_name = name
        u_name, v_name = id_to_name[u], id_to_name[v]
        t.add_edge(u_name, v_name)
        if weight != -1.0:
            dm.set(u_name, v_name, weight)
    if internal_name is not None:
        t.set_start_node(internal_name)
    return t
