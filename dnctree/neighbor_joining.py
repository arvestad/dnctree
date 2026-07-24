import itertools
import textwrap

from dnctree.tree import Tree


def nj_selection_function(dm, current_leaves):
    def sum_distances(a):
        sum_for_a = 0
        for b in current_leaves:
            if a != b:
                sum_for_a += dm.get(a, b)
        return sum_for_a

    row_sum = {x: sum_distances(x) for x in current_leaves}

    n_minus_2 = len(current_leaves) - 2
    Q_ij_min = n_minus_2 * dm.get(
        current_leaves[0], current_leaves[1]
    )  # Arbitrary but safe starting value that won't be chosen.
    best_so_far = current_leaves[0], current_leaves[1]
    for x, y in itertools.combinations(current_leaves, 2):
        Q_ij = n_minus_2 * dm.get(x, y) - row_sum[x] - row_sum[y]
        if Q_ij < Q_ij_min:
            Q_ij_min = Q_ij
            best_so_far = x, y

    return best_so_far


def dnc_neighborjoining(dm, taxa):
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
        return t
    else:
        current_leaves = taxa[:]  # Slicing out a copy
        while len(current_leaves) > 3:
            x, y = nj_selection_function(dm, current_leaves)
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
    s = textwrap.fill(
        ", ".join(map(lambda t: t.decode() if isinstance(t, bytes) else t, l)),
        width=100 - indent,
    )
    return textwrap.indent(s, " " * indent)
