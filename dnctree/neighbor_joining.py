import textwrap

import fastphylo

from dnctree.tree import Tree
from dnctree.fastphylo_bridge import TaxaDistanceView, tree_from_fastphylo


def dnc_neighborjoining(dm, taxa):
    """
    An implementation of NJ meant to compute base cases in dnctree.

    dm:       a PartialDistanceMatrix instance
    taxa:     a list of identifiers (part of dm) to infer the tree for
    """
    if len(taxa) == 1:
        t = Tree()
        t.add_vertex(taxa[0])
        return t
    elif len(taxa) == 2:
        t = Tree()
        t.add_edge(taxa[0], taxa[1])
        return t
    else:
        fp_tree = fastphylo.fnj(TaxaDistanceView(dm, taxa))
        return tree_from_fastphylo(fp_tree, dm)


def _taxa_string_helper(l, indent):
    s = textwrap.fill(
        ", ".join(map(lambda t: t.decode() if isinstance(t, bytes) else t, l)),
        width=100 - indent,
    )
    return textwrap.indent(s, " " * indent)
