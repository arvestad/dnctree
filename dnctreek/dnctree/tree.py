import sys

class Tree:
    '''
    Simple class to store unrooted trees by recording their edges.
    There is no support for detecting cycles.
    '''
    def __init__(self):
        self.adjacencies = dict()
        self.start = None        # Not a root, but decent starting point for conversion to string.

    def add_vertex(self, v):
        '''
        Add a vertex without connections to other nodes (yet).
        If you have an edge to add, use add_edge(self, e) instead.
        '''
        d = self.adjacencies
        if v not in d:
            d[v] = list()


    def add_edge(self, v, w):
        d = self.adjacencies
        if not v in d:
            d[v] = list()
        if not w in d:
            d[w] = list()
        d[v].append(w)
        d[w].append(v)


    def set_start_node(self, v):
        '''
        Recomment this node as a starting point for recursion
        when converting tree to a Newick string.
        '''
        self.start = v


    def merge(self, other, connecting_node):
        '''
        Add edges from other tree into this tree.

        It is assumed that the "connecting_node" is the only node present in both trees.
        Its neighborlists gets special treatment.

        There will be bugs if node ids (other than connecting_node) are not unique.
        '''
        self.adjacencies[connecting_node].extend(other.adjacencies[connecting_node])
        self.adjacencies = {**other.adjacencies, **self.adjacencies} # Order is important, because self.adjacencies[connecting_node] is the one  to keep!


    def __str__(self):
        vertices = list(self.adjacencies.keys())
        if len(vertices) == 0:
            raise Exception('Bug: there are no vertices in the tree.')
        else:
            if self.start is None:
                start_node = vertices[0]
            else:
                start_node = self.start
            s = _tree_to_string_helper(self.adjacencies, start_node, None)
            return s + ';'        # Ensure proper Newick format with terminating semi-colon

def _tree_to_string_helper(d, v, source):
    '''
    Convert the tree to a string by starting in the given node v.
    The source argument is pointing out where the recursion started,
    and we should not continue in that direction.
    '''
    n=0
    subtrees = list()
    for neighbor in d[v]:
        if neighbor != source:  # We did not come from this neighbor
            n += 1
            subtrees.append(_tree_to_string_helper(d, neighbor, v))
    if n == 0:
        return v.decode() if isinstance(v, bytes) else v
    else:
        return '(' + ','.join(subtrees) + ')'
