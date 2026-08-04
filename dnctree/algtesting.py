import ete3
import itertools
import logging
import random
import sys

from dnctree.partialdistancematrix import PartialDistanceMatrix
from dnctree.simple_algorithm import dnc_tree_simple_centroid
from dnctree.dnctree_multi import dnc_tree_multi


def testing_divide_n_conquer(
        model,
        algorithm='multi',
        max_n_attempts=100,
        max_clade_size=0.5,
        base_case_size=100,
        core_size=None,
        base_case_method='bionj',
        partition_method='tree'
):
    """
    The purpose with this function is to be able to test the algorithm and
    see how sensitive it is to errors in the distance estimation.

    Input and parameters:

    model_tree:   A tree defining the distances.
    epsilon:      Decides the error added to distances.

    Other parameters are the same as for "divide_n_conquer".
    """
    resolved_core_size = core_size if core_size is not None else base_case_size

    if algorithm == 'multi':
        t = dnc_tree_multi(
            model,
            model.taxa,
            base_case_size=base_case_size,
            core_size=resolved_core_size,
            base_case_method=base_case_method,
            partition_method=partition_method,
        )
    elif algorithm == 'simple':
        t = dnc_tree_simple_centroid(
            model,
            model.taxa,
            max_n_attempts,
            max_clade_size,
            base_case_size,
        )
    else:
        sys.exit(f'Programming error: algorithm {algorithm} is not implemented.')
    return t


def run_alg_testing(args):
    dt = TestingPartialDistanceMatrix(args.infile, args.alg_testing, args.alg_testing_err_distribution)
    rf_dnc_list = []
    anchor_dnc_list = []
    max_reuse_list = []
    max_generation_list = []
    for base_case_size in map(int, args.alg_testing_base_case_sizes.split(',')):
        dt.reset_anchor_counter()
        t_dnc = testing_divide_n_conquer(dt,
                                         algorithm=args.alg_testing_algorithm,
                                         max_n_attempts=args.max_n_attempts,
                                         max_clade_size=args.max_clade_size,
                                         base_case_size=base_case_size,
                                         core_size=args.core_tree_size,
                                         base_case_method=args.base_case_method,
                                         partition_method=args.partition)
        dnc_estimated_tree = ete3.Tree(str(t_dnc))

        rf_dnc, rf_dnc_max, *x = dnc_estimated_tree.robinson_foulds(dt._dt, unrooted_trees=True)
        rf_dnc_list.append(f'{rf_dnc/rf_dnc_max:.5}')

        anchor_dnc_list.append(f'{dt.n_anchor_is_virtual()}')
        max_reuse_list.append(f'{dt.max_anchor_reuse()}')
        max_generation_list.append(f'{dt.max_generation()}')

    rf_string = '# Base case sizes: ' + args.alg_testing_base_case_sizes + '\n# RF distance\nRF '
    rf_string += '  '.join(rf_dnc_list)

    anchor_string = '# number of times a virtual (center) vertex was picked as one of the 3 anchors:\nANCHOR '
    anchor_string += '  '.join(anchor_dnc_list)

    max_reuse_string = '# largest number of times any single virtual vertex was picked as an anchor:\nMAXREUSE '
    max_reuse_string += '  '.join(max_reuse_list)

    max_generation_string = ('# longest chain of anchor reuse (vertex built from a vertex built from '
                              'a vertex...):\nMAXGEN ')
    max_generation_string += '  '.join(max_generation_list)

    # NJ
    if args.alg_testing_nj:
        t_nj = testing_divide_n_conquer(dt,
                                        max_n_attempts=args.max_n_attempts,
                                        max_clade_size=args.max_clade_size,
                                        base_case_size=dt.n_taxa, # Note: always do the basecase!
                                        )
        nj_estimated_tree = ete3.Tree(str(t_nj))
        rf_nj, rf_nj_max, *y = nj_estimated_tree.robinson_foulds(dt._dt, unrooted_trees=True)

        print(f'{rf_string}   {rf_nj/rf_nj_max:5}')
        print(anchor_string)
        print(max_reuse_string)
        print(max_generation_string)
    else:
        print(rf_string)
        print(anchor_string)
        print(max_reuse_string)
        print(max_generation_string)


class TestingPartialDistanceMatrix(PartialDistanceMatrix):
    def __init__(self, treefile, error_parameter, error_distribution):
        self._dt = ete3.Tree(treefile)
        self._dm = dict()
        self._true_dm = dict()
        self._seen = dict()
        self._epsilon = error_parameter

        self.taxa = [taxon.name for taxon in self._dt]
        self.n_taxa = len(self.taxa)
        self._internal_vertex_counter = 0
        self._virtual_vertices = set()
        self._n_anchor_is_virtual = 0
        self._anchor_reuse_count = dict()
        self._generation = dict()
        self._last_progress = 0
        self._n_distances_computed = 0

        for s in self.taxa:
            self._dm[s] = dict()
            self._dm[s][s] = 0.0
            self._true_dm[s] = dict()
            self._seen[s] = dict()

        # get_distance() only does an O(n) tree-wide search when given a
        # *name*; given actual node objects it just walks ancestors, so
        # caching the name->node lookup here (instead of "self._dt&s" per
        # pair) turns this loop from O(n^3) into O(n^2 * tree depth).
        name_to_node = {leaf.name: leaf for leaf in self._dt}
        for s, t in itertools.combinations(self.taxa, 2):
            s_node = name_to_node[s]
            t_node = name_to_node[t]
            d = s_node.get_distance(t_node)
            self._true_dm[t][s] = d
            self._true_dm[s][t] = d
            self._seen[s][t] = False
            self._seen[t][s] = False

        self.resample_noise(error_distribution)


    def _draw_delta(self, error_distribution):
        if error_distribution == 'uniform':
            return random.uniform(- self._epsilon, self._epsilon)
        elif error_distribution == 'normal':
            return random.gauss(0, self._epsilon)
        elif error_distribution == 'evil':
            return random.choice([- 0.99*self._epsilon, 0.99*self._epsilon])
        else:
            raise Exception(f'Bad programming: distribution "{error_distribution}" not implemented')


    def resample_noise(self, error_distribution):
        '''
        Re-draw the random error term for every pair from the cached true
        distances (self._true_dm, computed once in __init__), and reset all
        recursion-derived state (virtual vertex ids, anchor-reuse counters).
        Lets the same instance be reused across repeated trials without
        redoing the O(n^2 * tree depth) true-distance computation, which
        dominates cost and doesn't depend on the noise draw.
        '''
        self._internal_vertex_counter = 0
        self._virtual_vertices = set()
        self._n_anchor_is_virtual = 0
        self._anchor_reuse_count = dict()
        self._generation = dict()
        self._n_distances_computed = 0
        self._last_progress = 0

        self._dm = dict()
        for s in self.taxa:
            self._dm[s] = dict()
            self._dm[s][s] = 0.0

        for s, t in itertools.combinations(self.taxa, 2):
            d = self._true_dm[s][t] + self._draw_delta(error_distribution)
            if d < 0:
                d = 0
            self._dm[s][t] = d
            self._dm[t][s] = d


    def get(self, t1, t2):
        return self._dm[t1][t2]

