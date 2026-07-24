import ete3
import itertools
import logging
import random
import sys
import warnings

from dnctree.partialdistancematrix import PartialDistanceMatrix
from dnctree.simple_algorithm import dnc_tree_simple
from dnctree.dnctree_k import dnc_tree_k

# tree_matching_distance's own __init__.py uses the deprecated pkg_resources
# API for its version lookup; harmless, but noisy every run.
warnings.filterwarnings("ignore", message=r"pkg_resources is deprecated.*", category=UserWarning)
from tree_matching_distance import distance as matching_distance


def testing_divide_n_conquer(
        model,
        algorithm='core-tree',
        max_n_attempts=100,
        max_clade_size=0.5,
        base_case_size=100
):
    """
    The purpose with this function is to be able to test the algorithm and
    see how sensitive it is to errors in the distance estimation.

    Input and parameters:

    model_tree:   A tree defining the distances.
    epsilon:      Decides the error added to distances.

    Other parameters are the same as for "divide_n_conquer".
    """

    if algorithm == 'core-tree':
        t = dnc_tree_k(
            model,
            model.taxa,
            base_case_size,
            base_case_size
        )
    elif algorithm == 'simple':
        t = dnc_tree_simple(
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
    tm_dnc_list = []
    for base_case_size in map(int, args.alg_testing_base_case_sizes.split(',')):
        t_dnc = testing_divide_n_conquer(dt,
                                         algorithm=args.alg_testing_algorithm,
                                         max_n_attempts=args.max_n_attempts,
                                         max_clade_size=args.max_clade_size,
                                         base_case_size=base_case_size)
        dnc_estimated_tree = ete3.Tree(str(t_dnc))

        rf_dnc, rf_dnc_max, *x = dnc_estimated_tree.robinson_foulds(dt._dt, unrooted_trees=True)
        rf_dnc_list.append(f'{rf_dnc/rf_dnc_max:.5}')

        tm_dnc = matching_distance(dnc_estimated_tree, dt._dt)
        tm_dnc_list.append(f'{tm_dnc}')

    rf_string = '# Base case sizes: ' + args.alg_testing_base_case_sizes + '\n# RF distance\nRF '
    rf_string += '  '.join(rf_dnc_list)

    tm_string = '# tree matching distance:\nTM '
    tm_string += '  '.join(tm_dnc_list)

    # NJ
    if args.alg_testing_nj:
        t_nj = testing_divide_n_conquer(dt,
                                        max_n_attempts=args.max_n_attempts,
                                        max_clade_size=args.max_clade_size,
                                        base_case_size=dt.n_taxa, # Note: always do the basecase!
                                        )
        nj_estimated_tree = ete3.Tree(str(t_nj))
        rf_nj, rf_nj_max, *y = nj_estimated_tree.robinson_foulds(dt._dt, unrooted_trees=True)
        tm_nj = matching_distance(dt._dt, nj_estimated_tree)

        print(f'{rf_string}   {rf_nj/rf_nj_max:5}')
        print(tm_string + f'  {tm_nj}')
    else:
        print(rf_string)
        print(tm_string)


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
        self._last_progress = 0
        self._n_distances_computed = 0

        for s in self.taxa:
            self._dm[s] = dict()
            self._dm[s][s] = 0.0
            self._true_dm[s] = dict()
            self._seen[s] = dict()

        for s, t in itertools.combinations(self.taxa, 2):
            s_node = self._dt&s
            d = s_node.get_distance(t)
            self._true_dm[t][s] = d
            self._true_dm[s][t] = d
            if error_distribution == 'uniform':
                delta = random.uniform(- self._epsilon, self._epsilon)
            elif error_distribution == 'normal':
                delta = random.gauss(0, self._epsilon)
            else:
                raise Exception(f'Bad programming: distribution "{args.error_model}" not implemented')

            d += delta          # Perturb
            if d < 0:
                d = 0
            self._dm[s][t] = d
            self._dm[t][s] = d
            self._seen[s][t] = False
            self._seen[t][s] = False


    def get(self, t1, t2):
        return self._dm[t1][t2]

