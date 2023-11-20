import ete3
import itertools
import random
import sys

from dnctree.partialdistancematrix import PartialDistanceMatrix
from dnctree import testing_divide_n_conquer
from tree_matching_distance import distance as matching_distance


def run_alg_testing(args, verbosity=[]):
    dt = TestingPartialDistanceMatrix(args.infile, args.alg_testing)
    rf_dnc_list = []
    tm_dnc_list = []
    for base_case_size in map(int, args.alg_testing_base_case_sizes.split(',')):
        t_dnc = testing_divide_n_conquer(dt,
                                         max_n_attempts=args.max_n_attempts,
                                         max_clade_size=args.max_clade_size,
                                         base_case_size=base_case_size,
                                         first_triple=args.first_triple,
                                         verbose=verbosity)
        dnc_estimated_tree = ete3.Tree(str(t_dnc));

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
                                        first_triple=args.first_triple,
                                        verbose=verbosity)
        nj_estimated_tree = ete3.Tree(str(t_nj))
        rf_nj, rf_nj_max, *y = nj_estimated_tree.robinson_foulds(dt._dt, unrooted_trees=True)
        tm_nj = matching_distance(dt._dt, nj_estimated_tree)

        print(f'{rf_string}   {rf_nj/rf_nj_max:5}')
        print(tm_string + f'  {tm_nj}')
    else:
        print(rf_string)
        print(tm_string)


class TestingPartialDistanceMatrix(PartialDistanceMatrix):
    def __init__(self, treefile, error_parameter, verbose=[]):
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
            self._true_dm[s] = dict()
            self._seen[s] = dict()

        for s, t in itertools.combinations(self.taxa, 2):
            s_node = self._dt&s
            d = s_node.get_distance(t)
            self._true_dm[t][s] = d
            self._true_dm[s][t] = d
            delta = random.uniform(- self._epsilon, self._epsilon)
#            delta = d * random.uniform(- self._epsilon, self._epsilon)
#            print(f'{s}  {t}\t {d:.2}\t {delta:.2}\t {abs(100*delta/d):.2f}%', file=sys.stderr)
            d += delta          # Perturb
            if d < 0:
                print(f'{s}  {t}\t {d:.2}\t {delta:.2}\t {abs(100*delta/d):.2f}%', file=sys.stderr)
            self._dm[s][t] = d
            self._dm[t][s] = d
            self._seen[s][t] = False
            self._seen[t][s] = False


    def get(self, t1, t2):
        return self._dm[t1][t2]


    # def quartet_test(self, l1, l2, l3, x):
    #     '''
    #     Return the leaf (l1, l2, or l3) of which leaf a quartet test suggests x belongs to.
    #     '''
    #     def q_diff(a, b, c, d):
    #         return self.get(a, c) + self.get(b, d) - self.get(a, b) - self.get(c, d)

    #     def true_q_diff(a, b, c, d):
    #         return self._true_dm[a][c] + self._true_dm[b][d] - self._true_dm[a][b] - self._true_dm[c][d]


    #     # First look at noisy distances
    #     options = [(l1, q_diff(x, l1, l2, l3)),
    #                (l2, q_diff(x, l2, l1, l3)),
    #                (l3, q_diff(x, l3, l1, l2)),]
    #     choice = max(options, key = lambda pair: pair[1])

        # Now use the ideal distances, but give up on internal vertices because not updating their distances
        # try:
        #     true_options = [(l1, true_q_diff(x, l1, l2, l3)),
        #                     (l2, true_q_diff(x, l2, l1, l3)),
        #                     (l3, true_q_diff(x, l3, l1, l2)),]
        #     true_choice = max(true_options, key = lambda pair: pair[1])

        #     if choice[0] != true_choice[0]:
        #         print(f'{x:4}\t     choice: {choice} options: {options}\n    \ttrue choice: {true_choice} from {true_options}\n', file=sys.stderr)
        #         for a, b in itertools.combinations([x, l1, l2, l3], 2):
        #             print(f'    \tD({a},{b})={self.get(a,b):.3}', file=sys.stderr)
        #         print(file=sys.stderr)
        # except:
        #     print(f'Not available for {x}')

        return choice[0]
