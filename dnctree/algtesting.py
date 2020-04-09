import ete3
import itertools
import random

from dnctree.partialdistancematrix import PartialDistanceMatrix
from dnctree.dnctree import testing_divide_n_conquer


def run_alg_testing(args):
    dt = TestingPartialDistanceMatrix(args.infile, args.alg_testing)
    t_dnc = testing_divide_n_conquer(dt,
                                 max_n_attempts=args.max_n_attempts,
                                 max_clade_size=args.max_clade_size,
                                 base_case_size=args.base_case_size,
                                 first_triple=args.first_triple,
                                 verbose=args.verbose)
    #    print(dt._dt.as_string(schema='newick', suppress_edge_lengths=True))
    # DivNConq algorithm:
    #print(t_dnc)
    dnc_estimated_tree = ete3.Tree(str(t_dnc));
    #print(dnc_estimated_tree)
    rf_dnc, rf_dnc_max, *x = dnc_estimated_tree.robinson_foulds(dt._dt, unrooted_trees=True)

    # NJ
    if args.alg_testing_nj:
        t_nj = testing_divide_n_conquer(dt,
                                        max_n_attempts=args.max_n_attempts,
                                        max_clade_size=args.max_clade_size,
                                        base_case_size=dt.n_taxa,
                                        first_triple=args.first_triple,
                                        verbose=args.verbose)
        nj_estimated_tree = ete3.Tree(str(t_nj))
        rf_nj, rf_nj_max, *y = nj_estimated_tree.robinson_foulds(dt._dt, unrooted_trees=True)
        print(rf_nj)
        print(f'{rf_dnc / rf_dnc_max:.5}\t{rf_nj / rf_nj_max:.5}')
    else:
        print(rf_dnc / rf_dnc_max)


class TestingPartialDistanceMatrix(PartialDistanceMatrix):
    def __init__(self, treefile, error_parameter, verbose=False):
        self._dt = ete3.Tree(treefile)
        self._dm = dict()
        self._seen = dict()
        self._epsilon = error_parameter

        self.taxa = [taxon.name for taxon in self._dt]
        self.n_taxa = len(self.taxa)
        self._internal_vertex_counter = 0
        self._last_progress = 0
        self._n_distances_computed = 0

        for s in self.taxa:
            self._dm[s] = dict()
            self._seen[s] = dict()

        for s, t in itertools.combinations(self.taxa, 2):
            s_node = self._dt&s
            d = s_node.get_distance(t)
            delta = d * random.uniform(- self._epsilon, self._epsilon)
#            print(d, delta)
            d += delta          # Perturb
            self._dm[s][t] = d
            self._dm[t][s] = d
            self._seen[s][t] = False
            self._seen[t][s] = False


    def get(self, t1, t2):
#        if not self._seen[t1][t2]:
#            self._n_distances_computed += 1
#            self._seen[t1][t2] = True
#            self._seen[t2][t1] = True
        return self._dm[t1][t2]
