import math
import sys
import dnctree.exceptions as dnc


class PartialDistanceMatrix:
    '''
    Represent distance matrices using a dictionary instead of a full (half) square matrix,
    because we are not going to compute all pairwise distances.

    This class pretends it is a distance matrix. Distances are computed on-demand, when
    they are asked for. If pahmm is set to True, then distances are computed using
    the pahmm library.
    '''
    def __init__(self, msa, distance_fcn, verbose=[]):
        '''
        The constructor initializes the distance matrix as a dict with dicts so that
        we can have distances as self.dm[t1][t2] elements. If no distance has been estimated to
        or from tx, then self.dm[tx] is empty.

        The "distance_fcn" is a function which takes a count matrix N,
        where element N_{r,c} counts the number of amino pairs (r, c) in an MSA.
        '''
        self.msa = msa
        self.taxa = list(msa.taxa())
        self.n_taxa = len(self.taxa)
        self._maximal_distance = 5

        self._distance_function = distance_fcn

        self._dm = dict()
        self._internal_vertex_counter = 0
        self._last_progress = 0
        for t in self.taxa:
            self._dm[t] = dict()
        self._n_distances_computed = 0
        if verbose:
            self.verbose = verbose
        else:
            self.verbose = []


    def get(self, t1, t2):
        '''
        Retrieve the distance between taxa t1 and t2. If it is not in the matrix, then it
        is computed using the given model (see constructor).
        '''
        if t2 not in self._dm[t1]:
            d = self._compute_distance(t1, t2)
            self.set(t1, t2, d)
            self._n_distances_computed += 1
            return d
        else:
            return self._dm[t1][t2]


    def set(self, t1, t2, d):
        '''
        Manipulate the distance matrix. Primary usage is to set up for unit testing.
        '''
#        print(f'set d({t1}, {t2}) = {d}')
        self._dm[t1][t2] = d
        self._dm[t2][t1] = d


    def select_clade(self, x, l1, l2, l3):
        '''
        Return the leaf (l1, l2, or l3) of which leaf a quartet test suggests x belongs to.
        '''

        def q_diff(a, b, c, d):
            '''
            Return the perceived length _m_ of the mid branch separating ab with cd.

            a     c
             \___/
             /   \
            b     d
            '''

            m1 = self.get(a, c) + self.get(b, d)
            m2 = self.get(a, d) + self.get(b, c)
            m = m1 + m2 - 2*(self.get(a, b) + self.get(c, d))
            return m

        options = [(l1, q_diff(x, l1, l2, l3)),
                   (l2, q_diff(x, l2, l1, l3)),
                   (l3, q_diff(x, l3, l1, l2)),]
        choice = max(options, key = lambda pair: pair[1])
        #    print(f'{x}, choice: {choice} options: {options}', file=sys.stderr)
        return choice[0]


    def create_unique_vertex_id(self):
        '''
        Return a string represening a unique identifier for a vertex
        '''
        v_unique_id = f'#{self._internal_vertex_counter}'
        self._internal_vertex_counter += 1
        return v_unique_id


    def create_central_vertex(self, three_taxa, three_clades):
        '''
        Given three taxa (three_taxa is a tuple of three taxa), create
        a vertex v that represents the center of those three. The distance from
        v to all other vertices is estimated using the given information.

        three_taxa   A 3-tuple of chosen taxa
        three_clades A 3-tuple of the taxa for which we are building clades.
        '''
        v_unique_id = self.create_unique_vertex_id()
        self._dm[v_unique_id] = dict()
        for i in range(3):
            # Estimate distances from v to the selected taxa
            other0 = three_taxa[i]
            other1 = three_taxa[(i - 1) % 3]
            other2 = three_taxa[(i + 1) % 3]
            d_12 = self.get(other1, other2)

            # print(f'other0: {other0}\tother1: {other1}\tother2: {other2}\t')
            d_sum = 0           # \sum_{u \in \{other1, other2\}} \sum_{v\in clade(opposite of u)} d(u, other0) + d(v, other0) - d(u, v)
            n_terms = 0
            for j in [-1, +1]:
                u = three_taxa[(i+j) % 3]
                d_u_other0 = self.get(u, other0) # First term in sum (see comment for d_sum): d(u_other0)
                for v in three_clades[(i-j) % 3]:
                    d_sum += d_u_other0 + self.get(v, other0) - self.get(u, v)
                    n_terms += 2 # Adding distance of v_unique to other0 in the sum, so compensating for that
            d = d_sum / n_terms
            self.set(v_unique_id, other0, d)
            # print(f'\td(new_vertex, {other0}) = {d:.3}')

            # Estimate distance from all taxa in clade to v
            for t in three_clades[i]:
                if t == other0:
                    continue    # We have already set the distance for this taxa
                d = 0.5 * (self.get(other1, t) + self.get(other2, t) - d_12)
                self.set(v_unique_id, t, d)
                # print(f'\tt: {t}, d(new_vertex, {t}) = {d:.4} \t 0.5*d({other1}, {t})={self.get(other1, t):.3} + d({other2}, {t})={self.get(other2, t):.3} - d({other1}, {other2})={d_12:.3}')


        return v_unique_id

    def create_representative_vertex(self, taxa1, taxa2, other_taxa):
        '''
        Perform the central NeighborJoining operation of replacing
        taxa1 and taxa2 by a new vertex in the distance matrix
        and update distances to listed vertices accordingly.

        Note that taxa1 and taxa2 are not removed from the PartialDistanceMatrix.
        '''
        v_new = self.create_unique_vertex_id()
        self._dm[v_new] = dict()

        d_12 = self.get(taxa1, taxa2)
        n = len(other_taxa)
        for w in other_taxa:
            d = 0.5 * (self.get(taxa1, w) + self.get(taxa2, w) - d_12)
            self.set(v_new, w, d)
        return v_new

    def _compute_distance(self, t1, t2):
        '''
        Return the estimated distance between the sequences of the two taxa.
        '''
        if self.msa.can_retrieve_distances():
            d =  self.msa.distance(t1, t2)
#            print(f'dnctree: {t1}, {t2}: {d}')
            return d

        supress_warnings = 'supress_warnings' in self.verbose

        try:
            N = self.msa.count_pairs(t1, t2)
        except dnc.NoSharedCharactersError:
            if not supress_warnings:
                print(f'No shared characters between {t1} and {t2}, so cannot estimate their distance. Defaulting to maximal distance.', file=sys.stderr)
            return self._maximal_distance
        except dnc.AllCharactersDifferentError:
            if not supress_warnings:
                print(f'All characters shared by {t1} and {t2}, are different. Defaulting to maximal distance.', file=sys.stderr)
            return self._maximal_distance

        distance = self._distance_function(N)
        return distance


    def _computational_savings(self):
        actual_work, fraction_work, n = self.estimate_computational_savings()
        normal_work = (n * (n-1)) // 2 # Full distance matrix
        return f'Computed {actual_work} distances for {n} taxa. A full distance matrix would contain {normal_work} pairs. Savings: {100 - 100 * fraction_work:.3} %'

    def estimate_computational_savings(self):
        n = len(self.taxa)
        actual_work = self._n_distances_computed
        normal_work = (n * (n-1)) // 2 # Full distance matrix
        fraction_work = actual_work / normal_work
        return actual_work, fraction_work, n

    def print_progress(self):
        '''
        Based on how many internal nodes we have defined, we can estimate how much of the final tree we have reconstructed.
        This function writes a simple statement about problem completion based on that.
        '''
        if self._internal_vertex_counter > self._last_progress: # Avoid redundant progress information
            self._last_progress = self._internal_vertex_counter
            progress = 100 * self._internal_vertex_counter / (self.n_taxa - 2)
            print(f'Completion: {progress:.1f} %', file=sys.stderr, flush=True)
