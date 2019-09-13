import math
import sys

def poisson_distance(t1, t2, s1, s2, supress_warnings=False):
    '''
    Default method for distance estimation in PartialDistanceMatrix.get().

    A better alternative, with the same signature, should be provided
    when instantiating PartialDistanceMatrix.
    '''
    n_chars = 0
    n_diffs = 0
    for c1, c2 in zip(s1, s2):
        if c1=='-' or c2=='-':
            continue
        if c1 != c2:
            n_diffs += 1
        n_chars += 1

    if n_chars == 0:
        if not supress_warnings:
            print(f'Warning: No shared characters between {t1} and {t2}, so cannot estimate their distance. Defaulting to distance=2.5.', file=sys.stderr)
        return 2.5
    if n_diffs == n_chars:
        if not supress_warnings:
            print(f'Warning: Degenerate sequence pair: {t1} and {t2}. All shared characters are different. Defaulting to distance=2.5.', file=sys.stderr)
        return 2.5

    distance = - math.log(1 - n_diffs/n_chars)
    return distance


class PartialDistanceMatrix:
    '''
    Represent distance matrices using a dictionary instead of a full (half) square matrix,
    because we are not going to compute all pairwise distances.

    This class pretends it is a distance matrix. Distances are computed on-demand, when
    they are asked for.
    '''
    def __init__(self, msa, verbose=False, distance_estimator=poisson_distance):
        '''
        The constructor initializes the distance matrix as a dict with dicts so that
        we can have distances as self.dm[t1][t2] elements. If no distance has been estimated to
        or from tx, then self.dm[tx] is empty.

        The "distance_estimator" is a function.
        '''
        self.msa = msa
        self.taxa = list(msa.taxa())
        self.n_taxa = len(self.taxa)

        self._distance_function = distance_estimator

        self._dm = dict()
        self._internal_vertex_counter = 0
        self._last_progress = 0
        for t in msa.taxa():
            self._dm[t] = dict()
        self._n_distances_computed = 0
        self.verbose = verbose

    def get(self, t1, t2):
        '''
        Retrieve the distanve between taxa t1 and t2. If it is not in the matrix, then it
        is computed using the given model (see constructor).
        '''
        if t2 not in self._dm[t1]:
            d = self._compute_distance(t1, t2)
            self.set(t1, t2, d)
            return d
        else:
            return self._dm[t1][t2]


    def set(self, t1, t2, d):
        self._dm[t1][t2] = d
        self._dm[t2][t1] = d


    def taxa(self):
        return self.msa.taxa()


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
            d = (self.get(other0, other1) + self.get(other0, other2)
                 - self.get(other1, other2)) / 2
            self.set(v_unique_id, other0, d)

            # Estimate distance from all taxa in clade to v
            for t in three_clades[i]:
                d = 0.5 * (self.get(other1, t) + self.get(other2, t))
                self.set(v_unique_id, t, d)

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
        As a first rough attempt, we just use Poisson distance.
        '''
        s1 = self.msa[t1].seq
        s2 = self.msa[t2].seq

        supress_warnings = 'supress_warnings'  in self.verbose
        distance = self._distance_function(t1, t2, s1, s2, supress_warnings)
        self._n_distances_computed += 1
        return distance


    def _print_computational_savings(self):
        n = len(self.taxa)
        normal_work = int(n * (n-1) / 2) # Full distance matrix
        actual_work = self._n_distances_computed
        print(f'[Computed {actual_work} distances for {n} taxa. A full distance matrix would contain {normal_work} pairs. Savings: {100 - 100 * actual_work/normal_work:.3} %]', file=sys.stderr)

    def print_progress(self):
        '''
        Based on how many internal nodes we have defined, we can estimate how much of the final tree we have reconstructed.
        This function writes a simple statement about problem completion based on that.
        '''
        if self._internal_vertex_counter > self._last_progress: # Avoid redundant progress information
            self._last_progress = self._internal_vertex_counter
            progress = 100 * self._internal_vertex_counter / (self.n_taxa - 2)
            print(f'Completion: {progress:.1f} %', file=sys.stderr, flush=True)
