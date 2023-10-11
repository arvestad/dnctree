import unittest

from modelmatcher.models import RateMatrix
import dnctree.distances as D
import numpy as np

class TestDistances(unittest.TestCase):

    def setUp(self):
        N = np.zeros((20,20))
        for i in range(20):
            N[i,i] = 2
        self.N_t_is_zero = N

    def test_kimura_distance(self):
        d = D.kimura_distance(self.N_t_is_zero)
        self.assertTrue(0 <= d)
        self.assertTrue(d < D._tolerance)

    def test_ml_distance(self):
        wag = RateMatrix.instantiate('WAG')
        d = D.ml_distance_estimate(wag, self.N_t_is_zero)
        self.assertTrue(0 <= d)
        self.assertTrue(d < D._tolerance)

        d = D.ml_distance_estimate(wag, self.N_t_is_zero, 1.0)
        self.assertTrue(0 <= d)
        self.assertTrue(d < D._tolerance)

        for t in np.arange(0.1,2.0,0.1):
            N = wag.sample_count_matrix(10000, [t])
            d = D.ml_distance_estimate(wag, N)
#            print(f'd = {d}, t={t}, abs(t-d) = {abs(t-d)}')
            self.assertTrue(abs(t - d) < 0.1) # Plausible mistake
