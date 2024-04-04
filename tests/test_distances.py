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


from dnctree.algtesting import TestingFullDistanceMatrix
class TestDistMatrixReading(unittest.TestCase):
    '''
    This class tests the class TestingFullDistanceMatrix,
    which is only used for algorithm evaluation when
    there is a distance matrix available (e.g.,
    using fastprot). 
    '''

    def test_read_matrix(self):
        '''
        Read matrix from file and verify it.
        '''
        dm = TestingFullDistanceMatrix('five.tree', 'five.distmatrix')
        self.assertAlmostEqual(dm.get('Alpha', 'Alpha'), 0.0,    4)

        self.assertAlmostEqual(dm.get('Alpha', 'Beta'),  0.2997, 4)
        self.assertAlmostEqual(dm.get('Alpha', 'Gamma'), 0.7820, 4)
        self.assertAlmostEqual(dm.get('Alpha', 'Delta'), 1.1716, 4)
        self.assertAlmostEqual(dm.get('Alpha', 'Epsilon'), 1.4617, 4)

        self.assertAlmostEqual(dm.get('Beta', 'Alpha'),  0.2997, 4)
        self.assertAlmostEqual(dm.get('Gamma', 'Alpha'), 0.7820, 4)
        self.assertAlmostEqual(dm.get('Delta', 'Alpha'), 1.1716, 4)
        self.assertAlmostEqual(dm.get('Epsilon', 'Alpha'), 1.4617, 4)        

        self.assertAlmostEqual(dm.get('Epsilon', 'Beta'), 0.5653, 4)        
        self.assertAlmostEqual(dm.get('Beta', 'Epsilon'), 0.5653, 4)        
        
        
        
    
