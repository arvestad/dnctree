import unittest

from dnctree.taxaset import TaxaSet

class TestTaxaSet(unittest.TestCase):
    def test_creation(self):
        ts = TaxaSet([])
        self.assertEqual(len(list(ts)), 0)

        ts = TaxaSet(['a'])
        self.assertEqual(len(list(ts)), 1)

        ts = TaxaSet(['a', 'b'])
        self.assertEqual(len(list(ts)), 2)

    def test_addition(self):
        ts = TaxaSet(['a', 'b'])
        ts.add_created_vertex('#internal1')
        self.assertEqual(len(list(ts)), 3)

        ts.add_created_vertex('#internal2')
        self.assertEqual(len(list(ts)), 4)
