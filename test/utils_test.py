import unittest

import afm.utils

class TestUtils(unittest.TestCase):

	def test_match_sequences1(self):

		seq1 = [1, 3, 1]
		seq2 = [2, 1, 2]

		matches = afm.utils.match_sequences(seq1, seq2)

		expected_matches = [[(0,0),1], 
							[(1,0),1], 
							[(1,1),1],
							[(1,2),1],
							[(2,2),1]]

		self.assertEqual(matches, expected_matches)


	def test_match_sequences2(self):
		"""
		Test match_sequences() can tolerate slight
		sum difference between two sequences, default tolerance
		is 1e-6.
		"""

		seq1 = [1, 3, 1-1e-6]
		seq2 = [2, 1, 2]

		matches = afm.utils.match_sequences(seq1, seq2)

		expected_matches = [[(0,0),1], 
							[(1,0),1], 
							[(1,1),1],
							[(1,2),1],
							[(2,2),1]]

		self.assertEqual(matches, expected_matches)

	def test_match_concentrations(self):

		conc1 = [('a', 1),
				('b', 3),
				('c', 1)]

		conc2 = [('x', 2),
				('y', 1),
				('z', 2)]

		matches = afm.utils.match_concentrations(conc1, conc2)

		expected_matches = [[('a','x'),1], 
							 [('b','x'),1], 
							 [('b','y'),1],
							 [('b','z'),1],
							 [('c','z'),1]]

		self.assertEqual(matches, expected_matches)

	def test_grind(self):

		conc = [('a', 1),
				('b', 3),
				('c', 1)]
		size = 0.6

		grinded_conc = afm.utils.grind(conc, size)
		expected_grinded_conc = [('a', 0.6),
								 ('a', 0.4),
								 ('b', 0.6),
								 ('b', 0.6),
								 ('b', 0.6),
								 ('b', 0.6),
								 ('b', 0.6),
								 ('c', 0.6),
								 ('c', 0.4)]
		self.assertEqual(grinded_conc, expected_grinded_conc)


