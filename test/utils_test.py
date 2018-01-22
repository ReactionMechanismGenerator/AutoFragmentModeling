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
		is 10e-6.
		"""

		seq1 = [1, 3, 1-10e-6]
		seq2 = [2, 1, 2]

		matches = afm.utils.match_sequences(seq1, seq2)

		expected_matches = [[(0,0),1], 
							[(1,0),1], 
							[(1,1),1],
							[(1,2),1],
							[(2,2),1]]

		self.assertEqual(matches, expected_matches)

