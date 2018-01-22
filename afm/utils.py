

def match_sequences(seq1, seq2, diff_tol=1e-6):

	"""
	Given two lists (each item is int or float):
	seq1 and seq2 with same sum, the method returns 
	matched indices and values.

	Example:

	seq1 = [1, 3, 1]
	seq2 = [2, 1, 2]

	return: [[(0,0),1], 
			 [(1,0),1], 
			 [(1,1),1],
			 [(1,2),1],
			 [(2,2),1]]
	"""

	# check if sums are close to same
	sum_diff = sum(seq2)-sum(seq1)
	assert abs(sum_diff/1.0/sum(seq1)) <= diff_tol, 'seq1 has different sum (diff={0}) than seq2.'.format(sum_diff)
	
	# force the sum to be same if the difference
	# is small enough
	if sum_diff >=0:
		seq1[-1] = seq1[-1] + sum_diff
	else:
		seq2[-1] = seq2[-1] - sum_diff
	
	# make cumulative sequences
	cum_seq1 = [seq1[0]]
	for item1 in seq1[1:]:
		cum_seq1.append(cum_seq1[-1] + item1)

	cum_seq2 = [seq2[0]]
	for item2 in seq2[1:]:
		cum_seq2.append(cum_seq2[-1] + item2)

	# add index tags two both cumulative seqs
	pin1 = 0
	pin2 = 0
	matched_indices = []
	matched_cum_values = []
	while pin1 < len(cum_seq1) and pin2 < len(cum_seq2):

		matched_indices.append((pin1, pin2))
		
		if cum_seq1[pin1] > cum_seq2[pin2]:
			matched_cum_values.append(cum_seq2[pin2])
			pin2 += 1
		elif cum_seq1[pin1] < cum_seq2[pin2]:
			matched_cum_values.append(cum_seq1[pin1])
			pin1 += 1
		else:
			matched_cum_values.append(cum_seq2[pin2])
			pin1 += 1
			pin2 += 1

	# get matches
	matches = []
	for i in range(len(matched_indices)):
		matched_index_tup = matched_indices[i]
		matched_cum_value = matched_cum_values[i]
		if i == 0:
			previous_cum_value = 0
		else:
			previous_cum_value = matched_cum_values[i-1]

		matches.append([matched_index_tup, matched_cum_value - previous_cum_value])

	return matches