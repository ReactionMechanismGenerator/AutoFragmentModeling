
def match_concentrations(conc1, conc2, diff_tol=1e-6):
	"""
	Given two lists with each item to be a tuple
	(species label, concentration)
	conc1 and conc2 with same total concentrations, 
	the method returns matched species labels and 
	concentrations.

	Example:

	conc1 = [('a', 1),
			('b', 3),
			('c', 1)]
	conc2 = [('x', 2),
			('y', 1),
			('z', 2)]

	return: [[('a','x'),1], 
			 [('b','x'),1], 
			 [('b','y'),1],
			 [('b','z'),1],
			 [('c','z'),1]]
	"""

	labels1 = [tup[0] for tup in conc1]
	labels2 = [tup[0] for tup in conc2]

	seq1 = [tup[1] for tup in conc1]
	seq2 = [tup[1] for tup in conc2]

	matches_seq = match_sequences(seq1, seq2, diff_tol)

	matches_conc = []
	for match_seq in matches_seq:
		matched_label_index1 = match_seq[0][0]
		matched_label_index2 = match_seq[0][1]
		matched_value = match_seq[1]

		matched_label1 = labels1[matched_label_index1]
		matched_label2 = labels2[matched_label_index2]
		match_conc = [(matched_label1, matched_label2), matched_value]
		matches_conc.append(match_conc)

	return matches_conc

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