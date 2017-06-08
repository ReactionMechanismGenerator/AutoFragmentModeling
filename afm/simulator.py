
# one iteration of the Monte Carlo simulation
# consists of the following steps
# s1: decide the time step to choose,
# need k, count of the related fragments
import rmgpy.constants

import afm.loader

class Simulator(object):

	def __init__(self, chemkin_path, dictionary_path):

		self.load_fragment_chemistry(chemkin_path, dictionary_path)

	def load_fragment_chemistry(self, chemkin_path, dictionary_path):

		fragment_dict, fragment_rxns = afm.loader.load_fragment_reactions_from_chemkin(chemkin_path,
                                        												dictionary_path)

		pseudo_fragrxns = afm.loader.load_pseudo_fragment_reactions(fragment_dict)

		self.fragment_reaction_list = fragment_rxns + pseudo_fragrxns
		self.fragment_dict = fragment_dict

class MonteCarloSimulator(Simulator):

	########################
	# Construction section #
	########################
	def __init__(self, 
				chemkin_path, 
				dictionary_path,
				input_fragment_count_dict,
				molecule_list,
				volume, 
				temperature):
		super(MonteCarloSimulator, self).__init__(chemkin_path, 
												dictionary_path)

		self.initialize_fragment_counts(input_fragment_count_dict)

		self.molecule_list = molecule_list
		self.volume = volume # unit: m^3
		self.temperature = temperature

	def initialize_fragment_counts(self, input_fragment_count_dict):

		# initialize the fragment_count
		self.fragment_count_dict = dict(zip(self.fragment_dict.keys(), 
											[0]*len(self.fragment_dict)))

		for frag in self.fragment_count_dict:
			if frag in input_fragment_count_dict:
				self.fragment_count_dict[frag] = input_fragment_count_dict[frag]

	

