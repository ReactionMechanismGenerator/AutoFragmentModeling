
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



