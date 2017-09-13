
# one iteration of the Monte Carlo simulation
# consists of the following steps
# s1: decide the time step to choose,
# need k, count of the related fragments
# s2: cast random number to decide which
# reaction rule to fire, need to create a
# dict which has <reaction, rate> pairs
# s3: realize the reaction rule on fragment level
# update fragment_count_dict
# s4: realize the reaction rule on molecule level
# need dict which has <fragment, [mol1, mol2, ...]> pairs
# update molecule list: some molecules to remove,
# some molecules to append.
import random
import numpy as np
import pandas as pd

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
				initial_molecules,
				volume, 
				temperature):
		super(MonteCarloSimulator, self).__init__(chemkin_path, 
												dictionary_path)

		self.initialize_fragment_counts(input_fragment_count_dict)

		self.initialize_molecule_fragment_df(initial_molecules)
		self.volume = volume # unit: m^3
		self.temperature = temperature
		self.reaction_flux_array = np.zeros(len(self.fragment_reaction_list))

	def initialize_fragment_counts(self, input_fragment_count_dict):

		# initialize the fragment_count
		self.fragment_count_dict = dict(zip(self.fragment_dict.keys(), 
											[0]*len(self.fragment_dict)))

		for frag in self.fragment_count_dict:
			if frag in input_fragment_count_dict:
				self.fragment_count_dict[frag] = input_fragment_count_dict[frag]

	def initialize_molecule_fragment_df(self, initial_molecules):

		self.molecule_fragment_df = pd.DataFrame(columns=self.fragment_count_dict.keys())

		for mol in initial_molecules:
			insert_row = {}
			for fragment_label  in self.fragment_count_dict:
				if fragment_label in mol:
					insert_row[fragment_label] = mol[fragment_label]
				else:
					insert_row[fragment_label] = 0

			if not self.molecule_fragment_df.index.empty:
				insert_index = 1 + max(self.molecule_fragment_df.index)
			else:
				insert_index = 0
			self.molecule_fragment_df.loc[insert_index] = insert_row


	######################
	# Simulation section #
	######################
	def calculate_unimolucular_rate(self, reaction):

		k_u = reaction.kinetics.getRateCoefficient(self.temperature)
		frag_label = reaction.reactants[0].label
		frag_count = self.fragment_count_dict[frag_label]
		rate_u = k_u * frag_count # unit: 1/s

		return rate_u

	def calculate_bimolucular_rate(self, reaction):

		Na = rmgpy.constants.Na
		k_b = reaction.kinetics.getRateCoefficient(self.temperature)/Na # unit: m^3/s
		frag_label1 = reaction.reactants[0].label
		frag_label2 = reaction.reactants[1].label

		frag_count1 = self.fragment_count_dict[frag_label1]
		frag_count2 = self.fragment_count_dict[frag_label2]
		rate_b = k_b * frag_count1 * frag_count2 / self.volume # unit: 1/s

		return rate_b

	def update_reaction_fluxes(self):

		for idx, frag_rxn in enumerate(self.fragment_reaction_list):
			# unimolecular reactions
			if len(frag_rxn.reactants) == 1:
				rate = self.calculate_unimolucular_rate(frag_rxn)

			# bimolecular reactions
			elif len(frag_rxn.reactants) == 2:
				rate = self.calculate_bimolucular_rate(frag_rxn)

			self.reaction_flux_array[idx] = rate

	def time_step(self):

		total_rate = np.sum(self.reaction_flux_array)

		return 1.0/total_rate # unit: s


	def random_reaction_selection(self, seed=None):

		# cast a random number
		random.seed(seed)
		rand_num = random.uniform(0, np.sum(self.reaction_flux_array))

		# pick the one fragment reaction
		# indicated by the random number
		current_fluxsum = 0
		for idx, flux in enumerate(self.reaction_flux_array):

			current_fluxsum += flux
			print current_fluxsum
			if current_fluxsum >= rand_num:
				return idx

	def update_fragment_counts(self, reaction_idx):

		reaction = self.fragment_reaction_list[reaction_idx]
		for reactant in reaction.reactants:
			self.fragment_count_dict[reactant.label] -= 1

		for product in reaction.products:
			self.fragment_count_dict[product.label] += 1

