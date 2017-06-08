import os
import unittest

import afm.simulator

class TestSimulator(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        chemkin_path = os.path.join(os.path.dirname(__file__), 
                                    'data', 
                                    'simulator_data',
                                    'chem.inp')

        dictionary_path = os.path.join(os.path.dirname(__file__), 
                                    'data', 
                                    'simulator_data',
                                    'species_dictionary.txt')

        self.simulator = afm.simulator.Simulator(chemkin_path, 
                                                dictionary_path)

    def test_fragment_chemistry(self):

        fragment_dict = self.simulator.fragment_dict
        fragment_reaction_list = self.simulator.fragment_reaction_list
        self.assertEqual(40, len(fragment_dict))
        self.assertEqual(5, len(fragment_reaction_list))

