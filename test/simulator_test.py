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

class TestMonteCarloSimulator(unittest.TestCase):

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
        input_fragment_count_dict = {
          "ArCCCCR": 1000,
          "ArC__C":10
        }

        volume = 4.8e-25 # unit: m^3
        temperature = 700 # unit: K
        self.mcs = afm.simulator.MonteCarloSimulator(chemkin_path, 
                                                     dictionary_path,
                                                     input_fragment_count_dict,
                                                     [],
                                                     volume, 
                                                     temperature)

    def test_calculate_unimolucular_rate(self):

        # select 157th reaction: ArCCCCR => ArC* + RCCC*
        uni_reaction = self.mcs.fragment_reaction_list[2]
        rate_u = self.mcs.calculate_unimolucular_rate(uni_reaction)

        expected_rate_u = 1.51e-4
        self.assertAlmostEqual(expected_rate_u, rate_u, 6)

    def test_calculate_bimolucular_rate(self):

        # select 182th reaction: ArC__C + ArCCCCR => ArC*CCCR + ArC*C
        bi_reaction = self.mcs.fragment_reaction_list[3]
        rate_b = self.mcs.calculate_bimolucular_rate(bi_reaction)

        expected_rate_b = 1.43e-3
        self.assertAlmostEqual(expected_rate_b, rate_b, 5)

    def test_time_step(self):

        time_step = self.mcs.time_step()

        expected_time_step = 633.0

        self.assertAlmostEqual(time_step, expected_time_step, 1)