import os
import unittest

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase

import afm.react
import afm.fragment

class TestReact(unittest.TestCase):

    def setUp(self):
        """
        A function run before each unit test in this class.
        """

        # load kinetics database
        db_path = settings['database.directory']
        self.database = RMGDatabase()

        # forbidden structure loading
        self.database.loadForbiddenStructures(os.path.join(db_path, 'forbiddenStructures.py'))
        # kinetics family loading
        self.database.loadKinetics(os.path.join(db_path, 'kinetics'),
                                   kineticsFamilies='default',
                                   reactionLibraries=[]
                                   )

    def test_react_fragments1(self):
        
        frag1 = afm.fragment.Fragment(label='frag1').from_SMILES_like_string('c1ccccc1CCCR')
        fragment_tuple = (frag1, )
        reactions = afm.react.react_fragments(self.database.kinetics, 
                                              fragment_tuple,
                                              only_families=['R_Recombination'],
                                              prod_resonance=False)

        self.assertTrue(len(reactions)==28) 

    def test_react_fragments2(self):
        
        frag1 = afm.fragment.Fragment(label='frag1').from_SMILES_like_string('c1ccccc1CCCR')
        frag2 = afm.fragment.Fragment(label='frag2').from_SMILES_like_string('[CH2]CR')
        fragment_tuple = (frag1, frag2)
        reactions = afm.react.react_fragments(self.database.kinetics, 
                                              fragment_tuple,
                                              only_families=['H_Abstraction'],
                                              prod_resonance=False)

        self.assertTrue(len(reactions)==11) 