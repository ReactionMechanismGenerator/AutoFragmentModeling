import os
import unittest

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase

import afm.react
import afm.fragment

TESTFAMILY = 'R_Recombination'

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
                                   kineticsFamilies=[TESTFAMILY],
                                   reactionLibraries=[]
                                   )

    def test_react_fragments1(self):
        
        frag1 = afm.fragment.Fragment(label='frag1').from_SMILES_like_string('c1ccccc1CCCR')
        fragment_tuple = (frag1, )
        reactions = afm.react.react_fragments(self.database.kinetics, 
                                              fragment_tuple,
                                              only_families=[TESTFAMILY],
                                              prod_resonance=False)

        self.assertTrue(len(reactions)==28) 