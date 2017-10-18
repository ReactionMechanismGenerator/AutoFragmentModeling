import os
import unittest

from rmgpy.molecule.molecule import Atom, Bond
from rmgpy.molecule.element import getElement

import afm.fragment



class TestFragment(unittest.TestCase):

    def setUp(self):
        """
        A function run before each unit test in this class.
        """

        # construct the first fragment
        atom_C1 = Atom(element=getElement('C'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        cutting_label_R1 = afm.fragment.CuttingLabel('R')
        cutting_label_L1 = afm.fragment.CuttingLabel('L')
        
        vertices = [
            atom_C1,
            cutting_label_R1,
            cutting_label_L1
        ]

        bonds = [
            Bond(atom_C1, cutting_label_R1),
            Bond(atom_C1, cutting_label_L1)
        ]
        
        self.fragment1 = afm.fragment.Fragment()
        for vertex in vertices: self.fragment1.addVertex(vertex)
        for bond in bonds: self.fragment1.addEdge(bond)

        # construct the first fragment
        atom_C2 = Atom(element=getElement('C'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        cutting_label_R2 = afm.fragment.CuttingLabel('R')
        cutting_label_L2 = afm.fragment.CuttingLabel('L')
        
        vertices = [
            atom_C2,
            cutting_label_R2,
            cutting_label_L2
        ]

        bonds = [
            Bond(atom_C2, cutting_label_R2),
            Bond(atom_C2, cutting_label_L2)
        ]
        
        self.fragment2 = afm.fragment.Fragment()
        for vertex in vertices: self.fragment2.addVertex(vertex)
        for bond in bonds: self.fragment2.addEdge(bond)

    def test_fragment_isomorphism(self):

        self.assertTrue(self.fragment1.isIsomorphic(self.fragment2))
        