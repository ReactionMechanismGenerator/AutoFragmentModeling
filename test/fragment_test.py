import os
import unittest

from rmgpy.species import Species
from rmgpy.molecule.molecule import Atom, Bond
from rmgpy.molecule.element import getElement
from rmgpy.molecule.atomtype import atomTypes

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

        # construct the second fragment
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

    def test_from_SMILES_like_string1(self):

        # generate fragment from SMILES like string
        # the atom type is also calculated
        smiles_like = 'C'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)
        
        # construct fragment manually
        atom_C = Atom(element=getElement('C'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H1 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H2 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H3 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H4 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_C.atomType=atomTypes['Cs']
        atom_H1.atomType=atomTypes['H']
        atom_H2.atomType=atomTypes['H']
        atom_H3.atomType=atomTypes['H']
        atom_H4.atomType=atomTypes['H']

        vertices = [
            atom_C,
            atom_H1,
            atom_H2,
            atom_H3,
            atom_H4
        ]

        bonds = [
            Bond(atom_C, atom_H1, 1),
            Bond(atom_C, atom_H2, 1),
            Bond(atom_C, atom_H3, 1),
            Bond(atom_C, atom_H4, 1)
        ]
        
        expected_fragment = afm.fragment.Fragment()
        for vertex in vertices: expected_fragment.addVertex(vertex)
        for bond in bonds: expected_fragment.addEdge(bond)

        self.assertTrue(expected_fragment.isIsomorphic(fragment))

    def test_from_SMILES_like_string2(self):

        # generate fragment from SMILES like string
        # the atom type is also calculated
        smiles_like = 'RCR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        atom_C = Atom(element=getElement('C'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H1 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        atom_H2 = Atom(element=getElement('H'), 
                    radicalElectrons=0, 
                    charge=0, 
                    lonePairs=0)

        # construct fragment manually
        atom_C.atomType=atomTypes['Cs']
        atom_H1.atomType=atomTypes['H']
        atom_H2.atomType=atomTypes['H']

        cutting_label_R1 = afm.fragment.CuttingLabel('R')
        cutting_label_R2 = afm.fragment.CuttingLabel('R')
        
        vertices = [
            atom_C,
            cutting_label_R1,
            cutting_label_R2,
            atom_H1,
            atom_H2
        ]

        bonds = [
            Bond(atom_C, cutting_label_R1, 1),
            Bond(atom_C, cutting_label_R2, 1),
            Bond(atom_C, atom_H1, 1),
            Bond(atom_C, atom_H2, 1)
        ]
        
        expected_fragment = afm.fragment.Fragment()
        for vertex in vertices: expected_fragment.addVertex(vertex)
        for bond in bonds: expected_fragment.addEdge(bond)

        self.assertTrue(expected_fragment.isIsomorphic(fragment))

    def test_assign_representative_species(self):

        smiles_like = 'RCR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        fragment.assign_representative_species()

        expected_repr_spec = Species().fromSMILES('CCCCC')

        self.assertTrue(expected_repr_spec.isIsomorphic(fragment.species_repr))
        