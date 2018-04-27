import os
import unittest

from rmgpy.species import Species
from rmgpy.molecule.molecule import Atom, Bond, Molecule
from rmgpy.molecule.element import getElement
from rmgpy.molecule.atomtype import atomTypes

import afm.fragment

class TestCuttingLabel(unittest.TestCase):

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.cutting_label_R = afm.fragment.CuttingLabel('R')

    def test_symbol(self):

        self.assertEqual('R', self.cutting_label_R.symbol)

    def test_copy(self):

        cutting_label_R_copy = self.cutting_label_R.copy()

        self.assertEqual('R', cutting_label_R_copy.name)
        self.assertEqual(self.cutting_label_R.label, 
                         cutting_label_R_copy.label)
        self.assertEqual(self.cutting_label_R.charge, 
                         cutting_label_R_copy.charge)
        self.assertEqual(self.cutting_label_R.radicalElectrons, 
                         cutting_label_R_copy.radicalElectrons)
        self.assertEqual(self.cutting_label_R.lonePairs, 
                         cutting_label_R_copy.lonePairs)
        self.assertEqual(self.cutting_label_R.isotope, 
                         cutting_label_R_copy.isotope)

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

    def test_from_SMILES_like_string3(self):

        # generate fragment from SMILES like string
        # the atom type is also calculated
        smiles_like = 'RCL'
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

        cutting_label_R = afm.fragment.CuttingLabel('R')
        cutting_label_L = afm.fragment.CuttingLabel('L')
        
        vertices = [
            atom_C,
            cutting_label_R,
            cutting_label_L,
            atom_H1,
            atom_H2
        ]

        bonds = [
            Bond(atom_C, cutting_label_R, 1),
            Bond(atom_C, cutting_label_L, 1),
            Bond(atom_C, atom_H1, 1),
            Bond(atom_C, atom_H2, 1)
        ]
        
        expected_fragment = afm.fragment.Fragment()
        for vertex in vertices: expected_fragment.addVertex(vertex)
        for bond in bonds: expected_fragment.addEdge(bond)

        self.assertTrue(expected_fragment.isIsomorphic(fragment))

    def test_isSubgraphIsomorphic1(self):

        from rmgpy.molecule.group import Group

        smiles_like = '[CH2]CR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        adj = """
1 * R u1
"""
        other = Group().fromAdjacencyList(adj)

        self.assertTrue(fragment.isSubgraphIsomorphic(other))

    def test_isSubgraphIsomorphic2(self):

        from rmgpy.molecule.group import Group

        smiles_like = '[CH2]CR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        adj = """
1 * Ct  u1 {2,T}
2   R!H u0 {1,T}
"""
        other = Group().fromAdjacencyList(adj)

        self.assertFalse(fragment.isSubgraphIsomorphic(other))

    def test_isSubgraphIsomorphic3(self):

        from rmgpy.molecule.group import Group

        smiles_like = '[CH2]CR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        adj = """
1 * R u1
"""
        other = Group().fromAdjacencyList(adj)

        # create initial map
        frag_atom_star = None
        for frag_atom in fragment.vertices:
            if isinstance(frag_atom, Atom) and frag_atom.radicalElectrons == 1:
                frag_atom_star = frag_atom
                break

        group_atom_star = other.vertices[0]
        initialMap = {frag_atom_star: group_atom_star}

        self.assertTrue(fragment.isSubgraphIsomorphic(other, initialMap=initialMap))

    def test_isSubgraphIsomorphic4(self):

        from rmgpy.molecule.group import Group

        smiles_like = '[CH2]CR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        fragment.assign_representative_molecule()

        adj = """
1 * Cs u1 {2,S}
2   N u0 {1,S}
"""
        other = Group().fromAdjacencyList(adj)

        # create initial map
        frag_atom_star = None
        for frag_atom in fragment.vertices:
            if isinstance(frag_atom, Atom) and frag_atom.radicalElectrons == 1:
                frag_atom_star = frag_atom
                break

        group_atom_star = other.vertices[0]
        initialMap = {frag_atom_star: group_atom_star}

        self.assertFalse(fragment.isSubgraphIsomorphic(other, initialMap=initialMap))

    def test_assign_representative_species(self):

        smiles_like = 'RCR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        fragment.assign_representative_species()

        expected_repr_spec = Species().fromSMILES('CCCCC')

        self.assertTrue(expected_repr_spec.isIsomorphic(fragment.species_repr))

    def test_assign_representative_molecule(self):

        smiles_like = 'RCR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)

        fragment.assign_representative_molecule()

        expected_repr_mol = Molecule().fromSMILES('CCCCC')

        self.assertTrue(expected_repr_mol.isIsomorphic(fragment.mol_repr))

    def test_get_fragmental_weight1(self):

        fragmental_weight = self.fragment1.get_fragmental_weight()
        self.assertAlmostEqual(fragmental_weight*1000, 12.01, 2)

    def test_get_fragmental_weight2(self):

        smiles_like = 'RCR'
        fragment = afm.fragment.Fragment().from_SMILES_like_string(smiles_like)
        fragmental_weight = fragment.get_fragmental_weight()
        self.assertAlmostEqual(fragmental_weight*1000, 14.03, 2)

    def test_updateAtomTypes(self):

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
        
        fragment = afm.fragment.Fragment()
        for vertex in vertices: fragment.addVertex(vertex)
        for bond in bonds: fragment.addEdge(bond)

        fragment.updateAtomTypes()

        for v in fragment.vertices:
            if isinstance(v, Atom) and v.isCarbon():
                break

        self.assertTrue(v.atomType == atomTypes['Cs'])

    def test_update(self):

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
        
        fragment = afm.fragment.Fragment()
        for vertex in vertices: fragment.addVertex(vertex)
        for bond in bonds: fragment.addEdge(bond)

        fragment.update()

        for v in fragment.vertices:
            if isinstance(v, Atom) and v.isCarbon():
                break

        self.assertTrue(v.atomType == atomTypes['Cs'])
        self.assertTrue(fragment.getNetCharge() == 0)
        self.assertTrue(fragment.multiplicity == 1)
