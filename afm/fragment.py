
from rmgpy.species import Species
from rmgpy.molecule.element import getElement
from rmgpy.molecule.graph import Graph, Vertex
from rmgpy.molecule.molecule import Atom, Bond, Molecule

class CuttingLabel(Vertex):

	def __init__(self, label=''):
		Vertex.__init__(self)
		self.label = label

	def __str__(self):
		"""
		Return a human-readable string representation of the object.
		"""
		return str(self.label)

	def __repr__(self):
		"""
		Return a representation that can be used to reconstruct the object.
		"""
		return "<CuttingLabel '{0}'>".format(str(self))

	def isSpecificCaseOf(self, other):
		"""
		Return ``True`` if `self` is a specific case of `other`, or ``False``
		otherwise. At this moment, this is the same as the :meth:`equivalent()`.
		"""
		return self.equivalent(other)

	def equivalent(self, other):
		"""
		Return ``True`` if `other` is indistinguishable from this CuttingLabel, or
		``False`` otherwise. If `other` is an :class:`CuttingLabel` object, then all
		attributes must match exactly. 
		"""
		if isinstance(other, CuttingLabel):
			return self.label == other.label
		else:
			return False

	def copy(self):
		"""
		Generate a deep copy of the current CuttingLabel. Modifying the
		attributes of the copy will not affect the original.
		"""

		c = CuttingLabel.__new__(CuttingLabel)
		c.edges = {}
		c.resetConnectivityValues()
		c.label = self.label
		return c

class Fragment(Graph):

	def __init__(self,
				label='',
				species_repr=None,
				vertices=None):
		self.index = -1
		self.label = label
		self.species_repr = species_repr
		Graph.__init__(self, vertices)

	def __str__(self):
		"""
		Return a string representation, in the form 'label(id)'.
		"""
		if self.index == -1: return self.label
		else: return '{0}({1:d})'.format(self.label, self.index)
	
	## need to rewrite later
	def _repr_png_(self):
		if len(self.species_repr.molecule) > 0:
			return self.species_repr.molecule[0]._repr_png_()
		else:
			return None

	def from_SMILES_like_string(self, SMILES_like_string):

		smiles_replace_dict = {
						'R': 'Cl', 
						'L': '[Si]'
		}

		atom_replace_dict = {
						'Cl': 'R', 
						'Si': 'L'
		}
		
		smiles =  SMILES_like_string
		for label_str, element in smiles_replace_dict.iteritems():
			smiles = smiles.replace(label_str, element)

		mol = Molecule().fromSMILES(smiles)
		# convert mol to fragment object
		final_vertices = []
		for atom in mol.atoms:
			element_symbol = atom.element.symbol
			if element_symbol in atom_replace_dict:
				
				label = atom_replace_dict[element_symbol]

				cutting_label = CuttingLabel(label=label)
				for bondedAtom, bond in atom.edges.iteritems():
					new_bond = Bond(bondedAtom, cutting_label, order=bond.order)
					
					bondedAtom.edges[cutting_label] = new_bond
					del bondedAtom.edges[atom]

					cutting_label.edges[bondedAtom] = new_bond
				final_vertices.append(cutting_label)
			
			else:
				final_vertices.append(atom)

		self.vertices = final_vertices
		return self

	def assign_representative_species(self):

		# create a molecule from fragment.vertices.copy
		fragment_copy = self.copy(deep=True)

		# replace CuttingLabel with CC
		atoms = []
		additional_atoms = []
		additional_bonds = []
		for vertex in fragment_copy.vertices:
			if isinstance(vertex, CuttingLabel):

				# replace cutting label with atom C
				atom_C1 = Atom(element=getElement('C'), 
							radicalElectrons=0, 
							charge=0, 
							lonePairs=0)

				for bondedAtom, bond in vertex.edges.iteritems():
					new_bond = Bond(bondedAtom, atom_C1, order=bond.order)
					
					bondedAtom.edges[atom_C1] = new_bond
					del bondedAtom.edges[vertex]

					atom_C1.edges[bondedAtom] = new_bond

				# add hydrogens and carbon to make it CC
				atom_H1 = Atom(element=getElement('H'), 
							radicalElectrons=0, 
							charge=0, 
							lonePairs=0)

				atom_H2 = Atom(element=getElement('H'), 
							radicalElectrons=0, 
							charge=0, 
							lonePairs=0)

				atom_C2 = Atom(element=getElement('C'), 
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

				atom_H5 = Atom(element=getElement('H'), 
							radicalElectrons=0, 
							charge=0, 
							lonePairs=0)

				atoms.append(atom_C1)

				additional_atoms.extend([atom_H1,
									atom_H2,
									atom_H3,
									atom_H4,
									atom_H5,
									atom_C2])

				additional_bonds.extend([Bond(atom_C1, atom_H1, 1),
									Bond(atom_C1, atom_H2, 1),
									Bond(atom_C2, atom_H3, 1),
									Bond(atom_C2, atom_H4, 1),
									Bond(atom_C2, atom_H5, 1),
									Bond(atom_C1, atom_C2, 1)])

			else:
				# assert isinstance(vertex, Atom), type(vertex)
				atoms.append(vertex)

		mol_repr = Molecule()
		mol_repr.atoms = atoms
		for atom in additional_atoms: mol_repr.addAtom(atom)
		for bond in additional_bonds: mol_repr.addBond(bond)
		# update connectivity
		mol_repr.update()

		# create a species object from molecule
		self.species_repr = Species(molecule=[mol_repr])

	def get_fragmental_weight(self):
		"""
		Return the fragmental weight of the fragment in kg/mol.
		"""
		mass = 0
		for vertex in self.vertices:
			if isinstance(vertex, Atom):
				mass += vertex.element.mass
		return mass

	def is_radical(self):
		"""
		Return ``True`` if the fragment contains at least one radical electron,
		or ``False`` otherwise.
		"""
		for vertex in self.vertices:
			if isinstance(vertex, Atom) and vertex.radicalElectrons > 0:
				return True
		return False