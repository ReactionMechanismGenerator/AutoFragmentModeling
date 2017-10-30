
from rmgpy.molecule.graph import Graph, Vertex
from rmgpy.molecule.molecule import Bond, Molecule

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

		smiles =  SMILES_like_string.replace('R', 'Cl')

		mol = Molecule().fromSMILES(smiles)
		# convert mol to fragment object
		final_vertices = []
		for atom in mol.atoms:
			if atom.element.symbol == 'Cl':
				
				R = CuttingLabel('R')
				for bondedAtom, bond in atom.edges.iteritems():
					new_bond = Bond(bondedAtom, R, order=bond.order)
					
					bondedAtom.edges[R] = new_bond
					del bondedAtom.edges[atom]

					R.edges[bondedAtom] = new_bond
				final_vertices.append(R)
			
			else:
				final_vertices.append(atom)

		self.vertices = final_vertices
		return self