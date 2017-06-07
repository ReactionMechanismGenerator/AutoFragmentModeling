
class Fragment(object):

	def __init__(self,
				label,
				species_repr):
		self.index = -1
		self.label = label
		self.species_repr = species_repr

	def __str__(self):
		"""
		Return a string representation, in the form 'label(id)'.
		"""
		if self.index == -1: return self.label
		else: return '{0}({1:d})'.format(self.label, self.index)
	
	def _repr_png_(self):
		if len(self.species_repr.molecule) > 0:
			return self.species_repr.molecule[0]._repr_png_()
		else:
			return None