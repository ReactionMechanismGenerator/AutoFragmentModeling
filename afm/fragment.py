import os

from rmgpy.species import Species
import rmgpy.molecule.group as gr
import rmgpy.molecule.element as elements
import rmgpy.molecule.converter as converter
from rmgpy.molecule.element import getElement
from rmgpy.molecule.graph import Graph, Vertex
from rmgpy.molecule.molecule import Atom, Bond, Molecule
from rmgpy.molecule.atomtype import getAtomType, AtomTypeError

class CuttingLabel(Vertex):

    def __init__(self, name='', label=''):
        Vertex.__init__(self)
        self.name = name # equivalent to Atom element symbol
        self.label = label # equivalent to Atom label attribute
        self.charge = 0
        self.radicalElectrons = 0
        self.lonePairs = 0
        self.isotope = -1

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return str(self.name)

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "<CuttingLabel '{0}'>".format(str(self))

    @property
    def symbol(self): return self.name

    @property
    def bonds(self): return self.edges

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
            return self.name == other.name
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
        c.name = self.name
        c.label = self.label
        c.charge = self.charge
        c.radicalElectrons = self.radicalElectrons
        c.lonePairs = self.lonePairs
        c.isotope = self.isotope
        return c

class Fragment(Graph):

    def __init__(self,
                label='',
                species_repr=None,
                vertices=None,
                multiplicity=-187,
                props=None):
        self.index = -1
        self.label = label
        self.species_repr = species_repr
        Graph.__init__(self, vertices)
        self._fingerprint = None
        self.props = props or {}
        self.multiplicity = multiplicity

    def __str__(self):
        """
        Return a string representation, in the form 'label(id)'.
        """
        if self.index == -1: return self.label
        else: return '{0}({1:d})'.format(self.label, self.index)

    def __getAtoms(self): return self.vertices
    def __setAtoms(self, atoms): self.vertices = atoms
    atoms = property(__getAtoms, __setAtoms)

    def _repr_png_(self):
        """
        Return a png picture of the molecule, useful for ipython-qtconsole.
        """
        from rmgpy.molecule.draw import MoleculeDrawer
        tempFileName = 'temp_molecule.png'
        MoleculeDrawer().draw(self, 'png', tempFileName)
        png = open(tempFileName,'rb').read()
        os.unlink(tempFileName)
        return png

    def copy(self, deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """
        g = Graph.copy(self, deep)
        other = Fragment(vertices=g.vertices)
        return other

    def clearLabeledAtoms(self):
        """
        Remove the labels from all atoms in the molecule.
        """
        for v in self.vertices:
            v.label = ''

    def containsLabeledAtom(self, label):
        """
        Return :data:`True` if the fragment contains an atom with the label
        `label` and :data:`False` otherwise.
        """
        for v in self.vertices:
            if v.label == label: return True
        return False

    def getLabeledAtom(self, label):
        """
        Return the atoms in the molecule that are labeled.
        """
        for v in self.vertices:
            if v.label == label: return v
        raise ValueError('No vertex in the fragment has the label "{0}".'.format(label))

    def getLabeledAtoms(self):
        """
        Return the labeled atoms as a ``dict`` with the keys being the labels
        and the values the atoms themselves. If two or more atoms have the
        same label, the value is converted to a list of these atoms.
        """
        labeled = {}
        for v in self.vertices:
            if v.label != '':
                if v.label in labeled:
                    if isinstance(labeled[v.label],list):
                        labeled[v.label].append(v)
                    else:
                        labeled[v.label] = [labeled[v.label]]
                        labeled[v.label].append(v)
                else:
                    labeled[v.label] = v
        return labeled

    def removeAtom(self, atom):
        """
        Remove `atom` and all bonds associated with it from the graph. Does
        not remove atoms that no longer have any bonds as a result of this
        removal.
        """
        self._fingerprint = None
        return self.removeVertex(atom)

    def hasBond(self, atom1, atom2):
        """
        Returns ``True`` if atoms `atom1` and `atom2` are connected
        by an bond, or ``False`` if not.
        """
        return self.hasEdge(atom1, atom2)

    def getBond(self, atom1, atom2):
        """
        Returns the bond connecting atoms `atom1` and `atom2`.
        """
        return self.getEdge(atom1, atom2)

    def addBond(self, bond):
        """
        Add a `bond` to the graph as an edge connecting the two atoms `atom1`
        and `atom2`.
        """
        self._fingerprint = None
        return self.addEdge(bond)

    def removeBond(self, bond):
        """
        Remove the bond between atoms `atom1` and `atom2` from the graph.
        Does not remove atoms that no longer have any bonds as a result of
        this removal.
        """
        self._fingerprint = None
        return self.removeEdge(bond)

    def getNetCharge(self):
        """
        Iterate through the atoms in the structure and calculate the net charge
        on the overall fragment.
        """
        charge = 0
        for v in self.vertices:
            charge += v.charge
        return charge

    def merge(self, other):
        """
        Merge two fragments so as to store them in a single :class:`Fragment`
        object. The merged :class:`Fragment` object is returned.
        """
        g = Graph.merge(self, other)
        fragment = Fragment(vertices=g.vertices)
        return fragment

    def split(self):
        """
        Convert a single :class:`Fragment` object containing two or more
        unconnected fragment into separate class:`Fragment` objects.
        """
        graphs = Graph.split(self)
        fragments = []
        for g in graphs:
            fragment = Fragment(vertices=g.vertices)
            fragments.append(fragment)
        return fragments

    def getSingletCarbeneCount(self):
        """
        Return the total number of singlet carbenes (lone pair on a carbon atom)
        in the fragment. Counts the number of carbon atoms with a lone pair.
        In the case of [C] with two lone pairs, this method will return 1.
        """
        carbenes = 0
        for v in self.vertices:
            if isinstance(v, Atom) and v.isCarbon() and v.lonePairs > 0:
                carbenes += 1
        return carbenes
    
    def getRadicalCount(self):
        """
        Return the total number of radical electrons on all atoms in the
        molecule. In this function, monoradical atoms count as one, biradicals
        count as two, etc.
        """
        radicals = 0
        for v in self.vertices:
            radicals += v.radicalElectrons
        return radicals
        
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
                
                cuttinglabel_name = atom_replace_dict[element_symbol]

                cuttinglabel = CuttingLabel(name=cuttinglabel_name)
                for bondedAtom, bond in atom.edges.iteritems():
                    new_bond = Bond(bondedAtom, cuttinglabel, order=bond.order)
                    
                    bondedAtom.edges[cuttinglabel] = new_bond
                    del bondedAtom.edges[atom]

                    cuttinglabel.edges[bondedAtom] = new_bond
                final_vertices.append(cuttinglabel)
            
            else:
                final_vertices.append(atom)

        self.vertices = final_vertices
        return self

    def isSubgraphIsomorphic(self, other, initialMap=None):
        """
        Fragment's subgraph isomorphism check is done by first creating 
        a representative molecule of fragment, and then following same procedure
        of subgraph isomorphism check of `Molecule` object aganist `Group` object
        """
        if not isinstance(other, gr.Group):
            raise TypeError('Got a {0} object for parameter "other", when a Molecule object is required.'.format(other.__class__))
        group = other

        mapping = self.assign_representative_molecule()

        # Check multiplicity
        if group.multiplicity:
            if self.mol_repr.multiplicity not in group.multiplicity: return False

        # Compare radical counts
        if self.mol_repr.getRadicalCount() < group.radicalCount:
            return False

        # Compare element counts
        element_count = self.mol_repr.get_element_count()
        for element, count in group.elementCount.iteritems():
            if element not in element_count:
                return False
            elif element_count[element] < count:
                return False

        # Do the isomorphism comparison
        new_initial_map = None
        if initialMap:
            new_initial_map = {}
            for fragment_vertex in initialMap:
                repr_mol_vertex = mapping[fragment_vertex]
                new_initial_map[repr_mol_vertex] = initialMap[fragment_vertex]

        result = Graph.isSubgraphIsomorphic(self.mol_repr, other, new_initial_map)
        return result

    def assign_representative_molecule(self):

        # create a molecule from fragment.vertices.copy
        mapping = self.copyAndMap()

        # replace CuttingLabel with CC
        atoms = []
        additional_atoms = []
        additional_bonds = []
        for vertex in self.vertices:

            mapped_vertex = mapping[vertex]
            if isinstance(mapped_vertex, CuttingLabel):

                # replace cutting label with atom C
                atom_C1 = Atom(element=getElement('C'), 
                            radicalElectrons=0, 
                            charge=0, 
                            lonePairs=0)

                for bondedAtom, bond in mapped_vertex.edges.iteritems():
                    new_bond = Bond(bondedAtom, atom_C1, order=bond.order)
                    
                    bondedAtom.edges[atom_C1] = new_bond
                    del bondedAtom.edges[mapped_vertex]

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
                atoms.append(mapped_vertex)

        mol_repr = Molecule()
        mol_repr.atoms = atoms
        for atom in additional_atoms: mol_repr.addAtom(atom)
        for bond in additional_bonds: mol_repr.addBond(bond)
        # update connectivity
        mol_repr.update()

        # create a species object from molecule
        self.mol_repr = mol_repr

        return mapping

    def assign_representative_species(self):

        self.assign_representative_molecule()
        self.species_repr = Species(molecule=[self.mol_repr])

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

    def update(self):

        for v in self.vertices:
            if isinstance(v, Atom):
                v.updateCharge()

        self.updateAtomTypes()
        self.updateMultiplicity()
        self.sortVertices()

    def updateAtomTypes(self, logSpecies=True, raiseException=True):
        """
        Iterate through the atoms in the structure, checking their atom types
        to ensure they are correct (i.e. accurately describe their local bond
        environment) and complete (i.e. are as detailed as possible).
        
        If `raiseException` is `False`, then the generic atomType 'R' will
        be prescribed to any atom when getAtomType fails. Currently used for
        resonance hybrid atom types.
        """
        #Because we use lonepairs to match atomtypes and default is -100 when unspecified,
        #we should update before getting the atomtype.
        self.updateLonePairs()

        for v in self.vertices:
            if not isinstance(v, Atom): continue
            try:
                v.atomType = getAtomType(v, v.edges)
            except AtomTypeError:
                if logSpecies:
                    logging.error("Could not update atomtypes for this fragment:\n{0}".format(self.toAdjacencyList()))
                if raiseException:
                    raise
                v.atomType = atomTypes['R']

    def updateMultiplicity(self):
        """
        Update the multiplicity of a newly formed molecule.
        """
        # Assume this is always true
        # There are cases where 2 radicalElectrons is a singlet, but
        # the triplet is often more stable, 
        self.multiplicity = self.getRadicalCount() + 1


    def updateLonePairs(self):
        """
        Iterate through the atoms in the structure and calculate the
        number of lone electron pairs, assuming a neutral molecule.
        """
        for v in self.vertices:
            if not isinstance(v, Atom): continue
            if not v.isHydrogen():
                order = v.getBondOrdersForAtom()
                v.lonePairs = (elements.PeriodicSystem.valence_electrons[v.symbol] - v.radicalElectrons - v.charge - int(order)) / 2.0
                if v.lonePairs % 1 > 0 or v.lonePairs > 4:
                    logging.error("Unable to determine the number of lone pairs for element {0} in {1}".format(v,self))
            else:
                v.lonePairs = 0

    def getFormula(self):
        """
        Return the molecular formula for the fragment.
        """
        
        # Count the number of each element in the molecule
        elements = {}
        cuttinglabels = {}
        for atom in self.vertices:
            symbol = atom.symbol
            elements[symbol] = elements.get(symbol, 0) + 1
        
        # Use the Hill system to generate the formula
        formula = ''
        
        # Carbon and hydrogen always come first if carbon is present
        if 'C' in elements.keys():
            count = elements['C']
            formula += 'C{0:d}'.format(count) if count > 1 else 'C'
            del elements['C']
            if 'H' in elements.keys():
                count = elements['H']
                formula += 'H{0:d}'.format(count) if count > 1 else 'H'
                del elements['H']

        # Other atoms are in alphabetical order
        # (This includes hydrogen if carbon is not present)
        keys = elements.keys()
        keys.sort()
        for key in keys:
            count = elements[key]
            formula += '{0}{1:d}'.format(key, count) if count > 1 else key
        
        return formula

    def toRDKitMol(self, removeHs=False, returnMapping=True):
        """
        Convert a molecular structure to a RDKit rdmol object.
        """
        if removeHs: 
            # because we're replacing
            # cutting labels with hydrogens
            # so do not allow removeHs to be True
            raise "Currently fragment toRDKitMol only allows keeping all the hydrogens."

        # create a molecule from fragment.vertices.copy
        mapping = self.copyAndMap()

        # replace CuttingLabel with CC
        atoms = []
        for vertex in self.vertices:

            mapped_vertex = mapping[vertex]
            if isinstance(mapped_vertex, CuttingLabel):

                # replace cutting label with atom H
                atom_H = Atom(element=getElement('H'), 
                            radicalElectrons=0, 
                            charge=0, 
                            lonePairs=0)

                for bondedAtom, bond in mapped_vertex.edges.iteritems():
                    new_bond = Bond(bondedAtom, atom_H, order=bond.order)
                    
                    bondedAtom.edges[atom_H] = new_bond
                    del bondedAtom.edges[mapped_vertex]

                    atom_H.edges[bondedAtom] = new_bond

                mapping[vertex] = atom_H
                atoms.append(atom_H)

            else:
                atoms.append(mapped_vertex)

        # Note: mapping is a dict with 
        # key: self.vertex and value: mol0.atom

        mol0 = Molecule()
        mol0.atoms = atoms

        rdmol, rdAtomIdx_mol0 = converter.toRDKitMol(mol0, removeHs=removeHs, 
                                                     returnMapping=returnMapping, 
                                                     sanitize=True)

        rdAtomIdx_frag = {}
        for frag_atom, mol0_atom in mapping.iteritems():
            rd_idx = rdAtomIdx_mol0[mol0_atom]
            rdAtomIdx_frag[frag_atom] = rd_idx

        # sync the order of fragment vertices with the order
        # of mol0.atoms since mol0.atoms is changed/sorted in 
        # converter.toRDKitMol().
        # Since the rdmol's atoms order is same as the order of mol0's atoms,
        # the synchronization between fragment.atoms order and mol0.atoms order
        # is necessary to make sure the order of fragment vertices
        # reflects the order of rdmol's atoms
        vertices_order = []
        for v in self.vertices:
            a = mapping[v]
            idx = mol0.atoms.index(a)
            vertices_order.append((v, idx))

        adapted_vertices = [tup[0] for tup in sorted(vertices_order, key=lambda tup: tup[1])]

        self.vertices = adapted_vertices

        return rdmol, rdAtomIdx_frag

    def toAdjacencyList(self, 
                        label='', 
                        removeH=False, 
                        removeLonePairs=False, 
                        oldStyle=False):
        """
        Convert the molecular structure to a string adjacency list.
        """
        from rmgpy.molecule.adjlist import toAdjacencyList
        result = toAdjacencyList(self.vertices, 
                                 self.multiplicity,  
                                 label=label, 
                                 group=False, 
                                 removeH=removeH, 
                                 removeLonePairs=removeLonePairs, 
                                 oldStyle=oldStyle)
        return result

