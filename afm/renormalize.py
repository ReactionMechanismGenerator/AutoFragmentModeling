import os
import re
from rmgpy.species import Species
from rmgpy.molecule.molecule import Molecule
from afm.fragment import Fragment
import numpy as np
import math
from rdkit import Chem
from rdkit.Chem import AllChem
import random

def merge_frag_to_frag(frag1, frag2, label): # label should match the desired merging label on frag2
    from rmgpy.molecule import Bond
    from afm.fragment import Fragment, CuttingLabel
    frag_spe1 = Fragment().from_smiles_like_string(frag1)
    frag_spe2 = Fragment().from_smiles_like_string(frag2)
    # find position of desired CuttingLabel
    # need to find CuttingLabel on frag2 first
    for vertex in frag_spe2.vertices:
        if isinstance(vertex, CuttingLabel):
            if vertex.symbol == label:
                cut2 = vertex
                atom2 = list(cut2.edges)[0]
                frag_spe2.remove_atom(cut2)
                break

    if cut2.symbol[0] == 'L':
        Ctl = cut2.symbol.replace('L', 'R')
    else: # that means this CuttingLabel is R something
        Ctl = cut2.symbol.replace('R', 'L')
        
    # merge to frag_spe1
    for vertex in frag_spe1.vertices:
        if isinstance(vertex, CuttingLabel):
            if vertex.symbol == Ctl:
                cut1 = vertex
                atom1 = list(cut1.edges)[0]
                frag_spe1.remove_atom(cut1)
                break
    
    # new merged fragment
    new_frag = frag_spe1.merge(frag_spe2)
    new_frag.add_bond(Bond(atom1=atom1, atom2=atom2, order=1))
    new_frag = new_frag.copy(deep=True)
    new_frag.update()
    return new_frag # return Fragment object

def random_pick_frag(target_dict, label):
    import random
    import re
    sum_dict = 0
    frag_dict_list = [] # keys in target_dict
    frag_dict_prob = [] # values in target_dict
    for frag, amt in target_dict.items():
        # only pick frag which has only 1R or 1L
        if len(re.findall(r'[LR]', frag.to_smiles())) == 1:
            sum_dict += amt
            frag_dict_list.append(frag)
    for frag in frag_dict_list:
        amt = target_dict[frag]
        frag_dict_prob.append(amt/sum_dict)
    
    x = random.uniform(0,1)
    cumulative_probability = 0.0
    for item, item_probability in zip(frag_dict_list, frag_dict_prob):
        cumulative_probability += item_probability
        if x < cumulative_probability:
            break
    return item

def pair_frag(amount, target_dict, label):
    additional_frag_list = []
    frag1 = random_pick_frag(target_dict, label)
    if target_dict[frag1] >= amount:
        target_dict[frag1] -= amount
        additional_frag_list.append((frag1, amount))
    else:
        # pick additional one to make sure no negative frag amount
        remain = amount - target_dict[frag1]
        additional_frag_list.append((frag1, amount))
        target_dict[frag1] = 0

        while remain > 0:
            frag1 = random_pick_frag(target_dict, label)
            if target_dict[frag1] >= remain:
                target_dict[frag1] -= remain
                additional_frag_list.append((frag1, remain))
                remain = 0
            else:
                frag_amt = target_dict[frag1]
                target_dict[frag1] = 0
                additional_frag_list.append((frag1, frag_amt))
                remain = remain - frag_amt
    return additional_frag_list

def cut_molecule(self, mol_to_cut):
    """
    For given input, output a list of cut fragments (either string or Fragment)
    """
    # if input is smiles string, output smiles
    if isinstance(mol_to_cut, str):
        mol = self.from_smiles_like_string(mol_to_cut)
        mol = mol.generate_resonance_structures()[0]
    # if input is fragment, output frafgment
    elif isinstance(mol_to_cut, Fragment):
        mol = mol_to_cut.generate_resonance_structures()[0]

    # slice mol
    frag_smiles_list = []
    arom_cut_frag = self.sliceitup_arom(mol.to_smiles())
    for frag in arom_cut_frag:
        aliph_cut_frag = self.sliceitup_aliph(frag)
        for ele in aliph_cut_frag:
            frag_smiles_list.append(ele)

    frag_list_new = []
    for frag_smiles in frag_smiles_list:
        n_frag_smiles = frag_smiles.replace('F', 'R')
        nn_frag_smiles = n_frag_smiles.replace('Cl', 'L')
        # cutting position near aromatic 
        nnn_frag_smiles = nn_frag_smiles.replace('Br', 'R')
        new_frag_smiles = nnn_frag_smiles.replace('I', 'L')

        frag_list_new.append(new_frag_smiles)

    if isinstance(mol_to_cut, str):
        return frag_list_new
    elif isinstance(mol_to_cut, Fragment):
        frag_list = []
        for frag in frag_list_new:
            frag = Fragment().from_smiles_like_string(frag)
            res_frag = frag.generate_resonance_structures()[0]
            frag_list.append(res_frag)
        return frag_list

def sliceitup_arom(molecule, smarts_rxn_list=None):
    """
    Several specified aromatic patterns
    """
    # if input is smiles string, output smiles
    if isinstance(molecule, str):
        molecule_smiles = molecule
    # if input is fragment, output frafgment
    elif isinstance(molecule, Fragment):
        mol = molecule.generate_resonance_structures()[0]
        molecule_smiles = mol.to_smiles()

    # if input is Fragment (with Cuttinglabel), need to transform to special Atom
    import re
    cutting_label_list = re.findall(r'([LR]\d?)', molecule_smiles)
    if cutting_label_list != []:
        n_frag_smiles = molecule_smiles.replace('R', '[Na]')
        nn_frag_smiles = n_frag_smiles.replace('L', '[Mg]')
        molecule_smiles = nn_frag_smiles

    molecule1 = Chem.MolFromSmiles(molecule_smiles)
    molecule_Br2 = Chem.MolFromSmiles('BrI')
    molecule_I2 = Chem.MolFromSmiles('BrIIBr')

    if smarts_rxn_list:
        arom_rxn_list = []
        for (sma,mol2) in smarts_rxn_list:
            rxn = AllChem.ReactionFromSmarts(sma)
            if mol2 == 1:
                arom_rxn_list.append((rxn,molecule_Br2))
            elif mol2 == 2:
                arom_rxn_list.append((rxn,molecule_I2))

    else: # use the default smarts rxn 
        rxn_near_arom_1 = AllChem.ReactionFromSmarts("[C:1]-[C:2]-[C:3]-[C:4]-[C:5](-[c:6]1:[c:7]:[c:8]:[c:9]:[c:10]:[c:11]:1)-[C:12](-[C:13]-[C:14]-[C:15]-[C:16])-[c:17]1:[c:18]:[c:19]:[c:20]:[c:21]:[c:22]:1.[Br:23]-[I:24]-[I:25]-[Br:26]>>\
        [Br:23]-[C:2]-[C:3]-[C:4]-[C:5](-[c:6]1:[c:7]:[c:8]:[c:9]:[c:10]:[c:11]:1)-[C:12](-[C:13]-[C:14]-[C:15]-[Br:26])-[c:17]1:[c:18]:[c:19]:[c:20]:[c:21]:[c:22]:1.[I:24]-[C:1].[I:25]-[C:16]")

        rxn_near_arom_2 = AllChem.ReactionFromSmarts("[c:1]:[c:2]:[c:3]:[c:4]:[c:5]:[c:6]-[C:7](-[C:8]-[C:9]-[C:10]-[C:11])-[C:12]-[C:13]-[C:14]-[C:15].[Br:16]-[I:17]-[I:18]-[Br:19]>>\
        [c:1]:[c:2]:[c:3]:[c:4]:[c:5]:[c:6]-[C:7](-[C:8]-[C:9]-[C:10]-[Br:19])-[C:12]-[C:13]-[C:14]-[Br:16].[I:17]-[C:15].[I:18]-[C:11]")

        rxn_near_arom_3 = AllChem.ReactionFromSmarts("[C:1]=[C:2](-[C:3]-[C:4]-[C:5]-[C:6])-[c:7]1:[c:8]:[c:9]:[c:10]:[c:11]:[c:12]:1.[Br:13]-[I:14]>>\
        [C:1]=[C:2](-[C:3]-[C:4]-[C:5]-[Br:13])-[c:7]1:[c:8]:[c:9]:[c:10]:[c:11]:[c:12]:1.[C:6]-[I:14]")

        rxn_near_arom_4 = AllChem.ReactionFromSmarts("[c:1]:[c:2]:[c:3]:[c:4]:[c:5]:[c:6]-[C:7]-[C:8]-[C:9]-[C:10]-[C:11].[Br:12]-[I:13]>>\
        [c:1]:[c:2]:[c:3]:[c:4]:[c:5]:[c:6]-[C:7]-[C:8]-[C:9]-[C:10]-[Br:12].[I:13]-[C:11]")

        arom_rxn_list = [(rxn_near_arom_1,molecule_I2), (rxn_near_arom_2,molecule_I2), (rxn_near_arom_3,molecule_Br2), (rxn_near_arom_4,molecule_Br2)]

    frag_list = []
    for i, (arom_rxn,pseudo_spe) in enumerate(arom_rxn_list):
        ps = arom_rxn.RunReactants((molecule1, pseudo_spe))
        if len(ps) != 0:
            for i in range(len(ps[0])):
                frag_list.append(Chem.MolToSmiles(ps[0][i]))
            break

    if frag_list == []:
        frag_list.append(molecule_smiles)
    else:
        if cutting_label_list != []:
            frag_list_replaced = []
            # replace metal atom back to Cuttinglabel
            for metal_frag in frag_list:
                n_frag_smiles = metal_frag.replace('[Na]', 'R')
                nn_frag_smiles = n_frag_smiles.replace('[Mg]', 'L')
                frag_list_replaced.append(nn_frag_smiles)
            frag_list = frag_list_replaced

    if isinstance(molecule, str):
        frag_list_new = []
        for frag in frag_list:
            n_frag_smiles = frag.replace('[Na]', 'R')
            nn_frag_smiles = n_frag_smiles.replace('[Mg]', 'L')
            nnn_frag_smiles = nn_frag_smiles.replace('Br', 'R')
            new_frag_smiles = nnn_frag_smiles.replace('I', 'L')
            frag_list_new.append(new_frag_smiles)
        return frag_list_new
    elif isinstance(molecule, Fragment):
        frag_list_new = []
        for frag in frag_list:
            n_frag_smiles = frag.replace('[Na]', 'R')
            nn_frag_smiles = n_frag_smiles.replace('[Mg]', 'L')
            nnn_frag_smiles = nn_frag_smiles.replace('Br', 'R')
            new_frag_smiles = nnn_frag_smiles.replace('I', 'L')

            frag = Fragment().from_smiles_like_string(new_frag_smiles)
            res_frag = frag.generate_resonance_structures()[0]
            frag_list_new.append(res_frag)
        return frag_list_new

def sliceitup_aliph(molecule, smarts_rxn_list=None, limit_length=True):
    """
    Several specified aliphatic patterns
    """
    # if input is smiles string, output smiles
    if isinstance(molecule, str):
        molecule_smiles = molecule
    # if input is fragment, output frafgment
    elif isinstance(molecule, Fragment):
        mol = molecule.generate_resonance_structures()[0]
        molecule_smiles = mol.to_smiles()

    # if input is Fragment (with Cuttinglabel), need to transform to special Atom
    import re
    cutting_label_list = re.findall(r'([LR]\d?)', molecule_smiles)
    if cutting_label_list != []:
        n_frag_smiles = molecule_smiles.replace('R', '[Na]')
        nn_frag_smiles = n_frag_smiles.replace('L', '[Mg]')
        molecule_smiles = nn_frag_smiles

    molecule1 = Chem.MolFromSmiles(molecule_smiles)
    molecule_F2 = Chem.MolFromSmiles('FCl')
    molecule_Br2 = Chem.MolFromSmiles('BrI')
    molecule_I2 = Chem.MolFromSmiles('BrIIBr')

    if smarts_rxn_list:
        aliph_rxn_list = []
        for (sma,mol2) in smarts_rxn_list:
            rxn = AllChem.ReactionFromSmarts(sma)
            if mol2 == 1:
                aliph_rxn_list.append((rxn,molecule_Br2))
            elif mol2 == 2:
                aliph_rxn_list.append((rxn,molecule_I2))

    else: # use the default smarts rxn 

        rxn_smallchain_1 = AllChem.ReactionFromSmarts("[C:1]-[C:2]-[!c;!R:3]=[!c;!R:4]-[C:5]-[C:6]-[C:7]-[C:8].[F:11]-[Cl:12]>>\
        [C:1]-[C:2]-[!c;!R:3]=[!c;!R:4]-[C:5]-[C:6]-[C:7][F:11].[Cl:12]-[C:8]") # CCC=CCCC

        rxn_smallchain_2 = AllChem.ReactionFromSmarts("[!c;!R:1]=[!c;!R:2]-[C:3]-[C:4]-[C:5]-[C:6].[F:7]-[Cl:8]>>\
        [!c;!R:1]=[!c;!R:2]-[C:3]-[C:4]-[C:5]-[F:7].[Cl:8]-[C:6]") # C=CCCC

        rxn_smallchain_3 = AllChem.ReactionFromSmarts("[C;!R:1]-[C;!R:2]-[C;!R:3]-[C;!R:4]-[C;!R:5]-[C;!R:6].[F:7]-[Cl:8]>>\
        [C;!R:1]-[C;!R:2]-[C;!R:3]-[F:7].[Cl:8]-[C;!R:4]-[C;!R:5]-[C;!R:6]") # CCCCCC

        aliph_rxn_list = [(rxn_smallchain_1,molecule_F2), (rxn_smallchain_2,molecule_F2), (rxn_smallchain_3,molecule_F2)]

    frag_list = []

    for i, (aliph_rxn,pseudo_spe) in enumerate(aliph_rxn_list):
        ps = aliph_rxn.RunReactants((molecule1, pseudo_spe))
        if len(ps) != 0:
            for loop in range(1000):
                ind = np.random.randint(0, len(ps))
                if limit_length == True:
                    if ps[ind][0].GetNumAtoms() > 4 and ps[ind][1].GetNumAtoms() > 5:
                        for i in range(len(ps[ind])):
                            frag_list.append(Chem.MolToSmiles(ps[ind][i]))
                        break
                else:
                    for i in range(len(ps[ind])):
                        frag_list.append(Chem.MolToSmiles(ps[ind][i]))
                    break
            if frag_list != []:
                break

    if frag_list == []:
        frag_list.append(molecule_smiles)
    else:
        if cutting_label_list != []:
            frag_list_replaced = []
            # replace metal atom back to Cuttinglabel
            for metal_frag in frag_list:
                n_frag_smiles = metal_frag.replace('[Na]', 'R')
                nn_frag_smiles = n_frag_smiles.replace('[Mg]', 'L')
                frag_list_replaced.append(nn_frag_smiles)
            frag_list = frag_list_replaced

    if isinstance(molecule, str):
        frag_list_new = []
        for frag in frag_list:
            n_frag_smiles = frag.replace('[Na]', 'R')
            nn_frag_smiles = n_frag_smiles.replace('[Mg]', 'L')
            nnn_frag_smiles = nn_frag_smiles.replace('Br', 'R')
            new_frag_smiles = nnn_frag_smiles.replace('I', 'L')
            frag_list_new.append(new_frag_smiles)
        return frag_list_new
    elif isinstance(molecule, Fragment):
        frag_list_new = []
        for frag in frag_list:
            n_frag_smiles = frag.replace('[Na]', 'R')
            nn_frag_smiles = n_frag_smiles.replace('[Mg]', 'L')
            nnn_frag_smiles = nn_frag_smiles.replace('Br', 'R')
            new_frag_smiles = nnn_frag_smiles.replace('I', 'L')

            frag = Fragment().from_smiles_like_string(new_frag_smiles)
            res_frag = frag.generate_resonance_structures()[0]
            frag_list_new.append(res_frag)
        return frag_list_new

# include the case of multiple cuttinglabel pattern smiles
def generate_smarts_rxn(pattern_smiles):
    # change fragment to pattern
    # need to replace the Cuttinglabel

    # if input is Fragment (with Cuttinglabel), need to transform to special Atom
    import re
    cutting_label_list = re.findall(r'([LR]\d?)', pattern_smiles)
    if cutting_label_list != []:
        n_frag_smiles = pattern_smiles.replace('R', '[Na]')
        nn_frag_smiles = n_frag_smiles.replace('L', '[Mg]')
        molecule_smiles = nn_frag_smiles
    
    molecule1 = Chem.MolFromSmiles(molecule_smiles)
    for i in range(molecule1.GetNumHeavyAtoms()):
        molecule1.GetAtomWithIdx(i).SetAtomMapNum(i+1)
    
    terminal_atom_index = []
    metal_atom_index = []
    for atom in molecule1.GetAtoms():
        if atom.GetSymbol() == 'Na' or atom.GetSymbol() == 'Mg':
            metal_atom_index.append(atom.GetIdx())
            atom_has_cuttinglabel = atom.GetNeighbors()[0]
            # print(atom_has_cuttinglabel.GetIdx())
            terminal_atom_index.append(atom_has_cuttinglabel.GetIdx())
    
    # create SMARTS pattern reactant 
    # remove the atom Cuttinglabel (metal atom) and the connected bond
    # add a new C atom and connected bond
    r1 = Chem.rdchem.RWMol(molecule1)
    C_map_num_list = []
    for i in range(len(terminal_atom_index)):
        new_C = Chem.rdchem.Atom('C')
        new_C.SetAtomMapNum(r1.GetNumHeavyAtoms()+i+1)
        r1.ReplaceAtom(metal_atom_index[i], new_C)
        C_map_num_list.append(r1.GetNumHeavyAtoms()+i+1)
    
    if len(terminal_atom_index) == 1:
        r2 = Chem.MolFromSmiles('BrI')
        mol2 = 1
    elif len(terminal_atom_index) == 2:
        r2 = Chem.MolFromSmiles('BrIIBr')
        mol2 = 2
    else:
        print('Currently need to extend this algorithm')

    for i in range(r2.GetNumAtoms()):
        r2.GetAtomWithIdx(i).SetAtomMapNum(C_map_num_list[-1]+i+1)
    
    p1 = Chem.rdchem.RWMol(molecule1)
    cutting_label_index = []
    for atom in r2.GetAtoms():
        if atom.GetSymbol() == 'Br':
            cutting_label_index.append(atom.GetAtomMapNum())

    for i in range(len(terminal_atom_index)):
        new_Br = Chem.rdchem.Atom('Br')
        new_Br.SetAtomMapNum(cutting_label_index[i])
        p1.ReplaceAtom(metal_atom_index[i], new_Br)
        
    p_other = []
    cutting_label_index = []
    for atom in r2.GetAtoms():
        if atom.GetSymbol() == 'I':
            cutting_label_index.append(atom.GetAtomMapNum())

    for i in range(len(terminal_atom_index)):
        p2 = Chem.MolFromSmiles('IC')
        p2.GetAtomWithIdx(0).SetAtomMapNum(cutting_label_index[i]) # 'I'
        p2.GetAtomWithIdx(1).SetAtomMapNum(C_map_num_list[i]) # 'C'
        p_other.append(p2)
    
    # create smarts reaction
    r = Chem.rdChemReactions.ChemicalReaction()
    r.AddReactantTemplate(r1)
    r.AddReactantTemplate(r2)
    r.AddProductTemplate(p1)
    for i in range(len(p_other)):
        r.AddProductTemplate(p_other[i])
    smarts_rxn = Chem.rdChemReactions.ReactionToSmarts(r)
    
    return_list = [(smarts_rxn,mol2)]
    
    if len(terminal_atom_index) >= 2:
        # need to generate additional smarts rxn 
        # so that it is available to handle the reactant to cut is fragment with cuttinglabel
        # ex: pattern smiles = RCCCCL
        # we should not only generate CCCCCC + RL -> CL + RCCCCL + RC
        # but also generate RCCCCC + RL -> RCCCCL + RC & CCCCCL + RL -> CL + RCCCCL
        # the following part is to generate the two smarts rxn
        
        for i in range(len(terminal_atom_index)):
            r1 = Chem.rdchem.RWMol(molecule1)
            C_map_num_list = []
            new_C = Chem.rdchem.Atom('C')
            new_C.SetAtomMapNum(r1.GetNumHeavyAtoms()+i+1)
            r1.ReplaceAtom(metal_atom_index[i], new_C)
            C_map_num_list.append(r1.GetNumHeavyAtoms()+i+1)

            # print(Chem.MolToSmarts(r1))

            r2 = Chem.MolFromSmiles('BrI')
            mol2 = 1
            
            for j in range(r2.GetNumAtoms()):
                r2.GetAtomWithIdx(j).SetAtomMapNum(C_map_num_list[-1]+j+1)
            # print(Chem.MolToSmarts(r2))

            p1 = Chem.rdchem.RWMol(molecule1)

            cutting_label_index = []
            for atom in r2.GetAtoms():
                if atom.GetSymbol() == 'Br':
                    cutting_label_index.append(atom.GetAtomMapNum())

            new_Br = Chem.rdchem.Atom('Br')
            new_Br.SetAtomMapNum(cutting_label_index[0])
            p1.ReplaceAtom(metal_atom_index[i], new_Br)
            # print(Chem.MolToSmarts(p1))

            p_other = []
            cutting_label_index = []
            for atom in r2.GetAtoms():
                if atom.GetSymbol() == 'I':
                    cutting_label_index.append(atom.GetAtomMapNum())

            p2 = Chem.MolFromSmiles('IC')
            p2.GetAtomWithIdx(0).SetAtomMapNum(cutting_label_index[0]) # 'I'
            p2.GetAtomWithIdx(1).SetAtomMapNum(C_map_num_list[0]) # 'C'
            p_other.append(p2)
            # print(Chem.MolToSmarts(p2))

            r = Chem.rdChemReactions.ChemicalReaction()
            r.AddReactantTemplate(r1)
            r.AddReactantTemplate(r2)
            r.AddProductTemplate(p1)
            for i in range(len(p_other)):
                r.AddProductTemplate(p_other[i])
            smarts_rxn = Chem.rdChemReactions.ReactionToSmarts(r)
            # print(smarts_rxn)
            
            return_list.append((smarts_rxn,mol2))
    
    return return_list

### main renormalization function
def frag_renormalize(core_species_concentrations, core_species, core_smiles):
    
    # replace core_species_concentrations with y if needed
    moles_dict = {} # key: species, value: moles
    for ind, spe in enumerate(core_species):
        moles_dict[spe] = core_species_concentrations[ind]
    
    r_m_sum = 0
    l_m_sum = 0
    moles_remain = []
    frag_list = []

    for spe, amt in moles_dict.items():
        import re
        if spe.molecule[0].is_radical() == False and amt >= 1e-6:
            if re.findall(r'R', spe.molecule[0].to_smiles()) == [] and re.findall(r'L', spe.molecule[0].to_smiles()) == []:
                moles_remain.append((spe, amt))
            else: # fragment
                frag_list.append((spe.molecule[0], amt))
                r_m_sum += len(re.findall(r'R', spe.molecule[0].to_smiles())) * amt
                l_m_sum += len(re.findall(r'L', spe.molecule[0].to_smiles())) * amt
    if abs(r_m_sum-l_m_sum) >= 1e-1:
        print('Cuttinglabel amount inbalance')

    # Categorize fragments
    general_R_list = []
    general_L_list = []
    multi_label_frag_3 = []
    multi_label_frag_4 = []
    r_l_moles = []
    rr_ll_list = []

    one_R_dict = {}
    one_L_dict = {}

    arom = 0

    for frag, amt in frag_list:
        smi = frag.to_smiles()
        if len(re.findall(r'R', frag.to_smiles())) == 1 and len(re.findall(r'L', frag.to_smiles())) == 0:
            one_R_dict[frag] = amt
            if 'c' in smi:
                num = smi.count('c')
                arom += amt*num/6
        elif len(re.findall(r'R', frag.to_smiles())) == 2 and len(re.findall(r'L', frag.to_smiles())) == 0:
            general_R_list.append((frag, amt * 2))
            rr_ll_list.append(frag)
            if 'c' in smi:
                num = smi.count('c')
                arom += amt*num/6
        elif len(re.findall(r'R', frag.to_smiles())) == 0 and len(re.findall(r'L', frag.to_smiles())) == 1:
            one_L_dict[frag] = amt
            if 'c' in smi:
                num = smi.count('c')
                arom += amt*num/6
        elif len(re.findall(r'R', frag.to_smiles())) == 0 and len(re.findall(r'L', frag.to_smiles())) == 2:
            general_L_list.append((frag, amt * 2))
            rr_ll_list.append(frag)
            if 'c' in smi:
                num = smi.count('c')
                arom += amt*num/6
        elif len(re.findall(r'R', frag.to_smiles())) == 1 and len(re.findall(r'L', frag.to_smiles())) == 1:
            r_l_moles.append((frag, amt))
            if 'c' in smi:
                num = smi.count('c')
                arom += amt*num/6
        else:
            # 3 or 4 labels on the fragment
            # 3 labels:
            if len(re.findall(r'[LR]', frag.to_smiles())) == 3:
                multi_label_frag_3.append((frag, amt)) # 2R1L, 1R2L, 3R, 3L
            # 4 labels:
            elif len(re.findall(r'[LR]', frag.to_smiles())) == 4:
                multi_label_frag_4.append((frag, amt))
        
    # Deal with 4-label and 3-label fragments
    multi_label_amount = 0

    for species, amount in multi_label_frag_4: # 4R, 3R1L, 2R2L, 1R3L, 4L
        print('enter')
        # 4R, reattach with 1L, then add to 3 label list
        if len(re.findall(r'R', species.to_smiles())) == 4 and len(re.findall(r'L', species.to_smiles())) == 0:
            paired_frag_list = pair_frag(amount, one_L_dict, 'L')
            multi_label_amount += amount
            # reattach
            for frag_amt in paired_frag_list:
                frag = frag_amt[0]
                amt = frag_amt[1]
                frag1 = frag # 1L
                frag2 = species # 4R
                frag_new = merge_frag_to_frag(frag1.to_smiles(), frag2.to_smiles(), 'R') # L,R,R -> 3R

                multi_label_frag_3.append((frag_new, amt))

        # 4L, reattach with 1R, then add to 3 label list
        elif len(re.findall(r'R', species.to_smiles())) == 0 and len(re.findall(r'L', species.to_smiles())) == 4:
            paired_frag_list = pair_frag(amount, one_R_dict, 'R')
            multi_label_amount += amount
            # reattach
            for frag_amt in paired_frag_list:
                frag = frag_amt[0]
                amt = frag_amt[1]
                frag1 = frag # 1R
                frag2 = species # 4L
                frag_new = merge_frag_to_frag(frag2.to_smiles(), frag1.to_smiles(), 'R') # L,R,R -> 3L

                multi_label_frag_3.append((frag_new, amt))

        # 2R2L, reattach with 1L
        elif len(re.findall(r'R', species.to_smiles())) == 2 and len(re.findall(r'L', species.to_smiles())) == 2:
            paired_frag_list = pair_frag(amount, one_L_dict, 'L')
            multi_label_amount += amount
            # reattach
            for frag_amt in paired_frag_list:
                frag = frag_amt[0]
                amt = frag_amt[1]
                frag1 = frag # 1L
                frag2 = species # 2R2L
                frag_new = merge_frag_to_frag(frag1.to_smiles(), frag2.to_smiles(), 'R') # L,R,R -> 1R2L

                multi_label_frag_3.append((frag_new, amt))

        # 3R1L, reattach with 1R
        elif len(re.findall(r'R', species.to_smiles())) == 3 and len(re.findall(r'L', species.to_smiles())) == 1:
            paired_frag_list = pair_frag(amount, one_R_dict, 'R')
            multi_label_amount += amount
            # reattach
            for frag_amt in paired_frag_list:
                frag = frag_amt[0]
                amt = frag_amt[1]
                frag1 = frag # 1R
                frag2 = species # 3R1L
                frag_new = merge_frag_to_frag(frag2.to_smiles(), frag1.to_smiles(), 'R') # L,R,R -> 3R

                multi_label_frag_3.append((frag_new, amt))

        # 1R3L, reattach with 1L
        elif len(re.findall(r'R', species.to_smiles())) == 1 and len(re.findall(r'L', species.to_smiles())) == 3:
            paired_frag_list = pair_frag(amount, one_L_dict, 'L')
            multi_label_amount += amount
            # reattach
            for frag_amt in paired_frag_list:
                frag = frag_amt[0]
                amt = frag_amt[1]
                frag1 = frag # 1L
                frag2 = species # 1R3L
                frag_new = merge_frag_to_frag(frag1.to_smiles(), frag2.to_smiles(), 'R') # L,R,R -> 3L

                multi_label_frag_3.append((frag_new, amt))

    for species, amount in multi_label_frag_3:
        print('here')
        # 2R1L
        if len(re.findall(r'R', species.to_smiles())) == 2 and len(re.findall(r'L', species.to_smiles())) == 1:
            paired_frag_list = pair_frag(amount, one_R_dict, 'R')
            multi_label_amount += amount
            # reattach
            for frag_amt in paired_frag_list:
                frag = frag_amt[0]
                amt = frag_amt[1]
                frag1 = frag # 1R
                frag2 = species # 2R1L
                frag_new = merge_frag_to_frag(frag2.to_smiles(), frag1.to_smiles(), 'R') # L,R,R

                general_R_list.append((frag_new, amt * 2))
                rr_ll_list.append(frag_new)

        # 1R2L
        elif len(re.findall(r'R', species.to_smiles())) == 1 and len(re.findall(r'L', species.to_smiles())) == 2:
            paired_frag_list = pair_frag(amount, one_L_dict, 'L')
            multi_label_amount += amount
            # reattach
            for frag_amt in paired_frag_list:
                frag = frag_amt[0]
                amt = frag_amt[1]
                frag1 = frag # 1L
                frag2 = species # 1R2L
                frag_new = merge_frag_to_frag(frag1.to_smiles(), frag2.to_smiles(), 'R') # L,R,R

                general_L_list.append((frag_new, amt * 2))
                rr_ll_list.append(frag_new)

        # 3R
        elif len(re.findall(r'R', species.to_smiles())) == 3 and len(re.findall(r'L', species.to_smiles())) == 0:
            paired_frag_list = pair_frag(amount, one_L_dict, 'L')
            multi_label_amount += amount
            # reattach
            for frag_amt in paired_frag_list:
                frag = frag_amt[0]
                amt = frag_amt[1]
                frag1 = frag # 1L
                frag2 = species # 3R
                frag_new = merge_frag_to_frag(frag1.to_smiles(), frag2.to_smiles(), 'R') # L,R,R

                general_R_list.append((frag_new, amt * 2))
                rr_ll_list.append(frag_new)

        # 3L
        elif len(re.findall(r'R', species.to_smiles())) == 0 and len(re.findall(r'L', species.to_smiles())) == 3:
            paired_frag_list = pair_frag(amount, one_R_dict, 'R')
            multi_label_amount += amount
            # reattach
            for frag_amt in paired_frag_list:
                frag = frag_amt[0]
                amt = frag_amt[1]
                frag1 = frag # 1R
                frag2 = species # 3L
                frag_new = merge_frag_to_frag(frag2.to_smiles(), frag1.to_smiles(), 'R') # L,R,R

                general_L_list.append((frag_new, amt * 2))
                rr_ll_list.append(frag_new)

    # print(multi_label_amount)

    for one_R_frag, amt in one_R_dict.items():
        general_R_list.append((one_R_frag, amt))

    for one_L_frag, amt in one_L_dict.items():
        general_L_list.append((one_L_frag, amt))

    # grind amounts and shuffle them
    import afm.utils
    grind_size = 1
    shuffle_seed = 0
    # cut large moles into smaller pieces
    grinded_r_moles = afm.utils.grind(general_R_list, grind_size)
    grinded_l_moles = afm.utils.grind(general_L_list, grind_size)

    # random shuffle
    r_moles_shuffle = afm.utils.shuffle(grinded_r_moles, shuffle_seed)
    l_moles_shuffle = afm.utils.shuffle(grinded_l_moles, shuffle_seed)

    # match concentrations for single-labeled fragments
    # including RCCCCR
    matches0_random = afm.utils.match_concentrations_with_same_sums(l_moles_shuffle, 
                                                                             r_moles_shuffle, 
                                                                             diff_tol=1e-3)

    # modified matches_resolve function 
    new_matches_random = []
    new_r_l_moles_random = []
    for match in matches0_random:
        pair = match[0]
        value = match[1]
        l_frag, r_frag = pair

        if l_frag not in rr_ll_list:
            if r_frag not in rr_ll_list:
                # cases like (L-Y, X-R)
                new_matches_random.append((pair, value))
            else:
                # cases like (L-Y, R-U1-R)
                new_matches_random.append(((l_frag, r_frag, l_frag), value/2.0))
        else:
            if r_frag not in rr_ll_list:
                # cases like (L-W1-L, X-R)
                new_matches_random.append(((r_frag, l_frag, r_frag), value/2.0))
            else:
                # cases like (L-W1-L, R-U1-R)
                new_r_l_moles_random.append((pair, value/2.0))

    r_l_moles.extend(new_r_l_moles_random)
    if r_l_moles:
        matches_random = afm.utils.match_concentrations_with_different_sums(new_matches_random, r_l_moles)
    else:
        matches_random = new_matches_random

    ## Calculate the MW of reattached fragmented molecules

    ## flatten the tuple set first
    flattened_matches_random = [(tuple(afm.utils.flatten(m[0])), m[1]) for m in matches_random]

    final_frags_moles_random = []

    for non_cut_mole, val in moles_remain:
        final_frags_moles_random.append(((non_cut_mole.molecule[0], ), val))

    final_frags_moles_random.extend(flattened_matches_random)

    # reattach fragment pairs to get molecules
    correct_reattach = []
    for frag_amt in final_frags_moles_random:
        frag_tup, amt = frag_amt
        if len(frag_tup) == 6: # 1L,2R,1L,1R1L,2L,2R
            frag1 = frag_tup[0].to_smiles()
            frag2 = frag_tup[1].to_smiles()
            frag3 = frag_tup[2].to_smiles()
            frag4 = frag_tup[3].to_smiles()
            frag5 = frag_tup[4].to_smiles()
            frag6 = frag_tup[5].to_smiles()
            #reattach frag 56 and then frag4 then frag3 then frag 2 then frag 1
            if 'L' in frag1:
                frag56 = merge_frag_to_frag(frag5, frag6, 'R') #LRR
                frag456 = merge_frag_to_frag(frag4, frag56.to_smiles(), 'R')
                frag3456 = merge_frag_to_frag(frag3, frag456.to_smiles(), 'R')
                frag23456 = merge_frag_to_frag(frag3456.to_smiles(), frag2, 'R')
                frag = merge_frag_to_frag(frag1, frag23456.to_smiles(), 'R')
            else: # 1R,2L,1R,1R1L,2L,2R
                frag56 = merge_frag_to_frag(frag5, frag6, 'R') #LRR
                frag456 = merge_frag_to_frag(frag4, frag56.to_smiles(), 'R')
                frag3456 = merge_frag_to_frag(frag3, frag456.to_smiles(), 'L') #RLL
                frag23456 = merge_frag_to_frag(frag3456.to_smiles(), frag2, 'L')
                frag = merge_frag_to_frag(frag1, frag23456.to_smiles(), 'L')            

        elif len(frag_tup) == 5: # 1R,2L,1R,2L,2R or 1L,2R,1L,2L,2R
            frag1 = frag_tup[0].to_smiles()
            frag2 = frag_tup[1].to_smiles()
            frag3 = frag_tup[2].to_smiles()
            frag4 = frag_tup[3].to_smiles()
            frag5 = frag_tup[4].to_smiles()
            #reattach frag 45 and then frag3 then frag 2 then frag 1
            if 'L' in frag1:
                frag45 = merge_frag_to_frag(frag4, frag5, 'R') #LRR
                frag345 = merge_frag_to_frag(frag3, frag45.to_smiles(), 'R')
                frag2345 = merge_frag_to_frag(frag345.to_smiles(), frag2, 'R')
                frag = merge_frag_to_frag(frag1, frag2345.to_smiles(), 'R')
            else:
                frag45 = merge_frag_to_frag(frag4, frag5, 'R') #LRR
                frag345 = merge_frag_to_frag(frag3, frag45.to_smiles(), 'L') #RLL
                frag2345 = merge_frag_to_frag(frag345.to_smiles(), frag2, 'L') #RLL
                frag = merge_frag_to_frag(frag1, frag2345.to_smiles(), 'L')

        elif len(frag_tup) == 4: # 1L,2R,1L,1R1L or 1R,2L,1R,1R1L
            frag1 = frag_tup[0].to_smiles()
            frag2 = frag_tup[1].to_smiles()
            frag3 = frag_tup[2].to_smiles()
            frag4 = frag_tup[3].to_smiles()
            #reattach frag 34 and then frag 2 then frag 1
            if 'L' in frag1:
                frag34 = merge_frag_to_frag(frag3, frag4, 'R')
                frag234 = merge_frag_to_frag(frag34.to_smiles(), frag2, 'R')
                frag = merge_frag_to_frag(frag1, frag234.to_smiles(), 'R')
            else:
                frag34 = merge_frag_to_frag(frag3, frag4, 'L')
                frag234 = merge_frag_to_frag(frag34.to_smiles(), frag2, 'L')
                frag = merge_frag_to_frag(frag1, frag234.to_smiles(), 'L')
        elif len(frag_tup) == 3: # 1L,1R,1R1L
            frag1 = frag_tup[0].to_smiles()
            frag2 = frag_tup[1].to_smiles()
            frag3 = frag_tup[2].to_smiles()
            if 'L' in frag1:
                frag23 = merge_frag_to_frag(frag2, frag3, 'L')
                frag = merge_frag_to_frag(frag1, frag23.to_smiles(), 'R')
            else:
                frag23 = merge_frag_to_frag(frag2, frag3, 'R') #LRR
                frag = merge_frag_to_frag(frag1, frag23.to_smiles(), 'L') #RLL
        elif len(frag_tup) == 2: # 1L,1R
            frag1 = frag_tup[0].to_smiles()
            frag2 = frag_tup[1].to_smiles()
            frag = merge_frag_to_frag(frag1, frag2, 'R')
        elif len(frag_tup) == 1: # molecule which does not need cut
            frag = frag_tup[0]
        else:
            print('sth wrong')
            frag_smi = [frag.to_smiles() for frag in frag_tup]
            print(frag_smi, amt)
        correct_reattach.append((frag,amt))

    # collect same frag combination tuples to speed up in frag merge & compare isomorphic
    molecule_dict = {}
    mole_frags_dict = {}
    for ind, mol_amt in enumerate(correct_reattach):
        molecule, amt = mol_amt
        check = 0
        for mole,_ in molecule_dict.items():
            if molecule.is_isomorphic(mole):
                molecule_dict[mole] += amt
    #             if mole_frags_dict == {}:
    #                 mole_frags_dict[molecule].append(ind)
                mole_frags_dict[mole].append(ind)
                check += 1
                break
        if check == 0:
            molecule_dict[molecule] = amt
            mole_frags_dict[molecule]=[ind]

    # build pattern list (SMARTS pattern rxn) and then cut
    non_rad_spe = []
    for spe in core_species:
        if spe.molecule[0].is_radical() != True and spe.molecule[0].is_isomorphic(Molecule().from_smiles('[H][H]')) == False:
            smi = spe.molecule[0].to_smiles()
            if '[' not in smi and '#' not in smi:
                non_rad_spe.append(spe)
                # display(spe)

    # sort the species based on C number
    non_rad_spe.sort(key= lambda spe: spe.molecule[0].get_element_count()['C'], reverse = True)

    # create pattern list
    pattern_list = []
    smiles_list = []

    for spe in non_rad_spe:
        smiles_list.append(spe.molecule[0].to_smiles())
        if isinstance(spe.molecule[0], Fragment) and re.findall(r'[LR]', spe.molecule[0].to_smiles()) != []:
            pattern_list.append(spe)

    # make smarts rxn for each element in pattern list
    pattern_arom_list = []
    pattern_aliph_list = []
    for pattern in pattern_list:

        pattern_smiles = pattern.molecule[0].to_smiles()
        if pattern.molecule[0].is_aromatic():
            pattern_arom_list.append(pattern_smiles)
        else:
            pattern_aliph_list.append(pattern_smiles)

    # create extended_smiles_dict to have easier comparison between same species smiles
    extended_smiles_dict = {}
    for smile in smiles_list:
        # change smiles (replace cutting label to metal element) to smiles in rdkit
        f = smile.replace('R', '[Na]')
        ff = f.replace('L', '[Mg]')
        frag = Chem.MolFromSmiles(ff)
        extended_smiles = []
        for i in range(1,20):
            new_smiles = Chem.MolToSmiles(frag, canonical=False, doRandom=True)
            # change back cutting labels
            new_f = new_smiles.replace('[Na]', 'R')
            new_frag = new_f.replace('[Mg]', 'L')
            if new_frag not in extended_smiles:
                extended_smiles.append(new_frag)
                # print(new_smiles)

        # create a dict to reference the extended smiles to the specific smiles in smiles_list
        # key: extended smiles of smile
        # value: smile
        for ext_smile in extended_smiles:
            extended_smiles_dict[ext_smile] = smile

    # generate_smarts_rxn
    smarts_arom_list = []
    for pattern_smiles in pattern_arom_list:
        return_list = generate_smarts_rxn(pattern_smiles)
        for (smarts_rxn, mol2) in return_list:
            smarts_arom_list.append((smarts_rxn, mol2))
    smarts_aliph_list = []
    for pattern_smiles in pattern_aliph_list:
        return_list = generate_smarts_rxn(pattern_smiles)
        for (smarts_rxn, mol2) in return_list:
            smarts_aliph_list.append((smarts_rxn, mol2))

    # perform cut on molecules
    import tqdm
    mole_new_frag_dict = {}
    for mol, amt in tqdm.tqdm(molecule_dict.items()):
        mol.reactive = True
        m = mol.generate_resonance_structures()[0]
        if m.is_aromatic():
            frag_smiles_list = []
            arom_cut_frag = sliceitup_arom(m.to_smiles(), smarts_arom_list)
            for frag in arom_cut_frag:
                if frag in smiles_list:
                    frag_smiles_list.append(frag)
                elif 'c' in frag:
                    # aromatic fragment, need to slice_arom
                    new_frag_cut = sliceitup_arom(frag, smarts_arom_list)
                    for new_frag in new_frag_cut:
                        arom_cut_frag.append(new_frag)
                else: #frag - aliphatic part will come to here
                    frag_aliph_list = []
                    sets = []
                    aliph_cut_frag = sliceitup_aliph(frag, smarts_aliph_list)
                    sets.append(aliph_cut_frag)
                    for ind, frag_sets in enumerate(sets):
                        if all(new_frag in smiles_list for new_frag in frag_sets):
                            for subfrag in frag_sets:
                                arom_cut_frag.append(subfrag)
                            break
                        else:
                            # check each frag to see whether there can be more cut on one frag
                            middle_smiles = []
                            for single_frag in frag_sets:
                                if single_frag in smiles_list:
                                    middle_smiles.append(single_frag)                        
                                else:
                                    # aliphatic fragment, need to slice_aliph
                                    new_frag_cut = sliceitup_aliph(single_frag, smarts_aliph_list)
                                    if len(new_frag_cut) == 1:
                                        pas = 0
                                        count = 0
                                        while pas == 0:
                                            shuffled_smarts_aliph = random.sample(smarts_aliph_list, len(smarts_aliph_list))
                                            new_frag_cut = sliceitup_aliph(single_frag, shuffled_smarts_aliph, limit_length=False)
                                            count += 1
                                            if all(new_frag in smiles_list for new_frag in new_frag_cut):
                                                pas = 1
                                            elif count == 1000:
                                                # probably wrong cut at the very beginning, need to shuffle smarts_aliph and cut from beginning
                                                break
                                    for new_frag in new_frag_cut:
                                        middle_smiles.append(new_frag)
                            if all(new_frag in smiles_list for new_frag in middle_smiles):
                                sets.append(middle_smiles)
                            else:
                                # probably wrong cut at the very beginning, need to shuffle smarts_aliph and cut from beginning  
                                shuffled_smarts_aliph = random.sample(smarts_aliph_list, len(smarts_aliph_list))
                                new_ali_cut = sliceitup_aliph(frag, shuffled_smarts_aliph, limit_length=False)
                                sets.append(new_ali_cut)

        else:
            frag_smiles_list = []
            sets = []
            aliph_cut_frag = sliceitup_aliph(m.to_smiles(), smarts_aliph_list, limit_length=False)
            sets.append(aliph_cut_frag)
            ccount = 0
            past_cut_set = []
            past_count = 0
            new_frag_cut = []
            for ind, frag_sets in enumerate(sets):
                if all(new_frag in smiles_list for new_frag in frag_sets):
                    for frag in frag_sets:
                        frag_smiles_list.append(frag)
                    past_cut_set = []
                    past_count = 0
                    break
                else:
                    # check each frag to see whether there can be more cut on one frag
                    middle_smiles = []
                    for frag in frag_sets:
                        if frag in smiles_list:
                            middle_smiles.append(frag)
                        # replace cutting label on each frag to see if they are in smiles_list
                        elif middle_smiles != [] and ccount >= 1000: 
                            switch_label = [frag]
                            for mid in middle_smiles:
                                switch_label.append(mid)
                            # switched label R with L & L with R
                            switched_label = []
                            for f in switch_label:
                                import re
                                p = re.compile("[RL]")
                                ind_cuttinglabel_list = []
                                replace_label = {'R':'L', 'L':'R'}
                                for mat in p.finditer(f):
                                    ind_cuttinglabel_list.append((mat.start(), mat.group()))

                                for (index, label) in ind_cuttinglabel_list:
                                    f = f[:index] + replace_label[label] + f[index + 1:]
                                switched_label.append(f)

                            # check if all frag are in smiles_list
                            if all(new_f in list(extended_smiles_dict.keys()) for new_f in switched_label):
                                # find matched frag smiles in smiles_list
                                switched_label[:] = [extended_smiles_dict[new_f] for new_f in switched_label]
                                middle_smiles = switched_label
                            ccount = 0
                        else:
                            # aliphatic fragment, need to slice_aliph
                            if new_frag_cut != []:
                                past_cut_set = new_frag_cut
                            new_frag_cut = sliceitup_aliph(frag, smarts_aliph_list, limit_length=False)
                            if new_frag_cut == past_cut_set:
                                past_count += 1
                            if len(new_frag_cut) == 1:
                                pas = 0
                                count = 0
                                while pas == 0:
                                    shuffled_smarts_aliph = random.sample(smarts_aliph_list, len(smarts_aliph_list))
                                    new_frag_cut = sliceitup_aliph(frag, shuffled_smarts_aliph, limit_length=False)
                                    count += 1
                                    if all(new_frag in list(extended_smiles_dict.keys()) for new_frag in new_frag_cut):
                                        # find matched frag smiles in smiles_list
                                        new_frag_cut[:] = [extended_smiles_dict[new_frag] for new_frag in new_frag_cut]

                                        pas = 1
                                    elif count == 1000:
                                        # probably wrong cut at the very beginning, need to shuffle smarts_aliph and cut from beginning
                                        break
                            elif new_frag_cut == past_cut_set and past_count >= 5:
                                # means it might be stuck in the middle smiles cutting
                                # should try with shuffles smarts aliph
                                shuffled_smarts_aliph = random.sample(smarts_aliph_list, len(smarts_aliph_list))
                                new_frag_cut = sliceitup_aliph(frag, shuffled_smarts_aliph, limit_length=False)

                            for new_frag in new_frag_cut:
                                middle_smiles.append(new_frag)
                    if all(new_frag in smiles_list for new_frag in middle_smiles):
                        sets.append(middle_smiles)
                    # try switch cutting labels
                    elif any(new_frag in smiles_list for new_frag in middle_smiles):
                        switch_label = middle_smiles
                        # switched label R with L & L with R
                        switched_label = []
                        for f in switch_label:
                            import re
                            p = re.compile("[RL]")
                            ind_cuttinglabel_list = []
                            replace_label = {'R':'L', 'L':'R'}
                            for mat in p.finditer(f):
                                ind_cuttinglabel_list.append((mat.start(), mat.group()))

                            for (index, label) in ind_cuttinglabel_list:
                                f = f[:index] + replace_label[label] + f[index + 1:]
                            switched_label.append(f)

                        # check if all frag are in smiles_list
                        if all(new_f in list(extended_smiles_dict.keys()) for new_f in switched_label):
                            # find matched frag smiles in smiles_list
                            switched_label[:] = [extended_smiles_dict[new_f] for new_f in switched_label]
                            sets.append(switched_label)
                        else:
                            shuffled_smarts_aliph = random.sample(smarts_aliph_list, len(smarts_aliph_list))
                            new_ali_cut = sliceitup_aliph(m.to_smiles(), shuffled_smarts_aliph, limit_length=False)
                            sets.append(new_ali_cut)

                    else:
                        # probably wrong cut at the very beginning, need to shuffle smarts_aliph and cut from beginning  
                        ccount += 1
                        shuffled_smarts_aliph = random.sample(smarts_aliph_list, len(smarts_aliph_list))
                        new_ali_cut = sliceitup_aliph(m.to_smiles(), shuffled_smarts_aliph, limit_length=False)
                        sets.append(new_ali_cut)



        mole_new_frag_dict[mol.to_smiles()] = (frag_smiles_list, amt)


    # check if there is atom inbalance
    for mol, amt in mole_new_frag_dict.items():
        before = mol.count('C')+mol.count('c')
        frag_smilist = amt[0]
        after = 0
        for frag in frag_smilist:
            after += frag.count('C') + frag.count('c')
        if before != after:
            print(mol)
            print('has atom inbalance after cutting')

    # calculate new total moles and get new self.y
    new_total_moles = 0
    for mol, frags_amt in mole_new_frag_dict.items():
        new_frags = frags_amt[0]
        amt = frags_amt[1]
        new_total_moles += len(new_frags) * amt

    mol_input_dict = {}

    for mol, frags_amt in mole_new_frag_dict.items():
        new_frags = frags_amt[0]
        amt = frags_amt[1]
    #     mf = amt
        mf = amt/new_total_moles
        for frag in new_frags:
            if frag in mol_input_dict.keys():
                mol_input_dict[frag]+= mf
            else:
                mol_input_dict[frag]= mf

    # check new input system molar fraction
    total_amt = 0
    for frags, amt in mol_input_dict.items():
        total_amt += amt
    # print(total_amt)
    if abs(total_amt-1.0) >= 1e-3:
        print('molar fraction sum is not 1.0')

    ### change the dict to np.array as new self.y ###
    new_y = np.zeros(len(core_species))
    for ind, spe_smiles in enumerate(core_smiles):
        if spe_smiles in mol_input_dict:
            new_y[ind] = mol_input_dict[spe_smiles]

    return new_y

