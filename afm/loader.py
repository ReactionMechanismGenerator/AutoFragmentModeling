import os

from rmgpy.kinetics import Arrhenius
from rmgpy.chemkin import loadChemkinFile

from afm.fragment import Fragment
from afm.reaction import FragmentReaction

def load_fragment_reactions_from_chemkin(chemkin_path,
                                        dictionary_path):
    """
    This method loads chemkin mechanism and 
    generate fragment reactions in irreversible
    format.
    """
    speciesList, reactionList = loadChemkinFile(chemkin_path, dictionary_path)

    species_dict = {}
    for spe in speciesList:
        label = spe.label
        if label not in species_dict:
            species_dict[label] = spe
        else:
            raise Exception('Duplicate species found with label {0}.'.format(label))
    


    fragments_dict = {}

    for label, spec in species_dict.iteritems():
        fragments_dict[label] = Fragment(label=label,species_repr=spec)

    orig_fragrxns = []
    for rxn0 in reactionList:
        fragrxts = [fragments_dict[spec.label] for spec in rxn0.reactants]
        
        fragprds = [fragments_dict[spec.label] for spec in rxn0.products]
        
        fragpairs = [(fragments_dict[rxt.label], fragments_dict[prod.label]) for rxt, prod in rxn0.pairs]
        fragrxn = FragmentReaction(index=-1,
                                    reactants=fragrxts,
                                    products=fragprds,
                                    kinetics=rxn0.kinetics,
                                   reversible=False,
                                    pairs=fragpairs,
                                    family=rxn0.family)
        orig_fragrxns.append(fragrxn)


    revs_fragrxns = []
    for rxn0 in reactionList:
        fragrxts = [fragments_dict[spec.label] for spec in rxn0.products]
        
        fragprds = [fragments_dict[spec.label] for spec in rxn0.reactants]
        
        fragpairs = [(fragments_dict[prod.label], fragments_dict[rxt.label]) for rxt, prod in rxn0.pairs]
        
        revs_kinetics = rxn0.generateReverseRateCoefficient()
        fragrxn = FragmentReaction(index=-1,
                                    reactants=fragrxts,
                                    products=fragprds,
                                    kinetics=revs_kinetics,
                                   reversible=False,
                                    pairs=fragpairs,
                                    family=rxn0.family)
        revs_fragrxns.append(fragrxn)

    return fragments_dict, orig_fragrxns + revs_fragrxns

def load_pseudo_fragment_reactions(fragments_dict):
    """
    Currently only returns a pseudo reaction. It can be
    extended to generate multiple ractions in the future.
    """

    pseudo_fragrxts = ['RC*C__C', 'RCCCCR']
    pseudo_fragprds = ['RCCCCC__CC*']
    pseudo_frag_pairs = [('RC*C__C', 'RCCCCC__CC*'), ('RCCCCR', 'RCCCCC__CC*')]

    fragrxts = [fragments_dict[label] for label in pseudo_fragrxts]
        
    fragprds = [fragments_dict[label] for label in pseudo_fragprds]

    fragpairs = [(fragments_dict[rxt_label], fragments_dict[prod_label]) for rxt_label, prod_label in pseudo_frag_pairs]

    pseudo_kinetics = Arrhenius(A=(2.000e+05, 'cm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kcal/mol'), T0=(1, 'K'))

    pseudo_fragrxn = FragmentReaction(index=-1,
                                reactants=fragrxts,
                                products=fragprds,
                                kinetics=pseudo_kinetics,
                                reversible=False,
                                pairs=fragpairs,
                                family='pseudo_rxn')

    return [pseudo_fragrxn]