"""
Overview
========

PySB implementation of the ErbB related MAPK and AKT signaling
pathways originally published in [Chen2009]_.

This file containst functions that implement the ErbB execution
pathway in three modules:

- Receptor layer events, taking into account all ErbB1-4 interactions with ligand.
- AKT pathway
- MAPK pathway

"""

from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components
#from egfr.shared import * # modified model aliases

# Receptor Layer

# Default rates?
KF = 1e-6
KR = 1e-3
KC = 1
KDIMF = 1e-6
KDIMR = 1e-3
KINTF = 1.0e-3
KINTR = 5.0e-5
KDEG = .1

# Monomer declarations
# ====================

def erbb_reclayer_monomers():
    """ Declares the ErbB receptor interactions.
    'bf' is the default site to be used for all binding/catalysis reactions.
    """
    Monomer('EGF', ['b']) # Epidermal Growth Factor ligand
    Monomer('HRG', ['b']) # Heregulin ligand
    Monomer('erbb', ['bl', 'bd', 'ba', 'ty', 'st', 'loc'], {'ty':['1','2','3','4'], 'st':['U','P'], 'loc':['C','E']}) # bl: lig, bd: dimer, ba: atp, ty: rec type, st: (U)n(P)hosphorylated, loc: (C)yto 'brane or (E)ndosome 'brane
    Monomer('DEP', ['b'])
    Monomer('ATP', ['b'])

def rec_events():
    """ TEXT HERE
    """

    # Parameter definitions
    # =====================


    # Alias model components for names in present namespace
    alias_model_components()
    
    # binding to receptors
    bind_table([[                    EGF,       HRG],
                [erbb(ty='1', loc='C'),  (1.0,1.0),      None],
                [erbb(ty='3', loc='C'),       None, (1.0,1.0)],
                [erbb(ty='4', loc='C'),       None, (1.0,1.0)]],
                'bl', 'b')
    
    # erbb dimerization
    bind_table([[                       erbb(ty='1', loc='C'), erbb(ty='2', loc='C'), erbb(ty='3', loc='C'), erbb(ty='4', loc='C')],
                [erbb(ty='1', loc='C'),        (KDIMF, KDIMR),                  None,                  None,                  None],
                [erbb(ty='2', loc='C'),        (KDIMF, KDIMR),        (KDIMF, KDIMR),                  None,                  None],
                [erbb(ty='3', loc='C'),        (KDIMF, KDIMR),        (KDIMF, KDIMR),        (KDIMF, KDIMR),                  None],
                [erbb(ty='4', loc='C'),        (KDIMF, KDIMR),        (KDIMF, KDIMR),        (KDIMF, KDIMR),        (KDIMF, KDIMR)]],
        'bd', 'bd')

    # ATP binding
    bind_table([[                            ATP],
                [erbb(ty='1', loc='C'), (KF, KR)],
                [erbb(ty='2', loc='C'), (KF, KR)],
                [erbb(ty='4', loc='C'), (KF, KR)]],
        'ba', 'b')
                
    # Receptor Cross Phosphorylation
    for i in ['1','2','4']:
        for j in ['1','2','4']:
            Rule("cross_phospho_"+i+"_"+j,
                erbb(ty=i, bd=1, ba=2) % erbb(ty=j, bd=1, st='U') >>
                erbb(ty=i, bd=1, ba=2) % erbb(ty=j, bd=1, st='P'),
                Parameter("kc"+i+j, KC))

    # Receptor Dephosphorylation
    Rule("dephospho",
         erbb(st='P') + DEP(


    # Receptor internalization
    # This internalizes all receptor combos 
    Rule("rec_intern",
         erbb(loc="C") <> erbb(loc="E"),
         Parameter("kintf", KINTF), Parameter("kintr", KINTR))

    # Receptor degradation
    # This degrades all receptor combis within an endosome
    degrade(erbb(loc="E"), Parameter("kdeg", KDEG))
