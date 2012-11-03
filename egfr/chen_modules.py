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


# Monomer declarations
# ====================

def erbb_reclayer_monomers():
    """ Declares the ErbB receptor interactions.
    'bf' is the default site to be used for all binding/catalysis reactions.
    """
    Monomer('EGF', ['b']) # Epidermal Growth Factor ligand
    Monomer('HRG', ['b']) # Heregulin ligand
    Monomer('erbb', ['bl', 'bd', 'ba', 'ty', 'st'], {'ty':['1','2','3','4'], 'st':['I','A']}) # bl: ligand, bd: dimerization, ba: atp, ty: receptor type, st: state
    Monomer('DEP', ['bf'])
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
                [erbb(ty='1'),  (1.0,1.0),      None],
                [erbb(ty='3'),       None, (1.0,1.0)],
                [erbb(ty='4'),       None, (1.0,1.0)]],
                'bl', 'b')
    
    # erbb dimerization
    bind_table([[           erbb(ty='1'),   erbb(ty='2'),   erbb(ty='3'),   erbb(ty='4')],
                [erbb(ty='1'), (KDIMF, KDIMR), (KDIMF, KDIMR), (KDIMF, KDIMR), (KDIMF, KDIMR)],
                [erbb(ty='2'), (KDIMF, KDIMR), (KDIMF, KDIMR), (KDIMF, KDIMR), (KDIMF, KDIMR)],
                [erbb(ty='3'), (KDIMF, KDIMR), (KDIMF, KDIMR),           None, (KDIMF, KDIMR)],
                [erbb(ty='4'), (KDIMF, KDIMR), (KDIMF, KDIMR), (KDIMF, KDIMR), (KDIMF, KDIMR)]],
        'bd', 'bd')

    # ATP binding
    bind_table([[                  ATP],
                [erbb(ty='1'), (KF, KR)],
                [erbb(ty='2'), (KF, KR)],
                [erbb(ty='4'), (KF, KR)]],
        'ba', 'b')
                
    # Receptor Cross Phosphorylation
    for i in ['1','2','4']:
        for j in ['1','2','4']:
            Rule("cross_phospho_"+i+"_"+j,
                erbb(ty=i, bd=1, ba=2) % erbb(ty=j, bd=1, st='I') >>
                erbb(ty=i, bd=1, ba=2) % erbb(ty=j, bd=1, st='A'),
                Parameter("kc"+i+j, KC))
            
    # Receptor internalization
    

            
