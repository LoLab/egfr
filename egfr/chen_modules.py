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
from pysb.util import alias_model_components
from egfr.shared import * # modified model aliases


# Receptor Layer

# Default rates?

# Monomer declarations
# ====================

def erbb_reclayer_monomers():
    """ Declares the ErbB receptor interactions.
    'bf' is the default site to be used for all binding/catalysis reactions.
    """
    Monomer('EGF', ['bf']) # Epidermal Growth Factor ligand
    Monomer('HRG', ['bf']) # Heregulin ligand
    Monomer('erbb1', ['bf'])
    Monomer('erbb2', ['bf'])
    Monomer('erbb3', ['bf'])
    Monomer('erbb4', ['bf'])
    Monomer('DEP1', ['bf'])
    Monomer('DEP2', ['bf'])
    Monomer('DEP3', ['bf'])
    Monomer('DEP4', ['bf'])


def rec_events():
    """ TEXT HERE
    """
    # binding to receptors
    bind_table([[              EGF,       HRG],
                [erbb1,  (1.0,1.0),      None],
                [erbb3,       None, (1.0,1.0)],
                [erbb4,       None, (1.0,1.0)],
                ])
    
    # erbb dimerization
    bind_table([[                erbb1,          erbb2,          erbb3,          erbb4],
                [erbb1, (kdimf, kdimr), (kdimf, kdimr), (kdimf, kdimr), (kdimf, kdimr)],
                [erbb2, (kdimf, kdimr), (kdimf, kdimr), (kdimf, kdimr), (kdimf, kdimr)],
                [erbb3, (kdimf, kdimr), (kdimf, kdimr),           None, (kdimf, kdimr)],
                [erbb4, (kdimf, kdimr), (kdimf, kdimr), (kdimf, kdimr), (kdimf, kdimr)]
                ])

    # ATP binding
    bind_table([[            ATP],
                [erbb1, (kf, kr)],
                [erbb2, (kf, kr)],
                [erbb4, (kf, kr)]
                ])
                
    # Receptor cross phosphorylation
    
