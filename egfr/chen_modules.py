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
    Monomer('ERBB1', ['bf'])
    Monomer('ERBB2', ['bf'])
    Monomer('ERBB3', ['bf'])
    Monomer('ERBB4', ['bf'])
    Monomer('DEP1', ['bf'])
    Monomer('DEP2', ['bf'])
    Monomer('DEP3', ['bf'])
    Monomer('DEP4', ['bf'])


def rec_events():
    """ TEXT HERE
    """
    
