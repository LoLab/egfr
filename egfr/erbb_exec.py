"""
Model M1a: Extrinsic apoptosis model with "embedded together" model of MOMP.
"""


from pysb import *
#from egfr import shared
from egfr import chen_modules

Model()

# Declare monomers
chen_modules.rec_monomers()
chen_modules.mapk_monomers()

# Generate the upstream and downstream sections
chen_modules.rec_events()
chen_modules.mapk_events()

# Initial protein concentrations
chen_modules.rec_initial()
chen_modules.mapk_initial()

# Declare observables
#shared.observables()

