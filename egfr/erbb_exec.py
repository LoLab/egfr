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
chen_modules.akt_monomers()
chen_modules.crosstalk_monomers()

# Generate the upstream and downstream sections
chen_modules.rec_events()
chen_modules.mapk_events()
chen_modules.akt_events()
chen_modules.crosstalk_events()

# Initial protein concentrations
chen_modules.rec_initial()
chen_modules.mapk_initial()
chen_modules.akt_initial()
chen_modules.crosstalk_initial()

# Declare observables
Observable('ErbB1_ErbB1', erbb(bd=1, ty='1', st='U', loc='C') % erbb(bd=1, ty='1', st='U', loc='C'))
Observable('EGF_any', EGF(b=ANY))
Observable('obsPIP3', PIP(bakt=None, both=None, bpi3k=None, S='PIP3'))
Observable('obsPTEN', PTEN(bpip3=None))
Observable('obsSHP', SHP(bpip3=None))
Observable('obsPIPPTEN', PIP(bakt=None, both=1, S='PIP3') % PTEN(bpip3=1))
Observable('obsPIPSHP', PIP(bakt=None, both=1, S='PIP3') % SHP(bpip3=1))
Observable('obsPDK1', PDK1(bakt=None, both=None))
Observable('obsAKT', AKT(bpip3=None, both=None, S='U'))
Observable('obsAKTPIP', AKT(bpip3=1, both=None, S='U') % PIP(bakt=1, both=None, S='PIP3'))
Observable('obsAKTPDK1PIP', AKT(bpip3=None, both=1, S='U') % PDK1(bakt=None, both=2) % PIP(bakt=None, both=1, S='PIP3'))
Observable('obsAKTP', AKT(bpip3=None, both=None, S='P'))
Observable('obsAKTPPIP', AKT(bpip3=1, both=None, S='P') % PIP(bakt=1, both=None, S='PIP3'))
Observable('obsPIP2', PIP(bakt=None, both=None, bpi3k=None, S='PIP2'))
Observable('obsPP2A_III', PP2A_III(bakt=None))
Observable('obsAKTPP', AKT(bpip3=None, both=None, S='PP'))
Observable('obsPP2A_IIIAKTPP', AKT(bpip3=None, both=1, S='PP') % PP2A_III(bakt=1))
Observable('obsAKTPPDK1PIP', AKT(bpip3=1, both=None, S='P') % PIP(bakt=1, both=2, S='PIP3') % PDK1(both=2))
Observable('obsGAB1_unbound', GAB1(bgrb2=None, bshp2=None, bpi3k=None, batp=None,bERKPP=None,bPase9t=None,S='U'))
Observable('obsGAB1_bound', GAB1(bgrb2=ANY, bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, S='U'))
Observable('obsGAB1P', GAB1(bgrb2=ANY, bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, S='P'))
Observable('obsATP', ATP(b=None))
Observable('obsGAPGRB2SOS', GRB2(b=None, bgap=ANY, bsos=ANY, bgab1=None))
Observable('obsGAPGRB2SOSRASGDP', GRB2(bgap=ANY, bsos=ANY) % SOS(bras=ANY, bgrb=ANY) % RAS(braf=None, bsos=ANY, st='GDP'))
Observable('obsRASGTP', RAS(braf=None, bsos=None, st='GTP'))
Observable('obsGAPSHCPGRB2SOSRASGDP', GRB2(b=ANY, bsos=ANY) % SOS(bras=ANY, bgrb=ANY) % RAS(braf=None, bsos=ANY, st='GDP'))
Observable('obsGAPGRB2bub', GRB2(bgap=ANY))
Observable('obsSHCPGRB2', GRB2(b=ANY))
Observable('obsGRB2ub', GRB2(bgap=None, b=None))
Observable('obsRAFPserP', RAF(st='P', ser295='P'))
Observable('obsGAB1PP', GAB1(S='PP'))
Observable('obsSOSP', SOS(st='P'))
Observable('obsRASGTPPI3K', RAS(bpi3k=ANY))
