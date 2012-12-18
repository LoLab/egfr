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

# Receptor Layer, 0

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

def rec_monomers():
    """ Declares the ErbB receptor interactions.
    'bf' is the default site to be used for all binding/catalysis reactions.
    """
    Monomer('EGF', ['b']) # Epidermal Growth Factor ligand
    Monomer('HRG', ['b']) # Heregulin ligand
    Monomer('erbb', ['bl', 'bd', 'bp', 'ty', 'st', 'loc'], {'ty':['1','2','3','4'], 'st':['U','P'], 'loc':['C','E']}) # bl: lig, bd: dimer, bp: atp, ty: rec type, st: (U)n(P)hosphorylated, loc: (C)yto 'brane or (E)ndosome 'brane
    Monomer('DEP', ['b'])
    Monomer('ATP', ['b'])

def rec_initial():
    """
    STATE WHERE PARAMS CAME FROM
    """
    Parameter('EGF_0',   1000)
    Parameter('HRG_0',   1000)
    Parameter('erbb1_0', 1000)
    Parameter('erbb2_0', 1000)
    Parameter('erbb3_0', 1000)
    Parameter('erbb4_0', 1000)
    Parameter('DEP_0',   1000)
    Parameter('ATP_0',   1000)

    alias_model_components()

    Initial(EGF(b=None), EGF_0)
    Initial(HRG(b=None), HRG_0)
    Initial(erbb(bl=None, bd=None, bp=None, ty='1', st='U', loc='C'), erbb1_0)
    Initial(erbb(bl=None, bd=None, bp=None, ty='2', st='U', loc='C'), erbb2_0)
    Initial(erbb(bl=None, bd=None, bp=None, ty='3', st='U', loc='C'), erbb3_0)
    Initial(erbb(bl=None, bd=None, bp=None, ty='4', st='U', loc='C'), erbb4_0)
    Initial(DEP(b=None), DEP_0)
    Initial(ATP(b=None), ATP_0)
            
def rec_events():
    """ Describe receptor-level events here. 
    """

    # Parameter definitions
    # =====================
    # Alias model components for names in present namespace
    alias_model_components()
    
    # binding to receptors
    bind_table([[                              EGF,       HRG],
                [erbb(ty='1', loc='C'),  (1.0,1.0),      None],
                [erbb(ty='3', loc='C'),       None, (1.0,1.0)],
                [erbb(ty='4', loc='C'),       None, (1.0,1.0)]],
                'bl', 'b')
    
    # make this more readable:
    
    # erbb dimerization
    bind_table([[                       erbb(ty='1', loc='C'), erbb(ty='2', loc='C'), erbb(ty='3', loc='C'), erbb(ty='4', loc='C')],
                [erbb(ty='1', loc='C'),        (KDIMF, KDIMR),                  None,                  None,                  None],
                [erbb(ty='2', loc='C'),        (KDIMF, KDIMR),        (KDIMF, KDIMR),                  None,                  None],
                [erbb(ty='3', loc='C'),        (KDIMF, KDIMR),        (KDIMF, KDIMR),                  None,                  None],
                [erbb(ty='4', loc='C'),        (KDIMF, KDIMR),        (KDIMF, KDIMR),                  None,                  None]],
        'bd', 'bd')

    # ATP binding
    bind_table([[                            ATP],
                [erbb(ty='1', loc='C'), (KF, KR)],
                [erbb(ty='2', loc='C'), (KF, KR)],
                [erbb(ty='4', loc='C'), (KF, KR)]],
        'bp', 'b')
                
    # Receptor Cross Phosphorylation
    for i in ['1','2','4']:
        for j in ['1','2','4']:
            Rule("cross_phospho_"+i+"_"+j,
                erbb(ty=i, bd=1, bp=2) % erbb(ty=j, bd=1, st='U') >>
                erbb(ty=i, bd=1, bp=2) % erbb(ty=j, bd=1, st='P'),
                Parameter("kcp"+i+j, KC))

    # Receptor Dephosphorylation
    # DEPHOSPHORYLATION: 
    #  * Density enhanced phosphatase1 (DEP1) dephosphorylates ERB1 (at the cell-membrane)
    #  * Protein Tyrosine Phosphatase1b (PTP1b) dephosphorylates all RTKs (at the endo-membrane)
    #  Bursett, TA, Hoier, EF, Hajnal, A: Genes Dev. 19:1328-1340 (2005)
    #  Haj, FG, Verver, PJ, Squire, A, Neel, BG, Bastiaens, PI: Science 295:1708-1711 (2002)
    for i in ['1','2','3','4']:
        Rule("dephospho_"+i+"_"+j,
             erbb(st='P', bp=None, ty=i) + DEP(b=None) <>
             erbb(st='P', bp=1, ty=i) % DEP(b=1),
             Parameter('kfdephos'+i+j, KF),Parameter('krdephos'+i+j, KR))
        Rule("dephosphoC_"+i+"_"+j,
             erbb(st='P', bp=1, ty=i) % DEP(b=1) >>
             erbb(st='U', bp=None, ty=i) + DEP(b=None),
             Parameter('kcdephos'+i+j, KC))
        
    # Receptor internalization
    # This internalizes all receptor combos 
    Rule("rec_intern",
         erbb(bd=1, loc="C") % erbb(bd=1, loc='C') <> erbb(bd=1, loc="E") % erbb(bd=1, loc="E"),
         Parameter("kintf", KINTF), Parameter("kintr", KINTR))

    # Receptor degradation
    # This degrades all receptor combos within an endosome
    degrade(erbb(bd=1, loc="E") % erbb(bd=1, loc="E"), Parameter("kdeg", KDEG))

def mapk_monomers():
    Monomer('SHC', [b])
    Monomer('SHCPase', [b])
    Monomer('GRB2', [b])
    Monomer('SOS', [b])
    Monomer('RAS', [b])
    Monomer('RAF', [b])
    Monomer('MEK', [b])
    Monomer('ERK', [b])
    

def mapk_events():

    # =====================
    # Alias model components for names in present namespace
    alias_model_components()

    
