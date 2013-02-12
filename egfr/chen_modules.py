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
# from egfr.shared import * # modified model aliases

# Receptor Layer, 0

# Rates obtained from Chen et al. Table 1 pg. 5
KF = 1e-5
KR = 1e-1
KCP = 1e-1
KCD = 1e-2
KDIMF = 1.6e-6
KDIMR = 1.6e-1
KINTF = 1.3e-3
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
    Monomer('erbb', ['bl', 'bd', 'b', 'ty', 'st', 'loc'], {'ty':['1','2','3','4'], 'st':['U','P'], 'loc':['C','E']}) # bl: lig, bd: dimer, b: binding, ty: rec type, st: (U)n(P)hosphorylated, loc: (C)yto 'brane or (E)ndosome 'brane

    Monomer('DEP', ['b'])
    Monomer('ATP', ['b'])
    Monomer('ADP')

def rec_initial():
    """
    STATE WHERE PARAMS CAME FROM
    """
    Parameter('EGF_0',   5e-9)
    Parameter('HRG_0',   5e-9)
    Parameter('erbb1_0', 5e4)
    Parameter('erbb2_0', 5e4)
    Parameter('erbb3_0', 5e4)
    Parameter('erbb4_0', 5e4)
    Parameter('DEP_0',   1e7)
    Parameter('ATP_0',   1.2e9)

    alias_model_components()

    Initial(EGF(b=None), EGF_0)
    Initial(HRG(b=None), HRG_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='1', st='U', loc='C'), erbb1_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='2', st='U', loc='C'), erbb2_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='3', st='U', loc='C'), erbb3_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='4', st='U', loc='C'), erbb4_0)
    Initial(DEP(b=None), DEP_0)
    Initial(ATP(b=None), ATP_0)
            
def rec_events():
    """ Describe receptor-level events here. 
    """

    # Parameter definitions
    # =====================
    # Alias model components for names in present namespace
    alias_model_components()
    
    # EGF / HRG binding to receptors
    # EGF / HRG receptor binding rates obtained from Chen et al Table 1 pg. 5
    bind_table([[                                                   EGF,             HRG],
                [erbb(ty='1', b=None, st='U', loc='C'),   (1e7, 3e-2),            None],
                [erbb(ty='3', b=None, st='U', loc='C'),            None,  (1e7, 7e-2)],
                [erbb(ty='4', b=None, st='U', loc='C'),            None,  (1e7, 7e-2)]],
                'bl', 'b')
    
    # ErbB dimerization
    # Dimerization rates obtained from Chen et al Table 1 pg. 5
    bind_table([[                       erbb(ty='1', b=None, st='U', loc='C'), erbb(ty='2', b=None, st='U', loc='C'), erbb(ty='3', b=None, st='U', loc='C'), erbb(ty='4', b=None, st='U', loc='C')],
                [erbb(ty='1', b=None, st='U', loc='C'),        (KDIMF, KDIMR),                          None,                 None,                  None],
                [erbb(ty='2', b=None, st='U', loc='C'),        (KDIMF, KDIMR),                (KDIMF, KDIMR),                 None,                  None],
                [erbb(ty='3', b=None, st='U', loc='C'),        (KDIMF, KDIMR),                (KDIMF, KDIMR),                 None,                  None],
                [erbb(ty='4', b=None, st='U', loc='C'),        (KDIMF, KDIMR),                (KDIMF, KDIMR),                 None,                  None]],
        'bd', 'bd')

    # ATP binding: ATP only binds to dimers
    # ATP binding rates obtained from Chen et al (k122, k_122)
    bind_table([[                                            ATP],
                [erbb(ty='1', bd=ANY, st='U', loc='C'), (1.87e-8, 1)],
                [erbb(ty='2', bd=ANY, st='U', loc='C'), (1.87e-8, 1)],
                [erbb(ty='4', bd=ANY, st='U', loc='C'), (1.87e-8, 1)]],
        'b', 'b')

    # This works b/c only erbb1, 2, and 4 have ATP, and they can cross-phosphorylate any other receptor
    # erbb2:erbb2 pairs only happen by dissociation of phosphorylated monomers
    # 
    for i in ['1','2','4']:
        for j in ['1','2','3','4']:
            Rule("cross_phospho_"+i+"_"+j,
                 ATP(b=1) % erbb(ty=i, b=1,    bd=2) % erbb(ty=j, bd=2, st='U') >>
                 ADP()    + erbb(ty=i, b=None, bd=2) % erbb(ty=j, bd=2, st='P'),
                 Parameter("kcp"+i+j, KCP))

    
    # Receptor Dephosphorylation
    # DEPHOSPHORYLATION: 
    #  * Density enhanced phosphatase1 (DEP1) dephosphorylates ERB1 (at the cell-membrane)
    #  * Protein Tyrosine Phosphatase1b (PTP1b) dephosphorylates all RTKs (at the endo-membrane)
    #  Bursett, TA, Hoier, EF, Hajnal, A: Genes Dev. 19:1328-1340 (2005)
    #  Haj, FG, Verver, PJ, Squire, A, Neel, BG, Bastiaens, PI: Science 295:1708-1711 (2002)
    #  FIXME: REPLACE W A NEW CATALYSIS TABLE????
    for i in ['1','2','3','4']:
        catalyze(DEP(), 'b', erbb(st='P', ty=i), 'b', erbb(st='U', ty=i),
                 (KF,KR,KCD))
        
    # Receptor internalization
    # This internalizes all receptor combos 
    Rule("rec_intern",
         erbb(bd=1, loc='C') % erbb(bd=1, loc='C') <> erbb(bd=1, loc='E') % erbb(bd=1, loc='E'),
         Parameter("kintf", KINTF), Parameter("kintr", KINTR))

    # Receptor degradation
    # This degrades all receptor combos within an endosome
    degrade(erbb(bd=1, loc='E') % erbb(bd=1, loc='E'), Parameter("kdeg", KDEG))

    # FIXME:need negative feedback from ERK and AKT. include that in the other modules?

def mapk_monomers():
    Monomer('GAP', ['bd', 'b'])
    Monomer('SHC', ['bgap', 'bgrb', 'batp', 'st'], {'st':['U','P']})
    # Monomer('SHCPase', ['b'])
    Monomer('GRB2', ['b', 'bsos'])
    Monomer('SOS', ['bgrb', 'bras'])
    Monomer('RAS', ['bsos', 'braf', 'st'], {'st':['GDP', 'GTP']})
    Monomer('RAF', ['b', 'st'], {'st':['U', 'P']})
    Monomer('PP1', ['b'])
    Monomer('PP2', ['b'])
    Monomer('PP3', ['b'])
    Monomer('MEK', ['b', 'st'], {'st':['U', 'P', 'PP']})
    Monomer('ERK', ['b', 'st'], {'st':['U', 'P', 'PP']})

def mapk_initial():

    # Initial conditions obtained from Schoeberl et al.
    Parameter('GAP_0', 1.2e4)
    Parameter('SHC_0', 1.01e6)
    # Parameter('SHCPase_0', 1000)
    Parameter('GRB2_0', 5.1e4)
    Parameter('SOS_0', 6.63e4)
    Parameter('RAS_0', 1.14e7)
    Parameter('RAF_0', 4e4)
    Parameter('PP1_0', 4e4)
    Parameter('PP2_0', 4e4)
    Parameter('PP3_0', 1e7)
    Parameter('MEK_0', 2.2e7)
    Parameter('ERK_0', 2.1e7)

    alias_model_components()

    Initial(GAP(bd=None, b=None), GAP_0)
    Initial(SHC(bgap=None, bgrb=None, batp=None, st='U'), SHC_0)
    Initial(GRB2(b=None, bsos=None), GRB2_0)
    Initial(SOS(bgrb=None, bras=None), SOS_0)
    Initial(RAS(bsos=None, braf=None, st='GDP'), RAS_0)
    Initial(RAF(b=None, st='U'), RAF_0)
    Initial(MEK(b=None, st='U'), MEK_0)
    Initial(ERK(b=None, st='U'), ERK_0)
    Initial(PP1(b=None), PP1_0)
    Initial(PP2(b=None), PP2_0)
    Initial(PP3(b=None), PP3_0)

    
def mapk_events():

    # =====================
    # Alias model components for names in present namespace
    alias_model_components()

    # GAP binds to phosphorylated dimers
    # in the present we use MatchOnce to insure correct representation of the binding
    # similar to Chen et al
    Rule("GAP_binding",
         MatchOnce(erbb(bd=1, b=None, st='P') % erbb(bd=1, b=None, st='P')) + GAP(bd=None, b=None) <>
         MatchOnce(erbb(bd=1, b=2,    st='P') % erbb(bd=1, b=None, st='P')  % GAP(bd=2, b=None)),
         Parameter("kerbb_dim_GAPf", KF), Parameter("kerbb_dim_GAPr", KR))
    
    # SHC binds to GAP-complex
    bind(GAP(bd=ANY), 'b', SHC(batp=None, st='U'), 'bgap', [KF, KR])

    # SHC phosphorylation
    Rule("Shc_bind_ATP",
         GAP(bd=ANY, b=1) % SHC(bgap=1, bgrb=None, batp=None, st='U') + ATP(b=None) <>
         GAP(bd=ANY, b=1) % SHC(bgap=1, bgrb=None, batp=2, st='U') % ATP(b=2),
         Parameter("ShcATPf",KF), Parameter("ShcATPr",KR))
    
    Rule("Shc_phos",
         GAP(bd=ANY, b=1) % SHC(bgap=1, bgrb=None, batp=2, st='U') % ATP(b=2) >>
         GAP(bd=ANY, b=1) % SHC(bgap=1, bgrb=None, batp=None, st='P') + ADP(),
         Parameter("ShcPhosc", KCP))

    # GRB2 binds to GAP-SHC:P
    bind(SHC(batp=None, st='P'), 'bgrb', GRB2(), 'b', [KF, KR])

    # SOS binds to GAP-SHC:P-GRB2 - rate obtained from Chen et al. pg 5. 
    bind(GRB2(b=ANY), 'bsos', SOS(bras=None), 'bgrb', [7.5e-6, 1.5])

    # RAS-GDP binds to GAP-SHC:P-GRB2-SOS
    bind(SOS(bgrb=ANY), 'bras', RAS(braf=None, st='GDP'), 'bsos',  [KF, KR])

    # RAS-GTP binds to GAP-SHC:P-GRB2-SOS
    catalyze(SOS(bgrb=ANY), 'bras', RAS(braf=None, st='GTP'), 'bsos', RAS(bsos=None, braf=None, st='GDP'),
             (KF,KR,KCD))

    # RAS-GDP dissociates from GAP complex to form RAS-GTPU
    Rule("RAS_GDP_to_RAS_GTP",
         SOS(bgrb=ANY, bras=4) % RAS(bsos=4, braf=None, st='GDP') <>
         SOS(bgrb=ANY, bras=None) + RAS(bsos=None, braf=None, st='GTP'),
         Parameter("RasGDP_GTPf",KF), Parameter("RasGDP_GTPr",KR))

    # Activation of RAF -> RAF:P by RAS-GTP 
    catalyze(RAS(bsos=None, st='GTP'), 'braf', RAF(st='U'), 'b', RAF(st='P'),
             (KF,KR,KCP))

    
    # Deactivation of RAF:P -> RAF by PP1
    catalyze(PP1(), 'b', RAF(st='P'), 'b', RAF(st='U'),
             (KF,KR,KCD))

    # Activation of MEK -> MEK:P by activated RAF
    catalyze(RAF(st='P'), 'b', MEK(st='U'), 'b', MEK(st='P'),
             (KF,KR,KCP))

    # Deactivation of MEK:P -> MEK by PP2
    catalyze(PP2(), 'b', MEK(st='P'), 'b', MEK(st='U'),
             (KF,KR,KCD))
    
    # Activation of MEK:P -> MEK:P:P by activated RAF
    catalyze(RAF(st='P'), 'b', MEK(st='P'), 'b', MEK(st='PP'),
             (KF,KR,KCP))

    # Deactivation of MEK:P:P -> MEK:P by PP2
    catalyze(PP2(), 'b', MEK(st='PP'), 'b', MEK(st='P'),
             (KF,KR,KCD))
    
    # Activation of ERK -> ERK:P by activated MEK:P:P
    catalyze(MEK(st='PP'), 'b', ERK(st='U'), 'b', ERK(st='P'),
             (KF,KR,KCP))

    # Deactivation of ERK:P -> ERK by PP3
    catalyze(PP3(), 'b', ERK(st='P'), 'b', ERK(st='U'),
             (KF,KR,KCD))

    # Activation of ERK:P -> ERK:P:P by activated MEK:P:P
    catalyze(MEK(st='PP'), 'b', ERK(st='P'), 'b', ERK(st='PP'),
             (KF,KR,KCD))

    # Deactivation of ERK:P:P -> ERK:P by PP3
    catalyze(PP3(), 'b', ERK(st='PP'), 'b', ERK(st='P'),
             (KF,KR,KCD))
