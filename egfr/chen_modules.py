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
    # Initial concentrations of ligands, receptor monomers, and ATP taken from Chen et al. 
    Parameter('EGF_0',      5e-9)
    Parameter('HRG_0',         0)
    Parameter('erbb1_0',  1.08e6)
    Parameter('erbb2_0',  4.62e5)
    Parameter('erbb3_0',  6.23e3)
    Parameter('erbb4_0',  7.94e2)
    Parameter('ATP_0',     1.2e9)
    Parameter('DEP_0',       7e4)
    

    alias_model_components()

    Initial(EGF(b=None), EGF_0)
    Initial(HRG(b=None), HRG_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='1', st='U', loc='C'), erbb1_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='2', st='U', loc='C'), erbb2_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='3', st='U', loc='C'), erbb3_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='4', st='U', loc='C'), erbb4_0)
    Initial(ATP(b=None), ATP_0)
    Initial(DEP(b=None), DEP_0)

            
def rec_events():
    """ Describe receptor-level events here. 
    """

    # Parameter definitions
    # =====================
    # Alias model components for names in present namespace
    alias_model_components()
    
    # EGF / HRG binding to receptors
    # EGF / HRG receptor binding rates obtained from Chen et al (Supplementary)
    bind_table([[                                                          EGF,             HRG],
                [erbb(ty='1', bd=None, b=None, st='U', loc='C'),   (1e7, 3.3e-2),          None],
                [erbb(ty='3', bd=None, b=None, st='U', loc='C'),            None,   (1e7, 7e-2)],
                [erbb(ty='4', bd=None, b=None, st='U', loc='C'),            None,   (1e7, 7e-2)]],
                'bl', 'b')
    
    # ErbB dimerization
    # Dimerization rates obtained from Chen et al (Supplementary)
    # erbb's is required to containe a ligand (except for erbb2 which cannot bind a ligand)
    erbb1Lig = erbb(ty='1', bl=ANY, b=None, st='U', loc='C')
    erbb2Lig = erbb(ty='2', bl=None, b=None, st='U', loc='C')
    erbb3Lig = erbb(ty='3',bl=ANY, b=None, st='U', loc='C')
    erbb4Lig = erbb(ty='4',bl=ANY, b=None, st='U', loc='C')
    bind_table([[                          erbb1Lig,            erbb2Lig, erbb3Lig, erbb4Lig],
                [erbb1Lig,        (7.45e-6, 1.6e-1),                None,     None,     None],
                [erbb2Lig,        (3.74e-8, 1.6e-2),  (1.67e-10, 1.6e-2),     None,     None],
                [erbb3Lig,        (3.74e-8, 1.6e-2),  (1.67e-10, 1.6e-2),     None,     None],
                [erbb4Lig,        (3.74e-8, 1.6e-2),  (1.67e-10, 1.6e-2),     None,     None]],
        'bd', 'bd')

    # ATP binding: ATP only binds to dimers
    # ATP binding rates obtained from Chen et al (Supplementary)
    # include DEP binding here since they both bind to the same site
    bind_table([[                                                ATP,  DEP],
                [erbb(ty='1', bd=ANY, st='U', loc='C'), (1.87e-8, 1), (5e-5, 1e-2)],
                [erbb(ty='2', bd=ANY, st='U', loc='C'), (1.87e-8, 1), (5e-5, 1e-2)],
                [erbb(ty='4', bd=ANY, st='U', loc='C'), (1.87e-8, 1), (5e-5, 1e-2)]],
        'b', 'b')

    # Cross phosphorylation: only erbb1, 2, and 4 have ATP, and they can cross-phosphorylate any other receptor
    # erbb2:erbb2 pairs only happen by dissociation of phosphorylated monomers
    # kcat phosphorylation obtained from Chen et al Table I pg. 5

    # Both dimers become phosphorylated/dephosphorylated to agree better w Chen/Schoeberl model
    # Receptor Dephosphorylation
    # DEPHOSPHORYLATION: 
    #  * Density enhanced phosphatase1 (DEP1) dephosphorylates ERB1 (at the cell-membrane)
    #  * Protein Tyrosine Phosphatase1b (PTP1b) dephosphorylates all RTKs (at the endo-membrane)
    #  Bursett, TA, Hoier, EF, Hajnal, A: Genes Dev. 19:1328-1340 (2005)
    #  Haj, FG, Verver, PJ, Squire, A, Neel, BG, Bastiaens, PI: Science 295:1708-1711 (2002)

    for i in ['1','2','4']:
        for j in ['1','2','3','4']:
            Rule("cross_phospho_"+i+"_"+j,
                 ATP(b=1) % erbb(ty=i, b=1,    bd=2, st='U') % erbb(ty=j, bd=2, st='U') >>
                 ADP()    + erbb(ty=i, b=None, bd=2, st='P') % erbb(ty=j, bd=2, st='P'),
                 Parameter("kcp"+i+j, 1e-1))
            Rule("cross_DEphospho_"+i+"_"+j,
                 DEP(b=1)   %  erbb(ty=i, b=1,    bd=2, st='P') % erbb(ty=j, bd=2, st='P') >>
                 DEP(b=None) + erbb(ty=i, b=None, bd=2, st='U') % erbb(ty=j, bd=2, st='U'),
                 Parameter("kcd"+i+j, 1e-1))

    # Receptor internalization
    # This internalizes all receptor combos
    # FIXME: this should just be a state transformation (i.e. the rule looks noisy)
    # Internalization rates taken from Chen et al Table I pg. 5
    Rule("rec_intern",
         erbb(bd=1, loc='C') % erbb(bd=1, loc='C') <> erbb(bd=1, loc='E') % erbb(bd=1, loc='E'),
         Parameter("kintf", 1.3e-3), Parameter("kintr", 5e-5))

    # Receptor degradation
    # This degrades all receptor combos within an endosome
    # Should this have a forward and reverse rate - k62b
    degrade(erbb(bd=1, loc='E') % erbb(bd=1, loc='E'), Parameter("kdeg", 4.16e-4))

    # FIXME: need negative feedback from ERK and AKT. include that in the other modules?

def mapk_monomers():
    Monomer('GAP', ['bd', 'b', 'bgrb2'])
    Monomer('SHC', ['bgap', 'bgrb', 'batp', 'st'], {'st':['U','P']})
    # Monomer('SHCPase', ['b'])
    Monomer('GRB2', ['b', 'bsos', 'bgap', 'bgab1'])
    Monomer('SOS', ['bgrb', 'bras', 'bERKPP', 'st'], {'st':['U', 'P']})
    Monomer('RAS', ['bsos', 'braf', 'bpi3k', 'st'], {'st':['GDP', 'GTP']})
    Monomer('RAF', ['b', 'st', 'ser295'], {'st':['U', 'P'], 'ser295':['U', 'P']})
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
    Parameter('MEK_0', 2.2e7)
    Parameter('ERK_0', 2.1e7)
    Parameter('PP1_0', 4e4)
    Parameter('PP2_0', 4e4)
    Parameter('PP3_0', 1e7)

    alias_model_components()

    Initial(GAP(bd=None, b=None, bgrb2=None), GAP_0)
    Initial(SHC(bgap=None, bgrb=None, batp=None, st='U'), SHC_0)
    Initial(GRB2(b=None, bsos=None, bgap=None, bgab1=None), GRB2_0)
    Initial(SOS(bgrb=None, bras=None, bERKPP=None, st='U'), SOS_0)
    Initial(RAS(bsos=None, braf=None, bpi3k=None, st='GDP'), RAS_0)
    Initial(RAF(b=None, st='U', ser295='U'), RAF_0)
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
    bind(GAP(bd=ANY, bgrb2=None), 'b', SHC(batp=None, st='U'), 'bgap', [KF, KR])

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
    bind(SHC(batp=None, st='P'), 'bgrb', GRB2(bgap=None, bgab1=None), 'b', [KF, KR])

    # SOS binds to GAP-SHC:P-GRB2 - rate obtained from Chen et al. pg 5. 
    bind(GRB2(b=ANY, bgap=None, bgab1=None), 'bsos', SOS(bras=None, st='U'), 'bgrb', [7.5e-6, 1.5])

    # SOS also binds GAP-GRB2
    Rule("GAP_GRB2_bind_SOS",
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=None) + SOS(bras=None, bgrb=None, st='U') <>
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1) % SOS(bras=None, bgrb=1, st='U'),
         Parameter("GAP_GRB2_bind_SOSf", 7.5e-6),
         Parameter("GAP_GRB2_bind_SOSr", 1.5))

    # GAP-GRB2-SOS and GAP-SHC:P-GRB2-SOS catalyze RAS-GDP->RAS-GTP:
    catalyze_state(SOS(bgrb=ANY), 'bras', RAS(braf=None), 'bsos', 'st', 'GDP', 'GTP', (KF, KR, KCD))

    # RAS-GDP binds to GAP-SHC:P-GRB2-SOS
    #bind(SOS(bgrb=ANY), 'bras', RAS(braf=None, st='GDP'), 'bsos',  [KF, KR])

    # RAS-GTP binds to GAP-SHC:P-GRB2-SOS
    #catalyze(SOS(bgrb=ANY), 'bras', RAS(braf=None, st='GTP'), 'bsos', RAS(bsos=None, braf=None, st='GDP'),
    #        (KF,KR,KCD))

    # RAS-GDP dissociates from GAP complex to form RAS-GTPU
    #Rule("RAS_GDP_to_RAS_GTP",
    #    SOS(bgrb=ANY, bras=4) % RAS(bsos=4, braf=None, st='GDP') <>
    #   SOS(bgrb=ANY, bras=None) + RAS(bsos=None, braf=None, st='GTP'),
    #   Parameter("RasGDP_GTPf",KF), Parameter("RasGDP_GTPr",KR))

    # Activation of RAF -> RAF:P by RAS-GTP 
    catalyze(RAS(bsos=None, st='GTP'), 'braf', RAF(st='U'), 'b', RAF(st='P'),
             (KF,KR,KCP))

    # Deactivation of RAF:P -> RAF by PP1
    catalyze(PP1(), 'b', RAF(st='P'), 'b', RAF(st='U'),
             (KF,KR,KCD))

    # Activation of MEK -> MEK:P by activated RAF
    catalyze(RAF(st='P', ser295='U'), 'b', MEK(st='U'), 'b', MEK(st='P'),
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

def akt_monomers():
    """ This is the akt part of the pathway from the Chen et al. 2009 paper.  Initial rules for all binding reactions were generated and then coded again using macros and higher order macros.  Initial parameters and conditions were taken from Chen et al. 2009 paper and supplementary, but were later modified in order to get the model working correctly.  This pathway follows AKT from its initial state to a phosphorylated and then double phosphorylated state before returning to unphosphorylated AKT.  The model works correctly, but parameters and rates may need to be modified in order to get best fit.  Parameters and rates included are from trial and error of what best fit the model.  The last unbinding reactions may not be needed because of the catalyze_state macros, but were left in just in case these are needed later.  
"""
    #This pathway coded by Tim O'Brien.
    Monomer('GAB1', ['bgrb2', 'bshp2', 'bpi3k', 'bpi3k2', 'bpi3k3', 'bpi3k4', 'bpi3k5', 'bpi3k6','batp','bERKPP','bPase9t','S'],{'S':['U','P','PP']})
    Monomer('PI3K',['bgab1','bpip', 'bras'])
    Monomer('SHP2',['bgab1'])
    Monomer('PIP', ['bakt', 'both', 'S', 'bpi3k'], {'S':['PIP2', 'PIP3']})
    Monomer('PTEN', ['bpip3', 'both'])
    Monomer('SHP', ['bpip3', 'both'])
    Monomer('AKT', ['bpip3', 'both', 'S'], {'S':['U', 'P', 'PP']})
    Monomer('PDK1', ['bakt', 'both'])
    Monomer('PP2A_III', ['bakt', 'both'])

def akt_initial():
    # Initial concentrations modified from Chen et al. 2009 Supplementary text values
    # Initial PIP_0 not from Chen et al.
    Parameter('GAB1_0', 94868.3)
    Parameter('PI3K_0', 3.55656e7)
    Parameter('SHP2_0', 1e6)
    Parameter('PIP_0',     9.00e3)
    Parameter('PTEN_0',    5.00e4)
    Parameter('SHP_0',     7.00e4)
    Parameter('AKT_0',     9.05e6)
    Parameter('PDK1_0',     9.5e6)
    Parameter('PP2A_III_0', 4.5e5)
    alias_model_components()
    
    # Initial conditions 
    Initial(GAB1(bgrb2=None, bshp2=None, bpi3k=None, bpi3k2=None, bpi3k3=None, bpi3k4=None, bpi3k5=None, bpi3k6=None, batp=None, bERKPP=None, bPase9t=None, S='U'), GAB1_0)
    Initial(PI3K(bgab1=None, bpip=None, bras=None), PI3K_0)
    Initial(SHP2(bgab1=None), SHP2_0)
    Initial(PIP(bakt=None, both=None, S='PIP2', bpi3k=None), PIP_0)
    Initial(PTEN(bpip3=None, both=None), PTEN_0)
    Initial(SHP(bpip3=None, both=None), SHP_0)
    Initial(AKT(bpip3=None, both=None, S='U'), AKT_0)
    Initial(PDK1(bakt=None, both=None), PDK1_0)
    Initial(PP2A_III(bakt=None, both=None), PP2A_III_0)
def akt_events():
    # parameter values taken from Chen et al. 2009
    Parameter('kf',  1e-5)
    Parameter('kr',  1e-3)
    Parameter('kcd', 1e-1)
    alias_model_components()
    #GRB2 binds GAP-complex (without requiring SHC bound to complex) k16, kd24
    bind(GRB2(b=None, bsos=None), 'bgap', GAP(bd=ANY, b=None), 'bgrb2', [1.67e-05, 5.5e-01])
    #GAB1 binds GAP-GRB2 k105, kd105
    bind(GRB2(b=None, bsos=None, bgap=ANY), 'bgab1', GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None), 'bgrb2', [6.67e-05, 1e-01])
    #GAP-GRB2-GAB1 phosphorylation - Rates from Table p. 5 Chen et al 2009
    Rule('GAB1_bind_ATP',
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='U') + ATP(b=None) <>
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=1, bERKPP=None, bPase9t=None, bgrb2=ANY, S='U') % ATP(b=1),
         Parameter('GAB1ATPf', 1e-5),
         Parameter('GAB1ATPr', 1e-1))

    Rule('GAB1_phos',
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=1, bERKPP=None, bPase9t=None, bgrb2=ANY, S='U') % ATP(b=1) >>
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + ADP(),
         Parameter('GAB1Phosc', 1e-1))

    #SHP2 can desphosphorylate GAB1-P
    catalyze_state(SHP2(), 'bgab1', GAB1(bgrb2=ANY, bpi3k=None, batp=None, bERKPP=None, bPase9t=None), 'bshp2', 'S', 'P', 'U', (1e-5, 1e-1, 1e-2))
   
    #After GAB1 phosphorylation, all receptor dimer combinations can bind a single PI3K
    Rule('GAB1_bind_PI3K_1',
         GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + PI3K(bpip=None, bgab1=None, bras=None) <>
         GAB1(bshp2=None, bpi3k=1, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=1, bras=None),
         Parameter('GAB1PI3Kf', 1e-5),
         Parameter('GAB1PI3Kr', 1e-1))

    #ErbB2-ErbB3 dimers bound to GAB1 contain 6 binding domains for PI3K.
    for i in range(1,6):
         pi3k_string_beg = '+ PI3K(bpip=None) + PIP(S="PIP2")'
         pi3k_string_end = '+ PI3K(bpip=None) + PIP(S="PIP3")'
         Rule('ErbB2_3bindPI3K'+str(i),
              erbb(bd=1, ty='2') % erbb(bd=1, ty='3') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') i*pi3k_string_beg<>
              erbb(bd=1, ty='2') % erbb(bd=1, ty='3') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') i*pi3k_string_end,
              Parameter('GAB1PI3Kf'+str(i), 1e-5),
              Parameter('GAB1PI3Kr'+str(i), 1e-1))
             
    #gsites= ['bpi3k2'] #, 'bpi3k3', 'bpi3k4', 'bpi3k5', 'bpi3k6']
                       #GAB1specs_none=[GAB1(bpi3k2=None, bshp2=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P')]
    #GAB1specs_b=[GAB1(bpi3k2=1, bshp2=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P')]
    
    #for gspec in gsites:
    # Rule('GAB1_bind_'+gspec,
    #erbb(bd=1, ty='2') % erbb(bd=1, ty='3') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) %  + PI3K(bgab1=None, bpip=None) <>
    #erbb(bd=1, ty='2') % erbb(bd=1, ty='3') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % gspec % PI3K(bgab1=1, bpip=None),
    #  Parameter('GAB1'+gspec+'f', 1e-5),
    #  Parameter('GAB1'+gspec+'r', 1e-1))

    #PI3K bound to complex catalyzes PIP2 -> PIP3
    catalyze_state(PI3K(bgab1=ANY), 'bpip', PIP(bakt=None, both=None), 'bpi3k', 'S', 'PIP2', 'PIP3', (1e-5, 1e-1, 1e-1))
             
     # Setting up the binding reactions necessary for AKT to be phosphorylated and move through the pathway
    bind_table([[                       AKT(S='U', both=None),       AKT(S='P', both=None),      AKT(S='PP', both=None)],
                [PIP(S='PIP3', both=None, bpi3k=None),       (1e-5, 1e-3),     (1e-5, 1e-3),    (1e-5, 1e-3)]],
                'bakt', 'bpip3')

    # AKT-PIP3 is phosphorylated by PDK1 to AKTP
    catalyze(PDK1, 'bakt', AKT(bpip3=ANY, S='U'), 'both', AKT(bpip3=ANY, S='P'), (kf, kr, kcd))

    # AKTP-PIP3 is phosphorylated by PDK1 to AKTPP
    catalyze(PDK1, 'bakt', AKT(bpip3=ANY, S='P'), 'both', AKT(bpip3=ANY, S='PP'), (kf, kr, kcd))

    # AKTP is dephosphorylated by PP2A-III back to AKT
    catalyze_state(PP2A_III, 'bakt', AKT(bpip3=None), 'both', 'S', 'P', 'U',(1e-5, 1e-3, 1e-3))
   
    # AKTPP is dephosphorylated by PP2A-III back to AKTP
    catalyze_state(PP2A_III, 'bakt', AKT(bpip3=None), 'both', 'S', 'PP', 'P',(1e-5, 1e-3, 1e-3))

    # PIP3 is dephosphorylated by PTEN to PIP2
    catalyze_state(PTEN, 'bpip3', PIP(bakt=None), 'both', 'S', 'PIP3', 'PIP2', (1e-5, 1e-3, 1e-3))

    # PIP3 is dephosphorylated by SHP to PIP2
    catalyze_state(SHP, 'bpip3', PIP(bakt=None), 'both', 'S', 'PIP3', 'PIP2', (1e-5, 1e-3, 1e-3))

    # Release of AKTP from PIP3 (this may not be needed because of binding table above)
    Rule('AKT_PDK1_undbind', AKT(bpip3=ANY, both=None, S='U') % PIP(bakt=1, both=2, S='PIP3') % PDK1(bakt=None, both=2) >> AKT(bpip3=None, both=None, S='P') + PDK1(bakt=None, both=1) % PIP(bakt=None, both=1, S='PIP3'), kf)

# PDK1-PIP3 dissociate after phosphorylation of AKT
    Rule('PDK1_PIP3_unbind', PDK1(bakt=None, both=1) % PIP(bakt=None, both=1, S='PIP3') >> PDK1(bakt=None, both=None) + PIP(bakt=None, both=None, S='PIP3'), kf)

# AKT-PDK1 Release of AKTPP from PIP3 (this may not be needed because of binding table above)
    Rule('AKT_PDK1_unbind', AKT(bpip3=1, both=None, S='P') % PIP(bakt=1, both=2, S='PIP3') % PDK1(bakt=None, both=2) >> AKT(bpip3=None, both=None, S='PP') + PDK1(bakt=None, both=1) % PIP(bakt=None, both=1, S='PIP3'), kf)

# AKT:P-PP2A-III dissociate (this may not be needed to due to catalyze_state macro)
    Rule('AKT_PP2A_III_unbind', AKT(bpip3=None, both=1, S='P') % PP2A_III(bakt=1) >> AKT(bpip3=None, both=None, S='U') + PP2A_III(bakt=None), kf)

# PIP3-PDK1 dissociate to PIP3 and PDK1 (this may not be needed due to catalyze_state macro)
    Rule('PIP3_PDK1_unbind', PIP(bakt=None, both=1, S='PIP3') % PDK1(bakt=None, both=1) >> PIP(bakt=None, both=None, S='PIP3') + PDK1(bakt=None, both=None), kf)

def crosstalk_monomers():
    Monomer('Pase9t', ['bgab1'])
    alias_model_components()
def crosstalk_initial():
    Parameter('Pase9t_0', 0)

def crosstalk_events():
    #ERK:P:P phosphorylates GAP-GRB2-GAB1:P (making it unable to bind PI3K)
    catalyze_state(ERK(st='PP'), 'b', GAB1(bgrb2=ANY, bpi3k=None, bpi3k2=None, bpi3k3=None, bpi3k4=None, bpi3k5=None, bpi3k6=None), 'bERKPP', 'S', 'P', 'PP', (1e-5, 1e-1, 1e-1))

    #GAP-GRB2-GAB1:P:P is dephosphorylated by Pase9t
    catalyze_state(Pase9t(), 'bgab1', GAB1(bgrb2=ANY), 'bPase9t', 'S', 'PP', 'P', (1e-5, 1e-1, 1e-2))

    #ERK:P:P phosphorylates GRB2-SOS, preventing RAS-GDP->RAS-GTP conversion
    catalyze_state(ERK(st='PP'), 'b', SOS(bgrb=ANY, bras=None), 'bERKPP', 'st', 'U', 'P', (1e-5, 1e-1, 1e-1))

    #AKT:P:P phosphorylates RAF:P at Ser295, preventing MEK phosphorylation.
    catalyze_state(AKT(S='PP', bpip3=None), 'both', RAF(st='P'), 'b', 'ser295', 'U', 'P', (1e-5, 1e-1, 1e-1))

    #RAS-GTP binds PI3K
    bind(RAS(bsos=None, braf=None, st='GTP'), 'bpi3k', PI3K(bgab1=None, bpip=None), 'bras', [1e-5, 1e-1])

    #RAS-GTP-PI3K binds PIP2
    bind(PI3K(bgab1=None, bpip=None, bras=ANY), 'bpip', PIP(bakt=None, both=None, bpi3k=None, S='PIP2'), 'bpi3k',
         [1e-5, 1e-1])

    #RAS-GTP-PI3K-PIP2 disassociates to give PIP3
    Rule('RASGTPPI3KcatPIP',
         RAS(bsos=None, braf=None, st='GTP', bpi3k=1) % PI3K(bgab1=None, bpip=2, bras=1) % PIP(bakt=None, both=None, bpi3k=2, S='PIP2') >>
         RAS(bsos=None, braf=None, st='GTP', bpi3k=1) % PI3K(bgab1=None, bpip=None, bras=1) + PIP(bakt=None, both=None, bpi3k=None, S='PIP3'),
         Parameter('RASGTPPI3K_kc', 1e-1))
    
