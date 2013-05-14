def crosstalk():
    """Defines crosstalk between MAPK and AKT pathways"""
    #AKT:P:P binds RAF:P; KF and KR from Chen et al Supplementary info
    bind(AKT(bpip3=None, both=None, bRAF=None, S='PP'), 'bRAF', RAF(b=None, st='P'), 'b', [2.5e-5, 1e-1])
    #ERK:P:P binds GAP-GRB2:SOS
    #ERK:P:P binds GAP-GRB2:GAB1
    #GAP:GRB2:GAB1:PI3K binds RAS:GDP and RAS:GTP and hydrolyzes GTP -> GDP
