(in-package :b-user)
(include b-user/ode-biochem)
(setf *max-monomers-per-complex* 25)
;;(include b-user/3d-ode-biochem) ;use when using reaction rates with units
;;;;
;;;; ------BIOLOGY HERE:
;;;;
;;;;
;;;  -----RATES:
(include @FOLDER/erb_rates)
;;;  A stereotyped cell which can
;;;  be used to define internalization and other complex reactions
;;;  in which species move between sublocations of the cell.
;;;
;; -------------------------------------------------------------------------;;;;
;; BEGIN DEFINE CELL STUFF
;; -------------------------------------------------------------------------;;;;
;; DEFINE CELL AND SPECIES IN COMPARTMENTS
;; CELL DEFINITION 
;; \wikidoc{CELLDEF}
;; set all volumes to 1.0 since the rates already account for volumes
(defcon cell-location (location)
  (&optional (id := *name*) 
	     &property
	     (cytoplasm compartment :#= [[compartment]{.size.value := 1.0}]) ; volume from birgit (pL) , volume of the cell
	     (dish      compartment :#= [[compartment]{.size.value := 1.0 }]) ; volume from birgit (pL),  1.0E9
	     (endosome compartment  :#= [[compartment]{.size.value := 1.0}]) ; endosome volume from birgit (pL), 4.2E-6
	     (cell-membrane membrane :#= (let ((d .dish)
					       (c .cytoplasm))
					; membrane volume based on a cell radius of 7.5um and a membrane thickness of .03um
					   [[membrane] :outer d :inner c {.size.value := 1.0}])) ; (pL) 1.7665E-3
	     (endo-membrane membrane :#= (let ((e .endosome)
					       (c .cytoplasm))
					   [[membrane] :inner e :outer c {.size.value := 1.0}])) ; pL based on a membrane thickness of 2nm 2.46E-7
	     (inverse-endo-membrane membrane :#= .endo-membrane.inverse))) ; invert inner and outer

;; -------------------------------------------------------------------------;;;;
;; BEGIN MONOMER DEFINITIONS
;; -------------------------------------------------------------------------;;;;
;; \wikidoc{LIGANDS}
(defmonomer egf R)
(defmonomer hrg R)
;; RECEPTOR:
;; \wikidoc{RTKs}
(with-substitution-table ($NAME erbb1 erbb2 erbb3 erbb4)
	(defmonomer ($NAME membrane)
	  (L :documentation "Ligand binding site")
	  (D :documentation "Dimerization site")
	  (psite :states (member :u :p) :default :u  :documentation "Phosphorylation site")
	  (CSHC :documentation "complexation site1, SH2 site in ERBB1")
	  (CGRB2 :documentation "complexation site, PTB site in ERBB1")
	  (CA :documentation "ATP-binding site")
	  (CPI3K :documentation "complexation site for PI3K")
	  (CPASE :documentation "phosphatase complexation site")
	  ))

;; EXPLICIT ATP
;; \wikidoc{atp}
(defmonomer ATP CA)
;(defmonomer ADP CA)
;; Only DEP1 is real. The others were added just as a fudge... (ask Paul about this)
;; \wikidoc{phosphatases}
(defmonomer DEP1            "Density enhanced phosphatase 1" CERB) 
(defmonomer DEP2            "Density enhanced phosphatase 2" CERB) 
(defmonomer DEP3            "Density enhanced phosphatase 3" CERB) 
(defmonomer DEP4            "Density enhanced phosphatase 4" CERB) 
;; \wikidoc{SHC}
(defmonomer SHC             "Shc complexation sites, CERB: erbb binding, CGRB2: grb2 binding"
  CERB CGRB2 (asite :states (member :active :inactive) :default :inactive) CPASE)
;; \wikidoc{SHCPASE}
(defmonomer SHCPASE         "Shc hoky phosphatase to deactivate it"
  CSHC)
;; \wikidoc{GRB2}
(defmonomer GRB2            "Grb2 complexation sites, CSOS: sos binding " 
   CERB CSHC CSOS CERK CGBSP)
;; \wikidoc{GBSP}
(defmonomer GBSP            "Combination Gab1 and Shp hoky protein to account for GAB1:SHP2"
  CGRB2 CGAB1 CSOS CPI3K)
;; \wikidoc{SOS}
(defmonomer SOS             "complexation site to GRB2/RAS" 
  CGRB2 CRAS CERK CGBSP)
;; \wikidoc{RAS}
(defmonomer (RAS membrane)  "RAS g-protein with GTP/GDP states" 
  CSOS CRAF CGAP CiRAFK (gtpsite :states (member :gdp :gtp) :default :gdp))
;; \wikidoc{GAP}
(defmonomer GAP             "GAP is a RAS deactivating protein" CRAS)
;; \wikidoc{RAF}
(defmonomer RAF             "RAF with ras binding site and active/inactive states"
  CRAS CMEK CPP2A CPI CiRAFK CAKT (psite :states (member :u :p) :default :u))
;; \wikidoc{iRAFK}
(defmonomer iRAFK           "Imaginary RAF kinase to fudge the RAF activation" 
  CRAF)
;; \wikidoc{MEK}
(defmonomer MEK             "MEK with binding site and two phospho sites" 
  CRAF CERK CPP2A
  (psite1 :states (member :u :p) :default :u) 
  (psite2 :states (member :u :p) :default :u))
;; \wikidoc{ERK}
(defmonomer ERK             "ERK with binding site and two phospho sites" 
  CMEK CMKP CSOS
  (psite1 :states (member :u :p) :default :u) 
  (psite2 :states (member :u :p) :default :u))
;; \wikidoc{PP2A}
(defmonomer PP2A            "pseudo-generic phosphatase" C) 
;; \wikidoc{MKP}
(defmonomer MKP             "MKP phosphatase for ERK" C) 
;; \wikidoc{PI3K}
(defmonomer PI3K            "PI3K with binding site for ERBB (what about p120cbl adaptor?!?!!!!!)" 
  CERB CPI CA CGBSP)
;; \wikidoc{PtdIns}
(defmonomer (PtdIns membrane) "Phosphatidylinositol in the membrane" 
  CPI3K CPTEN CAKT CPDK
  (pstate :states (member :p2_4-5 :p3_3-4-5) :default :p2_4-5 :documentation "4,5 or 3,4,5 phospho state"))
;; \wikidoc{PTEN}
(defmonomer PTEN            "Dephosphorylate PtdIns (generic for now, perhaps add others later" 
  CPI)
;; \wikidoc{AKT}
(defmonomer AKT             "AKT activates lots of stuff downstream " 
  CPI3K CPDK CPP2A CRAF
  (psite1 :states (member :u :p) :default :u)
  (psite2 :states (member :u :p) :default :u))
;; \wikidoc{PDK}
(defmonomer PDK            "PDK necessary for activation of AKT" 
  CPI3K CAKT CA)
;;; ----------------------------------------------------------------------------
;;; RECEPTOR-LEVEL EVENTS 
;;; ----------------------------------------------------------------------------
;; ERB MONOMERS BIND TO LIGAND:  R + L <-> R-L
;;
;; ERB + LIG -> ERB:LIG
;;
;; \wikidoc{erbb_lig_bind}
;;
(with-data-table 
 (:rows $R1 :cols $R2 :cells ($Kf $Kr) :ignore _)
 ((        egf                 hrg              )
  (erbb1  (KLbindF KLbindR)    _                )
  (erbb2   _                   _                )
  (erbb3   _                  (KHbindF KHbindR) )
  (erbb4    _                 (KHbindF KHbindR) )) 
 ;; FORWARD (BINDING)
[{[$R1 L._ CA.* __] + [$R2 R._] @ :outer <<->> [[$R1 L.1 CA.* __][$R2 R.1]]}
  (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)])

;; \wikidoc{erbb_dimerization}
;; notice erbb3 homodimers are not included here. we include the possiblity of erbb2
;; homodimers in the next reaction
;; notice we only need one instance of each pair since the rules are otherwise redundant
;; and simplify to the same rule in the end if they are instantiated.
(with-data-table 
    (:rows ($R1 $L1) :cols ($R2 $L2) :cells ($Kf $Kr) :ignore _)
    ((             (erbb1 *)     (erbb2 nil)   (erbb3 *)     (erbb4 *)     )
     ((erbb1   *)  (KdimF KdimR) (KdimF KdimR) (KdimF KdimR) (KdimF KdimR) )
     ((erbb2 NIL)  (KdimF KdimR) (KdimF KdimR) (KdimF KdimR) (KdimF KdimR) )
     ((erbb3   *)  (KdimF KdimR) (KdimF KdimR)       _       (KdimF KdimR) )
     ((erbb4   *)  (KdimF KdimR) (KdimF KdimR) (KdimF KdimR) (KdimF KdimR) ))
    ;;
  [{[[$R1 L.1 CA.* D._ psite.*][$L1 *.1]] + [[$R2 L.2 CA.* D._ psite.*][$L2 *.2]] <<->> 
    [[$R1 L.1 CA.* D.3 psite.*][$L1 *.1][$R2 L.2 CA.* D.3 psite.*][$L2 *.2]]}
   (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)])

;; ATP - RECEPTOR BINDING: R + ATP <-> R-ATP 
;;  erbb3 has no ATP binding domain
;;
;; \wikidoc{erbb_phosphorylation}
;; 
(with-substitution-table 
    (($R    $Kf $Kr) 
     (erbb1 KATPbindF   KATPbindR)
     (erbb2 KATPbindF   KATPbindR)
     (erbb4 KATPbindF   KATPbindR))
  ;; ATP BINDS/UNBINDS RECEPTOR:
  [{[$R L.* D.* CA._ __] + [ATP CA._ ] @ :inner <<->> [[$R L.* D.* CA.1 __][ATP CA.1]]} ; ATP at inner 'brane
   (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)])

;; RECEPTOR CROSS PHOSPHORYLATION:  ATP-R1-R2(U) -> R1-R2(P) + ADP
;;
;; When receptors dimerize, one kinase containing ATP can
;; phosphorylate the partner monomer.  notice that active dimers can
;; still undergo dissociation based on the reverse dimerization rule
;; (as desired)
(with-data-table 
    (:rows ($R1 $L1) :cols ($R2 $L2) :cells $K1 :ignore _)
    ((             (erbb1 *) (erbb2 nil) (erbb3 *) (erbb4 *))
     ((erbb1   *)   KphoC     KphoC       KphoC     KphoC   )
     ((erbb2 NIL)   KphoC     KphoC       KphoC     KphoC   )
     ((erbb3   *)   _         _           _          _      ) ;; erbb3 has no kinase activity
     ((erbb4   *)   KphoC     KphoC       KphoC      KphoC  ))
  [{[[$R1 L.1 D.3 CA.4 psite.* __][$L1 *.1][ATP CA.4][$R2 L.2 D.3 CA.* psite.u __][$L2 *.2]] ->>
    [[$R1 L.1 D.3 CA._ psite.* __][$L1 *.1][$R2 L.2 D.3 CA.* psite.p __][$L2 *.2]] + [ATP CA._] @ :inner}
   (.set-rate-function 'mass-action $K1)])

;; RECEPTOR INTERNALIZATION:  ERBB1-R @ cell-membrane -> ERBB1-R @ endo-membrane
;; Endocytosis can take 3 paths:
;; * degradation at the lysosome
;; * organization centers for signaling pathways (dealt with later)
;; * delivery to nucleus (ignored for all receptors)
;; ***NOTICE: not reversible
;; receptors can get internalized, but only the erbb1 receptors get degraded. 
;; all other receptors get sent back to the surface.
;; \wikidoc{erbb_internalization}
(with-data-table 
 (:rows $R1 :cols $R2 :cells $K1 :ignore _)
 ((          erbb1       erbb2       erbb3       erbb4  )
  (erbb1     KintF        _           _           _     )
  (erbb2      _           _           _           _     )
  (erbb3      _           _           _           _     )
  (erbb4      _           _           _           _     ))
 [{{[[$R1 D.1 psite.p **][$R2 D.1 psite.p **]] @ :cell-membrane} @ cell-location ->>
   [[$R1 D.1 psite.p **][$R2 D.1 psite.p **]] @ :inverse-endo-membrane}
  (.set-rate-function 'mass-action  $K1)])


;; RECEPTOR RECYCLING: ERBB1-R @ inverse-endo-membrane -> [nil]
;; * receptors can go back to the cellular membrane
;; ***NOTICE: not reversible
(with-data-table (:rows $R1 :cols $R2 :cells $K1 :ignore _)
 ((          erbb1       erbb2       erbb3       erbb4  )
  (erbb1  KintR     _           _           _     )
  (erbb2      _           _           _           _     )
  (erbb3      _           _           _           _     )
  (erbb4      _           _           _           _     ))
  [{{[[$R1 D.1 psite.p **][$R2 D.1 psite.p **]] @ :inverse-endo-membrane} @ cell-location ->>
    [[$R1 D.1 psite.p **][$R2 D.1 psite.p **]] @ :cell-membrane}
   (.set-rate-function 'mass-action  $K1)])

;; RECEPTOR DEGRADATION: ERBB1-R @ inverse-endo-membrane -> [nil]
;;***NOTICE: not reversible
(with-data-table (:rows $R1 :cols $R2 :cells $K1 :ignore _)
 ((          erbb1       erbb2       erbb3       erbb4  )
  (erbb1  KintF           _           _           _     )
  (erbb2      _           _           _           _     )
  (erbb3      _           _           _           _     )
  (erbb4      _           _           _           _     ))
  ;;destruction reaction
  [{{[[$R1 D.1 psite.p **][$R2 D.1 psite.p **]] @ :inverse-endo-membrane} @ cell-location ->> }
   (.set-rate-function 'mass-action $K1)])

;; DEPHOSPHORYLATION: 
;;  * Density enhanced phosphatase1 (DEP1) dephosphorylates ERB1 (at the cell-membrane)
;;  * Protein Tyrosine Phosphatase1b (PTP1b) dephosphorylates all RTKs (at the endo-membrane)
;;
;;  Bursett, TA, Hoier, EF, Hajnal, A: Genes Dev. 19:1328-1340 (2005)
;;  Haj, FG, Verver, PJ, Squire, A, Neel, BG, Bastiaens, PI: Science 295:1708-1711 (2002)
;;  
;; \wikidoc{erbb_dephosphorylation} CHOSE one phosphatase per erbb
(with-data-table (:rows $R :cols $P :cells $K2f :ignore _)
  ((        dep1      dep2      dep3      dep4         )
   (erbb1   KDepho1C  _         _         _            )
   (erbb2   _         KDepho2C  _         _            )
   (erbb3   _         _         KDepho3C  _            )
   (erbb4   _         _         _         KDepho4C     ))
  ;; Important: note here that dephosphorylation can take place at any membrane
  ;; (in this case the endo or the cell membrane). In this case @inner is the right
  ;; side for both membranes since we are using inverse-endo-membrane
  ;;
  [{[$R psite.p CPASE._ **] + [$P CERB._] @ :inner ->> [$R psite.u CPASE._ **] + [$P CERB._] @ :inner}
   (.set-rate-function 'mass-action $K2f)])

;; -------------------------------------------------------------------------;;;;
;; MAPK PATHWAY 
;; -------------------------------------------------------------------------;;;;

;; BINDING TO RECEPTORS:
;; Shc, Grb2, can bind to ERBB dimers at respective sites
;; Use "with-data-table" to avoid writing MANY rxns explicitly
;; Note that the proteins bind to only ONE monomer
;;

;; \wikidoc{shc_binds_erbb} 
(with-data-table 
 (:rows ($R1 $L1)  :cols ($R2 $L2) :cells ($Kf $Kr $Kc) :ignore _ )
 ((            (erbb1   *     )          (erbb2   NIL   )          (erbb3   *     )          (erbb4   *     )         )
  ((erbb1   *) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC))
  ((erbb2 NIL) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC))
  ((erbb3   *)          _                 _                         _                         _                       )
  ((erbb4   *)          _                 _                         _                         _                       ))
 ;; dimer binds SHC
 [{[[$R1 L.1 D.3 psite.p CSHC._ CGRB2.* ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2]] + [SHC CERB._] @ :inner <<->>
   [[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.* ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4]]}
  (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)] 
 ;; SHC gets activated ; suggested by Will/Laura
 [{[[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.* ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 asite.inactive]] ->>
   [[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.* ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 asite.active]]}
  (.set-rate-function 'mass-action $Kc)])
;; \wikidoc{shc_deactivation}
;; we use a hoky SHC phosphatase for SHC deactivation.
(with-data-table 
 (:rows ($R1 $L1)  :cols ($R2 $L2) :cells ($Kf $Kr $Kc) :ignore _ )
 ((            (erbb1  *     ) (erbb2  NIL     ) (erbb3  *     ) (erbb4  *     ))
  ((erbb1   *) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC))
  ((erbb2 NIL) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC))
  ((erbb3   *)        _               _        _               _       )
  ((erbb4   *)        _               _        _               _       ))
 [{[[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.* ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 asite.active CPASE._]] + [SHCPASE CSHC._] @ :inner <<->>
   [[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.* ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 asite.active CPASE.5][SHCPASE CSHC.5]]}
  (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)]
 [{[[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.* ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 asite.active CPASE.5][SHCPASE CSHC.5]] ->>
   [[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.* ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 asite.inactive CPASE._]] + [SHCPASE CSHC._] @ :inner }
  (.set-rate-function 'mass-action $Kc)])

;; \wikidoc{grb2_binds_erbb}
(with-data-table 
 (:rows ($R1 $L1) :cols ($R2 $L2) :cells ($Kf $Kr) :ignore _)
 ((            (erbb1  *)        (erbb2  NIL)      (erbb3  *)        (erbb4  *)       )
  ((erbb1   *) (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR))
  ((erbb2 NIL) (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR))
  ((erbb3   *) (KpbindF KpbindR) (KpbindF KpbindR)        _               _           )
  ((erbb4   *) (KpbindF KpbindR) (KpbindF KpbindR)        _               _           ))
 [{[[$R1 L.1 D.3 psite.p CSHC.* CGRB2._ ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2]] + [GRB2 CERB._ CSHC._ CSOS._ CERK._ ] @ :inner <<->>
   [[$R1 L.1 D.3 psite.p CSHC.* CGRB2.4 ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][GRB2 CERB.4 CSHC._ CSOS._ CERK._ ]]}
   (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)])

;; Grb2 binds to complex with bound Shc
;; \wikidoc{grb2_binds_erbbcomplex}
(with-data-table 
 (:rows ($R1 $L1) :cols ($R2 $L2) :cells ($Kf $Kr) :ignore _)
 ((            (erbb1  *     ) (erbb2  NIL     ) (erbb3  *     ) (erbb4  *     ))
  ((erbb1   *) (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR)   )
  ((erbb2 NIL) (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR)   )
  ((erbb3   *) (KpbindF KpbindR) (KpbindF KpbindR)        _               _          )
  ((erbb4   *) (KpbindF KpbindR) (KpbindF KpbindR)        _               _          ))
 [{[[$R1 L.1 D.3 psite.p CSHC.4 CGRB2._ ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2._]] + [GRB2 CERB._ CSHC._ CERK._ ] @ :inner <<->> 
   [[$R1 L.1 D.3 psite.p CSHC.4 CGRB2._ ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CERB._ CSHC.5 CERK._ ]]}
  (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)])

;; SOS binds to bound Grb2, whether Grb2 is bound to Shc or directly to the ErbB monomer. We account for both cases here.
;; \wikidoc{sos_binds_erbbcomplex}
(with-data-table 
 (:rows ($R1 $L1) :cols ($R2 $L2) :cells ($Kf1 $Kr1 $Kf2 $Kr2) :ignore _)
 ((             (erbb1  *)                    (erbb2  NIL)                    (erbb3  *)                    (erbb4  *)                   )
  ((erbb1   *)  (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR))
  ((erbb2 NIL)  (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR))
  ((erbb3   *)  (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR)               _                             _              )
  ((erbb4   *)  (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR)               _                             _              ))
 ;; SOS binds to [STUFF][ERBB:GRB2]
 [{[[$R1 L.1 D.3 psite.p CSHC.* CGRB2.4][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][GRB2 CERB.4 CSOS._]] + [SOS CGRB2._ CRAS._] @ :inner <<->> 
   [[$R1 L.1 D.3 psite.p CSHC.* CGRB2.4][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][GRB2 CERB.4 CSOS.5][SOS CGRB2.5 CRAS._]]}
  (.set-rate-function 'mass-action :fwd $Kf1 :rev $Kr1)]
 ;; SOS binds to [STUFF][ERBB:SHC:GRB2]
 [{[[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.*][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CSHC.5 CSOS._]] + [SOS CGRB2._ CRAS._] @ :inner <<->> 
   [[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.*][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CSHC.5 CSOS.6][SOS CGRB2.6 CRAS._]]}
  (.set-rate-function 'mass-action :fwd $Kf2 :rev $Kr2)])

;; The [dimer][grb2][sos] complex activates RAS
;; \wikidoc{ras_binds_erbbgrb2complex}
(with-data-table 
 (:rows ($R1 $L1) :cols ($R2 $L2) :cells ($Kf1 $Kr1 $Kf2 $Kr2) :ignore _)
 ((             (erbb1  *)                    (erbb2  NIL)                    (erbb3  *)                    (erbb4  *)                   )
  ((erbb1   *)  (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR))
  ((erbb2 NIL)  (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR))
  ((erbb3   *)                _                             _                             _                             _              )
  ((erbb4   *)                _                             _                             _                             _              ))
 ;; E + S <-> [E:S]
 [{[[$R1 L.1 D.3 psite.p CSHC.* CGRB2.4][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][GRB2 CERB.4 CSOS.5][SOS CGRB2.5 CRAS._]] + [RAS CSOS._  gtpsite.gdp __]  <<->>
   [[$R1 L.1 D.3 psite.p CSHC.* CGRB2.4][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][GRB2 CERB.4 CSOS.5][SOS CGRB2.5 CRAS.6][RAS CSOS.6 gtpsite.gdp]]}
  (.set-rate-function 'mass-action :fwd $Kf1 :rev $Kr1)]
 ;; [E:S] -> E + P
 [{[[$R1 L.1 D.3 psite.p CSHC.* CGRB2.4][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][GRB2 CERB.4 CSOS.5][SOS CGRB2.5 CRAS.6][RAS CSOS.6 gtpsite.gdp]]  <<->>
   [[$R1 L.1 D.3 psite.p CSHC.* CGRB2.4][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][GRB2 CERB.4 CSOS.5][SOS CGRB2.5 CRAS._]] + [RAS CSOS._ gtpsite.gtp __]}
  (.set-rate-function 'mass-action :fwd $Kf2 :rev $Kr2)])

;; The [dimer][shc][grb2][sos] complex activates RAS
;; \wikidoc{ras_binds_erbbshcgrb2complex}
(with-data-table 
 (:rows ($R1 $L1) :cols ($R2 $L2) :cells ($Kf1 $Kr1 $Kf2 $Kr2) :ignore _)
 ((             (erbb1  *)                    (erbb2  NIL)                    (erbb3  *)                    (erbb4  *)                   )
  ((erbb1   *)  (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR))
  ((erbb2 NIL)  (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR))
  ((erbb3   *)  (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR)               _                             _              )
  ((erbb4   *)  (KpbindF KpbindR KpbindF KpbindR) (KpbindF KpbindR KpbindF KpbindR)               _                             _              ))
 ;; E + S <-> [E:S]
 [{[[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.* ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CSHC.5 CSOS.6][SOS CGRB2.6 CRAS._]] + [RAS CSOS._ gtpsite.gdp __]  <<->>
   [[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.* ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CSHC.5 CSOS.6][SOS CGRB2.6 CRAS.7][RAS CSOS.7 gtpsite.gdp]]}
  (.set-rate-function 'mass-action :fwd $Kf1 :rev $Kr1)]
 ;; [E:S] -> E + P
 [{[[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.* ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CSHC.5 CSOS.6][SOS CGRB2.6 CRAS.7][RAS CSOS.7 gtpsite.gdp]] <<->>
   [[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.* ][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CSHC.5 CSOS.6][SOS CGRB2.6 CRAS._]] + [RAS CSOS._ gtpsite.gtp __]}
  (.set-rate-function 'mass-action :fwd $Kf2 :rev $Kr2)])

;; \wikidoc{ras_deactivation}
;; RAS Deactivation :
;; LAURA: ADD GAP  Ras deactivation
[{[RAS CGAP._ gtpsite.gtp]  + [GAP CRAS._] @ :inner <<->> [[RAS CGAP.1 gtpsite.gtp][GAP CRAS.1]]}
 (.set-rate-function 'mass-action :fwd 1 :rev 1)]
[{[[RAS CGAP.1 gtpsite.gtp][GAP CRAS.1]]  <<->> [RAS CGAP._ gtpsite.gdp]  + [GAP CRAS._] @ :inner}
 (.set-rate-function 'mass-action :fwd 1 :rev 1)]

;; Active RAS binds RAF and activates it
;;
;; RAF activation is VERY complex. We fudge it by a two step process:
;; First RAS activates a site in RAF (emulating the biochemical modification that RAS has on RAF)
;; Second RAF gets phosphorylated by an imaginary "RAF Kinase" which in this binds to the RAS:RAF complex.
;; E + S <-> [E:S]
;; \wikidoc{raf_binds_activeras}
;; First Step:
[{[RAS CSOS._ CRAF._ gtpsite.gtp] + [RAF CRAS._ CMEK._ psite.u] @ :inner <<->> [[RAS CSOS._ CRAF.1 gtpsite.gtp][RAF CRAS.1 CMEK._ psite.u]]}
 (.set-rate-function 'mass-action :fwd KpbindF :rev KpbindR)]
;; [E:S] + E <<-> [E:S:E]
[{[[RAS CSOS._ CRAF.1 gtpsite.gtp][RAF CRAS.1 CiRAFK._ psite.u]] + [iRAFK CRAF._] @ :inner <<->> 
  [[RAS CSOS._ CRAF.1 gtpsite.gtp][RAF CRAS.1 CiRAFK.2 psite.u][iRAFK CRAF.2]]}
 (.set-rate-function 'mass-action  :fwd KpbindF :rev KpbindR)]
;; Second Step:
[{[[RAS CSOS._ CRAF.1 gtpsite.gtp][RAF CRAS.1 CiRAFK.2 psite.u][iRAFK CRAF.2]] ->> 
  [RAS CSOS._ CRAF._ gtpsite.gtp] + [RAF CRAS._ CiRAFK._ psite.p] @ :inner + [iRAFK CRAF._] @ :inner}
 (.set-rate-function 'mass-action  KpbindF)]

;; RAF Deactivation:
;; \wikidoc{raf_deactivation}
[{[RAF CRAS._ CMEK._ CPP2A._ psite.p] + [PP2A C._] <<->> [[RAF CRAS._ CMEK._ CPP2A.1 psite.p][PP2A C.1]]}
 (.set-rate-function 'mass-action :fwd KpbindF :rev KpbindR)]
[{[[RAF CRAS._ CMEK._ CPP2A.1 psite.p][PP2A C.1]] ->> [RAF CRAS._ CMEK._ CPP2A._ psite.u] + [PP2A C._]}
 (.set-rate-function 'mass-action  KpbindF )]

;; Active RAF binds MEK and activates it in two phosphorylation steps
;; Again, this should probably be scaffolded... 
;;
;; \wikidoc{activeraf_binds_mek}
;; E + S <-> [E:S]
[{[RAF CRAS._ CMEK._ CPP2A._ psite.p] + [MEK CRAF._] <<->> [[RAF CRAS._ CMEK.1 CPP2A._ psite.p][MEK CRAF.1]]}
 (.set-rate-function 'mass-action :fwd KpbindF :rev KpbindR)]
;; [E:S] -> E + P
(with-data-table (:rows ($R1 $R2) :cols ($P1 $P2) :cells ($Kf $Kr) :ignore _)
  ((                         (psite1.p psite2.u)     (psite1.u psite2.p)     (psite1.p psite2.p)  )
   ((psite1.u psite2.u)    (KpbindF KpbindR) (KpbindF KpbindR)            _          )
   ((psite1.u psite2.p)               _                       _            (KpbindF KpbindR))
   ((psite1.p psite2.u)               _                       _            (KpbindF KpbindR)))
  [{[[RAF CRAS.1 CMEK._ CPP2A._ psite.p][MEK CRAF.1 CERK._ CPP2A._ $R1 $R2]] ->> [RAF CRAS._ CMEK._ CPP2A._ psite.p] + [MEK CRAF._ CERK._ CPP2A._ $P1 $P2]}
   (.set-rate-function 'mass-action $Kf)])
  
;; MEK deactivation (two step dephosphorylation)
;;
;; \wikidoc{mek_deactivation}
;; E + S <-> [ES]
(with-substitution-table 
    ((  $R1      $R2    $Kf $Kr) 
     (psite1.p psite2.u  KpbindF KpbindR)
     (psite1.u psite2.p  KpbindF KpbindR)
     (psite1.p psite2.p  KpbindF KpbindR))
  [{[MEK CRAF._ CERK._ CPP2A._ $R1 $R2] + [PP2A C._] <<->> [[MEK CRAF._ CERK._ CPP2A.1 $R1 $R2][PP2A C.1]]}
   (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)])
;; [ES] <-> E + P -- dephosphorylation
(with-data-table (:rows ($R1 $R2) :cols ($P1 $P2) :cells $Kc :ignore _)
  ((                       (psite1.u psite2.p)   (psite1.p psite2.u)   (psite1.u psite2.u))
   ((psite1.p psite2.u)     _                     _             Kdephospho )
   ((psite1.u psite2.p)     _                     _             Kdephospho )
   ((psite1.p psite2.p)     Kdephospho            Kdephospho    _          ))
  [{[[MEK CRAF._ CERK._ CPP2A.1 $R1 $R2][PP2A C.1]] ->> [MEK CRAF._ CERK._ CPP2A._ $P1 $P2] + [PP2A C._]}
   (.set-rate-function 'mass-action $Kc)])

;; Active MEK binds and activates ERK in a two phosphorylation steps
;; This should probably also be using some sort of scaffold.
;; \wikidoc{activemek_binds_erk}
;;
;; E + S <-> [E:S]
[{[MEK CRAF._ CERK._ psite1.p psite2.p] + [ERK CMEK._ CMKP._ ] <<->> [[MEK CRAF._ CERK.1 psite1.p psite2.p][ERK CMEK.1 CMKP._ ]]}
 (.set-rate-function 'mass-action :fwd 1 :rev 1)]
;; [E:S] <-> E + P -- phosphorylation (two steps)
(with-data-table (:rows ($R1 $R2) :cols ($P1 $P2) :cells $Kc :ignore _)
  ((                       (psite1.p psite2.u)    (psite1.u psite2.p)   (psite1.p psite2.p))
   ((psite1.u psite2.u)     KphoC                  KphoC                 _                )
   ((psite1.u psite2.p)             _                    _               KphoC            )
   ((psite1.p psite2.u)             _                    _               KphoC            ))
[{[[MEK CRAF._ CERK.1 psite1.p psite2.p][ERK CMEK.1 CMKP._ $R1 $R2]] ->> [MEK CRAF._ CERK._ psite1.p psite2.p] + [ERK CMEK._ CMKP._ $P1 $P2]}
 (.set-rate-function 'mass-action $Kc)])

;; ERK binds to SOS to deactivate it:
;; NOTICE this works b/c the default state is __
;; the previous rules assume that CSOS when not explicitly defined is in the __ state
;; Therefore, having something bound to the CSOS site effectively deactivates SOS
[{[ERK CSOS._ psite1.p psite2.p] + [[SOS CERK._ CRAS._ CGRB2.1][* *.1]] <<->>
  [[ERK CSOS.2 psite1.p psite2.p][SOS CERK.2 CRAS._ CGRB2.1][* *.1]]}
  (.set-rate-function 'mass-action :fwd KpbindF :rev KpbindR)]

;; ERK two step dephosphorylation
;;
;; \wikidoc{erk_deactivation}
;; E + S <-> [ES]
(with-substitution-table 
    ((  $R1      $R2    $Kf $Kr) 
     (psite1.p psite2.u  KpbindF KpbindR)
     (psite1.u psite2.p  KpbindF KpbindR)
     (psite1.p psite2.p  KpbindF KpbindR))
  [{[ERK CMEK._ CMKP._ $R1 $R2] + [MKP C._] <<->> [[ERK CMEK._ CMKP._ $R1 $R2][MKP C._]]}
   (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)])
;; [E:S] -> dephosphorylation (two steps)
(with-data-table (:rows ($R1 $R2) :cols ($P1 $P2) :cells $Kc :ignore _)
  ((                     (psite1.u psite2.p)   (psite1.p psite2.u)  (psite1.u psite2.u))
   ((psite1.p psite2.u)           _                     _            Kdephospho)
   ((psite1.u psite2.p)           _                     _            Kdephospho)
   ((psite1.p psite2.p)   Kdephospho            Kdephospho           _          ))
  [{[[ERK CMEK._ CMKP._ $R1 $R2][MKP C._]] ->> [ERK CMEK._ CMKP._ $P1 $P2] + [MKP C._]}
   (.set-rate-function 'mass-action  $Kc)])

;;; -------------------------------------------------------------------------;;;;
;;; END MAPK PATHWAY
;;; -------------------------------------------------------------------------;;;;

;;; -------------------------------------------------------------------------;;;;
;;; AKT PATHWAY
;;; -------------------------------------------------------------------------;;;;

;; BINDING OF GAB1 to ERBBs
;; The Gab1 Shp2 complex binds to the DIMER : [SHC] : GRB2 : SOS complex
;; This in turn activates PI3K. ******** IS GBSP always active??
(with-data-table 
 (:rows ($R1 $L1) :cols ($R2 $L2) :cells ($Kf1 $Kr1) :ignore _)
 ((             (erbb1  *)      (erbb2  NIL)    (erbb3  *)      (erbb4  *)     )
  ((erbb1   *)  (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR))
  ((erbb2 NIL)  (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR))
  ((erbb3   *)         _               _               _               _       )
  ((erbb4   *)         _               _               _               _       ))
;; DIMER:GRB2:SOS case
 [{[[$R1 L.1 D.3 psite.p CSHC.* CGRB2.4][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][GRB2 CERB.4 CSOS.5 CGBSP._][SOS CGRB2.5 CGBSP._]] + [GBSP CGRB2._ CSOS._] @ :inner <<->>
   [[$R1 L.1 D.3 psite.p CSHC.* CGRB2.4][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][GRB2 CERB.4 CSOS.5 CGBSP.6][SOS CGRB2.5 CGBSP.7][GBSP CGRB2.6 CSOS.7]]}
  (.set-rate-function 'mass-action :fwd $Kf1 :rev $Kr1)]
;; DIMER:SHC:GRB2:SOS case
[{[[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.*][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CSHC.5 CSOS.6 CGBSP._][SOS CGRB2.6 CGBSP._]] + [GBSP CGRB2._ CSOS._] @ :inner <<->>
  [[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.*][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CSHC.5 CSOS.6 CGBSP.7][SOS CGRB2.6 CGBSP.8][GBSP CGRB2.7 CSOS.8]]}
  (.set-rate-function 'mass-action :fwd $Kf1 :rev $Kr1)])

;; BINDING OF PI3K to ERBBs through GAB1
(with-data-table 
 (:rows ($R1 $L1) :cols ($R2 $L2) :cells ($Kf1 $Kr1) :ignore _)
 ((             (erbb1  *)      (erbb2  NIL)    (erbb3  *)      (erbb4  *)     )
  ((erbb1   *)  (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR))
  ((erbb2 NIL)  (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR) (KpbindF KpbindR))
  ((erbb3   *)         _               _               _               _       )
  ((erbb4   *)         _               _               _               _       ))
 [{[[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.*][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CSHC.5 CSOS.6 CGBSP.7][SOS CGRB2.6 CGBSP.8][GBSP CGRB2.7 CSOS.8 CPI3K._]] + [PI3K CGBSP._] @ :inner <<->>
   [[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.*][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CSHC.5 CSOS.6 CGBSP.7][SOS CGRB2.6 CGBSP.8][GBSP CGRB2.7 CSOS.8 CPI3K.9][PI3K CGBSP.9]]}
  (.set-rate-function 'mass-action :fwd $Kf1 :rev $Kr1)])

;; activation of PtdIns by PI3K bound indirectly to ERBB through GAB1
(with-data-table 
 (:rows ($R1 $L1) :cols ($R2 $L2) :cells ($Kf1 $Kr1 $Kc) :ignore _)
 ((             (erbb1  *)      (erbb2  NIL)    (erbb3  *)      (erbb4  *)     )
  ((erbb1   *)  (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC))
  ((erbb2 NIL)  (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC) (KpbindF KpbindR KpbindC))
  ((erbb3   *)         _               _               _               _       )
  ((erbb4   *)         _               _               _               _       ))
 [{[[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.*][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CSHC.5 CSOS.6 CGBSP.7][SOS CGRB2.6 CGBSP.8][GBSP CGRB2.7 CSOS.8 CPI3K.9][PI3K CGBSP.9 CPI._]] + [PtdIns CPI3K._ pstate.p2_4-5] @ :inner <<->>
   [[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.*][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CSHC.5 CSOS.6 CGBSP.7][SOS CGRB2.6 CGBSP.8][GBSP CGRB2.7 CSOS.8 CPI3K.9][PI3K CGBSP.9 CPI.10][PtdIns CPI3K.10 pstate.p2_4-5]]}
  (.set-rate-function 'mass-action :fwd $Kf1 :rev $Kr1)]
 [{[[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.*][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CSHC.5 CSOS.6 CGBSP.7][SOS CGRB2.6 CGBSP.8][GBSP CGRB2.7 CSOS.8 CPI3K.9][PI3K CGBSP.9 CPI.10][PtdIns CPI3K.10 pstate.p2_4-5]] ->>
   [[$R1 L.1 D.3 psite.p CSHC.4 CGRB2.*][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][SHC CERB.4 CGRB2.5][GRB2 CSHC.5 CSOS.6 CGBSP.7][SOS CGRB2.6 CGBSP.8][GBSP CGRB2.7 CSOS.8 CPI3K.9][PI3K CGBSP.9 CPI._]] + [PtdIns CPI3K._ pstate.p3_3-4-5]}
  (.set-rate-function 'mass-action $Kc)])

;; DIRECT BINDING OF PI3K to ERBBs
;; PI3K binds directly to ERBB3 and to some ERBB4 dimers
;; PI3K binds through an adaptor to ERBB1 and ERBB2 dimers (******** HOW TO ADDRESS THIS?*****)
;;
;; \wikidoc{pi3k_binds_erbb}
;; ONLY ERBB3 and ERBB4 have binding sites for PI3K
(with-data-table 
 (:rows ($R1 $L1)  :cols ($R2 $L2) :cells ($Kf $Kr) :ignore _)
 ((                (erbb1  *)        (erbb2  NIL)          (erbb3  *)            (erbb4  *)           )
  ((erbb1   *)      _                 _                     _                     _                   )
  ((erbb2 NIL)      _                 _                     _                     _                   )
  ((erbb3   *) (grb2sosF grb2sosR) (grb2sosF grb2sosR) (grb2sosF grb2sosR) (grb2sosF grb2sosR))
  ((erbb4   *) (grb2sosF grb2sosR) (grb2sosF grb2sosR) (grb2sosF grb2sosR) (grb2sosF grb2sosR)))
 [{[[$R1 L.1 D.3 psite.p CPI3K._][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2]] + [PI3K CERB._ ] @ :inner <<->>
   [[$R1 L.1 D.3 psite.p CPI3K.4][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][PI3K CERB.4 ]]}
  (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)])

;;  activation of PtdIns by PI3K directly bound to ERBB
(with-data-table
 (:rows ($R1 $L1)  :cols ($R2 $L2) :cells ($Kf $Kr) :ignore _)
 ((                (erbb1  *)        (erbb2  NIL)          (erbb3  *)            (erbb4  *)           )
  ((erbb1   *)      _                 _                     _                     _                   )
  ((erbb2 NIL)      _                 _                     _                     _                   )
  ((erbb3   *) (grb2sosF grb2sosR) (grb2sosF grb2sosR) (grb2sosF grb2sosR) (grb2sosF grb2sosR))
  ((erbb4   *) (grb2sosF grb2sosR) (grb2sosF grb2sosR) (grb2sosF grb2sosR) (grb2sosF grb2sosR)))
 ;; binding of ptdins to complex
 [{[[$R1 L.1 D.3 psite.p CPI3K.4][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][PI3K CERB.4 CPI._]] + [PtdIns CPI3K._ pstate.p2_4-5] @ :inner <<->>
   [[$R1 L.1 D.3 psite.p CPI3K.4][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][PI3K CERB.4 CPI.5][PtdIns CPI3K.5 pstate.p2_4-5]]}
  (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)]
 ;; phosphorylation of ptdins 
 [{[[$R1 L.1 D.3 psite.p CPI3K.4][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][PI3K CERB.4 CPI.5][PtdIns CPI3K.5 pstate.p2_4-5]] ->>
   [[$R1 L.1 D.3 psite.p CPI3K.4][$L1 *.1][$R2 L.2 D.3 CPASE._ **][$L2 *.2][PI3K CERB.4 CPI._ ]] + [PtdIns CPI3K._ pstate.p3_3-4-5]}
  (.set-rate-function 'mass-action $Kf)])

;; PtdIns:p3 to PtdIns:p2 by PTEN
;; COMMENT ABOUT SHIP1/2 converting PtdIns:p3 to PtdIns:p2_3-4 and why we ignore it
;; \wikidoc{pip3_dephosphorylation}
[{[PTEN CPI._ ] @ :inner + [PtdIns CPTEN._ pstate.p3_3-4-5 ]  <<->> [[PTEN CPI.1 ][PtdIns CPTEN.1 pstate.p3_3-4-5 ]]}
 (.set-rate-function 'mass-action :fwd KpbindF :rev Kpbindr )]
[{[[PTEN CPI.1 ][PtdIns CPTEN.1 pstate.p3_3-4-5 ]]  ->> [PTEN CPI._ ] @ :inner + [PtdIns CPTEN._ pstate.p2_4-5 ]}
 (.set-rate-function 'mass-action KphoC)]

;; inactive AKT binds to PIP3
;; 
;; \wikidoc{akt_binds_pip3}
[{[PtdIns CAKT._ pstate.p3_3-4-5] + [AKT CPI3K._ psite1.u psite2.u ] @ :inner <<->> 
  [[PtdIns CAKT.1 pstate.p3_3-4-5][AKT CPI3K.1 psite1.u psite2.u ]]}
 (.set-rate-function 'mass-action :fwd KpbindF :rev KpbindR)]

;; PDK binds to PIP3
;;
;; \wikidoc{pdk1_binds_pip3}
;; we decided to use only PDK as a generic protein for PDK and PDK2
[{[PtdIns CPDK._ pstate.p3_3-4-5] + [PDK CPI3K._ ] @ :inner <<->> [[PtdIns CPDK.1 pstate.p3_3-4-5 ][PDK CPI3K.1 ]]}
 (.set-rate-function 'mass-action :fwd KpbindF :rev KpbindR)]

;; inactive or partially active  AKT:PIP3 activated by PDK:PIP3 (****IS PDK IN ACTIVE STATE ALWAYS???****)
;; 
;; \wikidoc{aktpip3_pdkpip3_complex}
(with-data-table (:rows ($RP1 $RP2) :cols ($PP1 $PP2) :cells ($Kf $Kr $Kc) :ignore _)
  ((                       (psite1.p psite2.u)     (psite1.u psite2.p)     (psite1.p psite2.p)    )
   ((psite1.u psite2.u)    (KpbindF KpbindR KphoC) (KpbindF KpbindR KphoC)          _             )
   ((psite1.u psite2.p)               _                       _            (KpbindF KpbindR KphoC))
   ((psite1.p psite2.u)               _                       _            (KpbindF KpbindR KphoC)))
  ;; partial phosphorylation
  [{[[PtdIns CPDK.1 pstate.p3_3-4-5][PDK CPI3K.1 CAKT._]]+[[PtdIns CAKT.2 pstate.p3_3-4-5][AKT CPI3K.2 CPDK._ $RP1 $RP2]] <<->>
    [[PtdIns CPDK.1 pstate.p3_3-4-5][PDK CPI3K.1 CAKT.3][PtdIns CAKT.2 pstate.p3_3-4-5][AKT CPI3K.2 CPDK.3 $RP1 $RP2]]}
   (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)]
  ;; total phosphorylation
  [{[[PtdIns CPDK.1 pstate.p3_3-4-5][PDK CPI3K.1 CAKT.3][PtdIns CAKT.2 pstate.p3_3-4-5][AKT CPI3K.2 CPDK.3 $RP1 $RP2]]->>
    [[PtdIns CPDK.1 pstate.p3_3-4-5][PDK CPI3K.1 CAKT._]] + [[PtdIns CAKT.2 pstate.p3_3-4-5][AKT CPI3K.2 CPDK._ $PP1 $PP2]]}
   (.set-rate-function 'mass-action $Kc)])

;; 
;; active AKT binds to unphosphorylated RAF and deactivates it (AKT-MAPK crosstalk)
[{[[PtdIns CAKT.2 pstate.p3_3-4-5][AKT CPI3K.2 CRAF._  psite1.p psite2.p]] + [RAF CAKT._ psite.u] @ :inner <<->>
  [[PtdIns CAKT.2 pstate.p3_3-4-5][AKT CPI3K.2 CRAF.3  psite1.p psite2.p][RAF CAKT.3 psite.u]]}
  (.set-rate-function 'mass-action :rev KpbindF :fwd KpbindR)]

;; deactivation of AKT by generic phosphatase... (!!! see refs for this in sawyers 2002)
;; \wikidoc{akt_dephosphorylation}
;;
(with-data-table (:rows ($RP1 $RP2) :cols ($PP1 $PP2) :cells ($Kf $Kr $Kc) :ignore _)
  ((                     (psite1.u psite2.p)      (psite1.p psite2.u)      (psite1.u psite2.u)    )
   ((psite1.p psite2.u)           _                     _                  (KpbindF KpbindR KphoC))
   ((psite1.u psite2.p)           _                     _                  (KpbindF KpbindR KphoC))
   ((psite1.p psite2.p)  (KpbindF KpbindR KphoC) (KpbindF KpbindR KphoC)             _            ))
  [{[[PtdIns CAKT.3 pstate.p3_3-4-5][AKT CPI3K.3 CPP2A._ $RP1 $RP2]] + [PP2A C._] @ :inner <<->>
    [[PtdIns CAKT.3 pstate.p3_3-4-5][AKT CPI3K.3 CPP2A.1 $RP1 $RP2][PP2A C.1]]}
   (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)]
  [{[[PtdIns CAKT.3 pstate.p3_3-4-5][AKT CPI3K.3 CPP2A.1 $RP1 $RP2][PP2A C.1]]  ->>
    [[PtdIns CAKT.3 pstate.p3_3-4-5][AKT CPI3K.3 CPP2A._ $PP1 $PP2]] + [PP2A C._] @ :inner}
   (.set-rate-function 'mass-action $Kc)])

;;; -------------------------------------------------------------------------;;;;
;;; END AKT PATHWAY
;;; -------------------------------------------------------------------------;;;;
;;
;; \wikidoc{cell_spec_compartments}
(define cell [cell-location])
cell.dish.(contains [egf] [hrg])
cell.cytoplasm.(contains [ATP]   [DEP1] [DEP2] [DEP3] [DEP4] [SHC]  [SHCPASE] [GRB2] [GBSP] [SOS] [GAP] [RAF] 
			 [iRAFK] [MEK]  [ERK]  [PP2A] [MKP]  [PI3K] [PTEN]    [AKT]  [PDK])
cell.cell-membrane.(contains [erbB1] [erbB2] [erbB3] [erbB4] [RAS] [PtdIns])
cell.endo-membrane.(contains [erbB1] [erbB2] [erbB3] [erbB4] [RAS] [PtdIns])
; There is no need for next line since they are transported there automagically by littleb
;cell.inverse-endo-membrane.(contains [erbB1] [erbB2] [erbB3] [erbB4])
;;
;; SET INITIAL CONDITIONS:
;; \wikidoc{init_conditions}
;;
{[egf].(in cell.dish).moles.t0 :=               1.0E12   } ; (#/pL), 1.0E12/cell, v_cell = 1E-12L
{[hrg].(in cell.dish).moles.t0 :=               1.0E8    } ; (#/pL), 1.0E8/cell
{[erbB1].(in cell.cell-membrane).moles.t0 :=   70.00E3   } ; (#/pL), 70000.0/cell
{[erbB2].(in cell.cell-membrane).moles.t0 :=  162.00E3   } ; (#/pL), 162000.0/cell
{[erbB3].(in cell.cell-membrane).moles.t0 :=    4.44E3   } ; (#/pL), 4440.0/cell
{[erbB4].(in cell.cell-membrane).moles.t0 :=    1.23E3   } ; (#/pL), 740.0/cell 
{[ATP].(in cell.cytoplasm).moles.t0 :=          1.20E9   } ; (#/pL), 1.2E9/dish, dish = 1mL
{[DEP1].(in cell.cytoplasm).conc.t0 :=        50.0E3    } ; (#/pL)
{[DEP2].(in cell.cytoplasm).conc.t0 :=        50.0E3    } ; (#/pL)
{[DEP3].(in cell.cytoplasm).conc.t0 :=        50.0E3    } ; (#/pL)
{[DEP4].(in cell.cytoplasm).conc.t0 :=        50.0E3    } ; (#/pL)
{[SHC].(in cell.cytoplasm).conc.t0 := .5}
{[SHCPASE].(in cell.cytoplasm).conc.t0 := .5}
{[GRB2].(in cell.cytoplasm).conc.t0 := .5}
{[GBSP].(in cell.cytoplasm).conc.t0 := .5}
{[SOS].(in cell.cytoplasm).conc.t0 := .5}
{[RAS].(in cell.cell-membrane).conc.t0 := .5}     ; membrane
{[GAP].(in cell.cytoplasm).conc.t0 := .5} 
{[RAF].(in cell.cytoplasm).conc.t0 := .5} 
{[iRAFK].(in cell.cytoplasm).conc.t0 := .5} 
{[MEK].(in cell.cytoplasm).conc.t0 := .5}
{[ERK].(in cell.cytoplasm).conc.t0 := .5} 
{[PP2A].(in cell.cytoplasm).conc.t0 := .5}
{[MKP].(in cell.cytoplasm).conc.t0 := .5}
{[PI3K].(in cell.cytoplasm).conc.t0 := .5}
{[PtdIns].(in cell.cell-membrane).conc.t0 := 1}   ; membrane
{[PTEN].(in cell.cytoplasm).conc.t0 := .5}
{[AKT].(in cell.cytoplasm).conc.t0 := .5}
{[PDK].(in cell.cytoplasm).conc.t0 := .5}


;; \wikidoc{create_model}
;(include b/matlab/ode-translation)
;(create-ode-model "erb4" :ode-comments nil :overwrite t)
;(format t "SETTING MAX MONOMERS PER COMPLEX TO 0 ~%")
;(setf *max-monomers-per-complex* NIL)
;(include b/numerica/ode-translation)
;(setf *NUMERICA-RATE-STRING-MAX-LENGTH* NIL)
;(create-numerica-model "input/erb6_test" :ode-comments nil :overwrite t :vars (query species.moles)
;		       :sim-options '("CSVOUTPUT := TRUE" "DYNAMIC_REPORTING_INTERVAL := 1.0") :sim-steps 1800)
