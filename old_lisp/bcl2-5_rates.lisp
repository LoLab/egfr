;; BCL2 MODEL RATES
(in-package :b-user)
(defmacro define-rates (&rest var-defs)
  `(progn ,@(loop for (name value) in var-defs
	       collect `(define ,name [[reference-var] :value, value]))))

(define-rates
(KBIDCCSP8F 	3.0E-8 ) ; (#/cell)^-1 sec^-1 								    
(KBIDCCSP8R 	1.0E-3 ) ; sec^-1										    	 
(KBIDCCSP8C 	1.0    ) ; sec^-1										    	 
(KBIDCBAXCF 	4.5E-8 ) ; (#/cell)^-1 sec^-1 								    
(KBIDCBAXCR 	1.0E-3 ) ; sec^-1										    	 
(KBIDCBAXCC 	1.0    ) ; sec^-1 										    	 
(KBIDMBAXCF 	2.2E-8 ) ; (#/cell)^-1 sec^-1 								    	 
(KBIDMBAXCR 	1.0E-3 ) ; sec^-1										    	 
(KBIDMBAXCC 	1.0    ) ; sec^-1 										    	 
(KBIDCBAKF 	2.0E-6 ) ; (#/cell)^-1 sec^-1           							    	 
(KBIDCBAKR 	1.0E-3 ) ; sec^-1           									    	 
(KBIDCBAKC 	1.0    ) ; sec^-1           									    	 
(KBIDMBAKF 	2.0E-6 ) ; (#/cell)^-1 sec^-1           							    	 
(KBIDMBAKR 	1.0E-3 ) ; sec^-1           									    	 
(KBIDMBAKC 	1.0    ) ; sec^-1           									    	 
(KBIMBAXCF 	0.0    ) ; (#/cell)^-1 sec^-1 2.4E-7	 							    
(KBIMBAXCR 	0.0    ) ; sec^-1 	       1.0E-3	 							    
(KBIMBAXCC 	0.0    ) ; sec^-1 	       1.0	 							    
(KBIMBAKF 	0.0    ) ; (#/cell)^-1 sec^-1 2.0E-6	            						    
(KBIMBAKR 	0.0    ) ; sec^-1 	       1.0E-3	 							    
(KBIMBAKC 	0.0    ) ; sec^-1             1.0	 							    
(KBAXAUTOACTF 	0.0    ) ; sec^-1 *****NOT SURE								    	 
(KBAXAUTOACTR 	0.0    ) ; sec^-1 *****NOT SURE								    	 
(KBAXMDIMF 	1.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KBAXMDIMR 	1.0e-3 ) ; sec^-1 										    
(KBAXMTETF 	1.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KBAXMTETR 	1.0e-3 ) ; sec^-1 										    
(KBAXMPORE      1.0    ) ; sec 										    
(KBAKAUTOACTF 	0.0    ) ; sec^-1 *****NOT SURE								    
(KBAKAUTOACTR 	0.0    ) ; sec^-1 *****NOT SURE								    
(KBAKDIMF 	2.0e-6 ) ; (#/cell)^-1 sec^-1            							    
(KBAKDIMR 	1.0e-3 ) ; sec^-1            									    
(KBAKTETF 	2.0e-6 ) ; (#/cell)^-1 sec^-1            							    
(KBAKTETR 	1.0e-3 ) ; sec^-1             								    
(KBAKPORE 	1.0    ) ; sec^-1       									    
(KBCLXLCBAXCF 	8.3e-10) ; (#/cell)^-1 sec^-1 								    
(KBCLXLCBAXCR 	1.0e-3 ) ; sec^-1										    
(KBCLXLMBAXCF 	8.3e-10) ; (#/cell)^-1 sec^-1 								    
(KBCLXLMBAXCR 	1.0e-3 ) ; sec^-1										    
(KBCLXLCBAXMF 	8.3e-10) ; (#/cell)^-1 sec^-1 								    
(KBCLXLCBAXMR 	1.0e-3 ) ; sec^-1										    
(KBCLXLMBAXMF 	8.3e-10) ; (#/cell)^-1 sec^-1 								    
(KBCLXLMBAXMR 	1.0e-3 ) ; sec^-1										    
(KMCL1CBAXCF 	1.0e-6 ) ; ******CHECK******									    
(KMCL1CBAXCR 	1.0e-3 ) ; ******CHECK******									    
(KMCL1CBAXMF 	1.0e-6 ) ; ******CHECK******									    
(KMCL1CBAXMR 	1.0e-3 ) ; ******CHECK******									    
(KMCL1MBAXCF 	1.0e-6 ) ; ******CHECK******									    
(KMCL1MBAXCR 	1.0e-3 ) ; ******CHECK******									    
(KMCL1MBAXMF 	1.0e-6 ) ; ******CHECK******									    
(KMCL1MBAXMR 	1.0e-3 ) ; ******CHECK******									    
(KBCL2BAXCF 	1.7e-9 ) ; (#/cell)^-1 sec^-1 								    
(KBCL2BAXCR 	1.0e-3 ) ; sec^-1										    
(KBCL2BAXMF 	1.0E-6 ) ; (#/cell)^-1 sec^-1 								    
(KBCL2BAXMR 	1.0e-3 ) ; sec^-1										    
(KBCLXLCBAKF 	3.3e-8 ) ; (#/cell)^-1 sec^-1									    
(KBCLXLCBAKR 	1.0e-3 ) ; sec^-1       									    
(KBCLXLMBAKF 	3.3e-8 ) ; (#/cell)^-1 sec^-1         							    
(KBCLXLMBAKR 	1.0e-3 ) ; sec^-1										    
(KMCL1CBAKF 	1.7e-6 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1CBAKR 	1.0e-3 ) ; sec^-1										    
(KMCL1MBAKF 	1.7e-6 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1MBAKR 	1.0e-3 ) ; sec^-1										    
(KBCL2BAKF 	1.7e-6 ) ; ******CHECK*******									    
(KBCL2BAKR 	1.0e-3 ) ; ******CHECK*******									    
(KBCLXLCBIDCF 	3.2e-7 ) ; (#/cell)^-1 sec^-1 								    
(KBCLXLCBIDCR 	5.5e-3 ) ; sec^-1										    
(KBCLXLMBIDCF 	3.2e-7 ) ; (#/cell)^-1 sec^-1 								    
(KBCLXLMBIDCR 	5.5e-3 ) ; sec^-1										    
(KBCLXLCBIDMF 	3.2e-7 ) ; (#/cell)^-1 sec^-1 								    
(KBCLXLCBIDMR 	5.5e-3 ) ; sec^-1										    
(KBCLXLMBIDMF 	3.2e-7 ) ; (#/cell)^-1 sec^-1 								    
(KBCLXLMBIDMR 	5.5e-3 ) ; sec^-1										    
(KMCL1CBIDCF 	1.7e-6 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1CBIDCR 	1.0e-3 ) ; sec^-1										    
(KMCL1CBIDMF 	1.7e-6 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1CBIDMR 	1.0e-3 ) ; sec^-1										    
(KMCL1MBIDCF 	1.7e-6 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1MBIDCR 	1.0e-3 ) ; sec^-1										    
(KMCL1MBIDMF 	1.7e-6 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1MBIDMR 	1.0e-3 ) ; sec^-1										    
(KBCL2BIDCF 	1.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KBCL2BIDCR 	1.0e-3 ) ; sec^-1										    
(KBCL2BIDMF 	9.1e-7 ) ; (#/cell)^-1 sec^-1 								    
(KBCL2BIDMR 	1.0e-3 ) ; sec^-1										    
(KBCLXLCBIMF 	5.3e-7 ) ; (#/cell)^-1 sec^-1 								    
(KBCLXLCBIMR 	2.3e-3 ) ; sec^-1										    
(KBCLXLMBIMF 	5.3e-7 ) ; (#/cell)^-1 sec^-1 								    
(KBCLXLMBIMR 	2.3e-3 ) ; sec^-1										    
(KBCL2BIMF 	5.0e-9 ) ; (#/cell)^-1 sec^-1 								    
(KBCL2BIMR 	1.4e-4 ) ; sec^-1										    
(KMCL1CBIMF 	1.3e-6 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1CBIMR 	2.2e-3 ) ; sec^-1										    
(KMCL1MBIMF 	1.3e-6 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1MBIMR 	2.2e-3 ) ; sec^-1										    
(KBCLXLCBADF 	6.8e-7 ) ; (#/cell)^-1 sec^-1 								    
(KBCLXLCBADR 	6.5e-4 ) ; sec^-1										    
(KBCLXLMBADF 	6.8e-7 ) ; (#/cell)^-1 sec^-1 								    
(KBCLXLMBADR 	6.5e-4 ) ; sec^-1										    
(KBCL2BADF 	1.3e-7 ) ; (#/cell)^-1 sec^-1 								    
(KBCL2BADR 	1.0e-3 ) ; sec^-1										    
(KMCL1CNOXAF 	1.3e-8 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1CNOXAR 	3.7    ) ; sec^-1										    
(KMCL1CNOXAC 	0.0001 ) ; sec^-1 ****NOT SURE								    
(KMCL1MNOXAF 	1.3e-8 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1MNOXAR 	3.7    ) ; sec^-1										    
(KMCL1MNOXAC 	0.0001 ) ; sec^-1 ****NOT SURE								    
(KMCL1CMULEF 	2.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1CMULER 	1.0e-3 ) ; sec^-1										    
(KMCL1CMULEC 	0.0001 ) ; sec^-1										    
(KMCL1MMULEF 	2.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1MMULER 	1.0e-3 ) ; sec^-1										    
(KMCL1MMULEC 	0.0001 ) ; sec^-1										    
(KCSP3CSP8F 	8.0e-8 ) ; (#/cell)^-1 sec^-1 								    
(KCSP3CSP8R 	1.0e-3 ) ; sec^-1										    
(KCSP3CSP8C 	1.0    ) ; sec^-1										    
(KMCL1CCSP8F 	2.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1CCSP8R 	1.0e-3 ) ; sec^-1										    
(KMCL1CCSP8C 	1.0    ) ; sec^-1										    
(KMCL1MCSP8F 	2.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1MCSP8R 	1.0e-3 ) ; sec^-1										    
(KMCL1MCSP8C 	1.0    ) ; sec^-1										    
(KMCL1CCSP3F 	2.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1CCSP3R 	1.0e-3 ) ; sec^-1										    
(KMCL1CCSP3C 	1.0    ) ; sec^-1										    
(KMCL1MCSP3F 	2.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KMCL1MCSP3R 	1.0e-3 ) ; sec^-1										    
(KMCL1MCSP3C 	1.0    ) ; sec^-1										    
(KCYTCBAXPF 	2.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KCYTCBAXPR 	1.0e-3 ) ; sec^-1										    
(KCYTCBAXPC 	10.0   ) ; sec^-1										    
(KSMACBAXPF 	2.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KSMACBAXPR 	1.0e-3 ) ; sec^-1										    
(KSMACBAXPC 	10.0   ) ; sec^-1										    
(KCYTCBAKPF 	2.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KCYTCBAKPR 	1.0e-3 ) ; sec^-1										    
(KCYTCBAKPC 	10.0   ) ; sec^-1										    
(KSMACBAKPF 	2.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KSMACBAKPR 	1.0e-3 ) ; sec^-1										    
(KSMACBAKPC 	10.0   ) ; sec^-1										    
(KBIDCBIDMF 	1.0e2  ) ; sec^-1										    
(KBIDCBIDMR 	0.1    ) ; sec^-1										    
(KBAXCBAXMF 	1.0e-2 ) ; sec^-1										    
(KBAXCBAXMR 	1.0e-2 ) ; sec^-1										    
(KBCLXLCBCLXLMF 1.0e-2 ) ; sec^-1										    
(KBCLXLCBCLXLMR 1.0e-3 ) ; sec^-1										    
(KMCL1CMCL1MF 	1.0e-2 ) ; sec^-1										    
(KMCL1CMCL1MR 	1.0e-3 ) ; sec^-1										    
(KCSP8ICSP8A 	7.7E-3 ) ; sec^-1 *****MADE UP assuming exponential growth -- Emiko****			    
(KMCL1CAMCL1CD 	1664.0 ) ; sec^-1 *****MADE UP assuming exponential decay with half life of 40 min ****	    
(KMCL1MAMCL1MD 	1664.0 ) ; sec^-1										    
(KTRAILRECF 	4.0e-7 ) ; (#/cell)^-1 sec^-1 								    
(KTRAILRECR 	1.0e-3 ) ; sec^-1 										    
(KTRAILRECC 	1.0e-5 ) ; sec^-1 										    
(KDRECCSP8F 	6.0E-8 ) ; (#/cell)^-1 sec^-1 								    
(KDRECCSP8R 	1.0e-3 ) ; sec^-1 										    
(KDRECCSP8C 	1.0    ) ; sec^-1 										    
(KCSP3CSP6F 	1.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KCSP3CSP6R 	1.0e-3 ) ; sec^-1 										    
(KCSP3CSP6C 	1.0    ) ; sec^-1 										    
(KCSP3XIAPF 	5.0e-6 ) ; (#/cell)^-1 sec^-1 								    
(KCSP3XIAPR 	1.0e-3 ) ; sec^-1 										    
(KCSP3XIAPC 	0.1    ) ; sec^-1 										    
(KCSP3PARPF 	2.0E-7 ) ; (#/cell)^-1 sec^-1 								    
(KCSP3PARPR 	1.0e-2 ) ; sec^-1 										    
(KCSP3PARPC 	1.0    ) ; sec^-1 										    
(KCSP6CSP8F 	3.0e-8 ) ; (#/cell)^-1 sec^-1 								    
(KCSP6CSP8R 	1.0e-3 ) ; sec^-1 										    
(KCSP6CSP8C 	1.0    ) ; sec^-1 										    
(KCYTCAPAFF 	5.0e-7 ) ; sec^-1 										    
(KCYTCAPAFR 	1.0e-3 ) ; sec^-1 										    
(KCYTCAPAFC 	1.0    ) ; sec^-1 										    
(KAPAFCSP9F 	5.0e-8 ) ; sec^-1 										    
(KAPAFCSP9R 	1.0e-3 ) ; sec^-1 										    
(KAPOPCSP3F 	1.0e-6 ) ; sec^-1 										    
(KAPOPCSP3R 	1.0e-3 ) ; sec^-1 										    
(KAPOPCSP3C 	1.0    ) ; sec^-1 										    
(KAPOPXIAPF 	1.0e-4 ) ; sec^-1 										    
(KAPOPXIAPR 	1.0e-3 ) ; sec^-1 										    
(KSMACXIAPF 	7.0e-6 ) ; sec^-1 										    
(KSMACXIAPR 	1.0e-3 ) ; sec^-1
)
;(defparameter *Na* 6.022E23)
;(defparameter molecule molar/*Na*) 


