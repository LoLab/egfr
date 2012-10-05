(in-package :b-user)
(defmacro define-rates (&rest var-defs)
  `(progn ,@(loop for (name value) in var-defs
  collect `(define ,name [[reference-var] :value, value]))))
(define-rates
(KLbindF     1.66E-14 )   ; (#/pL)^-1 s^-1   | ligand binding forward  1.66E-5 was changed to this to match the EC50
(KLbindR     3.0E-2   )   ; s^-1             | ligand binding reverse
(KHbindF     1.66E-14 )   ; (#/pL)^-1 s^-1   | ligand binding forward  1.66E-5 was changed to this to match the EC50
(KHbindR     3.0E-2   )   ; s^-1             | ligand binding reverse
(KATPbindF   1.66E-8  )   ; (#/pL)^-1 s^-1   | ATP binding is like ligand binding 
(KATPbindR   3.8E-3   )   ; s^-1             | ATP unbinding is a bit different
(KdimF       1.6E-6   )   ; (#/pL)^-1 s^-1   | receptor dimerization forward
(KdimR       1.6E-1   )   ; s^-1             | receptor dimerization reverse
(KphoC       1.0E-1   )   ; s^-1             | k_cat for phosphorylation
(KpbindF     1.0E-5   )   ; (#/pL)^-1 s^-1   | adaptor/protein binding forward
(KpbindR     1.0E-3   )  ; s^-1             | adaptor/protein binding reverse
(KDepho1C    3.0E-7   )   ; s^-1             | k_cat for dephosphorylation
(KDepho2C    3.0E-7   )   ; s^-1             | k_cat for dephosphorylation
(KDepho3C    6.0E-7   )   ; s^-1             | k_cat for dephosphorylation
(KDepho4C    4.0E-7   )   ; s^-1             | k_cat for dephosphorylation
(KintF       1.3E-2   )   ; s^-1             | receptor internalization forward
(KintR       5.0E-5   ))   ; s^-1             | receptor internalization reverse
;(KADPATP     1.0E6    )  ; s^-1             | ADP converts to ATP

;(KpbindC    1.0E-3 )      ; adaptor/protein binding reverse
;(grb2sosF   1.67E-5)      ; Grb2-sos binding forward
;(grb2sosR   6.0E-12))     ; Grb2-sos binding reverse
;(defparameter *Na* 6.022E23)
;(defparameter molecule molar/*Na*) 
;(define-rates
;($ligbindF   {1.0E7   / Molar / second})     ; ligand binding forward
;($ligbindR   {3.0E-2  / second})             ; ligand binding reverse
;($recdimF    {1.6E-6  / molecule / second})  ; receptor dimerization forward
;($recdimR    {1.6E-1  / second})             ; receptor dimerization reverse
;($recinternF {1.3E-3  / second})             ; receptor internalization forward
;($recinternR {5.0E-5  / second})             ; receptor internalization reverse
;($phospho    {1.0E-1  / second})             ; phosphorylation (no reverse)
;($dephospho  {1.0     / second})             ; dephosphorylation (no reverse)
;($protbindF  {1.0E-5  / molecule / second})  ; adaptor/protein binding forward
;($protbindR  {1.0E-3  / second})             ; adaptor/protein binding reverse
;($grb2sosF   {1.67E-5 / molecule / second})  ; Grb2-sos binding forward
;($grb2sosR   {6.0E-12 / molecule / second})) ; Grb2-sos binding reverse
