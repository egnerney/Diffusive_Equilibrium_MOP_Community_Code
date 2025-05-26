;========================================================
;  HYP2F1(a, b, c, z)
;
;  Returns the value of the Gauss hypergeometric function
;      2F1(a,b;c;z)
;  using the standard Euler-Type integral representation:
;
;    2F1(a,b;c;z) = [Gamma(c)/(Gamma(b)*Gamma(c-b))] *
;                   ∫[ x^(b-1)*(1 - x)^(c-b-1)*(1 - z*x)^(-a) ] dx
;                         from x=0 to x=1
;                         
;   For dummy integration variable x (not to be confused with x defined elsewehere)
;   See: https://en.wikipedia.org/wiki/Hypergeometric_function#Integral_formulas
;   Good for Real(c)>Real(b)>0, and Real(z) < 1
;   which we have for kappa_perp > 1, kappa_par>1/2, B/B0>1, \Delta PE/(kappa_par*T_par) >-1
;   
;
;  The integration is performed using QPINT1D in IDL 
;  adaptively calculates an approximation result to a given definite integral
;  QPINT1D is based on the QUADPACK fortran package 
;
; 
;========================================================
function HYP2F1, a, b, c, z, $
  REL_ERROR=rel_error, STAT=stat

  compile_opt idl2
  on_error, 2  ; Return to caller if error

  ;--
  ; Build the IDL expression for the Euler type integrand:
  ;   For dummy integration variable x (not to be confused with x defined elsewehere)
  ;   integrand(x) = x^(b-1)*(1 - x)^(c-b-1)*(1 - z*x)^(-a)
  ;   
  ;
  ; We'll pass parameters a,b,c,z via PRIVATE array P = [a,b,c,z]
  ; hence we do:
  ;
  ;   a = P(0), b = P(1), c = P(2), z = P(3)
  ;   X is the integration variable
  ;--
  integrand_expr = $
    '( (X^(P(1) - 1.0D)) * ((1D - X)^(P(2) - P(1) - 1.0D)) * ((1D - P(3)*X)^(-P(0))) )'

  ;--
  ; Call QPINT1D on interval [0, 1]
  ;   - /EXPRESSION => integrand_expr is compiled as an expression
  ;   - PRIVATE => pass the array [a,b,c,z] to be accessible as P(...)
  ;--
  double_a  = double(a)
  double_b  = double(b)
  double_c  = double(c)
  double_z  = double(z)
  P = [ double_a, double_b, double_c, double_z ]

  result = qpint1d(integrand_expr, 0D, 1D, P, $
    /EXPRESSION, $
    EPSABS = 0, $ ; Thus value ignored and only use relative error tolerance below
    EPSREL = 1d-8, $
    LIMIT = 500, $ ; Max Sub intervals to achieve tolerance
    ERROR = error, $
    STATUS = status)

  ;--
  ; Multiply by the front gamma-factor:
  ;    alpha = Gamma(c)/(Gamma(b)*Gamma(c - b))
  ;--
  ;alpha = gamma(double_c) / ( gamma(double_b)*gamma(double_c - double_b) )
  ; exponential of log of gamma functions better behaved for large c, b, and c-b 
  ; So for large kappa values
  alpha = exp(LNGAMMA(double_c) - LNGAMMA(double_b) - LNGAMMA(double_c - double_b ))
  hyper_val = alpha * result

  ;-- Return optional error & status
  if n_elements(error)  ne 0 then rel_error  = error/result
  if n_elements(status) ne 0 then stat = status

  return, hyper_val
end