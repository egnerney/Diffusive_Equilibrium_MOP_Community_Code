;========================================================
;  HYPERU(a, b, z)
;
;  Returns the value of the Gauss hypergeometric function
;      U(a,b,z)
;  using the standard Laplace-Type integral representation:
;
;   U(a,b,z) = [1/(Gamma(a))] *
;                   ∫[ e^(-z*t)*t^(a-1)*(1 + t)^(b-a-1) ] dt
;                         from t=0 to t=Infinity
;      
;   See: https://en.wikipedia.org/wiki/Confluent_hypergeometric_function#Integral_representations
;   Good for Real(a)>0, and Real(z) > 0
;   which we have for kappa_perp = kappa > 1, B>B0, A0>0
;   Using Limiting form instead for B/B0=1 => x =A0*(B/B0 - 1)= 0 
;
;  The integration is performed using QPINT1D in IDL 
;  adaptively calculates an approximation result to a given definite integral
;  QPINT1D is based on the QUADPACK fortran package 
;
; 
;========================================================
function HYPERU, a, b, z, $
  REL_ERROR=rel_error, STAT=stat

  compile_opt idl2
  on_error, 2  ; Return to caller if error

  ;--
  ; Build the IDL expression for the Laplace type integrand:
  ;   
  ;   integrand(t) = exp(-z*t)*t^(a-1)*(1 + t)^(b-a-1)
  ;   
  ;
  ; We'll pass parameters a,b,z via PRIVATE array P = [a,b,z]
  ; hence we do:
  ;
  ;   a = P(0), b = P(1), z = P(2)
  ;   X is the integration variable (used t above to prevent confusion but QPINT1D needs it as X)
  ;   For dummy integration variable X (not to be confused with x = A0*(B/B0 - 1) defined elsewehere)
  ;-- exp(-z*t)*t^(a-1)*(1 + t)^(b-a-1)
  
  integrand_expr = $
    '(  exp(-P(2)*X)* (X^(P(0) - 1.0D)) * ((1D + X)^(P(1) - P(0) - 1.0D)) )'

  ;--
  ; Call QPINT1D on interval [0, Infinity]
  ;   - /EXPRESSION => integrand_expr is compiled as an expression
  ;   - PRIVATE => pass the array [a,b,z] to be accessible as P(...)
  ;--
  double_a  = double(a)
  double_b  = double(b)
  double_z  = double(z)
  P = [ double_a, double_b, double_z ]

  result = qpint1d(integrand_expr, 0D, !VALUES.D_INFINITY, P, $
    /EXPRESSION, $
    EPSABS = 0, $ ; Thus value ignored and only use relative error tolerance below
    EPSREL = 1d-8, $ ; relative error t
    LIMIT = 500, $ ; Max Sub intervals to achieve tolerance
    ERROR = error, $
    STATUS = status)

  ;--
  ; Multiply by the front gamma-factor:
  ;    alpha = 1/Gamma(a)
  ;--
  ; alpha = 1D/  gamma(double_a) 
  ; exponential of log of gamma functions better behaved for large a
  ; So for large kappa_0 = k0 values for U_q = U[k0+1,k0+q,z]
  alpha = exp(- LNGAMMA(double_a ))
  hyperU_val = alpha * result

  ;-- Return optional error & status
  if n_elements(error)  ne 0 then rel_error  = error/result
  if n_elements(status) ne 0 then stat = status

  return, hyperU_val
end