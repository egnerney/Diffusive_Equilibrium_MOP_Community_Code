;  diffusive_equilibrium.pro
;
;  @author: Edward (Eddie) G. Nerney
;
;Nerney (2025)
;
;V1.01 May 2025
;
;Diffusive Equilibrium Code for finding steady state plasma Densities
;at extrapolated location s given plasma parameters at s along same
;magnetic field line trace. s is arc length along field line.
;In assumed non inertial rotating frame with planet at \Omega in rad/s,
;or with fraction of corotation fcor given
;along fixed magnetic field lines given field line traces,
;including gravitational potential energy (only important close to planet),
;Centrifugal potential energy, and electrostatic potenial energy.
;
;Assumes dynamic accesibility to
;all locations along field line in steady state,
;plasma frozen to field line,
;So No cross‐field drift (beyond small gyration around B),
;No significant scattering to transport particles across field lines,
;equivalently for MHD Alfvén's theorem frozen-in flux theorem or large magnetic reynolds number,
;no collisions or sources or losses of particles (no loss cone),
;no drift or transport (though boundary condition can be informed by),
;the dist. function same at reference location s0 as s,
;energy conserved (no sources or losses of energy or dissipation),
;magnetic moment conserved (adiabatic) or timescales all >>> gyroperiod,
;Assumes neutrality spatial scales >>> debye length,
;and assumes plasma frame is the same as the rotating frame.
;
;T parallel or Tpar is given T parameter
;with T perpendicular or Tperp given by A*T.
;or A = Tperp/Tpar (Anistropy Factor)
;For Kappa distributions Temps are characterstic temps or T_c not "real" temps. T
;Where real temps are defined as as Tperp = pperp/n and Tpar = Tpar/n
;For standard kappas the relation is
;Tpar = Tpar_c/(1-3/(2*kappa)), Tperp = Tperp_c/(1-3/(2*kappa))
;For product kappas there can be different kappa_par and kappa_perp
;Tpar = Tpar_c/(1-1/(2*kappa_par)), Tperp = Tperp_c/(1 - 1/(kappa_perp))
;kappa parameter given is assumed to be parralel kappa or kappa_par
;kappa_par = kappa
;then perp. kappa value is kappa_perp = lambda * kappa
;lambda = kappa_perp/kappa_par
;

;***********************************************************
;**  1) Define a simple "Species" structure               **
;***********************************************************
FUNCTION SPECIES_INIT, name, m, q, T, A, kappa, lam, n, dist_type
  ;+
  ; Initialize a Species "object" as an IDL structure:
  ;   name       - string species name, e.g. 'O+'
  ;   m          - mass (in same units used consistently, e.g. kg or eV/(c^2), etc.)
  ;   q          - charge (in same consistent units, e.g. Coulombs if SI)
  ;   T          - characteristic temperature
  ;   A          - anisotropy factor (Tperp / Tpar)
  ;   kappa      - kappa parameter (parallel kappa by convention)
  ;   lam        - ratio (kappa_perp / kappa_par)
  ;   n          - number density at reference location
  ;   dist_type  - string indicating which distribution model to use
  ;-
  s = { $
    NAME:        name, $
    M:           m, $
    Q:           q, $
    T:           T, $
    A:           A, $
    KAPPA:       kappa, $
    LAM:         lam, $
    N:           n, $
    TYPE:        dist_type $
  }
  RETURN, s
END


;***********************************************************
;**  2) Define a helper to get planet parameters           **
;***********************************************************
FUNCTION DEFINE_PLANET, planet_
  ;+
  ; Returns planet-specific parameters for:
  ;   - RP:    planet radius (m)
  ;   - OMEGA: rotation rate (rad/s)
  ;   - GM:    gravitational parameter (m^3 / s^2)
  ;
  ; Example:
  ;    result = DEFINE_PLANET('jupiter')
  ;    rp = result[0]  &  omega = result[1]  &  gm = result[2]
  ;-
  planetlc = STRLOWCASE(planet_)
  CASE planetlc OF
    'earth': BEGIN
      RP = 3.3781d6
      OMEGA = 7.29210d-5
      GM = 3.9860d14
    END
    'jupiter': BEGIN
      RP = 7.1492d7
      OMEGA = 1.7585d-4
      GM = 1.26687d17
    END
    'saturn': BEGIN
      RP = 6.0268d7
      OMEGA = 1.6379d-4
      GM = 3.7931d16
    END
    'uranus': BEGIN
      RP = 2.5559d7
      OMEGA = 1.012d-4
      GM = 5.7940d15
    END
    'neptune': BEGIN
      RP = 2.4764d7
      OMEGA = 1.083d-4
      GM = 6.8351d15
    END
    ELSE: BEGIN
      PRINT, 'Planet ' + planet_ + ' not supported in DEFINE_PLANET.'
      stop
      RETURN, [-1,-1,-1]
    END
  ENDCASE

  RETURN, [RP, OMEGA, GM]
END


;***********************************************************
;**  3) DENSITY FUNCTIONS for different distributions      **
;***********************************************************
FUNCTION MAXWELLIAN_DENSITY, sp, deltaU, phi, Bratio

  ;  Isotropic Maxwellian distribution: n(s) = n0*exp[-(m*deltaU + q*Phi)/T ].
  ;
  ;  B-field ratio (Bratio) is not used in the isotropic Maxwellian formula,
  ;  but is included for uniform interface.
  ;
  ;  References:
  ;  -----------
  ;  Bagenal & Sullivan (1981), eq. for isotropic Maxwellian
  ;
  ;  Parameters:
  ;  -----------
  ;  species : Species
  ;  Contains n, m, q, T (SI or consistent units).
  ;  deltaU : float
  ;  Gravitational + centrifugal potential energy difference per unit mass.
  ;  Phi : float
  ;  Electrostatic potential difference (voltage).
  ;  Bratio : float
  ;  B(s) / B0. (Unused here but included for interface consistency.)
  ;
  ;  Returns:
  ;  --------
  ;  n_local : double
  ;  Number density at location s.
    
  RETURN, sp.n * EXP( - (sp.m*deltaU + sp.q*phi)/sp.T )
END

FUNCTION ANISO_MAXWELLIAN_DENSITY, sp, deltaU, phi, Bratio
  ;  Anisotropic Maxwellian distribution (standard form), where T_par is constant
  ;  along the field line, but T_perp(s) ~ 1 / [ A + (1 - A)*B0/B ].
  ;
  ;  n(s) = n0 / [ A0 + (1 - A0)*(B0/B) ] * exp[ - (m*deltaU + q*Phi)/T_par ]
  ;
  ;  References:
  ;  -----------
  ;  Huang & Birmingham (1992)
  ;
  ;  Parameters:
  ;  -----------
  ;  species : Species
  ;  Must have .A = T_perp0 / T_par0 at reference location.
  ;  deltaU : double
  ;  Potential energy difference per unit mass (grav + cent).
  ;  Phi : double
  ;  Electrostatic potential difference.
  ;  Bratio : double
  ;  B(s)/B0.
  ;
  ;  Returns:
  ;  --------
  ;  n_local : double
  ;  Computed number density at s.
  b0_over_b = 1d / Bratio
  denom = sp.A + (1d - sp.A)*b0_over_b
  IF denom LE 0d THEN denom = 1d-15  ; Just clamp if negative

  fac_exp = EXP( - (sp.m*deltaU + sp.q*phi)/sp.T )
  RETURN, (sp.n / denom)*fac_exp
END

FUNCTION ANISO_MAXWELLIAN_CONST_TPERP_DENSITY, sp, deltaU, phi, Bratio
  ;
  ;    Alternate anisotropic Maxwellian variant where T_perp is also assumed constant
  ;    (along with T_par). This leads to a different expression:
  ;
  ;    n(s) ~ (B(s)/B0)^(1 - A0) * exp[-(m*deltaU + q*Phi)/T_par ].
  ;
  ;    References:
  ;    -----------
  ;    Bagenal & Sullivan (1981), Bagenal (1994), Mei, Thorne, Bagenal (1995)
  ;
  ;    Parameters:
  ;    -----------
  ;    species : Species
  ;        .A is T_perp/T_par (constant).
  ;    deltaU : double
  ;        Potential energy difference per unit mass.
  ;    Phi : double
  ;        Electrostatic potential difference.
  ;    Bratio : double
  ;        B(s)/B0.
  ;
  ;    Returns:
  ;    --------
  ;    n_local : double
  ;        Number density at s.
  fac = Bratio^(1d - sp.A)
  expterm = EXP( - (sp.m*deltaU + sp.q*phi)/sp.T )
  RETURN, sp.n * fac * expterm
END

FUNCTION ANISO_KAPPA_DENSITY, sp, deltaU, phi, Bratio
  ;
  ;    Standard anisotropic kappa distribution with single kappa (same parallel & perp):
  ;      n(s) = [ n0 / ( A0 + (1-A0)(B0/B) ) ] * [ 1 + (m deltaU + q Phi) / (kappa T_par) ]^(0.5 - kappa)
  ;
  ;    References:
  ;    -----------
  ;    Meyer-Vernet et al. (1995); Moncuquet et al. (2002)
  ;
  ;    Parameters:
  ;    -----------
  ;    species : Species
  ;        Has .kappa, .A, .T, .m, .q
  ;    deltaU : float
  ;        Potential difference per unit mass (grav + cent).
  ;    Phi : double
  ;        Electrostatic potential difference.
  ;    Bratio : double
  ;        B(s)/B0.
  ;
  ;    Returns:
  ;    --------
  ;    n_local : double
  ;        Number density at s.
  b0_over_b = 1d / Bratio
  denom = sp.A + (1d - sp.A)*b0_over_b
  IF denom LE 0d THEN denom = 1d-15

  PEratio = (sp.m*deltaU + sp.q*phi)/(sp.kappa*sp.T)
  IF PEratio LE -1d THEN PEratio = -1d + 1d-15

  exponent_factor = (1d + PEratio)^(0.5d - sp.kappa)
  RETURN, (sp.n / denom)*exponent_factor
END

FUNCTION ANISO_PRODUCT_KAPPA_DENSITY, sp, deltaU, phi, Bratio
  
  ;  Product kappa distribution with possibly different kappa_par and kappa_perp:
  ;  f ~ [1 + (m v_par^2)/(2 kappa_par T_par)]^(-kappa_par - a) *
  ;  [1 + (m v_perp^2)/(2 kappa_perp T_perp)]^(-kappa_perp - b)
  ;  The integral solution typically involves hypergeometric functions.
  ;
  ;  This function includes a rough check for B/B0 >= 1. If B < B0, we set Bratio = 1 + eps
  ;  as integral diverges for B/B0<1
  ;
  ;  References:
  ;  -----------
  ;  Nerney (2025)
  ;
  ;  Parameters:
  ;  -----------
  ;  species : Species
  ;  Has .kappa (interpreted as kappa_par), .lam ( = kappa_perp / kappa_par ), etc.
  ;  deltaU : double
  ;  Phi : float
  ;  Bratio : double
  ;  B(s)/B0.
  ;
  ;  Returns:
  ;  --------
  ;  n_local : double
  ;  Number density at s (via hypergeometric for small kappa; a simpler expression for large).

  k0 = double(sp.kappa)
  kperp0 = double(k0*sp.lam)
  ktot = double(k0 + kperp0)
  n0 = double(sp.n)
  Tpar0 = double(sp.T)
  A0 = double(sp.A)

  ; Clamp B ratio
  IF Bratio LT 1d THEN Bratio = 1d + 1d-15

  beta_ = double((sp.m*deltaU + sp.q*phi)/(k0*Tpar0))
  IF beta_ LE -1d THEN beta_ = -1d + 1d-15

  x = double(A0*(Bratio - 1d))

  IF (Bratio NE 1d) THEN BEGIN
    ; Decide small vs large kappa:
    IF (k0 LT 17d) AND (kperp0 LT 17d) THEN BEGIN
      delta_val = (sp.lam*sp.A)*((Bratio - 1d)/(1d + beta_))
      ; We call HYP2F1(1., k_perp+1, k_tot+1.5, 1 - 1/delta_val)
      hf = HYP2F1(1d, kperp0 + 1d, ktot + 1.5d, 1d - 1d/double(delta_val))

      factor = k0/(ktot + 0.5d)
      outer = (n0/x)*((1d + beta_)^(0.5d - k0))
      RETURN, Bratio*outer*factor*hf
    ENDIF ELSE BEGIN
      ; Large kappa fallback
      exponent_factor = Bratio*((1d + beta_)^(0.5d - k0))
      hypfac = 1d/(x + 1d + beta_)
      RETURN, n0*exponent_factor*hypfac
    ENDELSE
  ENDIF ELSE BEGIN
    ; If B/B0=1
    RETURN, n0*((1d + beta_)^(0.5d - k0))
  ENDELSE
END

FUNCTION FRIED_EGG_DENSITY, sp, deltaU, phi, Bratio
  ;  "Fried-Egg" distribution: Maxwellian in parallel (kappa_par -> inf),
  ;  kappa in the perpendicular direction. The velocity-space integral leads
  ;  to exponential integrals E_n(x).
  ;
  ;  n(s) ~ [some prefactor]* exp[-(m*deltaU + q*Phi)/T_par + ... ] * E_{kappa+1}(kappa * x)
  ;  with x = A0*(B/B0 - 1).
  ;
  ;  For physically consistent solutions in the derived formula, we
  ;  require B >= B0.
  ;
  ;  References:
  ;  -----------
  ;  Nerney (2025)
  ;
  ;  Parameters:
  ;  -----------
  ;  species : Species
  ;  Has .kappa (the perp kappa), .A = T_perp0/T_par0, .T = T_par in eV or J, etc.
  ;  deltaU : double
  ;  Potential difference per unit mass
  ;  Phi : double
  ;  Electrostatic potential difference
  ;  Bratio : double
  ;  B(s)/B0
  ;
  ;  Returns:
  ;  --------
  ;  n_local : double
  ;  Number density at s, using the "fried-egg" limit.
  IF Bratio LT 1d THEN Bratio = 1d + 1d-15

  k0 = double(sp.kappa)   ; perp kappa
  A0 = double(sp.A)
  x = double(A0*(Bratio - 1d))
  z = double(k0*x)

  exponent = double(- (sp.m*deltaU + sp.q*phi)/sp.T)
  ; clamp exponent to avoid floating over/underflow
  IF ABS(exponent) GT 700d THEN exponent = (exponent>0d) ? 700d : -700d
  exp_factor = EXP(exponent)
  
  ; T assumed parallel T here
  ; e^z * E_{k0+1}[z] = U[1, 1 - k0, z] 

  IF (Bratio NE 1d) THEN BEGIN
    ; Small x or Bratio expansion to second order in x
    ; good to  0.9%. Much better away from k0 = 1
    ; but near k0 = 1 this expansion doesn't do well so keep to within x< 0.001
    ; Prevents nans from small Bratio>1 and x>0 in confluent hypergeometric call
    ; Some approximate expansions for large k0 or large x, else call HYPERU
    IF x LT 0.0001d THEN BEGIN
      factor = 1d
      ; if kappa=kappa_perp>30 and x>30. 0th order Approximation good to 0.1% that is same as Aniso-Maxwellian
      ; But easier and better approximation to just use first order approximation for k0>15.
      RETURN, sp.n*Bratio*factor*exp_factor
    ENDIF ELSE IF k0 GE 15d THEN BEGIN
      ; large k0 expansion
      ; First order in 1/k0 expansion
      ; k0 e^(k0*x)*E_{k0+1}[k0*x] ~1 / (1 + x) - x/(k0*(1+ x)^3)
      ; Good to 0.035% for k0>15. and all x values >0 => For Bratio>1
      factor = 1d/(1d + x) - x/(k0*(1d + x)^3)
      RETURN, sp.n*Bratio*factor*exp_factor
    ENDIF ELSE IF ((k0 LT 15d) AND (x lT 30d) AND (x gt 15d)) THEN BEGIN
      ; Third order in 1/k0 expansion
      ; k0 e^(k0*x)*E_{k0+1}[k0*x] ~1 / (1 + x) - x/(k0*(1+ x)^3) + x(2x-1)/(k0^2*(1+ x)^5)-  (6x^3 - 8x^2 +x)/(k0^3*(1+ x)^7)
      ; Didn't expand in x but 1/k0 but still works but for large x that x^7 gets huge so use below for x>30
      ; if x> 15. then for any k0 (though refer to above approximation if k0>15 for any x)
      ; Good to 0.02% for k0=1.001, better for larger k0 and x
      factor = 1d/(1d + x) - x / (k0*((x + 1d)^3d)) + x*(2.d*x - 1.)/((k0^2d)*((1.d + x)^5d)) -  (6.d*(x^3d) - 8.d*(x^2d) + x)/((k0^3d)*((1.d + x)^7d))
      RETURN, sp.n*Bratio*factor*exp_factor
    ENDIF ELSE IF ((k0 LT 15d) AND (x GE 30d)) THEN BEGIN
      ; fourth order in 1/x expansion
      ; For x>30 use this so that we don't have to calculate the x^7 above which gets large
      ; if x> 30. then for any k0 (though refer to above approximation if k0>15)
      ; Good to 0.01%
      numer = -6.d + k0 *(-11.d + 2.d*x - (6.d + (-3.d + x)*x) * k0 + (-1.d + x)*(1d + x^2d)*(k0^2d))
      denom = (x^4d)*(k0^3d)
      factor = numer/denom
      RETURN, sp.n*Bratio*factor*exp_factor
    ENDIF ELSE BEGIN
      func_value = k0 * hyperu(1d, 1d - k0, z)
      RETURN, sp.n*Bratio*func_value*exp_factor
    ENDELSE
  ENDIF ELSE BEGIN
    ; If Bratio = 1, x= 0, and we are at reference location s0
    ; but just in case \Delta PE  does not equal 0 where Bratio = 1
    ; then we use limiting form for non zero \delta PE same as Isotropic Maxwellian
    RETURN, sp.n*exp_factor
  ENDELSE
END


;*************************************************************
;**  4) TEMPERATURE FUNCTIONS for each distribution          **
;*************************************************************
FUNCTION MAXWELLIAN_TEMPERATURE, sp, deltaU, phi, Bratio
  ;  Temperature function for isotropic Maxwellian distribution.
  ;  Temperature is constant along field line.
  ;
  ;  Parameters:
  ;  -----------
  ;  species : Species
  ;  Contains T (SI or consistent units).
  ;  deltaU, Phi, Bratio : float
  ;  Not used but included for interface consistency.
  ;
  ;  Returns:
  ;  --------
  ;  T_par_local : double
  ;  Parallel temperature at location s.
  ;  T_perp_local : double
  ;  Perpendicular temperature at location s.
  ;  Assumed isotropic and therefor no variation either along field line 
  RETURN, [sp.T, sp.T]  ; [T_par, T_perp]
END

FUNCTION ANISO_MAXWELLIAN_TEMPERATURE, sp, deltaU, phi, Bratio
  ;  Temperature function for anisotropic Maxwellian distribution.
  ;  T_parallel is constant, T_perpendicular varies with B.
  ;
  ;  Parameters:
  ;  -----------
  ;  species : Species
  ;  Contains T and A.
  ;  deltaU, Phi, Bratio : double
  ;  Bratio = B(s)/B0.
  ;
  ;  Returns:
  ;  --------
  ;  T_par_local : double
  ;  Parallel temperature at location s.
  ;  T_perp_local : double
  ;  Perpendicular temperature at location s.
  T_par_local = sp.T
  b0_over_b = 1d/Bratio
  denom = sp.A + (1d - sp.A)*b0_over_b
  IF denom LE 0d THEN denom = 1d-15
  T_perp_local = sp.A*sp.T/denom
  RETURN, [T_par_local, T_perp_local]
END

FUNCTION ANISO_MAXWELLIAN_CONST_TPERP_TEMPERATURE, sp, deltaU, phi, Bratio
  ;
  ;  Temperature function for anisotropic Maxwellian with constant T_perp.
  ;  Both T_parallel and T_perpendicular are constant.
  ;
  ;  Parameters:
  ;  -----------
  ;  species : Species
  ;  Contains T and A.
  ;  deltaU, Phi, Bratio : float
  ;  Not used but included for interface consistency.
  ;
  ;  Returns:
  ;  --------
  ;  T_par_local : double
  ;  Parallel temperature at location s.
  ;  T_perp_local : double
  ;  Perpendicular temperature at location s.
  ; Assumed constant along field line... Inconsistent but a used formulation so kept
  T_par_local = sp.T
  T_perp_local = sp.A*sp.T
  RETURN, [T_par_local, T_perp_local]
END

FUNCTION ANISO_KAPPA_TEMPERATURE, sp, deltaU, phi, Bratio
  ;  Temperature function for standard anisotropic kappa distribution.
  ;
  ;  Parameters:
  ;  -----------
  ;  species : Species
  ;  Contains T, A, kappa, m, q.
  ;  deltaU : double
  ;  Gravitational + centrifugal potential energy difference per unit mass.
  ;  Phi : double
  ;  Electrostatic potential difference (voltage).
  ;  Bratio : double
  ;  B(s) / B0.
  ;
  ;  Returns:
  ;  --------
  ;  T_par_local : double
  ;  Parallel temperature at location s.
  ;  T_perp_local : double
  ;  Perpendicular temperature at location s.
  beta_ = (sp.m*deltaU + sp.q*phi)/(sp.kappa*sp.T)
  IF beta_ LE -1d THEN beta_ = -1d + 1d-15

  T_par_local = sp.T*(1d + beta_)

  b0_over_b = 1d/Bratio
  denom = sp.A + (1d - sp.A)*b0_over_b
  IF denom LE 0d THEN denom=1d-15

  T_perp_local = sp.A*sp.T*(1d + beta_)/denom
  RETURN, [T_par_local, T_perp_local]
END

FUNCTION ANISO_PRODUCT_KAPPA_TEMPERATURE, sp, deltaU, phi, Bratio
  ;
  ;    Temperature function for anisotropic product kappa distribution.
  ;    
  ;    Parameters:
  ;    -----------
  ;    species : Species
  ;        Contains T, A, kappa (interpreted as kappa_par), lam (kappa_perp/kappa_par).
  ;    deltaU : double
  ;        Gravitational + centrifugal potential energy difference per unit mass.
  ;    Phi : double
  ;        Electrostatic potential difference (voltage).
  ;    Bratio : double
  ;        B(s) / B0.
  ;        
  ;    Returns:
  ;    --------
  ;    T_par_local : double
  ;        Parallel temperature at location s.
  ;    T_perp_local : double
  ;        Perpendicular temperature at location s.
  IF Bratio LT 1d THEN Bratio=1d + 1d-15

  beta_ = (sp.m*deltaU + sp.q*phi)/(sp.kappa*sp.T)
  IF beta_ LE -1d THEN beta_=-1d + 1d-15

  kpar0 = sp.kappa
  kperp0 = kpar0*sp.lam
  A0=sp.A
  Tpar0=sp.T
  Tperp0=A0*Tpar0
  x = A0*(Bratio - 1d)

  g = 1d + beta_
  IF (Bratio NE 1d) THEN BEGIN
    ; For larger kappa values, use simpler approximations
    IF ((kpar0 GE 17d) AND (kperp0 GE 17d)) THEN BEGIN
      ; Large kappa approx
      T_par_local = Tpar0*g
      T_perp_local = Bratio*Tperp0*g/(x + g)
    ENDIF ELSE BEGIN
      ; Full calculation with hypergeometric functions
      delta_val = (sp.lam * sp.A) * ((Bratio - 1d) / (1d + beta_))
      kappa_tot = kpar0 + kperp0

      ; Calculate F_q values for Tperp
      ; F_q = 2_F_1(q,kperp0 + 0.5, kappa_tot + 1.5, 1.0 - 1./delta_val)
      F1 = hyp2f1(1.0d, kperp0 + 1.0d, kappa_tot + 1.5d, 1.0d - 1.d/delta_val)
      F2 = hyp2f1(2.0d, kperp0 + 1.0d, kappa_tot + 1.5d, 1.0d - 1.d/delta_val)
      F3q = hyp2f1(0.75d, kperp0 + 1.0d, kappa_tot + 1.5d, 1.0d - 1.d/delta_val)
      F7q = hyp2f1(1.75d, kperp0 + 1.0d, kappa_tot + 1.5d, 1.0d - 1.d/delta_val)


      ; Calculate regularized 2_F_1 values for Tpar
      ; H_q = 2_F_1(1.0, kperp0 + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
      H1h = hyp2f1(1.0d, kperp0 + 1.0d, kappa_tot + 0.5d, 1.0d - 1.d/delta_val)/gamma(kappa_tot + 0.5d)
      H3q = hyp2f1(1.0d, kperp0 + 1.0d, kappa_tot + 0.75d, 1.0d - 1.d/delta_val)/gamma(kappa_tot + 0.75d)
      H3h = F1/gamma(kappa_tot + 1.5d)
      H7q = hyp2f1(1.0d, kperp0 + 1.0d, kappa_tot + 1.75d, 1.0d - 1.d/delta_val)/gamma(kappa_tot + 1.75d)

      ; Calculate constants

      C1 = 4.0d * kappa_tot - 1.0d
      C2 = 2.0d * kappa_tot - 1.0d
      C3 = 4.0d * kpar0 - 1.0d
      C4 = 2.0d * kpar0 - 1.0d

      ; Calculate temperatures
      T_perp_num = (beta_ + 1.0d)*kpar0*Tperp0*(A0 + x)/(3.d * A0 * x)

      Dperp1 =  C1*F3q / (3.d*F7q)
      Dperp2 = C2*F1 / (2.d*F2)

      T_perp_denom = Dperp1 - Dperp2

      ; Value of Perp Temperature at s
      T_perp_local = T_perp_num / T_perp_denom


      T_par_num = 8.0d*(beta_ + 1.0d)*kpar0*Tpar0

      Dpar1 =  C3*C1*H7q / H3q
      Dpar2 = 2.d*C4*C2*H3h / H1h
      T_par_denom = Dpar1 - Dpar2

      ; Value of Parallel Temperature at s
      T_par_local = T_par_num / T_par_denom
    ENDELSE
  ENDIF ELSE BEGIN
    ; If B/B0 = 1 then we are at reference location s0
    T_par_local = Tpar0
    T_perp_local = Tperp0
  ENDELSE

  RETURN, [T_par_local, T_perp_local]
END

FUNCTION FRIED_EGG_TEMPERATURE, sp, deltaU, phi, Bratio
 
  ;  Temperature function for "Fried-Egg" distribution.
  ;
  ;  Parameters:
  ;  -----------
  ;  species : Species
  ;  Contains T, A, kappa (perp kappa).
  ;  deltaU : double
  ;  Gravitational + centrifugal potential energy difference per unit mass.
  ;  Phi : double
  ;  Electrostatic potential difference (voltage).
  ;  Bratio : float
  ;  B(s) / B0.
  ;
  ;  Returns:
  ;  --------
  ;  T_par_local : double
  ;  Parallel temperature at location s.
  ;  T_perp_local : double
  ;  Perpendicular temperature at location s.
  ; Check if B/B0 >= 1
  IF Bratio LT 1d THEN Bratio=1d + 1d-15


  
  
  
  ; T_parallel is constant for Fried-Egg
  T_par_local = sp.T
  T_perp0 = sp.T*sp.A
  A0 = sp.A
  k0 = sp.kappa
  ; For T_perpendicular, we use the derived formula with x and z
  x = A0 * (Bratio - 1.0d)
  z = k0 * x



  if (Bratio ne 1d) then begin
    if (k0 gt 100d) then begin
      ; 0th order approximation is just
      ; Tperp result is same as Aniso Maxwellian (Huang & Birmingham)
      ; factor = 1./(1. + x) as Bratio/(1. + x) = 1/(A0+(1-A0)B0/B)
      ; Good to 0.25 % for k0>100 for all x>0
      ; Better for larger k0 and larger x
      factor = 1.d/(1.d + x)
      T_perp_local = Bratio*T_perp0*factor
      
    endif else if (k0 ge 15d) then begin
      ; Good to 0.05% for all x>0
      ; Better for larger kappa and larger x
      factor = 1.d/(1.d + x) - x/(((1.d + x)^3.d)*k0)
      T_perp_local =  Bratio*T_perp0*factor
    endif else if ((x ge 7d) and (k0 lt 15d)) then begin
      ; For k0> 1, Bratio>1, and x>7. this 5th order approximation
      ; Good to better than 0.22% for all all k0>1
      ; better for larger x or k0
      C1 = (3.d + k0)*(4.d + k0 * (7.d + k0))
      C2 = -x*((1.d + k0)^2.d)*(4.d + k0)
      C3 = -( x^3.d) * (k0^2.d) *(1.d + k0) + ( x^2.d) * k0 *(1.d + k0) *(2.d + k0)
      C4 = ( x^4.d) * (k0^3.d) + C3 + C2 + C1
      numerator =  -4.d + k0 *C4
      denominator = ( x^5.d) * (k0^4.d)
      factor = numerator/denominator
      T_perp_local = Bratio*T_perp0*factor
    endif else begin
      ;for x<7 and k0 <15
      U5q = hyperu(k0 + 1.d, k0 + 5.d/4.d, z)
      U1q = hyperu(k0 + 1.d, k0 + 1.d/4.d, z)
      U1 = hyperu(k0 + 1.d, k0 + 1.d, z)
      U0 = hyperu(k0 + 1.d, k0, z)
      denominator = 4.d* (U5q/U1q) - 3.d *(U1/U0)
      factor = (1.d/(x*denominator))
      T_perp_local = Bratio * T_perp0 * factor
    endelse
     
  endif else begin
    T_perp_local = T_perp0
  endelse
  
  

  RETURN, [T_par_local, T_perp_local]
END


;*************************************************************
;**  5) KAPPA FUNCTIONS for each distribution               **
;*************************************************************
FUNCTION MAXWELLIAN_KAPPA, sp, deltaU, phi, Bratio
  ;
  ;  Kappa function for Maxwellian distribution.
  ;  Not applicable but included for interface consistency.
  ;
  ;  Returns:
  ;  --------
  ;  NaN, NaN : No kappa parameters for Maxwellian.
  ;  Though actually effectivly Infinite as Maxwellian is infinite kappa limit
  RETURN, [!VALUES.D_NAN, !VALUES.D_NAN]
END

FUNCTION ANISO_MAXWELLIAN_KAPPA, sp, deltaU, phi, Bratio
  ;
  ;  Kappa function for Aniso Maxwellian distribution.
  ;  Not applicable but included for interface consistency.
  ;
  ;  Returns:
  ;  --------
  ;  NaN, NaN : No kappa parameters for Maxwellian.
  ;  Though actually effectivly Infinite as Maxwellian is infinite kappa limit
  RETURN, [!VALUES.D_NAN, !VALUES.D_NAN]
END

FUNCTION ANISO_MAXWELLIAN_CONST_TPERP_KAPPA, sp, deltaU, phi, Bratio
  ;  Kappa function for Aniso Maxwellian constant Tperp distribution.
  ;  Not applicable but included for interface consistency.
  ;
  ;  Returns:
  ;  --------
  ;  NaN, NaN : No kappa parameters for Maxwellian.
  ;  Though actually effectivly Infinite as Maxwellian is infinite kappa limit
  RETURN, [!VALUES.D_NAN, !VALUES.D_NAN]
END

FUNCTION ANISO_KAPPA_KAPPA, sp, deltaU, phi, Bratio
  ;  Kappa function for standard anisotropic kappa distribution.
  ;  Kappa parameter is constant along field line.
  ;
  ;  Parameters:
  ;  -----------
  ;  species : Species
  ;  Contains kappa.
  ;  deltaU, Phi, Bratio : double
  ;  Not used but included for interface consistency.
  ;
  ;  Returns:
  ;  --------
  ;  kappa_par_local : double
  ;  Parallel kappa parameter at location s.
  ;  kappa_perp_local : double
  ;  Perpendicular kappa parameter at location s.
  ;  there is no par or perp in standard kappa only kappa
  ;  but include both for interface consistency.
  RETURN, [sp.kappa, sp.kappa]
END

FUNCTION ANISO_PRODUCT_KAPPA_KAPPA, sp, deltaU, phi, Bratio
  ;  Kappa function for anisotropic product kappa distribution.
  ;
  ;  Parameters:
  ;  -----------
  ;  species : Species
  ;  Contains kappa (kappa_par), lam (kappa_perp/kappa_par).
  ;  deltaU : double
  ;  Gravitational + centrifugal potential energy difference per unit mass.
  ;  Phi : double
  ;  Electrostatic potential difference (voltage).
  ;  Bratio : double
  ;  B(s) / B0.
  ;
  ;  Returns:
  ;  --------
  ;  kappa_par_local : double
  ;  Parallel kappa parameter at location s.
  ;  kappa_perp_local : double
  ;  Perpendicular kappa parameter at location s.
  IF Bratio LT 1d THEN Bratio=1d + 1d-15

  kpar0=sp.kappa
  kperp0 = kpar0*sp.lam
  
  beta_ = (sp.m * deltaU + sp.q * phi) / (sp.kappa * sp.T)
  if (beta_ le -1.) then beta_ = -1.0 + 1d-15

  IF (Bratio NE 1d) THEN BEGIN
    IF ((kpar0 GE 30d) AND (kperp0 GE 30d)) THEN BEGIN
      ; For larger kappa values, return constant kappa
      ; (approximation for large kappa matches standard kappa behavior)
      RETURN, [kpar0, kperp0]
    ENDIF ELSE BEGIN
      ; Full calculation with hypergeometric functions
      delta_val = (sp.lam * sp.A) * ((Bratio - 1.d) / (1.d + beta_))
      kappa_tot = kpar0 + kperp0

      ; Calculate F_q values for Tperp
      ; F_q = 2_F_1(q,kappa_perp + 0.5, kappa_tot + 1.5, 1.0 - 1./delta_val)
      F1 = hyp2f1(1.0d, kperp0 + 1.0d, kappa_tot + 1.5d, 1.0d - 1.d/delta_val)
      F2 = hyp2f1(2.0d, kperp0 + 1.0d, kappa_tot + 1.5d, 1.0d - 1.d/delta_val)
      F3q = hyp2f1(0.75d, kperp0 + 1.0d, kappa_tot + 1.5d, 1.0d - 1.d/delta_val)
      F7q = hyp2f1(1.75d, kperp0 + 1.0d, kappa_tot + 1.5d, 1.0d - 1.d/delta_val)


      ; Calculate regularized 2_F_1 values for Tpar
      ; H_q = 2_F_1(1.0, kperp0 + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
      H1h = hyp2f1(1.0d, kperp0 + 1.0d, kappa_tot + 0.5d, 1.0d - 1.d/delta_val)/gamma(kappa_tot + 0.5d)
      H3q = hyp2f1(1.0d, kperp0 + 1.0d, kappa_tot + 0.75d, 1.0d - 1.d/delta_val)/gamma(kappa_tot + 0.75d)
      H3h = F1/gamma(kappa_tot + 1.5d)
      H7q = hyp2f1(1.0d, kperp0 + 1.0d, kappa_tot + 1.75d, 1.0d - 1.d/delta_val)/gamma(kappa_tot + 1.75d)

      ; Calculate constants

      C1 = 4.0d * kappa_tot - 1.0d
      C2 = 2.0d * kappa_tot - 1.0d
      C3 = 4.0d * kpar0 - 1.0d
      C4 = 2.0d * kpar0 - 1.0d

      num_kappa_perp = 2.0d * C1* F3q * F2 / (C2* F1 * F7q) - 4.0d
      ; Calculate kappa_perp at s
      kappa_perp_local = 1.0d/num_kappa_perp + 1.0d



      ; Calculate kappa_par at s
      fac1 = C3*C1/(C4*C2)
      fac2 = fac1*H1h*H7q/(H3q*H3h)
      kappa_par_local = (fac2 -2.0d) / (2d*fac2 - 8.0d )
      RETURN, [kappa_par_local, kappa_perp_local]  
    ENDELSE
  ENDIF ELSE BEGIN
    ; If Bratio = 1 then at reference location s0
    RETURN, [kpar0, kperp0]
  ENDELSE
END


FUNCTION FRIED_EGG_KAPPA, sp, deltaU, phi, Bratio
  ;
  ;  Kappa function for "Fried-Egg" distribution.
  ;
  ;  Parameters:
  ;  -----------
  ;  species : Species
  ;  Contains kappa (perp kappa).
  ;  deltaU : double
  ;  Gravitational + centrifugal potential energy difference per unit mass.
  ;  Phi : double
  ;  Electrostatic potential difference (voltage).
  ;  Bratio : double
  ;  B(s) / B0.
  ;
  ;  Returns:
  ;  --------
  ;  None : Since there's no parallel kappa
  ;  (Maxwellian in parallel direction so parallel kappa effectivly infinite).
  ;  kappa_perp_local : double
  ;  Perpendicular kappa parameter at location s.
  
  IF (Bratio LT 1d) THEN Bratio=1d + 1d-15
  
  ; No parallel kappa parameter (Maxwellian in parallel direction) effectivly infinite
  kappa_par_local = !VALUES.D_INFINITY

  ;k0 = kappa_perp0
  k0 = sp.kappa
  A0 = sp.A
  ; For kappa_perpendicular, use the derived formula with x and z
  x = A0 * (Bratio - 1.0d)
  z = k0 * x
  
  if (Bratio ne 1d) then begin
    if (k0 ge 50d) then begin
      kappa_perp_local = k0

    endif else if (x ge 15d) then begin
      ;x>=15.
      ; Series expansion about x = infinity to 2nd order
      ; Good to 0.4% for k0 >1 and x>15
      term1 = ((x^2.d)* (k0^2.d))/(1.d + k0)
      term2 =  (289.d* (2.d + k0))/(16.d * x*k0 *(1.d + k0))
      term3 = (x*k0 *(27.d + 8.d*k0))/(4.d* (1.d + k0))
      kappa_perp_local = term1 + term2 + term3
    endif else begin
      ;x<15. and kappa_0 <50.0
      ; Via numerical quadrature integration defined confluent hypergeometric function U(a,b,z)
      U0 = hyperu(k0 + 1.d, k0, z)
      U1 = hyperu(k0 + 1.d, k0 + 1., z)
      U1q = hyperu(k0 + 1.d, k0 + 1.d/4.d, z)
      U5q = hyperu(k0 + 1.d, k0 + 5.d/4.d, z)

      fac = U0*U5q/(U1q*U1)
      kappa_perp_local = (0.75d - fac)/(1.0d - fac)
    endelse

  endif else begin
  ; If Bratio = 1 we are at reference location
    kappa_perp_local = k0
  endelse

  RETURN, [kappa_par_local, kappa_perp_local]
END


;*************************************************************
;**  6) Helper to compute deltaU (grav+cent)                 **
;*************************************************************
FUNCTION CALC_DELTAU, r0_in, lat0_in, r_in, lat_in, planet_, fcor
  ; r0_in, r_in in planet radii
  ; lat0_in, lat_in in degrees
  ; planet is string
  ; fcor is fraction of corotation
  

  p = DEFINE_PLANET(planet_) ; Set planet specific variables
  IF p[0] LT 0 THEN RETURN, 0d  ; invalid planet
  rp = p[0]
  omega0 = p[1]
  gm = p[2]

  ; scale rotation by fraction
  omega = omega0*fcor

  r0 = r0_in*rp
  lat0 = lat0_in*(!DPI/180d)
  r = r_in*rp
  lat = lat_in*(!DPI/180d)

  U0 = - gm/r0 - 0.5d*(r0^2d) * ((COS(lat0))^2d) * (omega^2d)
  U  = - gm/r  - 0.5d*(r^2d)  * ((COS(lat))^2d)  * (omega^2d)
  ; difference in gravitational + centrifdugal potential energy per mass
  ; along field line as a function of s
  RETURN, (U - U0)  ; \Delta U (Joules/kg)
END


;*************************************************************
;**  7) Find local maxima in a 1D array (deltaU vs lat)      **
;*************************************************************
FUNCTION FIND_LOCAL_MAXIMA, deltaUarr, latarr
  n = N_ELEMENTS(deltaUarr)
  IF n NE N_ELEMENTS(latarr) THEN RETURN, -1

  srt = SORT(latarr)
  dS  = deltaUarr[srt]
  lS  = latarr[srt]

  maxidx = -1L  ; sentinel

  FOR i = 1, n-2 DO BEGIN
    IF (dS[i] GT dS[i-1]) AND (dS[i] GT dS[i+1]) THEN BEGIN
      IF (N_ELEMENTS(maxidx) EQ 1) AND (maxidx[0] EQ -1) THEN $
        maxidx = i $
      ELSE $
        maxidx = [maxidx, i]
    ENDIF
  ENDFOR

  ; ------ updated test -------------
  IF (N_ELEMENTS(maxidx) EQ 1) AND (maxidx[0] EQ -1) THEN RETURN, -1
  ; ---------------------------------

  orig = srt[maxidx]
  RETURN, {idx:orig, lat:latarr[orig], val:deltaUarr[orig]}
END




;*************************************************************
;**  8) Solve for net neutrality => find phi via bisection   **
;*************************************************************
FUNCTION NET_CHARGE_DENSITY, species0, deltaU, phi, Bratio

  ; Compute total net charge from all species
  nspec = N_ELEMENTS(species0)
  nvals = DBLARR(nspec)

  FOR i=0, nspec-1 DO BEGIN
    disttype = STRLOWCASE(species0[i].type)
    nvals[i] = call_function( $
      disttype + '_density',species0[i], deltaU, phi, Bratio)
  ENDFOR

  
  ; so sum over all species
  totalCharge = 0d
  FOR i=0, nspec-1 DO totalCharge += species0[i].q * nvals[i]

  RETURN, totalCharge
END

  FUNCTION local_relative_error,  species0, deltaU, thisPhi, br
  ; netCharge:
  
  nQ = NET_CHARGE_DENSITY( species0, deltaU, thisPhi, br)
  nspec = N_ELEMENTS(species0)
  ; cold electron density is species0[nspec-1]
  distFunc_e = STRLOWCASE(species0[nspec-1].type) + '_density'
  ne_val = CALL_FUNCTION(distFunc_e, species0[nspec-1], deltaU, thisPhi, br)

  ; e- charge
  qe = species0[nspec-1].q

  ; If ne_val ~0 => avoid divide by zero => return huge
  IF (ABS(ne_val) LT 1d-30) THEN RETURN, 1d30

  RETURN, ABS(nQ / (qe * ne_val))
END



FUNCTION DIFF_EQ, pos_in, pos0_in, bratio, species0, deltaU, planet_, fcor
  ;+
  ; Solve for local densities "n" and potential "phi" to ensure net charge=0,
  ; given an array of species, a position s, reference s0, etc.
  ; Currently have the Following assumptions:
  ;   species0[nspec-1] = cold e- which is neutral with all species
  ; returns a structure with tags:
  ;   N:   array of densities [nspec]
  ;   PHI: potential
  ;   DU:  deltaU
  ;-
  nspec = N_ELEMENTS(species0)
  ;-------------------------------------------------
  ; 1) Force neutrality at reference => set cold-e
  ;    density to cancel the total ion charge.
  ;-------------------------------------------------
  total_ion_charge = 0d
  FOR i=0, nspec-2 DO total_ion_charge += (species0[i].q * species0[i].n)
  ; The last species must be cold e- => solve n_e = - (ion charge)/q_e
  species0[nspec-1].n = - total_ion_charge / species0[nspec-1].q

  
  ;-------------------------------------------------
  ; 2) Bracket the root in phi
  ;-------------------------------------------------
  phi1 = 0d
  phi2 = 0d
  dphi = 1d
 
  nq1 = NET_CHARGE_DENSITY(species0, deltaU, phi1, Bratio)
  nq2 = NET_CHARGE_DENSITY(species0, deltaU, phi2, Bratio)

  maxiters = 10000
  iter = 0
  ; Expand outward until we find nq1 & nq2 with opposite signs:
  WHILE (nq1*nq2 GT 0d AND iter LT maxiters) DO BEGIN
    phi1 -= dphi
    phi2 += dphi
    nq1 = NET_CHARGE_DENSITY(species0, deltaU, phi1, Bratio)
    nq2 = NET_CHARGE_DENSITY(species0, deltaU, phi2, Bratio)
    iter += 1
  ENDWHILE

  IF (iter GE maxiters) THEN BEGIN
    PRINT, 'Failed to bracket root in DIFF_EQ.'
    RETURN, { n: DBLARR(nspec), PHI: 0d, DU: deltaU }
  ENDIF

  ;-------------------------------------------------
  ; 3) Bisection until net charge ~ 0 or maxiters
  ;-------------------------------------------------
  tolerance = 1d-8
  it2 = 0
  phiRoot = 0.5d*(phi1 + phi2)
  rel1 = local_relative_error(species0, deltaU, phi1, Bratio)
  rel2 = local_relative_error(species0, deltaU, phi2, Bratio)
  relMid = local_relative_error(species0, deltaU, phiRoot, Bratio)


  WHILE (relMid GT tolerance AND it2 LT maxiters) DO BEGIN

    ; Evaluate sign of netCharge at midpoint
    nqMid = NET_CHARGE_DENSITY(species0, deltaU, phiRoot, Bratio)

    IF (nq1 * nqMid LT 0d) THEN BEGIN
      phi2 = phiRoot
      nq2  = nqMid
    ENDIF ELSE BEGIN
      phi1 = phiRoot
      nq1  = nqMid
    ENDELSE

    phiRoot = 0.5d*(phi1 + phi2)
    relMid = local_relative_error(species0, deltaU, phiRoot, Bratio)
    it2 += 1
  ENDWHILE

  ;-------------------------------------------------
  ; 4) Compute final densities at phiRoot
  ;-------------------------------------------------
  nvals = DBLARR(nspec)
  FOR i=0, nspec-1 DO BEGIN
    distFuncName = STRLOWCASE(species0[i].type) + '_density'
    nvals[i] = CALL_FUNCTION(distFuncName, species0[i], deltaU, phiRoot, Bratio)
  ENDFOR

  RETURN, { n: nvals, PHI: phiRoot, DU: deltaU }
END


;*************************************************************
;**  9) Compute Temperatures for each species at location    **
;*************************************************************
PRO CALC_TEMPS, species0, deltaU, phi, Bratio, Tpar, Tperp
  nspec = N_ELEMENTS(species0)
  Tpar = DBLARR(nspec)
  Tperp = DBLARR(nspec)

  FOR i=0, nspec-1 DO BEGIN
    disttype = STRLOWCASE(species0[i].type)
    tmp = call_function( $
      disttype + '_temperature',species0[i], deltaU, phi, Bratio)
    Tpar[i] = tmp[0]
    Tperp[i] = tmp[1]
  ENDFOR
END


;*************************************************************
;**  10) Compute Kappa Values for each species at location   **
;*************************************************************
PRO CALC_KAPPA_VALS, species0, deltaU, phi, Bratio, kpar, kperp
  nspec = N_ELEMENTS(species0)
  kpar = DBLARR(nspec)
  kperp = DBLARR(nspec)

  FOR i=0, nspec-1 DO BEGIN
    disttype = STRLOWCASE(species0[i].type)
    tmp = call_function( $
      disttype + '_kappa',species0[i], deltaU, phi, Bratio)
    kpar[i] = tmp[0]
    kperp[i] = tmp[1]
  ENDFOR
END


;*************************************************************
;**  11) Main dummy procedure                               **
;*************************************************************
PRO DIFFUSIVE_EQUILIBRIUM
  ;COMPILE_OPT idl2
  
  PRINT, 'diffusive_equilibrium.pro loaded.'
  PRINT, 'You can now call: '
  PRINT, '   species_s = SPECIES_INIT(...)'
  PRINT, '   result    = DIFF_EQ(..., species_s, deltaU, planet_ , fcor)'
  PRINT, '   etc.'
END
