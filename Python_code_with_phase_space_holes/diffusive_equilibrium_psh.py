#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@author: Edward (Eddie) G. Nerney

Nerney (2025)

Extended Diffusive Equilibrium Code with Phase Space Hole (PSH) support
for finding steady state plasma Densities at extrapolated location s 
given plasma parameters at s0 along same magnetic field line trace.

This version includes PSH density functions for cases where reference 
location s0 and target location s are on opposite sides of a local 
maximum in deltaU (choke point).
"""

__all__ = ['Species', 'define_planet', 'maxwellian_density', 'aniso_maxwellian_density', 'aniso_maxwellian_const_Tperp_density', 
           'aniso_kappa_density', 'aniso_product_kappa_density', 'fried_egg_density', 
           'maxwellian_density_psh', 'aniso_maxwellian_density_psh', 'aniso_maxwellian_const_Tperp_density_psh',
           'aniso_kappa_density_psh', 'aniso_product_kappa_density_psh', 'fried_egg_density_psh',
           'maxwellian_temperature', 'aniso_maxwellian_temperature',
           'aniso_maxwellian_const_Tperp_temperature', 'aniso_kappa_temperature', 'aniso_product_kappa_temperature', 'fried_egg_temperature',
           'maxwellian_kappa', 'aniso_maxwellian_kappa','aniso_maxwellian_const_Tperp_kappa', 'aniso_kappa_kappa', 'aniso_product_kappa_kappa', 'fried_egg_kappa',
           'density_functions', 'density_functions_psh', 'temperature_functions', 'kappa_functions', 
           'diff_eq','diff_eq_psh','find_local_maxima','calc_deltaU', 'calc_temps','calc_kappa_vals' ]

import numpy as np
import sys

# Numerical libraries
from scipy.optimize import bisect
from scipy.special import hyp2f1, gamma, hyperu, erfc, betainc
from scipy.integrate import quad

# CODATA RECOMMENDED VALUES OF THE FUNDAMENTAL PHYSICAL CONSTANTS: 2022, physics.nist.gov/constants
# Unnecessary precision
kg_per_u = 1.66053906892e-27  # Atomic mass unit (u or daltons, defined as 1/12 carbon 12 mass) in kg/u
u_per_kg = 1./kg_per_u        # u/kg
m_e = 5.485799090441e-4       # electron mass in u ~ 1.0/1822.9 u 

#Standard Periodic table values taking into account isotopic weighting and differences due to nuclear binding energy
m_H = 1.0078                  # Hydrogen mass in u
m_O = 15.999                  # Oxygen mass in u
m_S = 32.065                  # Sulfur mass in u
m_Na = 22.990                 # Sodium mass in u
# Even more Unwarranted Unnecessary precision (ignoring electron binding energy  decrease too)
m_Hp = m_H - m_e              # Singly (Fully) Ionized Hydrogen or H^{+} or proton mass in u
m_Op = m_O - m_e              # Singly Ionized Oxygen or O^{+} mass in u
m_O2p = m_O - 2.*m_e          # Doubly Ionized Oxygen or O^{++} mass in u
m_Sp = m_S - m_e              # Singly Ionized Sulfur or S^{+} mass in u
m_S2p = m_S - 2.*m_e          # Doubly Ionized Oxygen or S^{++} mass in u
m_S3p = m_S - 3.*m_e          # Triply Ionized Oxygen or S^{+++} mass in u
m_Nap = m_Na - m_e            # Singly Ionized Sodium or Na^{+} mass in u

ELEMENTARY_CHARGE = 1.602176634e-19  # Elementary charge in Coulombs
EV_TO_JOULE = ELEMENTARY_CHARGE   # Conversion from eV to Joules
JOULE_TO_EV = 1./ELEMENTARY_CHARGE 


###############################################################################
#                               CLASS DEFINITIONS
###############################################################################
class Species:
    """
    Container for plasma species parameters.

    Attributes:
    -----------
    name : str
        Species name (e.g., 'O+', 'e-', 'S++', etc.).
    m : float
        Mass in AMU (if prior to conversion) or in kg (if SI).
    q : float
        Charge in elementary charges (if prior to conversion) or in Coulombs (if SI).
    T : float
        Characteristic temperature in eV (if prior to conversion) or Joules (if SI).
    A : float
        Anisotropy factor A0 = T_perp / T_par at reference location, if relevant.
    kappa : float
        Kappa parameter. For Maxwellian, kappa -> infinity.
    lam : float
        λ = kappa_perp / kappa_par (only needed in product-kappa or more exotic models).
    n : float
        Number density at reference location (in cm^-3 or other consistent units).
    type : str
        Distribution type. Must match a key in density_functions dict (e.g., 'Maxwellian',
        'Aniso_Maxwellian', 'Aniso_kappa', 'Fried_Egg', etc.).
    """
    def __init__(self, name, m, q, T, A, kappa, lam, n, dist_type):
        self.name = name
        self.m = m
        self.q = q
        self.T = T
        self.A = A
        self.kappa = kappa
        self.lam = lam
        self.n = n
        self.type = dist_type


###############################################################################
#                            HELPER FUNCTIONS
###############################################################################
def define_planet(planet):
    """
    Returns planet-specific parameters: radius (m), rotation rate (rad/s),
    and GM (m^3/s^2).

    Parameters
    ----------
    planet : str
        Name of the planet

    Returns
    -------
    RP : float
        Planetary Equatorial radius (1 bar level) in meters.
    Omega : float
        Rotation rate in rad/s. 
    GM : float
        Gravitational parameter G*M_planet in m^3/s^2.
    """
    if planet.lower() == 'earth':      # From NASA Earth Fact Sheet
        RP = 3.3781e6                  # Equatorial Radius (1 bar level) in meters
        Omega = 7.29210e-5             # Rotation rate in rad/s 2.*np.pi/(23.9345*3600) using sidereal rotation period
        GM = 3.9860e14                 # G*M (m^3/s^2) or 0.39860*10^6 km^3/s^2 
        return RP, Omega, GM
    
    if planet.lower() == 'jupiter':    # From NASA Jupiter Fact Sheet
        RP = 7.1492e7                  # Equatorial Radius (1 bar level) in meters
        Omega = 1.7585e-4              # Rotation rate in rad/s 2.*np.pi/(9.9250*3600) using sidereal rotation period
        GM = 1.26687e17                # G*M (m^3/s^2) or 126.687*10^6 km^3/s^2
        return RP, Omega, GM

    elif planet.lower() == 'saturn':   # From NASA Saturn Fact Sheet
        RP = 6.0268e7                  # Equatorial Radius (1 bar level) in meters
        Omega = 1.6379e-4              # Rotation rate in rad/s 2.*np.pi/(10.656*3600) using sidereal rotation period
        GM = 3.7931e16                 # G*M (m^3/s^2) or 37.931*10^6 km^3/s^2
        return RP, Omega, GM
    
    if planet.lower() == 'uranus':     # From NASA Uranus Fact Sheet
        RP = 2.5559e7                  # Equatorial Radius (1 bar level) in meters
        Omega = 1.012e-4               # Rotation rate in rad/s 2.*np.pi/(17.24*3600) using sidereal rotation period
        GM = 5.7940e15                 # G*M (m^3/s^2) or 5.7940*10^6 km^3/s^2
        return RP, Omega, GM
    
    elif planet.lower() == 'neptune':  # From NASA Neptune Fact Sheet
        RP = 2.4764e7                  # Equatorial Radius (1 bar level) in meters
        Omega = 1.083e-4               # Rotation rate in rad/s 2.*np.pi/(16.11*3600) using sidereal rotation period
        GM = 6.8351e15                 # G*M (m^3/s^2) or 6.8351*10^6 km^3/s^2
        return RP, Omega, GM
    
    else:
        print(f"Planet {planet} is not supported in define_planet(). Exiting.")
        sys.exit(1)

###############################################################################
#                  DENSITY FUNCTIONS FOR DIFFERENT DISTRIBUTIONS
###############################################################################

# Original density functions
def maxwellian_density(species, deltaU, Phi, Bratio):
    """
    Isotropic Maxwellian distribution: n(s) ~ exp[-(m*deltaU + q*Phi)/T ].

    B-field ratio (Bratio) is not used in the isotropic Maxwellian formula,
    but is included for uniform interface.

    References:
    -----------
    Bagenal & Sullivan (1981), eq. for isotropic Maxwellian

    Parameters:
    -----------
    species : Species
        Contains n, m, q, T (SI or consistent units).
    deltaU : float
        Gravitational + centrifugal potential energy difference per unit mass.
    Phi : float
        Electrostatic potential difference (voltage).
    Bratio : float
        B(s) / B0. (Unused here but included for interface consistency.)

    Returns:
    --------
    n_local : float
        Number density at location s.
    """
    n_local = species.n * np.exp(-(species.m * deltaU + species.q * Phi) / species.T)
    return n_local


def aniso_maxwellian_density(species, deltaU, Phi, Bratio):
    """
    Anisotropic Maxwellian distribution (standard form), where T_par is constant
    along the field line, but T_perp(s) ~ 1 / [ A + (1 - A)*B0/B ].

    n(s) = n0 / [ A0 + (1 - A0)*(B0/B) ] * exp[ - (m*deltaU + q*Phi)/T_par ]

    References:
    -----------
    Huang & Birmingham (1992)

    Parameters:
    -----------
    species : Species
        Must have .A = T_perp0 / T_par0 at reference location.
    deltaU : float
        Potential energy difference per unit mass (grav + cent).
    Phi : float
        Electrostatic potential difference.
    Bratio : float
        B(s)/B0.

    Returns:
    --------
    n_local : float
        Computed number density at s.
    """
    B0_over_B = 1.0 / Bratio
    denom = species.A + (1.0 - species.A) * B0_over_B

    if denom <= 0:
        print("Warning: Non-physical denom in aniso_maxwellian_density. Setting denom -> eps.")
        denom = np.finfo(float).eps

    factor_exp = np.exp(-(species.m * deltaU + species.q * Phi) / species.T)
    n_local = (species.n / denom) * factor_exp
    return n_local


def aniso_maxwellian_const_Tperp_density(species, deltaU, Phi, Bratio):
    """
    Alternate anisotropic Maxwellian variant where T_perp is also assumed constant
    (along with T_par). This leads to a different expression:

    n(s) ~ (B(s)/B0)^(1 - A0) * exp[-(m*deltaU + q*Phi)/T_par ].

    References:
    -----------
    Bagenal & Sullivan (1981), Bagenal (1994), Mei, Thorne, Bagenal (1995)

    Parameters:
    -----------
    species : Species
        .A is T_perp/T_par (constant).
    deltaU : float
        Potential energy difference per unit mass.
    Phi : float
        Electrostatic potential difference.
    Bratio : float
        B(s)/B0.

    Returns:
    --------
    n_local : float
        Number density at s.
    """
    factor = Bratio**(1.0 - species.A)
    exponent_term = np.exp(-(species.m * deltaU + species.q * Phi) / species.T)
    n_local = species.n * factor * exponent_term
    return n_local


def aniso_kappa_density(species, deltaU, Phi, Bratio):
    """
    Standard anisotropic kappa distribution with single kappa (same parallel & perp):
      n(s) = [ n0 / ( A0 + (1-A0)(B0/B) ) ] * [ 1 + (m deltaU + q Phi) / (kappa T_par) ]^(0.5 - kappa)

    References:
    -----------
    Meyer-Vernet et al. (1995); Moncuquet et al. (2002)

    Parameters:
    -----------
    species : Species
        Has .kappa, .A, .T, .m, .q
    deltaU : float
        Potential difference per unit mass (grav + cent).
    Phi : float
        Electrostatic potential difference.
    Bratio : float
        B(s)/B0.

    Returns:
    --------
    n_local : float
        Number density at s.
    """
    B0_over_B = 1.0 / Bratio
    denom = species.A + (1.0 - species.A) * B0_over_B
    if denom <= 0.:
        print("Warning: Negative denom in aniso_kappa_density. Setting denom -> eps.")
        denom = np.finfo(float).eps

    PEratio = (species.m * deltaU + species.q * Phi) / (species.kappa * species.T)
    # For physical consistency in root-finding, clamp if PEratio <= -1
    if PEratio <= -1.:
        PEratio = -1.0 + np.finfo(float).eps
    
    exponent_factor = (1.0 + PEratio)**(0.5 - species.kappa)
    n_local = (species.n / denom) * exponent_factor
    
    return n_local


def aniso_product_kappa_density(species, deltaU, Phi, Bratio):
    """
    Product kappa distribution with possibly different kappa_par and kappa_perp:
       f ~ [1 + (m v_par^2)/(2 kappa_par T_par)]^(-kappa_par - a) *
            [1 + (m v_perp^2)/(2 kappa_perp T_perp)]^(-kappa_perp - b)
    The integral solution typically involves hypergeometric functions.

    This function includes a rough check for B/B0 >= 1. If B < B0, we set Bratio = 1 + eps
    as integral diverges for B/B0<1

    References:
    -----------
    Nerney (2025) 

    Parameters:
    -----------
    species : Species
        Has .kappa (interpreted as kappa_par), .lam ( = kappa_perp / kappa_par ), etc.
    deltaU : float
    Phi : float
    Bratio : float
        B(s)/B0.

    Returns:
    --------
    n_local : float
        Number density at s (via hypergeometric for small kappa; a simpler expression for large).
    """
    
    k0 = species.kappa
    k_perp0 = k0 * species.lam
    k_tot0 = k0 + k_perp0
    n0 = species.n
    T_par0 = species.T
    A0 = species.A
    #T_perp0 = T_par0*A0
    l0 = species.lam
    
    
    # For convergence of integral and positive density (physical)
    if Bratio < 1.0:
        Bratio = 1.0 + np.finfo(float).eps

    
    # For convergence of integral and positive density (physical)
    beta = (species.m * deltaU + species.q * Phi) / (k0 * T_par0 )
    if beta <= -1.:
        beta = -1.0 + np.finfo(float).eps


    
    #alpha = beta*k0 
    x = A0*(Bratio - 1.0)

    if Bratio != 1.0:
        
        # For smaller kappa (< ~17), attempt hypergeometric approach. Otherwise fallback.
        if (k0 < 17.) and (k_perp0 < 17.):
            delta = (l0 * A0) * ((Bratio - 1.) / (1. + beta))
            # Using 2F1 from scipy.special
 
            hf = hyp2f1(1., k_perp0 + 1.0, k_tot0  + 1.5, 1. - 1./delta)
    
            factor = k0 / (k_tot0  + 0.5)
            outer = (n0/ x) * ((1.0 + beta)**(0.5 - k0))
            n_local = Bratio *outer * factor * hf
    
        else:
            # Large kappa fallback: simpler expression ignoring the hypergeometric detail
            # Use approximation for large kappa_perp +k_par
            # k0*hf/(k0*(1+lambda0)+1/2) ~ 1 / (x + 1 + beta)
            # and good to a couple % for kappa_perp and kappa_par >17 (better for larger kappa)
            exponent_factor = Bratio * ((1. + beta)**(0.5 - k0))
            hypfac = 1. / (x + 1. + beta)
            n_local = n0 * exponent_factor * hypfac
    else:
        n_local = n0*((1. + beta)**(0.5 - k0))

    return n_local


def fried_egg_density(species, deltaU, Phi, Bratio):
    """
    "Fried-Egg" distribution: Maxwellian in parallel (kappa_par -> inf),
    kappa in the perpendicular direction. The velocity-space integral leads
    to exponential integrals E_n(x).

    n(s) ~ [some prefactor]* exp[-(m*deltaU + q*Phi)/T_par + ... ] * E_{kappa+1}(kappa * x)
    with x = A0*(B/B0 - 1).

    For physically consistent solutions in the derived formula, we
    require B >= B0.

    References:
    -----------
    Nerney (2025)

    Parameters:
    -----------
    species : Species
        Has .kappa (the perp kappa), .A = T_perp0/T_par0, .T = T_par in eV or J, etc.
    deltaU : float
        Potential difference per unit mass
    Phi : float
        Electrostatic potential difference
    Bratio : float
        B(s)/B0

    Returns:
    --------
    n_local : float
        Number density at s, using the "fried-egg" limit.
    """
    if Bratio < 1.0:
        print("Warning: B/B0 < 1 for fried_egg_density => non-convergent integral. Setting Bratio=1.")
        Bratio = 1.0 + np.finfo(float).eps
        
    k0 = species.kappa # k0 = kappa_0 = kappa_perp0 as there is no kappa_par0 for Fried Egg
    A0 = species.A # A0 = Tperp0/Tpar0
    x = A0 * (Bratio - 1.0)
    z = k0*x
    
    # exponent = \Delta PE / Tpar0
    exponent = -(species.m * deltaU + species.q * Phi) / species.T
    if np.abs(exponent) > 700:
        exponent = np.sign(exponent)*700.
    exp_factor = np.exp(exponent)
    
    # T assumed parallel T here
    # Though we can use mpmath our approximations are highly accurate and faster
    # In addition scipy only allows positive integer n for E_n(x) 
    # and scipy special incomplete upper gamma function  only allows positive arguments so that respresntation won't work
    # Scipy special allows for real valued a, b, and x for hyperu(a,b,x) and checks against mpmath at dps=100 show good aggremment for kappa0 and x<30 
    # e^z * E_{k0+1}[z] = U[1, 1 - k0, z] 
    
    
    if Bratio != 1.0:
        # Small x or Bratio expansion to second order in x 
        # good to  0.9%. Much better away from k0 = 1
        # but near k0 = 1 this expansion doesn't do well so keep to within x< 0.001
        # Prevents nans from small Bratio>1 and x>0 in confluent hypergeometric scipy special call
        if x< 0.0001:
            factor = 1.
            n_local = species.n * Bratio * factor * exp_factor
        # if kappa=kappa_perp>30 and x>30. 0th order Approximation good to 0.1% that is same as Aniso-Maxwellian
        # But easier and better approximation to just use first order approximation for k0>15.
        elif k0 >= 15.:
            # First order in 1/k0 expansion
            # k0 e^(k0*x)*E_{k0+1}[k0*x] ~1 / (1 + x) - x/(k0*(1+ x)^3) 
            # Good to 0.035% for k0>15. and all x values >0 => For Bratio>1
            factor = 1./(1. + x) - x / (k0*((x + 1.)**3.)) 
            n_local = species.n * Bratio * factor * exp_factor 
        elif (k0 < 15.) and (15.<x<30.):
            # Third order in 1/k0 expansion
            # k0 e^(k0*x)*E_{k0+1}[k0*x] ~1 / (1 + x) - x/(k0*(1+ x)^3) + x(2x-1)/(k0^2*(1+ x)^5)-  (6x^3 - 8x^2 +x)/(k0^3*(1+ x)^7)
            # Didn't expand in x but 1/k0 but still works but for large x that x^7 gets huge so use below for x>30
            # if x> 15. then for any k0 (though refer to above approximation if k0>15 for any x) 
            # Good to 0.02% for k0=1.001, better for larger k0 and x
            factor = 1./(1. + x) - x / (k0*((x + 1.)**3.)) + x*(2.*x - 1.)/((k0**2.)*((1.+ x)**5.)) -  (6.*(x**3.) - 8.*(x**2.) + x)/((k0**3.)*((1. + x)**7.))
            n_local = species.n * Bratio * factor * exp_factor
        elif (k0 < 15.) and (x>=30.):
            # fourth order in 1/x expansion
            # For x>30 use this so that we don't have to calculate the x^7 above which gets large
            # if x> 30. then for any k0 (though refer to above approximation if k0>15) 
            # Good to 0.01% 
            numer = -6. + k0 *(-11. + 2.*x - (6. + (-3. + x)*x) * k0 + (-1. + x)*(1 + x**2.)*(k0**2.))
            denom = (x**4.)*(k0**3.)
            factor = numer/denom
            n_local = species.n * Bratio * factor * exp_factor
        else:
            func_value = k0 * hyperu(1., 1. - k0, z)
            n_local = species.n * Bratio * func_value * exp_factor
    else:
        # If Bratio = 1, x= 0, and we are at reference location s0
        # but just in case \Delta PE  does not equal 0 where Bratio = 1 
        # then we use limiting form for non zero \delta PE same as Isotropic Maxwellian
        n_local = species.n*exp_factor
        
    return n_local


###############################################################################
#        PHASE SPACE HOLE (PSH) DENSITY FUNCTIONS
###############################################################################

def maxwellian_density_psh(species, deltaU, Phi, Bratio, deltaU_c, Phi_c):
    """
    Phase space hole version of isotropic Maxwellian distribution.
    
    n_psh = n0 * exp(-Delta_PE/T) * erfc(sqrt(Delta_PE_c/T))
    
    This implementation ensures proper handling of potential calculations to avoid
    numerical issues near choke points.
    
    Parameters:
    -----------
    species : Species
        Contains n, m, q, T (SI or consistent units).
    deltaU : float
        Gravitational + centrifugal potential energy difference per unit mass from reference.
    Phi : float
        Electrostatic potential difference (voltage) from reference.
    Bratio : float
        B(s) / B0.
    deltaU_c : float
        Gravitational + centrifugal potential energy difference at choke point from reference.
    Phi_c : float
        Electrostatic potential at choke point from reference.
        
    Returns:
    --------
    n_local : float
        Number density at location s.
    """
    # Calculate Delta_PE = m*(U - U0) + q*Phi
    Delta_PE = species.m * deltaU + species.q * Phi
    
    # Calculate Delta_PE_c = m*(U_c - U) + q*(Phi_c - Phi)
    Delta_PE_c = species.m * (deltaU_c - deltaU) + species.q * (Phi_c - Phi)
    
    # Calculate the PSH density
    exp_term = np.exp(-Delta_PE / species.T)
    
    erfc_term = erfc(np.sqrt(max(0.0, Delta_PE_c / species.T)))
    
    n_local = species.n * exp_term * erfc_term
    return n_local


def aniso_maxwellian_density_psh(species, deltaU, Phi, Bratio, deltaU_c, Phi_c):
    """
    Phase space hole version of anisotropic Maxwellian distribution.
    
    n_psh = (n0 * exp(-Delta_PE/T_par0) * erfc(sqrt(Delta_PE_c/T_par0))) / (A0 + (1-A0)*B0/B)
    
    where Delta_PE_c = m*(U_c - U) + q*(Phi_c - Phi)
    
    This implementation ensures proper handling of potential calculations to avoid
    numerical issues near choke points.
    
    Parameters:
    -----------
    species : Species
        Must have .A = T_perp0 / T_par0 at reference location.
    deltaU : float
        Potential energy difference per unit mass (grav + cent) from reference: U - U0
    Phi : float
        Electrostatic potential difference from reference.
    Bratio : float
        B(s)/B0.
    deltaU_c : float
        Potential energy difference at choke point from reference: U_c - U0
    Phi_c : float
        Electrostatic potential at choke point from reference.
        
    Returns:
    --------
    n_local : float
        Computed PSH number density at s.
    """
    # Calculate Delta_PE = m*(U - U0) + q*Phi
    Delta_PE = species.m * deltaU + species.q * Phi
    
    # Calculate Delta_PE_c = m*(U_c - U) + q*(Phi_c - Phi)
    # Since deltaU = U - U0 and deltaU_c = U_c - U0
    # Then U_c - U = deltaU_c - deltaU
    Delta_PE_c = species.m * (deltaU_c - deltaU) + species.q * (Phi_c - Phi)
    
    B0_over_B = 1.0 / Bratio
    denom = species.A + (1.0 - species.A) * B0_over_B
    
    if denom <= 0:
        print("Warning: Non-physical denom in aniso_maxwellian_density_psh. Setting denom -> eps.")
        denom = np.finfo(float).eps
    
    # Calculate the PSH density
    exp_term = np.exp(-Delta_PE / species.T)
  
    erfc_term = erfc(np.sqrt(max(0.0, Delta_PE_c / species.T)))
    
    n_local = (species.n * exp_term * erfc_term) / denom
    return n_local


def aniso_maxwellian_const_Tperp_density_psh(species, deltaU, Phi, Bratio, deltaU_c, Phi_c):
    """
    Phase space hole version of anisotropic Maxwellian with constant T_perp.
    
    For this distribution, we apply the same erfc correction as standard aniso-Maxwellian.
    
    This implementation ensures proper handling of potential calculations to avoid
    numerical issues near choke points.
    
    Parameters:
    -----------
    species : Species
        .A is T_perp/T_par (constant).
    deltaU : float
        Potential energy difference per unit mass from reference.
    Phi : float
        Electrostatic potential difference from reference.
    Bratio : float
        B(s)/B0.
    deltaU_c : float
        Potential energy difference at choke point from reference.
    Phi_c : float
        Electrostatic potential at choke point from reference.
        
    Returns:
    --------
    n_local : float
        Number density at s.
    """
    # Calculate Delta_PE = m*(U - U0) + q*Phi
    Delta_PE = species.m * deltaU + species.q * Phi
    
    # Calculate Delta_PE_c = m*(U_c - U) + q*(Phi_c - Phi)
    Delta_PE_c = species.m * (deltaU_c - deltaU) + species.q * (Phi_c - Phi)
    
    factor = Bratio**(1.0 - species.A)
    exp_term = np.exp(-Delta_PE / species.T)
    
    erfc_term = erfc(np.sqrt(max(0.0, Delta_PE_c / species.T)))
    n_local = species.n * factor * exp_term * erfc_term
    return n_local


def aniso_kappa_density_psh(species, deltaU, Phi, Bratio, deltaU_c, Phi_c):
    """
    Phase space hole version of anisotropic kappa distribution.
    
    n_psh = (n0 * (1+beta)^(0.5-kappa) * I_{(beta+1)/(beta_c+beta+1)}(kappa-0.5, 0.5)) / (A0 + (1-A0)*B0/B)
    
    where:
    - beta = Delta_PE / (kappa * T_par0)
    - beta_c = Delta_PE_c / (kappa * T_par0)
    - I_x(a,b) is the regularized incomplete beta function
    
    This implementation ensures proper handling of potential calculations to avoid
    numerical issues near choke points.
    
    Parameters:
    -----------
    species : Species
        Has .kappa, .A, .T, .m, .q
    deltaU : float
        Potential difference per unit mass (grav + cent) from reference.
    Phi : float
        Electrostatic potential difference from reference.
    Bratio : float
        B(s)/B0.
    deltaU_c : float
        Potential difference at choke point from reference.
    Phi_c : float
        Electrostatic potential at choke point from reference.
        
    Returns:
    --------
    n_local : float
        Number density at s.
    """
    B0_over_B = 1.0 / Bratio
    denom = species.A + (1.0 - species.A) * B0_over_B
    if denom <= 0.:
        print("Warning: Negative denom in aniso_kappa_density_psh. Setting denom -> eps.")
        denom = np.finfo(float).eps
    
    # Calculate beta and beta_c
    Delta_PE = species.m * deltaU + species.q * Phi
    Delta_PE_c = species.m * (deltaU_c - deltaU) + species.q * (Phi_c - Phi)
    
    beta = Delta_PE / (species.kappa * species.T)
    
    # Clamp beta if needed for physical consistency
    if beta <= -1.:
        beta = -1.0 + np.finfo(float).eps
    
    beta_c = Delta_PE_c / (species.kappa * species.T)
    
    # Calculate the regularized incomplete beta function argument
    x = (beta + 1.0) / (beta_c + beta + 1.0)
    
    # Ensure x is in valid range [0,1]
    x = max(0.0, min(1.0, x))
    
    # Calculate I_x(kappa - 0.5, 0.5) using scipy's betainc
    # betainc(a, b, x) computes the regularized incomplete beta function
    inc_beta = betainc(species.kappa - 0.5, 0.5, x)
    
    power_term = (1.0 + beta)**(0.5 - species.kappa)
    
    n_local = (species.n * power_term * inc_beta) / denom
    
    return n_local


def aniso_product_kappa_density_psh(species, deltaU, Phi, Bratio, deltaU_c, Phi_c):
    """
    Phase space hole version of anisotropic product kappa distribution.
    
    n_psh = n_full - n2_psh
    
    where n_full is the regular aniso_product_kappa_density and n2_psh involves
    an integral with hypergeometric functions.
    
    This implementation ensures proper handling of potential calculations to avoid
    numerical issues near choke points.
    
    Parameters:
    -----------
    species : Species
        Has .kappa (kappa_par), .lam (kappa_perp/kappa_par), etc.
    deltaU : float
        Potential difference per unit mass from reference.
    Phi : float
        Electrostatic potential difference from reference.
    Bratio : float
        B(s)/B0.
    deltaU_c : float
        Potential difference at choke point from reference.
    Phi_c : float
        Electrostatic potential at choke point from reference.
        
    Returns:
    --------
    n_local : float
        PSH number density at s.
    """
    # First calculate the full density
    n_full = aniso_product_kappa_density(species, deltaU, Phi, Bratio)
    
    # Calculate parameters
    k0 = species.kappa  # kappa_par
    k_perp0 = k0 * species.lam
    k_tot0 = k0 + k_perp0
    A0 = species.A
    T_par0 = species.T
    T_perp0 = T_par0 * A0
    
    # For convergence
    if Bratio < 1.0:
        Bratio = 1.0 + np.finfo(float).eps
    
    # Calculate beta and beta_c
    Delta_PE = species.m * deltaU + species.q * Phi
    Delta_PE_c = species.m * (deltaU_c - deltaU) + species.q * (Phi_c - Phi)
    
    beta = Delta_PE / (k0 * T_par0)
    
    if beta <= -1.:
        beta = -1.0 + np.finfo(float).eps
    
    
    beta_c = Delta_PE_c / (k0 * T_par0)
    beta_c = max(0.0, beta_c)
    
    # Calculate n2_psh
 
    # Prefactor for n2_psh
    sqrt_beta_c = np.sqrt(beta_c) 
    prefactor = (species.n * sqrt_beta_c * k0 * T_par0 * (beta + 1.0)**(-k0) * gamma(k0 + 1.0)) / \
                (np.sqrt(np.pi) * (1.0 - 1.0/Bratio) * T_perp0 * (k_tot0 + 1.0) * gamma(k0 + 0.5))
    
    # Define the integrand for the integral I
    def integrand(t):
        if t <= 0 or t >= 1:
            return 0.0
        
        term1 = t**(-0.5)
        term2 = (1.0 + (beta_c/(1.0 + beta)) * t)**(-k0)
        
        # Calculate the hypergeometric function argument
        z = 1.0 - (1.0 + beta) / (A0 * species.lam * (Bratio - 1.0) * (1.0 + (beta_c/(1.0 + beta)) * t))
        
        hf = hyp2f1(1.0, k_perp0 + 1.0, k_tot0 + 2.0, z)

        return term1 * term2 * hf
    
    # Perform the integration
    I, _ = quad(integrand, 0., 1., limit=500)
    
    n2_psh = prefactor * I

    
    # Calculate final PSH density
    n_local = n_full - n2_psh
    
    # Ensure non-negative density
    n_local = max(0.0, n_local)
    
    return n_local


def fried_egg_density_psh(species, deltaU, Phi, Bratio, deltaU_c, Phi_c):
    """
    Phase space hole version of "Fried-Egg" distribution.
    
    n_psh = (B/B0) * n0 * kappa * erfc(sqrt(Delta_PE_c/T_par0)) * U[1, 1-kappa, (B/B0-1)*A0*kappa] * exp(-Delta_PE/T_par0)
    
    This implementation ensures proper handling of potential calculations to avoid
    numerical issues near choke points.
    
    Parameters:
    -----------
    species : Species
        Has .kappa (perp kappa), .A = T_perp0/T_par0, .T = T_par, etc.
    deltaU : float
        Potential difference per unit mass from reference.
    Phi : float
        Electrostatic potential difference from reference.
    Bratio : float
        B(s)/B0
    deltaU_c : float
        Potential difference at choke point from reference.
    Phi_c : float
        Electrostatic potential at choke point from reference.
        
    Returns:
    --------
    n_local : float
        PSH number density at s.
    """

    # CalculateDelta_PE_c
    Delta_PE_c = species.m * (deltaU_c - deltaU) + species.q * (Phi_c - Phi)
    

    erfc_term = erfc(np.sqrt(max(0.0, Delta_PE_c / species.T)))
    

    
    # Combine all terms
    n_local =  erfc_term  * fried_egg_density(species, deltaU, Phi, Bratio)
    
    return n_local


###############################################################################
#                  TEMPERATURE FUNCTIONS FOR DIFFERENT DISTRIBUTIONS
###############################################################################

def maxwellian_temperature(species, deltaU, Phi, Bratio):
    """
    Temperature function for isotropic Maxwellian distribution.
    Temperature is constant along field line.
    
    Parameters:
    -----------
    species : Species
        Contains T (SI or consistent units).
    deltaU, Phi, Bratio : float
        Not used but included for interface consistency.
        
    Returns:
    --------
    T_par_local : float
        Parallel temperature at location s.
    T_perp_local : float
        Perpendicular temperature at location s.
    """
    T_par_local = species.T
    T_perp_local = species.T
    return T_par_local, T_perp_local


def aniso_maxwellian_temperature(species, deltaU, Phi, Bratio):
    """
    Temperature function for anisotropic Maxwellian distribution.
    T_parallel is constant, T_perpendicular varies with B.
    
    Parameters:
    -----------
    species : Species
        Contains T and A.
    deltaU, Phi, Bratio : float
        Bratio = B(s)/B0.
        
    Returns:
    --------
    T_par_local : float
        Parallel temperature at location s.
    T_perp_local : float
        Perpendicular temperature at location s.
    """
    T_par_local = species.T
    
    B0_over_B = 1.0 / Bratio
    denom = species.A + (1.0 - species.A) * B0_over_B
    
    if denom <= 0:
        print("Warning: Non-physical denom in aniso_maxwellian_temperature. Setting denom -> eps.")
        denom = np.finfo(float).eps
    
    T_perp_local = species.A * species.T / denom
    return T_par_local, T_perp_local

def aniso_maxwellian_const_Tperp_temperature(species, deltaU, Phi, Bratio):
    """
    Temperature function for anisotropic Maxwellian with constant T_perp.
    Both T_parallel and T_perpendicular are constant.
    
    Parameters:
    -----------
    species : Species
        Contains T and A.
    deltaU, Phi, Bratio : float
        Not used but included for interface consistency.
        
    Returns:
    --------
    T_par_local : float
        Parallel temperature at location s.
    T_perp_local : float
        Perpendicular temperature at location s.
    """
    T_par_local = species.T
    T_perp_local = species.A * species.T
    return T_par_local, T_perp_local


def aniso_kappa_temperature(species, deltaU, Phi, Bratio):
    """
    Temperature function for standard anisotropic kappa distribution.
    
    Parameters:
    -----------
    species : Species
        Contains T, A, kappa, m, q.
    deltaU : float
        Gravitational + centrifugal potential energy difference per unit mass.
    Phi : float
        Electrostatic potential difference (voltage).
    Bratio : float
        B(s) / B0.
        
    Returns:
    --------
    T_par_local : float
        Parallel temperature at location s.
    T_perp_local : float
        Perpendicular temperature at location s.
    """
    PEratio = (species.m * deltaU + species.q * Phi) / (species.kappa * species.T)
    
    # For physical consistency, clamp if PEratio <= -1
    if PEratio <= -1.:
        PEratio = -1.0 + np.finfo(float).eps
    
    beta = PEratio
    T_par_local = species.T * (1.0 + beta)
    
    B0_over_B = 1.0 / Bratio
    denom = species.A + (1.0 - species.A) * B0_over_B
    
    if denom <= 0:
        print("Warning: Non-physical denom in aniso_kappa_temperature. Setting denom -> eps.")
        denom = np.finfo(float).eps
    
    T_perp_local = species.A * species.T * (1.0 + beta) / denom
    return T_par_local, T_perp_local


def aniso_product_kappa_temperature(species, deltaU, Phi, Bratio):
    """
    Temperature function for anisotropic product kappa distribution.
    
    Parameters:
    -----------
    species : Species
        Contains T, A, kappa (interpreted as kappa_par), lam (kappa_perp/kappa_par).
    deltaU : float
        Gravitational + centrifugal potential energy difference per unit mass.
    Phi : float
        Electrostatic potential difference (voltage).
    Bratio : float
        B(s) / B0.
        
    Returns:
    --------
    T_par_local : float
        Parallel temperature at location s.
    T_perp_local : float
        Perpendicular temperature at location s.
    """
    # Check if B/B0 >= 1
    if Bratio < 1.0:
        print("Warning: B/B0 < 1 in aniso_product_kappa_temperature. Setting Bratio -> 1 + eps.")
        Bratio = 1.0 + np.finfo(float).eps
    
    beta = (species.m * deltaU + species.q * Phi) / (species.kappa * species.T)
    if beta <= -1.:
        beta = -1.0 + np.finfo(float).eps
    
    # Reference Values at s0
    kappa_par = species.kappa
    kappa_perp = kappa_par * species.lam
    A0 = species.A 
    T_par =  species.T
    T_perp =  A0*T_par
    x = A0 *(Bratio - 1.)
    g = (1.0 + beta)
    if Bratio != 1.0:
        # For larger kappa values, use simpler approximations
        if (kappa_par >= 17.) and (kappa_perp >= 17.):
            # Approximation for large kappa
            T_par_local = T_par * g 
            T_perp_local = Bratio* T_perp * g /(x + g )
        else:
            # Full calculation with hypergeometric functions
            delta_val = (species.lam * species.A) * ((Bratio - 1.) / (1. + beta))
            kappa_tot = kappa_par + kappa_perp
            
            # Calculate F_q values for Tperp
            # F_q = 2_F_1(q,kappa_perp + 0.5, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F1 = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F2 = hyp2f1(2.0, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F3q = hyp2f1(0.75, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F7q = hyp2f1(1.75, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            
            
            # Calculate regularized 2_F_1 values for Tpar
            # H_q = 2_F_1(1.0, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            H1h = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 0.5, 1.0 - 1./delta_val)/gamma(kappa_tot + 0.5)
            H3q = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 0.75, 1.0 - 1./delta_val)/gamma(kappa_tot + 0.75)
            H3h = F1/gamma(kappa_tot + 1.5)
            H7q = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 1.75, 1.0 - 1./delta_val)/gamma(kappa_tot + 1.75)
            
            # Calculate constants
    
            C1 = 4.0 * kappa_tot - 1.0
            C2 = 2.0 * kappa_tot - 1.0
            C3 = 4.0 * kappa_par - 1.0
            C4 = 2.0 * kappa_par - 1.0
    
            # Calculate temperatures
            T_perp_num = (beta + 1.0)*kappa_par*T_perp*(A0 + x)/(3. * A0 * x)
            
            Dperp1 =  C1*F3q / (3.*F7q)
            Dperp2 = C2*F1 / (2.*F2)
            
            T_perp_denom = Dperp1 - Dperp2
            
            # Value of Perp Temperature at s 
            T_perp_local = T_perp_num / T_perp_denom
    
    
            T_par_num = 8.0*(beta + 1.0)*kappa_par*T_par
            
            Dpar1 =  C3*C1*H7q / H3q
            Dpar2 = 2.*C4*C2*H3h / H1h
            T_par_denom = Dpar1 - Dpar2
            
            # Value of Parallel Temperature at s 
            T_par_local = T_par_num / T_par_denom
    else:
        # If B/B0 = 1 then we are at reference location s0
        T_perp_local = T_perp 
        T_par_local = T_par
    
    return T_par_local, T_perp_local



def fried_egg_temperature(species, deltaU, Phi, Bratio):
    """
    Temperature function for "Fried-Egg" distribution.
    
    Parameters:
    -----------
    species : Species
        Contains T, A, kappa (perp kappa).
    deltaU : float
        Gravitational + centrifugal potential energy difference per unit mass.
    Phi : float
        Electrostatic potential difference (voltage).
    Bratio : float
        B(s) / B0.
        
    Returns:
    --------
    T_par_local : float
        Parallel temperature at location s.
    T_perp_local : float
        Perpendicular temperature at location s.
    """
    # Check if B/B0 >= 1
    if Bratio <1.0:
        print("Warning: B/B0 < 1 in fried_egg_temperature. Setting Bratio -> 1 + eps.")
        Bratio = 1.0 + np.finfo(float).eps
    
    # T_parallel is constant for Fried-Egg
    T_par_local = species.T
    T_perp0 =  species.T *species.A
    
    k0 = species.kappa
    # For T_perpendicular, we use the derived formula with x and z
    x = species.A * (Bratio - 1.0)
    z = k0 * x
    
    
    
    if Bratio != 1.0:
        if k0 > 100.:
            # 0th order approximation is just
            # Tperp result is same as Aniso Maxwellian (Huang & Birmingham)
            # factor = 1./(1. + x) as Bratio/(1. + x) = 1/(A0+(1-A0)B0/B)
            # Good to 0.25 % for k0>100 for all x>0
            # Better for larger k0 and larger x 
            factor = 1./(1. + x) 
            T_perp_local = Bratio*T_perp0*factor   
        elif k0 >= 15.:
            # Good to 0.05% for all x>0 
            # Better for larger kappa and larger x 
            factor = 1./(1. + x) - x/(((1. + x)**3.)*k0)
            T_perp_local =  Bratio*T_perp0*factor
        elif (x>=7.) and (k0<15.):
            # For k0> 1, Bratio>1, and x>7. this 5th order approximation
            # Good to better than 0.22% for all all k0>1 
            # better for larger x or k0
            C1 = (3. + k0)*(4. + k0 * (7. + k0))
            C2 = -x*((1. + k0) ** 2.)*(4. + k0)
            C3 = -( x**3.) * (k0**2.) *(1. + k0) + ( x**2.) * k0 *(1. + k0) *(2. + k0)
            C4 = ( x**4.) * (k0**3.) + C3 + C2 + C1
            numerator =  -4. + k0 *C4
            denominator = ( x**5.) * (k0**4.)
            factor = numerator/denominator
            T_perp_local = Bratio*T_perp0*factor
        else:
            #for x<7 and k0 <15
            U5q = hyperu(k0 + 1., k0 + 5./4., z)
            U1q = hyperu(k0 + 1., k0 + 1./4., z)
            U1 = hyperu(k0 + 1., k0 + 1., z)
            U0 = hyperu(k0 + 1., k0, z)
            denominator = 4.* (U5q/U1q) - 3. *(U1/U0)
            factor = (1./(x*denominator))
            T_perp_local = Bratio * T_perp0 * factor
 
    else:
        # If B/B0 = 1 then we are at reference location s0
        T_perp_local = T_perp0 
    
    return T_par_local, T_perp_local




###############################################################################
#                  KAPPA FUNCTIONS FOR DIFFERENT DISTRIBUTIONS
###############################################################################

def maxwellian_kappa(species, deltaU, Phi, Bratio):
    """
    Kappa function for Maxwellian distribution.
    Not applicable but included for interface consistency.
    
    Returns:
    --------
    None, None : No kappa parameters for Maxwellian.
    Though actually effectivly Infinite as Maxwellian is infinite kappa limit
    """
    return None, None


def aniso_maxwellian_kappa(species, deltaU, Phi, Bratio):
    """
    Kappa function for anisotropic Maxwellian distribution.
    Not applicable but included for interface consistency.
    
    Returns:
    --------
    None, None : No kappa parameters for anisotropic Maxwellian.
    Though actually effectivly Infinite as Maxwellian is infinite kappa limit
    """
    return None, None


def aniso_maxwellian_const_Tperp_kappa(species, deltaU, Phi, Bratio):
    """
    Kappa function for anisotropic Maxwellian with constant T_perp.
    Not applicable but included for interface consistency.
    
    Returns:
    --------
    None, None : No kappa parameters for this distribution.
    Though actually effectivly Infinite as Maxwellian is infinite kappa limit
    """
    return None, None


def aniso_kappa_kappa(species, deltaU, Phi, Bratio):
    """
    Kappa function for standard anisotropic kappa distribution.
    Kappa parameter is constant along field line.
    
    Parameters:
    -----------
    species : Species
        Contains kappa.
    deltaU, Phi, Bratio : float
        Not used but included for interface consistency.
        
    Returns:
    --------
    kappa_par_local : float
        Parallel kappa parameter at location s.
    kappa_perp_local : float
        Perpendicular kappa parameter at location s.
    there is no par or perp in standard kappa only kappa 
    but include both for interface consistency.
    """
    kappa_par_local = species.kappa
    kappa_perp_local = species.kappa  # Same kappa for parallel and perpendicular
    # but include both for interface consistency.
    return kappa_par_local, kappa_perp_local


def aniso_product_kappa_kappa(species, deltaU, Phi, Bratio):
    """
    Kappa function for anisotropic product kappa distribution.
    
    Parameters:
    -----------
    species : Species
        Contains kappa (kappa_par), lam (kappa_perp/kappa_par).
    deltaU : float
        Gravitational + centrifugal potential energy difference per unit mass.
    Phi : float
        Electrostatic potential difference (voltage).
    Bratio : float
        B(s) / B0.
        
    Returns:
    --------
    kappa_par_local : float
        Parallel kappa parameter at location s.
    kappa_perp_local : float
        Perpendicular kappa parameter at location s.
    """
    # Check if B/B0 >= 1
    if Bratio < 1.0:
        print("Warning: B/B0 < 1 in aniso_product_kappa_kappa. Setting Bratio -> 1 + eps.")
        Bratio = 1.0 + np.finfo(float).eps
    
    beta = (species.m * deltaU + species.q * Phi) / (species.kappa * species.T)
    if beta <= -1.:
        beta = -1.0 + np.finfo(float).eps
    
    # kappa values at reference location s0
    # kappa_par0
    # kappa_perp0
    kappa_par = species.kappa
    kappa_perp = kappa_par * species.lam
    
    if Bratio != 1.0:
        # For larger kappa values, return constant kappa 
        # (approximation for large kappa matches standard kappa behavior)
        if (kappa_par >= 30.) and (kappa_perp >= 30.):
            kappa_par_local = kappa_par
            kappa_perp_local = kappa_perp
        else:
            # Full calculation with hypergeometric functions
            delta_val = (species.lam * species.A) * ((Bratio - 1.) / (1. + beta))
            kappa_tot = kappa_par + kappa_perp
            
            # Calculate F_q values for Tperp
            # F_q = 2_F_1(q,kappa_perp + 0.5, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F1 = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F2 = hyp2f1(2.0, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F3q = hyp2f1(0.75, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F7q = hyp2f1(1.75, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            
            
            # Calculate regularized 2_F_1 values for Tpar
            # H_q = 2_F_1(1.0, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            H1h = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 0.5, 1.0 - 1./delta_val)/gamma(kappa_tot + 0.5)
            H3q = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 0.75, 1.0 - 1./delta_val)/gamma(kappa_tot + 0.75)
            H3h = F1/gamma(kappa_tot + 1.5)
            H7q = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 1.75, 1.0 - 1./delta_val)/gamma(kappa_tot + 1.75)
            
            # Calculate constants

            C1 = 4.0 * kappa_tot - 1.0
            C2 = 2.0 * kappa_tot - 1.0
            C3 = 4.0 * kappa_par - 1.0
            C4 = 2.0 * kappa_par - 1.0
             
            num_kappa_perp = 2.0 * C1* F3q * F2 / (C2* F1 * F7q) - 4.0 
            # Calculate kappa_perp at s
            kappa_perp_local = 1.0/num_kappa_perp + 1.0
            
            
            
            # Calculate kappa_par at s
            fac1 = C3*C1/(C4*C2)
            fac2 = fac1*H1h*H7q/(H3q*H3h)
            kappa_par_local = (fac2 -2.0) / (2*fac2 - 8.0 )
    else:
        # If Bratio = 1 then at reference location s0
        kappa_par_local, kappa_perp_local = kappa_par, kappa_perp
    
    return kappa_par_local, kappa_perp_local


def fried_egg_kappa(species, deltaU, Phi, Bratio):
    """
    Kappa function for "Fried-Egg" distribution.
    
    Parameters:
    -----------
    species : Species
        Contains kappa (perp kappa).
    deltaU : float
        Gravitational + centrifugal potential energy difference per unit mass.
    Phi : float
        Electrostatic potential difference (voltage).
    Bratio : float
        B(s) / B0.
        
    Returns:
    --------
    None : Since there's no parallel kappa 
    (Maxwellian in parallel direction so parallel kappa effectivly infinite).
    kappa_perp_local : float
        Perpendicular kappa parameter at location s.
    """
    # Check if B/B0 >= 1
    if Bratio <1.0:
        print("Warning: B/B0 < 1 in fried_egg_kappa. Setting Bratio -> 1 + eps.")
        Bratio = 1.0 + np.finfo(float).eps
    
    # No parallel kappa parameter (Maxwellian in parallel direction)
    kappa_par_local = None
    
    #k0 = kappa_perp0
    k0 = species.kappa 
    # For kappa_perpendicular, use the derived formula with x and z
    x = species.A * (Bratio - 1.0)
    z = k0 * x
 
    if Bratio != 1.0:
        if  (k0 >= 50.0):
            kappa_perp_local = k0
        elif (x >=15.0):
            #x>=15.
            # Series expansion about x = infinity to 2nd order
            # Good to 0.4% for k0 >1 and x>15
            term1 = ((x**2.)* (k0**2.))/(1. + k0) 
            term2 =  (289.* (2. + k0))/(16. * x*k0 *(1. + k0))
            term3 = (x*k0 *(27. + 8.*k0))/(4.* (1. + k0))
            kappa_perp_local = term1 + term2 + term3
        else: 
            # Via scipy special direct calculation using confluent hypergeometric function U(a,b,z)
            U0 = hyperu(k0 + 1., k0, z)
            U1 = hyperu(k0 + 1., k0 + 1., z)
            U1q = hyperu(k0 + 1., k0 + 1./4., z)
            U5q = hyperu(k0 + 1., k0 + 5./4., z)
            
            fac = U0*U5q/(U1q*U1)
            kappa_perp_local = (0.75 - fac)/(1.0 - fac)
    else:
        # If Bratio = 1 we are at reference location
        kappa_perp_local = species.kappa
    
    return kappa_par_local, kappa_perp_local


###############################################################################
#         DICTIONARIES OF DENSITY, TEMPERATURE, AND KAPPA FUNCTIONS
###############################################################################

density_functions = {
    'Maxwellian': maxwellian_density,
    'Aniso_Maxwellian': aniso_maxwellian_density,
    'Aniso_Maxwellian_const_Tperp': aniso_maxwellian_const_Tperp_density,
    'Aniso_kappa': aniso_kappa_density,
    'Aniso_product_kappa': aniso_product_kappa_density,
    'Fried_Egg': fried_egg_density
}

density_functions_psh = {
    'Maxwellian': maxwellian_density_psh,
    'Aniso_Maxwellian': aniso_maxwellian_density_psh,
    'Aniso_Maxwellian_const_Tperp': aniso_maxwellian_const_Tperp_density_psh,
    'Aniso_kappa': aniso_kappa_density_psh,
    'Aniso_product_kappa': aniso_product_kappa_density_psh,
    'Fried_Egg': fried_egg_density_psh
}

temperature_functions = {
    'Maxwellian': maxwellian_temperature, # Constant but included for consistent structures
    'Aniso_Maxwellian': aniso_maxwellian_temperature,
    'Aniso_Maxwellian_const_Tperp': aniso_maxwellian_const_Tperp_temperature, # Constant but included for consistent structures
    'Aniso_kappa': aniso_kappa_temperature,
    'Aniso_product_kappa': aniso_product_kappa_temperature,
    'Fried_Egg': fried_egg_temperature
}

kappa_functions = {
    'Maxwellian': maxwellian_kappa, #NA But included for consistency accross structures
    'Aniso_Maxwellian': aniso_maxwellian_kappa, #NA But included for consistency accross structures
    'Aniso_Maxwellian_const_Tperp': aniso_maxwellian_const_Tperp_kappa, #NA But included for consistency accross structures
    'Aniso_kappa': aniso_kappa_kappa, #Constant But included for consistency accross structures
    'Aniso_product_kappa': aniso_product_kappa_kappa, 
    'Fried_Egg': fried_egg_kappa
}


def calc_deltaU(r0_in, lat0_in, r_in, lat_in, planet, fcor):
    # Set planet-specific variables
    RP, Omega, GM = define_planet(planet)
    Omega *= fcor
    # Calculate gravitational-centrifugal potentials per mass at s0 and s along field line
    r0 = r0_in * RP    # Convert radius to meters
    lat0 = np.deg2rad(lat0_in) # convert to radians
    r = r_in * RP    # Convert radius to meters
    lat = np.deg2rad(lat_in) # convert to radians
    
    U0 = -GM / r0 - 0.5 * r0 ** 2. * np.cos(lat0) ** 2. * Omega ** 2.
    U = -GM / r - 0.5 *r ** 2. * np.cos(lat) ** 2. * Omega ** 2.
    
    # difference in gravitational + centrifdugal potential energy per mass
    # along field line as a function of s
    dU = U - U0
    return dU

def find_local_maxima(deltaU, lat):
    """
    Find local maxima in deltaU as a function of lat.
    
    Parameters:
    -----------
    deltaU : 1D numpy array
        The values to find maxima in
    lat : 1D numpy array
        The corresponding latitude values in degrees
    
    Returns:
    --------
    max_indices : 1D numpy array
        Indices of the local maxima in the original arrays
    max_lats : 1D numpy array
        Latitude values at local maxima
    max_deltaU : 1D numpy array
        Values of deltaU at local maxima
    """
    # Ensure inputs are numpy arrays
    deltaU = np.array(deltaU)
    lat = np.array(lat)
    
    # Check if arrays have the same length
    if len(deltaU) != len(lat):
        raise ValueError("deltaU and lat arrays must have the same length")
    
    # Sort by latitude (x-axis)
    sort_indices = np.argsort(lat)
    #lat_sorted = lat[sort_indices]
    deltaU_sorted = deltaU[sort_indices]
    
    # Find local maxima
    max_sorted_indices = []
    
    # We need at least 3 points to find local maxima
    if len(deltaU_sorted) < 3:
        return np.array([]), np.array([]), np.array([])
    
    # Check each point (except first and last)
    for i in range(1, len(deltaU_sorted) - 1):
        if deltaU_sorted[i] > deltaU_sorted[i-1] and deltaU_sorted[i] > deltaU_sorted[i+1]:
            max_sorted_indices.append(i)
    
    # No maxima found
    if not max_sorted_indices:
        return np.array([]), np.array([]), np.array([])
    
    # Convert sorted indices to original indices
    max_indices = sort_indices[max_sorted_indices]
    
    # Get values at maxima
    max_lats = lat[max_indices]
    max_deltaU = deltaU[max_indices]
    
    return max_indices, max_lats, max_deltaU

    

# Function to compute densities and find phi from neutrality condition
def diff_eq(Pos_in, Pos0_in, Bratio, species0, deltaU, planet, fcor):

    nspec = len(species0)

    # Ensure the plasma is charge neutral at Pos0
    # Adjust electron density
    # Assuming species0[nspec-1] is electrons
    total_ion_charge_density = np.sum([species0[i].q * species0[i].n for i in range(nspec - 1)])
    # Adjust cold electron  to make sure initially neutral at s0 reference location
    species0[nspec - 1].n = -total_ion_charge_density / species0[nspec - 1].q  


    def net_charge_density(Phi):
        n_local = np.zeros(nspec)
        for i in range(nspec):
            density_func = density_functions.get(species0[i].type, maxwellian_density)
            n_local[i] = density_func(species0[i], deltaU, Phi, Bratio)
        #charge neutrality calculation
        nq_temp = [species0[i].q * n_local[i] for i in range(nspec)]
        nq = np.sum(nq_temp)
        return nq

    # Initialize variables for root-finding
    Phi = 0.0
    Phi2 = 0.0
    dPhi = 1.0

    # Compute initial net charge densities
    nq = net_charge_density(Phi)
    nq2 = net_charge_density(Phi2)

    # Adjust Phi and Phi2 to bracket the root
    max_iterations = 10000
    iteration = 0
    while nq * nq2 > 0 and iteration < max_iterations:
        Phi -= dPhi
        Phi2 += dPhi
        nq = net_charge_density(Phi)
        nq2 = net_charge_density(Phi2)
        iteration += 1
    if iteration >= max_iterations:
        print("Failed to bracket the root.")
        return None

    # Use bisection method to find Phi where net charge density is zero
    # default tolerances, typically gives \sum z_i n_i/(e*n_e) ~ 1e-8 relative neutrality to electrons
    # Could be sped up a bit by relaxing this
    Phi_root = bisect(net_charge_density, Phi, Phi2, xtol=2e-12, rtol=8.881784197001252e-16) 
    # Compute densities at Phi_root
    n = np.zeros(nspec)
    for i in range(nspec):
        density_func = density_functions.get(species0[i].type, maxwellian_density)
        n[i] = density_func(species0[i], deltaU, Phi_root, Bratio)
    phi = Phi_root

    # Return densities and phi
    return n, phi, deltaU


def diff_eq_psh(Pos_in, Pos0_in, lat_in, lat0_in, Bratio, species0, deltaU, deltaU_array, lat_array, B_array, Bratio_array, planet, fcor):
    """
    Extended version of diff_eq that handles phase space holes.
    
    Determines whether to use regular or PSH density functions based on whether
    reference location s0 and target location s are on opposite sides of a local
    maximum in deltaU.
    
    This function properly handles phase space holes following the logic from
    the IDL version, ensuring correct potential calculations at choke points.
    
    Parameters:
    -----------
    Pos_in : array
        Position vector at target location s
    Pos0_in : array
        Position vector at reference location s0
    lat_in : float
        Latitude at target location s (degrees)
    lat0_in : float
        Latitude at reference location s0 (degrees)
    Bratio : float
        B(s)/B0
    species0 : list of Species
        Species parameters at reference location
    deltaU : float
        Potential energy difference at target location (U - U0)
    deltaU_array : array
        Array of deltaU values along field line
    lat_array : array
        Array of latitude values along field line
    B_array : array
        Array of B values along field line
    Bratio_array : array
        Array of B/B0 values along field line
    planet : str
        Planet name
    fcor : float
        Fraction of corotation
        
    Returns:
    --------
    n : array
        Number densities at target location
    phi : float
        Electrostatic potential at target location
    deltaU : float
        Potential energy difference
    """
    nspec = len(species0)

    # Ensure the plasma is charge neutral at Pos0
    # Adjust electron density
    # Assuming species0[nspec-1] is electrons
    total_ion_charge_density = np.sum([species0[i].q * species0[i].n for i in range(nspec - 1)])
    # Adjust cold electron to make sure initially neutral at s0 reference location
    species0[nspec - 1].n = -total_ion_charge_density / species0[nspec - 1].q  

    # Find local maxima in deltaU
    max_indices, max_lats, max_deltaU = find_local_maxima(deltaU_array, lat_array)
    
    # Determine if we need to use PSH functions
    # Following IDL logic: torus_side = 1 means we're between choke points or PSH not used
    torus_side = 1  # Default: between choke points or no PSH
    Phi_c = None
    deltaU_c = None
    Bratio_c = None
    
    if len(max_lats) > 0:
        # Check if s0 and s are on opposite sides of any maximum
        # We assume s0 is on the equatorward side of choke points
        for i, max_lat in enumerate(max_lats):
            # Check if the choke point is between s0 and s
            if min(lat0_in, lat_in) < max_lat < max(lat0_in, lat_in):
                # We're crossing a choke point
                torus_side = 0  # Planetward of choke point
                # Get values at choke point
                choke_idx = max_indices[i]
                deltaU_c = max_deltaU[i]
                Bratio_c = Bratio_array[choke_idx]
                
                # Calculate potential at choke point using regular density functions
                # This matches the IDL code's approach to finding Phi_c
                def net_charge_density_choke(Phi):
                    n_local = np.zeros(nspec)
                    for j in range(nspec):
                        density_func = density_functions.get(species0[j].type, maxwellian_density)
                        n_local[j] = density_func(species0[j], deltaU_c, Phi, Bratio_c)
                    nq_temp = [species0[j].q * n_local[j] for j in range(nspec)]
                    nq = np.sum(nq_temp)
                    return nq
                
                # Find Phi_c by bracketing and bisection
                Phi = 0.0
                Phi2 = 0.0
                dPhi = 1.0
                nq = net_charge_density_choke(Phi)
                nq2 = net_charge_density_choke(Phi2)
                
                # Bracket the root
                max_iterations = 10000
                iteration = 0
                while nq * nq2 > 0 and iteration < max_iterations:
                    Phi -= dPhi
                    Phi2 += dPhi
                    nq = net_charge_density_choke(Phi)
                    nq2 = net_charge_density_choke(Phi2)
                    iteration += 1
                
                if iteration < max_iterations:
                    # Use bisection to find Phi_c
                    Phi_c = bisect(net_charge_density_choke, Phi, Phi2, xtol=2e-12, rtol=8.881784197001252e-16)
                else:
                    print("Warning: Failed to find Phi_c at choke point - using 0")
                    Phi_c = 0.0
                
                break

    # Now solve for densities and potential at the target location
    if torus_side == 1:
        # Torus side: use regular density functions
        def net_charge_density(Phi):
            n_local = np.zeros(nspec)
            for i in range(nspec):
                density_func = density_functions.get(species0[i].type, maxwellian_density)
                n_local[i] = density_func(species0[i], deltaU, Phi, Bratio)
            nq_temp = [species0[i].q * n_local[i] for i in range(nspec)]
            nq = np.sum(nq_temp)
            return nq
        
        # Initialize variables for root-finding
        Phi = 0.0
        Phi2 = 0.0
        dPhi = 1.0

        # Compute initial net charge densities
        nq = net_charge_density(Phi)
        nq2 = net_charge_density(Phi2)

        # Adjust Phi and Phi2 to bracket the root
        max_iterations = 10000
        iteration = 0
        while nq * nq2 > 0 and iteration < max_iterations:
            Phi -= dPhi
            Phi2 += dPhi
            nq = net_charge_density(Phi)
            nq2 = net_charge_density(Phi2)
            iteration += 1
        
        if iteration >= max_iterations:
            print("Failed to bracket the root.")
            return None

        # Use bisection method to find Phi where net charge density is zero
        Phi_root = bisect(net_charge_density, Phi, Phi2, xtol=2e-12, rtol=8.881784197001252e-16)
        
        # Compute densities at Phi_root
        n = np.zeros(nspec)
        for i in range(nspec):
            density_func = density_functions.get(species0[i].type, maxwellian_density)
            n[i] = density_func(species0[i], deltaU, Phi_root, Bratio)
            
    else:
        # Planetward side: use PSH density functions
        # Following IDL logic for setting Phi bounds when using PSH
        
        # Set Phi and Phi2 to theoretical minimum and maximum values
        # Phi must ensure all species have positive argument in erfc
        # For electrons (last species), we need the most restrictive condition
        Phi_min = (species0[nspec-1].m/species0[nspec-1].q) * (deltaU_c - deltaU) + Phi_c
        
        # For positive ions, find the minimum allowed Phi
        Phi_max_list = []
        for i in range(nspec):
            if species0[i].q > 0:
                Phi_max_val = (species0[i].m/species0[i].q) * (deltaU_c - deltaU) + Phi_c
                Phi_max_list.append(Phi_max_val)
        Phi_max = min(Phi_max_list)
        
        # Define PSH net charge density function
        def net_charge_density_psh(Phi):
            n_local = np.zeros(nspec)
            for i in range(nspec):
                density_func = density_functions_psh.get(species0[i].type, maxwellian_density_psh)
                n_local[i] = density_func(species0[i], deltaU, Phi, Bratio, deltaU_c, Phi_c)
            nq_temp = [species0[i].q * n_local[i] for i in range(nspec)]
            nq = np.sum(nq_temp)
            return nq
        
        # The root should be between Phi_min and Phi_max
        Phi = Phi_min
        Phi2 = Phi_max
        
        # Evaluate charge density at bounds
        nq = net_charge_density_psh(Phi)
        nq2 = net_charge_density_psh(Phi2)
        
        # Bisect until charge density is under 0.1% of electron density
        max_iterations = 10000
        iteration = 0
        
        while abs(nq/species0[nspec-1].n/species0[nspec-1].q) > 1e-3 and iteration < max_iterations:
            # Evaluate midpoint charge density
            Phim = 0.5 * (Phi + Phi2)
            nqm = net_charge_density_psh(Phim)
            
            # Move Phi or Phi2 to Phim, so that the root is between the new Phi and Phi2
            if nq * nqm > 0:
                Phi = Phim
                nq = nqm
            else:
                Phi2 = Phim
                nq2 = nqm
            
            iteration += 1
        
        if iteration >= max_iterations:
            print("Warning: PSH bisection did not converge to desired tolerance")
        
        Phi_root = Phi
        
        # Compute densities at Phi_root using PSH functions
        n = np.zeros(nspec)
        for i in range(nspec):
            density_func = density_functions_psh.get(species0[i].type, maxwellian_density_psh)
            n[i] = density_func(species0[i], deltaU, Phi_root, Bratio, deltaU_c, Phi_c)
    
    phi = Phi_root

    # Return densities and phi
    return n, phi, deltaU


def calc_temps(species0, deltaU, phi, Bratio):
    nspec = len(species0)
    Tpar = np.zeros(nspec)
    Tperp = np.zeros(nspec)
    for i in range(nspec):
        #print(species0[i].type)
        temp_func = temperature_functions.get(species0[i].type, maxwellian_temperature)
        Tpar[i], Tperp[i] = temp_func(species0[i], deltaU, phi, Bratio)
    return Tpar, Tperp

def calc_kappa_vals(species0, deltaU, phi, Bratio):
    nspec = len(species0)
    kappa_par = np.zeros(nspec)
    kappa_perp = np.zeros(nspec)
    for i in range(nspec):
        kappa_func = kappa_functions.get(species0[i].type, maxwellian_kappa)
        kappa_par[i], kappa_perp[i] = kappa_func(species0[i], deltaU, phi, Bratio)
    return kappa_par, kappa_perp