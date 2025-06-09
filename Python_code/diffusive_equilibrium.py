#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Feb 20 17:10:10 2025

@author: Edward (Eddie) G. Nerney

Nerney (2025)

Diffusive Equilibrium Code for finding steady state plasma Densities
at extrapolated location s given plasma parameters at s along same 
magnetic field line trace. s is arc length along field line.
In assumed non inertial rotating frame with planet at \Omega in rad/s,
or with fraction of corotation fcor given
along fixed magnetic field lines given field line traces,
including gravitational potential energy (only important close to planet), 
Centrifugal potential energy, and electrostatic potenial energy.

Assumes dynamic accesibility to 
all locations along field line in steady state,
plasma frozen to field line,
So No cross‐field drift (beyond small gyration around B),
No significant scattering to transport particles across field lines,
equivalently for MHD Alfvén's theorem frozen-in flux theorem or large magnetic reynolds number,
no collisions or sources or losses of particles (no loss cone), 
no drift or transport (though boundary condition can be informed by),
the dist. function same at reference location s0 as s,
energy conserved (no sources or losses of energy or dissipation),
magnetic moment conserved (adiabatic) or timescales all >>> gyroperiod,
Assumes neutrality spatial scales >>> debye length,
and assumes plasma frame is the same as the rotating frame.

T parallel or Tpar is given T parameter 
with T perpendicular or Tperp given by A*T.
or A = Tperp/Tpar (Anistropy Factor)
For Kappa distributions Temps are characterstic temps or T_c not "real" temps. T
Where real temps are defined as as Tperp = pperp/n and Tpar = Tpar/n
For standard kappas the relation is
Tpar = Tpar_c/(1-3/(2*kappa)), Tperp = Tperp_c/(1-3/(2*kappa))  
For product kappas there can be different kappa_par and kappa_perp
Tpar = Tpar_c/(1-1/(2*kappa_par)), Tperp = Tperp_c/(1 - 1/(kappa_perp))  
kappa parameter given is assumed to be parralel kappa or kappa_par
kappa_par = kappa
then perp. kappa value is kappa_perp = lambda * kappa
lambda = kappa_perp/kappa_par
"""


__all__ = ['Species', 'define_planet', 'maxwellian_density', 'aniso_maxwellian_density', 'aniso_maxwellian_const_Tperp_density', 
           'aniso_kappa_density', 'aniso_product_kappa_density', 'fried_egg_density', 'maxwellian_temperature', 'aniso_maxwellian_temperature',
           'aniso_maxwellian_const_Tperp_temperature', 'aniso_kappa_temperature', 'aniso_product_kappa_temperature', 'fried_egg_temperature',
           'maxwellian_kappa', 'aniso_maxwellian_kappa','aniso_maxwellian_const_Tperp_kappa', 'aniso_kappa_kappa', 'aniso_product_kappa_kappa', 'fried_egg_kappa',
           'density_functions', 'temperature_functions', 'kappa_functions', 'diff_eq','find_local_maxima','calc_deltaU', 'calc_temps','calc_kappa_vals' ]



import numpy as np
import sys

# Numerical libraries
from scipy.optimize import bisect
from scipy.special import hyp2f1, gamma, hyperu

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
    



