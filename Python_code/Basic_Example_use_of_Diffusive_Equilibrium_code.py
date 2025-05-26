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

Test Example use of Diffusive Equilibrium MOP community code
"""


from diffusive_equilibrium import *

"""
__all__ = ['Species', 'define_planet', 'maxwellian_density', 'aniso_maxwellian_density', 'aniso_maxwellian_const_Tperp_density', 
           'aniso_kappa_density', 'aniso_product_kappa_density', 'fried_egg_density', 'maxwellian_temperature', 'aniso_maxwellian_temperature',
           'aniso_maxwellian_const_Tperp_temperature', 'aniso_kappa_temperature', 'aniso_product_kappa_temperature', 'fried_egg_temperature',
           'maxwellian_kappa', 'aniso_maxwellian_kappa','aniso_maxwellian_const_Tperp_kappa', 'aniso_kappa_kappa', 'aniso_product_kappa_kappa', 'fried_egg_kappa',
           'density_functions', 'temperature_functions', 'kappa_functions', 'diff_eq','find_local_maxima','calc_deltaU', 'calc_temps','calc_kappa_vals' ]
"""

import numpy as np

import matplotlib.pyplot as plt


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
m_Hp = m_H - m_e              # Sinlgy (Fully) Ionized Hydrogen or H^{+} or proton mass in u
m_Op = m_O - m_e              # Sinlgy Ionized Oxygen or O^{+} mass in u
m_O2p = m_O - 2.*m_e          # Doubly Ionized Oxygen or O^{++} mass in u
m_Sp = m_S - m_e              # Sinlgy Ionized Sulfur or S^{+} mass in u
m_S2p = m_S - 2.*m_e          # Doubly Ionized Oxygen or S^{++} mass in u
m_S3p = m_S - 3.*m_e          # Triply Ionized Oxygen or S^{+++} mass in u
m_Nap = m_Na - m_e            # Sinlgy Ionized Sodium or Na^{+} mass in u

ELEMENTARY_CHARGE = 1.602176634e-19  # Elementary charge in Coulombs
EV_TO_JOULE = ELEMENTARY_CHARGE   # Conversion from eV to Joules
JOULE_TO_EV = 1./ELEMENTARY_CHARGE 

if __name__ == '__main__':
    # Define species
    planet = 'Jupiter'
    fcor = 1.0  # azimuthal velocity relative to corotation
    # perfectly corotating with planet if fcor =1 
    # define planetary equatorial radius at 1 bar level R_P (m)
    # Sidereal rotation rate of planet Omega_P (rad/s)
    # G*M_P  Gravitational constant times Mass of planet (m^3/s^2)
    RP, Omega, GM = define_planet(planet) 
    Omega *= fcor
    
    
    nrfls = 601
    
    iofl_idx = 191
    
    # Load input data files s in RJ (s=0 corresponds to southern ionosphere at min latitude), x(s) in RJ, y(s) in RJ, z(s) in RJ, 
    # rho(s) in RJ, r(s) in RJ, lat(s) in Degrees at 0.1 degree resoluton, east longitude(s) in Degrees, w
    # west Longitude(s) in Degrees, and B(s) in nT for Io Aligned Field line
    # Field line trace was interpolated to this resolution using scipy.interpolate.PchipInterpolator
    # going through aligned equator for JRM 33 + Con2020 field line traces and magnetic field magnitude B 
    # Using MOP Community Code Wilson et al. (2023)
    # 13 Spherical Harmonics for internal JRM33 Field model
    # Integral version for Con2020 extrenal jovian magnetodisc/current sheet field model
    # Load the first file
    with np.load('Io_aligned_field_line_trace_JRM33+con2020_dlat=0.1_degree.npz') as data:
        s = data['s']
        x = data['x']
        y = data['y']
        z = data['z']
        rho = data['rho']
        r = data['r']
        lat = data['lat']
        elong = data['elong']
        wlong = data['wlong']
        B = data['B']
    npointsfl = len(s) # 1365
    
    
    # centrifugal equator cylindrical radius reference values for 1D radial model
    rho_c0 = 0.01 * np.arange(601) + 4.0
    # For this aligned case all equators are same so these work for reference values at latitiude = 0 degrees
    # but for other longitudes to use this model one should first find find centrifugal equator latitude
    # Centrifugal equator latitude is where rho_III is max along a given field line
    # then at that latitude along a given field line rho_ceq = rho_III
    # Then model radials reference and densities as a function of centrifugal rho or rho_ceq
    # can be interpolated to the rho_III value at that location
    
    # Load 1D radial centrifugal equator model reference values of densities, temperatures, kappa values, 
    # and temps for kappas matching sum of hot and cold maxwellians see Nerney (2025)
    with np.load('Nerney2025_1D_radial_ceq_reference_model_4-10RJ_601pts_dr=0.01RJ.npz', allow_pickle=True) as data:
        # Each loaded item is a 0-d array containing your dictionary.
        n_0_loaded = data['n_0'].item()
        T_0_loaded = data['T_0'].item()
        kappa_0_loaded = data['kappa_0'].item()
        kappa_temps_0_loaded = data['kappa_temps_0'].item()



    # Species names
    # For Singly and doubly ionized Oxygen, Single, doubly, and triply ionized Sulfur
    # protons or H^+, Na^+, the proxy hot O^+ (representing both O^{+} S^{++} to match Dougherty (2017) analysis)
    # then a hot electron and core electron population
    # Here hots are treated as their own individual species of distribution type given below (typically Maxwellians)
    # "Proper" Treatment would be to have hots included with cold core as a single kappa, See Nerney (2025)
    species_names = ['O+', 'O++', 'S+', 'S++', 'S+++', 'H+', 'Na+', 'O+(hot)', 'eh-', 'e-']
    

    nspec = len(species_names) # Number of species = 10
    
    #Mass of each species in u or AMU
    species_m = [m_Op, m_O2p, m_Sp, m_S2p, m_S3p, m_Hp, m_Nap, m_Op, m_e , m_e]
    
    #Charge of each species in e or elementary charge units
    species_q = [1.0, 2.0, 1.0, 2.0, 3.0, 1.0, 1.0, 1.0, -1.0, -1.0]
    
    
    
    
    # Input plasma paramters at s_0 to be extrapolated to arbitrary s value along field line
    # Create empty lists
    species_n = []
    species_T = []
    species_kappa_vals = []
    species_kappa_temps = []
    for key in species_names:
        # n_0 number densities in cm^{-3} at s_0 but can be any desired unit 
        # as their value is just a multiplicitive constant final output will be in whatver unit these are
        species_n.append(n_0_loaded[key][iofl_idx].item()) 
        # T_par0 or Parallel Temperatures at s_0 in eV. As we work with kB = 1
        species_T.append(T_0_loaded[key][iofl_idx].item())
        
        # kappa fits to istropic double Maxwellians for kappa values for standard kappa
        # See Nerney (2025)
        species_kappa_vals.append(kappa_0_loaded[key][iofl_idx].item())
        # kappa fits to istropic double Maxwellians for T_par0 or Parallel Temperatures at s_0 in eV. As we work with kB = 1
        # See Nerney (2025), These are characterstic temperatures not the ones exactly propto <0.5 mv^2> 
        species_kappa_temps.append(kappa_temps_0_loaded[key][iofl_idx].item())
        

    
    
    # Anisotropy A_0 = T_perp0/T_par0 => T_perp0 = A_0 * T_par0
    # Set all to 2 here so Anisotropic 
    species_A = [2.0] * nspec 

    # lambda_0 = kappa_perp0/kappa_par0 for product kappa only but including for consistency in structures
    species_lam = [1.] * nspec 
    
    # Set all 10 species to distribution type to Anisotropic Maxwellians, They need not all be the same for each species
    #species_type = ['Maxwellian'] * nspec  # Bagenal & Sullivan (1981) etc. Basic default Assumed isotropic always, ignores anisotropy inputs
    #species_type = ['Aniso_Maxwellian'] * nspec  # Huang & Birmingham (1992) Anisotropic Maxwellian
    #species_type = ['Aniso_Maxwellian_const_Tperp'] * nspec  # Bagenal & Sullivan (1981), Bagenal (1994), Mei, Thorne, Bagenal (1995) Assumed Constant Anisotrpy along field line, inconsistent but good enough often 
    #species_type = ['Aniso_kappa'] * nspec  # Meyer-Vernet et al. (1995); Moncuquet et al. (2002) Anisotropic Standard Kappa distribution
    #species_type = ['Aniso_product_kappa'] * nspec  # Nerney (2025) Anisotropic Product Kappa Distribution 
    species_type = ['Fried_Egg'] * nspec  # Nerney (2025) "Fried-Egg" Distribution Function (Maxwellian limit in parallel of aniso product kappa and still kappa in perp direction) 
    # Create empty species list to be filled
    species_0 = []
    # Fill Species class with plasma parameters at s_0 defined above for use with diff_eq code
    # Also set to m, q, and T to SI units for calculations. T in Joules as kB = 1
    for i in range(nspec):
        sp = Species(
            name=species_names[i],
            m=species_m[i]* kg_per_u,             # Convert mass from AMU (u or daltons) to kg
            q=species_q[i]* ELEMENTARY_CHARGE,    # Convert charge from elementary charges to Coulombs
            T=species_T[i]* ELEMENTARY_CHARGE,    # Convert temperature from eV to Joules
            A=species_A[i],
            kappa=species_kappa_vals[i],
            lam=species_lam[i],
            n=species_n[i],
            dist_type=species_type[i]
        )
        species_0.append(sp)
        
        
        
        

    # Define the output containers as dictionaries of numpy arrays for plasma parameters along field line
    n = {key: np.zeros(npointsfl) for key in species_names} #n(s)
    T_par= {key: np.zeros(npointsfl) for key in species_names} #T_par(s)
    T_perp= {key: np.zeros(npointsfl) for key in species_names} #T_perp(s)
    kappa_par= {key: np.zeros(npointsfl) for key in species_names} #kappa_par(s)
    kappa_perp = {key: np.zeros(npointsfl) for key in species_names} #kappa_perp(s)
    
    
    
    # Centrigual equator defined as where rho = max(rho) or farthest point from spin axis along field line. 
    ceq_idx = np.where(rho == max(rho))
    # Magnetic field ratio
    Bratio = B / B[ceq_idx]
    
    Pos = np.array([r, lat, wlong]).T
    Pos0 = [r[ceq_idx], lat[ceq_idx], wlong[ceq_idx]]
    
    
    
    deltaU = calc_deltaU(r[ceq_idx], lat[ceq_idx], r, lat, planet, fcor)
    
    max_indices, max_lats, max_deltaU = find_local_maxima(deltaU, lat)
    plt.figure(figsize=(10, 6))
    plt.plot(lat, deltaU)
    # Plot maxima points
    #plt.plot(max_lats, max_deltaU, label='Local maxima')
    # Add vertical dotted lines at maxima
    for i in range(len(max_lats)):
        plt.axvline(x=max_lats[i], color='C1', linestyle=':', alpha=0.7)
    plt.xlabel(r'Latitude ($\degree$)')
    plt.ylabel(r'$\Delta$U (Joules/kg)')
    plt.title('Centrifugal + Gravitational Potential Energy/mass vs Latitude with Local Maxima')
    #plt.legend()
    plt.grid(True)
    plt.show()
    # Prepare arrays for computation
    n_out = np.zeros((nspec, npointsfl))
    Tpar_out = np.zeros((nspec, npointsfl))
    Tperp_out = np.zeros((nspec, npointsfl))
    kappapar_out = np.zeros((nspec, npointsfl))
    kappaperp_out = np.zeros((nspec, npointsfl))
    phi = np.zeros(npointsfl)
    #phi = None
    posc = None
    nchoke = None
    phic = None

    # Loop over each point to compute densities along field line
    for i in range(npointsfl):
        if i % 100 == 0:
            print(i)
        result = diff_eq(
            Pos[i, :],
            Pos0,
            Bratio[i],
            species_0,
            deltaU[i],
            planet,
            fcor)
        
        
        
        if result is not None:
            n_out[:, i] = result[0]
            phi[i] = result[1]
            Tpar_out[:, i], Tperp_out[:, i]= calc_temps(species_0, deltaU[i], phi[i], Bratio[i])
            kappapar_out[:, i], kappaperp_out[:, i] = calc_kappa_vals(species_0, deltaU[i], phi[i], Bratio[i])
        else:
            n_out[:, i] = np.zeros(nspec)
            Tpar_out[:, i] = np.zeros(nspec)
            Tperp_out[:, i] = np.zeros(nspec)
            kappapar_out[:, i] = np.zeros(nspec)
            kappaperp_out[:, i] = np.zeros(nspec)
            phi[i] = np.nan

    species_labels = [r'O$^{+}$', r'O$^{++}$', r'S$^{+}$', r'S$^{++}$', r'S$^{+++}$', r'H$^{+}$', r'Na$^{+}$', r'O$^{+}$(hot)', r'eh$^{-}$', r'e$^{-}$']
    
    for i in range(nspec):
        key = species_names[i]
        lbl = species_labels[i]
        n[key] = n_out[i,:]
        plt.plot(s,n[key],label=lbl)
        T_par[key] = Tpar_out[i,:]/ELEMENTARY_CHARGE
        T_perp[key] = Tperp_out[i,:]/ELEMENTARY_CHARGE
        kappa_par[key] = kappapar_out[i,:]
        kappa_perp[key] = kappaperp_out[i,:]
        
    plt.yscale('log')
    plt.title(species_type[0].replace('_', ' ') + " Diffusive Equilibrium, Aligned Io Field Line")
    plt.xlabel(r"s (R$_J$)")
    plt.ylabel(r"n (cm$^{-3}$)")
    plt.legend(fontsize='x-small', ncol=2)
    plt.grid()
    plt.savefig(species_type[0] + "_diff_eq_densities_vs_s_log_aligned_Io_fl.png", dpi=600, bbox_inches="tight")
    plt.show()
    
    
    for i in range(nspec):
        key = species_names[i]
        lbl = species_labels[i]
        plt.plot(lat,n[key],label=lbl)
    plt.yscale('log')
    plt.title(species_type[0].replace('_', ' ') + " Diffusive Equilibrium, Aligned Io Field Line")
    plt.xlabel(r"Latitude$_{\text{III}}$ ($\degree$)")
    plt.ylabel(r"n (cm$^{-3}$)")
    plt.legend(fontsize='x-small', ncol=2)
    plt.grid()
    plt.savefig(species_type[0] + "_diff_all_eq_densities_vs_lat_log_aligned_Io_fl.png", dpi=600, bbox_inches="tight")
    plt.show()
    



    for i in range(nspec):
        key = species_names[i]
        lbl = species_labels[i]
        plt.plot(lat,T_par[key],label=lbl)
    plt.yscale('log')
    plt.title(species_type[0].replace('_', ' ') + " Diffusive Equilibrium, Aligned Io Field Line")
    plt.xlabel(r"Latitude$_{\text{III}}$ ($\degree$)")
    plt.ylabel(r"T$_{\parallel}$ (eV)")
    plt.legend(fontsize='x-small', ncol=2)
    plt.grid()
    plt.savefig(species_type[0] + "_all_diff_eq_Tpar_vs_lat_log_aligned_Io_fl.png", dpi=600, bbox_inches="tight")
    plt.show()
    
    
    for i in range(nspec):
        key = species_names[i]
        lbl = species_labels[i]
        plt.plot(lat,T_perp[key],label=lbl)
    plt.yscale('log')
    plt.title(species_type[0].replace('_', ' ') + " Diffusive Equilibrium, Aligned Io Field Line")
    plt.xlabel(r"Latitude$_{\text{III}}$ ($\degree$)")
    plt.ylabel(r"T$_{\perp}$ (eV)")
    plt.legend(fontsize='x-small', ncol=2)
    plt.grid()
    plt.savefig(species_type[0] + "_all_diff_eq_Tperp_vs_lat_log_aligned_Io_fl.png", dpi=600, bbox_inches="tight")
    plt.show()

    for i in range(nspec):
        key = species_names[i]
        lbl = species_labels[i]
        plt.plot(lat,kappa_par[key],label=lbl)
    plt.yscale('log')
    plt.title(species_type[0].replace('_', ' ') + " Diffusive Equilibrium, Aligned Io Field Line")
    plt.xlabel(r"Latitude$_{\text{III}}$ ($\degree$)")
    plt.ylabel(r"$\kappa_{\parallel}$ (eV)")
    plt.legend(fontsize='x-small', ncol=2)
    plt.grid()
    plt.savefig(species_type[0] + "_all_diff_eq_kappa_par_vs_lat_log_aligned_Io_fl.png", dpi=600, bbox_inches="tight")
    plt.show()
    
    for i in range(nspec):
        key = species_names[i]
        lbl = species_labels[i]
        plt.plot(lat,kappa_perp[key],label=lbl)
    plt.yscale('log')
    plt.title(species_type[0].replace('_', ' ') + " Diffusive Equilibrium, Aligned Io Field Line")
    plt.xlabel(r"Latitude$_{\text{III}}$ ($\degree$)")
    plt.ylabel(r"$\kappa_{\perp}$ (eV)")
    plt.legend(fontsize='x-small', ncol=2)
    plt.grid()
    plt.savefig(species_type[0] + "_all_diff_eq_kappa_perp_vs_lat_log_aligned_Io_fl.png", dpi=600, bbox_inches="tight")
    plt.show()







    """
    # Save dictionary of numpy arrays of plasma paramter outputs
    # number densities, temperatures, and kappa values
    # in addition to electrostatic potential phi in volts 
    # and centrifugal + gravitational potential difference per unit mass deltaU
    np.savez_compressed(species_type[0] + "Diff_eq_plasma_parameter_along_Io_aligned_fl_JRM33+Con2020.npz',
                    n=n, T_par=T_par, T_perp = T_perp, kappa_par = kappa_par, kappa_perp = kappa_perp, phi = phi, deltaU = deltaU)
    """
    """
    # Would be read back in like so:
    with np.load('Diff_eq_plasma_parameter_along_Io_aligned_fl_JRM33+Con2020.npz', allow_pickle=True) as data:
        n = data['n'].item()
        T_par = data['T_par'].item()
        T_perp= data['T_perp'].item()
        kappa_par = data['kappa_par'].item()
        kappa_perp= data['kappa_perp'].item()
        
        #Then can reference  electron density values along field line for example as: nec =  n['e-'], etc. see species_names above
    """



