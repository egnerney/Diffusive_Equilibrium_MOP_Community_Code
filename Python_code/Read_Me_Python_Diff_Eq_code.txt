Read Me For python Diffisuve Equilbrium Code:
Tested on Python 3.12.10 but should work on any modern python.

Requires numpy and scipy libraries. 

Test Example For Io aligned Field line is in basic_example_use_of_diffusive_equilibrium_code.py

Which will load helper functions/procedures from diffusive_equilibrium.py via from diffusive_equilibrium import * at top of code. 

and run a Single Io Aligned field line at Jupiter example taking the user through setup and then saving plots of output.


ASSUMES Last species is always core/cold electrons given by neutrality with all other species

Reads in field line trace at Io (JRM33+ Con2020), reference densities, temperatures, and kappa values/kappa temps. Gives reference at Centrifugal equator from Nerney (2025) 1D radial Centrifugal equator emperical model to then be extrapolated along a field line.

Available distribution function options set in species_type =  [string] in example are:



'Maxwellian' # Bagenal & Sullivan (1981) etc. Basic default Assumed isotropic always, ignores anisotropy
'Aniso_Maxwellian' # Huang & Birmingham (1992) Anisotropic Maxwellian
'Aniso_Maxwellian_Const_Tperp' # Bagenal & Sullivan (1981), Bagenal (1994), Mei, Thorne, Bagenal (1995) Assumed Constant Anisotrpy along field line, inconsistent but good enough often 
'Aniso_Kappa' # Meyer-Vernet et al. (1995); Moncuquet et al. (2002) Anisotropic Standard Kappa distribution
'Aniso_Product_Kappa' # Nerney (2025) Anisotropic Product Kappa Distribution 
'Fried_Egg' # Nerney (2025) "Fried-Egg" Distribution Function (Maxwellian limit in parallel of aniso product kappa and still kappa in perp direction) 

They need not all be the same for each species

Output For running the default basic_example_use_of_diffusive_equilibrium_code.py with only changing type_ is shown in subdirectory Example_Test_Output

Output from paper used a different fit kappa and kappa temp for fried egg and product kappa and standard kappa this is just the one for standard kappa and only included hot electrons with double maxwellian model and didn't include hot electrons in neutrality but for the basic example this functionality is not shown but can easily be implimented or found in other github code.
