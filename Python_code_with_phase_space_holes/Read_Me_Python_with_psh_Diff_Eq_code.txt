Cite **Nerney (2025) and http://doi.org/10.5281/zenodo.15623974 ** if you use the solver in published work or code.
We give permission to incorporate these routines into a user program if the copyright is acknowledged and Nerney (2025) is cited.

Read Me For Python Diffusive Equilibrium Code including psh:
Tested on Python 3.12.10 but should work on any modern Python. Though if using numpy 1.x instead of 2.x, uncomment lines at the top of the test script for loading .npz files due to pickling incompatibility issue. 

Requires numpy and scipy libraries for computation and matplotlib for plotting in basic use test script. Tested with numpy version '2.2.6', scipy version '1.15.2', and matplotlib version '3.10.3'. With uncommented lines at top of test script it also was tested and works with python 3.8.20, numpy '1.24.4', scipy '1.10.1', and matplotlib '3.7.3'. 

These can be installed via conda install numpy scipy matplotlib 
or if you prefer pip intsall numpy scipy matplotlib

Test Example For Io aligned Field line is in basic_example_use_of_diffusive_equilibrium_code.py

Which will load helper functions/procedures from diffusive_equilibrium.py via from diffusive_equilibrium import * at top of code
and run a single Io-aligned field line at Jupiter, for example, taking the user through setup and then saving plots of the output. Uses phase space hole approiximation if above local maxima along field line at s_c. 


ASSUMES the Last species is always core/cold electrons, given by neutrality with all other species

Reads in field line trace at Io (JRM33+ Con2020), reference densities, temperatures, and kappa values/kappa temps. 
Gives reference to the Centrifugal equator from Nerney (2025), a 1D radial Centrifugal equator empirical model to then be extrapolated along a field line.

Available distribution function options set in species_type =  [string] in the example are:



'Maxwellian' # Bagenal & Sullivan (1981) etc. Basic default, Assumed isotropic always, ignores anisotropy
'Aniso_Maxwellian' # Huang & Birmingham (1992) Anisotropic Maxwellian
'Aniso_Maxwellian_Const_Tperp' # Bagenal & Sullivan (1981), Bagenal (1994), Mei, Thorne, Bagenal (1995) Assumed Constant Anisotrpy along field line, inconsistent but good enough often 
'Aniso_Kappa' # Meyer-Vernet et al. (1995); Moncuquet et al. (2002) Anisotropic Standard Kappa distribution
'Aniso_Product_Kappa' # Nerney (2025) Anisotropic Product Kappa Distribution 
'Fried_Egg' # Nerney (2025) "Fried-Egg" Distribution Function (Maxwellian limit in parallel of aniso product kappa and still kappa in perp direction) 

They need not all be the same for each species

Output For running the default basic_example_use_of_diffusive_equilibrium_code.py with only changing type_ is shown in subdirectory Example_Test_Output

Output from paper used a different fit kappa and kappa temp for fried egg and product kappa and standard kappa this is just the one for standard kappa and only included hot electrons with double maxwellian model and didn't include hot electrons in neutrality but for the basic example this functionality is not shown but can easily be implimented or found in other github code for the paper.


## License

MIT License

Copyright (c) 2025 Edward (Eddie) G. Nerney

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


Cite **Nerney (2025) http://doi.org/10.5281/zenodo.15623974 ** if you use the solver in published work or code.
We give permission to incorporate this routine into a user program if the copyright is acknowledged and Nerney (2025) is cited.
