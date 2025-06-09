Read Me For MATLAB Diffusive Equilibrium Code:
Tested with MATLAB_R2024b but should work with any modern version. Only basic MATLAB with no external toolkits required. 

Cite **Nerney (2025)** if you use the solver in published work or code.
We give permission to incorporate these routines into a user program if the copyright is acknowledged and Nerney (2025) is cited.

Test Example For Io aligned Field line can be run from Basic_Example_use_of_Diffusive_Equilibrium.m

Which will load helper modules/functions from plasma.m and Species.m

and run a Single Io Aligned field line at Jupiter example taking the user through setup and then saving plots of output

ASSUMES Last species is always core/cold electrons given by neutrality with all other species

Reads in field line trace at Io (JRM33+ Con2020), reference densities, temperatures, and kappa values/kappa temps. Gives reference at Centrifugal equator from Nerney (2025) 1D radial Centrifugal equator emperical model to then be extrapolated along a field line.

%| **Swap distribution function** | Change `dist_type` before the species loop; each species can differ but in this test case they are all the same.
%Avaiable Distribution Types:
%dist_type = repmat({'Maxwellian'}, 1, nspec); % Bagenal & Sullivan (1981) etc. Basic default Assumed isotropic always, ignores anisotropy
%dist_type = repmat({'Aniso_Maxwellian'}, 1, nspec); % Huang & Birmingham (1992) Anisotropic Maxwellian
%dist_type = repmat({'Aniso_Maxwellian_const_Tperp'}, 1, nspec); % Bagenal & Sullivan (1981), Bagenal (1994), Mei, Thorne, Bagenal (1995) Assumed Constant Anisotrpy along field line, inconsistent but good enough often 
%dist_type = repmat({'Aniso_kappa'}, 1, nspec); % Meyer-Vernet et al. (1995), Moncuquet et al. (2002) Anisotropic Standard Kappa distribution
%dist_type = repmat({'Aniso_product_kappa'}, 1, nspec); % Nerney (2025) Anisotropic Product Kappa Distribution 
dist_type = repmat({'Fried_Egg'}, 1, nspec); % Nerney (2025)
%"Fried-Egg" Distribution Function (Maxwellian limit in parallel of
%aniso product kappa and still kappa in perp direction)

Output from paper used a different fit kappa and kappa temp for fried
egg and product kappa and only included hot electrons with double
maxwellian model and didn't include hot electrons in neutrality but
for the basic example this functionality is not shown but can easily
be implimented or found in other github code for the paper. 
 


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


Cite **Nerney (2025)** if you use the solver in published work or code.
We give permission to incorporate this routine into a user program if the copyright is acknowledged and Nerney (2025) is cited.


