
Read Me For IDL Diffusive Equilibrium Code:

Tested on IDL 8.8 but should work on any modern IDL. 

Cite **Nerney (2025)** if you use the solver in published work or code.
We give permission to incorporate these routines into a user program if the copyright is acknowledged and Nerney (2025) is cited.


Test Example For Io aligned Field line is in basic_example_use_of_diffusive_equilibrium_code.pro

Which will compile helper functions/procedures from diffusive_equilibrium.pro

and run a Single Io Aligned field line at Jupiter, for example, taking the user through setup and then saving plots of output.


ASSUMES Last species is always core/cold electrons given by neutrality with all other species

Reads in field line trace at Io (JRM33+ Con2020), reference densities, temperatures, and kappa values/kappa temps. Gives reference to the Centrifugal equator from Nerney (2025), a 1D radial Centrifugal equator empirical model to then be extrapolated along a field line.

Available distribution function options set in type_ string in the example are:

'Maxwellian' ; Bagenal & Sullivan (1981) etc. Basic default assumes isotropic always, ignores anisotropy
'Aniso_Maxwellian' ; Huang & Birmingham (1992) Anisotropic Maxwellian
'Aniso_Maxwellian_Const_Tperp' ; Bagenal & Sullivan (1981), Bagenal (1994), Mei, Thorne, Bagenal (1995) Assumed Constant Anisotrpy along field line, inconsistent but good enough often 
'Aniso_Kappa' ; Meyer-Vernet et al. (1995); Moncuquet et al. (2002) Anisotropic Standard Kappa distribution
'Aniso_Product_Kappa' ; Nerney (2025) Anisotropic Product Kappa Distribution 
'Fried_Egg' ; Nerney (2025) "Fried-Egg" Distribution Function (Maxwellian limit in parallel of aniso product kappa and still kappa in perp direction) 
  
For Aniso_Product_Kappa and Fried_Egg, the special functions used (not available in default IDL) are given by hyp2f1. pro and hyperu.pro, which are given by 1D numerical quad integrals given by qpint1d. pro, which are included in the IDL code directory.
These are much slower as a result, and if possible, parallelized code in Python and Fortran should be used instead if running a lot of field lines. 

Output For running the default basic_example_use_of_diffusive_equilibrium_code.pro with only changing type_ is shown in the subdirectory Example_Test_Output

Output from paper used a different fit kappa and kappa temp for fried egg and product kappa and only included hot electrons with double maxwellian model and didn't include hot electrons in neutrality but for the basic example this functionality is not shown but can easily be implimented or found in other github code for the paper. 


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

