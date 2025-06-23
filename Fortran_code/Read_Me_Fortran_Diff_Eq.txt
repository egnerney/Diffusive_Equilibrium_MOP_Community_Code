Read Me For Fortran Diffusive Equilibrium Code:
Tested with gfortran 15.1, but should work on  GNU Fortran ≥ 10 | Any compiler that supports F2003–F2018 should work.
gnuplot is required for plotting only in the example. Tested on gnuplot 6.0 

Cite **Nerney (2025) and code at http://doi.org/10.5281/zenodo.15623974** if you use the solver in published work or code.
We give permission to incorporate these routines into a user program if the copyright is acknowledged and Nerney (2025) is cited.

Test Example For Io aligned Field line is in basic_example_use_of_diffusive_equilibrium_code.f90, or if uncommented in Makefile, MPI parallelized version basic_example_use_of_diffusive_equilibrium_code_mpi.f90, replacing compilation of basic_example_use_of_diffusive_equilibrium_code.f90 run via Makefile via: make all

Tested on FC = gfortran for serial version and FC = mpifort for MPI parallelized version, both installed in Mac terminal via brew install gcc and brew install open-mpi (install open-mpi after gcc to have them properly linked), but can be managed via conda as well or using your preferred compiler or package manager. 

Which will compile helper modules/functions/subroutines from diffusive_equilibrium_mod.f90, then run via ./run_diff_eq or if using the MPI parallelized version, mpirun -np numprocs ./run_diff_eq with numprocs = number of processors you wish to or have available to run the code in parallel on. 

and run a Single Io Aligned field line at Jupiter example taking the user through setup and then saving plots of output (gnuplot required for plots). The parallelized version is implemented via static block decomposition, where the field line is chunked into approximately equal pieces, given the number of available processors specified in the MPI call. 

ASSUMES the last species is always core/cold electrons, given by neutrality with all other species

Reads in field line trace at Io (JRM33+ Con2020), reference densities, temperatures, and kappa values/kappa temps. Gives reference to the Centrifugal equator from Nerney (2025), a 1D radial Centrifugal equator empirical model to then be extrapolated along a field line.

| **Swap distribution function** | Change `dist_tag` before the species loop; each species can differ, but in this test case, they are all the same.
!Avaiable Distribution Types:
   !character(len=32)              :: dist_tag = 'Maxwellian' ! Bagenal & Sullivan (1981) etc. Basic default, assumed isotropic always, ignores anisotropy
   !character(len=32)              :: dist_tag = 'Aniso_Maxwellian' ! Huang & Birmingham (1992) Anisotropic Maxwellian
   !character(len=32)              :: dist_tag = 'Aniso_Maxwellian_Const_Tperp' ! Bagenal & Sullivan (1981), Bagenal (1994), Mei, Thorne, Bagenal (1995) Assumed Constant Anisotrpy along field line, inconsistent but good enough, often 
   !character(len=32)              :: dist_tag = 'Aniso_Kappa' ! Meyer-Vernet et al. (1995), Moncuquet et al. (2002) Anisotropic Standard Kappa distribution
   !character(len=32)              :: dist_tag = 'Aniso_Product_Kappa' ! Nerney (2025) Anisotropic Product Kappa Distribution 
   !character(len=32)              :: dist_tag = 'Fried_Egg' ! Nerney (2025) "Fried-Egg" Distribution Function (Maxwellian limit in parallel of aniso product kappa and still kappa in perp direction) 
 
| File                               | Purpose                                                             |
|------------------------------------|---------------------------------------------------------------------|
| `diffusive_equilibrium_mod.f90`    | Physics core (densities, temperatures, κ, root-finder, etc.)        |
| `special_functions.f90`            | Confluent U implementation (from Zhang & Jin)                       |
| `hyper_2f1.f90`                    | Gauss hyper-geometric ₂F₁ (Michel & Stoitsov)                        |
| `npy_reader.f90`                   | Minimal NPY reader (little-endian, 1-D / 2-D double arrays)  to read in input/reference .npy files        |
| `basic_example_use_of_diffusive_equilibrium_code.f90` | Io field line demo/example use of code                     |


---

## Quick start

```bash
git clone https://github.com/egnerney/Diffusive_Equilibrium_MOP_Community_Code
cd /Diffusive_Equilibrium_MOP_Community_Code/Fortran_code

# Compile & run in one go

Change the compiler and flags in the Makefile to match your use case/system.


make all (with FC=gfortran for serial version or FC=mpifort for mpi and basic_example_use_of_diffusive_equilibrium_code.f90 replaced with basic_example_use_of_diffusive_equilibrium_code_mpi.f90 in Makefile
./run_diff_eq      # For Serial version    # prints progress and writes PNGs in ./output
````

or for mpi parallelized version use:

mpirun -np numprocs ./run_diff_eq with numprocs = number of processors you wish to or have available to run the code in parallel on. 

> **Tip:** On Windows use *WSL* or MSYS2; pass `FC=ifx`, `FC=ifort`, etc.,
> depending on your compiler.

---


* Prints progress (`point 0 / 1340 …`).
* Solves the electrostatic potential at each point via bisection.
* Writes six figures into **`output/`**:

  * `deltaU_vs_lat.png`
  * `<dist>_n_vs_s.png`
  * `<dist>_n_vs_lat.png`
  * `<dist>_Tpar_vs_lat.png`
  * `<dist>_Tperp_vs_lat.png`
  * `<dist>_kappa_par_vs_lat.png`
  * `<dist>_kappa_perp_vs_lat.png`

(Here `<dist>` is, e.g., `Aniso_Maxwellian`, etc....)

---

## Cleaning

```bash
make clean        # removes *.o, *.mod, run_diff_eq
```
___________


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


External special-function subroutines carry their original licenses (see headers).

Cite **Nerney (2025) and code at http://doi.org/10.5281/zenodo.15623974** if you use the solver in published work or code.
We give permission to incorporate this routine into a user program if the copyright is acknowledged and Nerney (2025) is cited.


