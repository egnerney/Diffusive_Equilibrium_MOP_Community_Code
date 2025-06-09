# ------------------------------------------------------------------
# diffusive-equilibrium Makefile
# ------------------------------------------------------------------
FC      = gfortran
FFLAGS  = -O3
LDFLAGS =
SRC     = npy_reader.f90 special_functions.f90 hyper_2f1.f90 \
          diffusive_equilibrium_mod.f90 basic_example_use_of_diffusive_equilibrium_code.f90
OBJ     = $(SRC:.f90=.o)

.PHONY: all run_diff_eq clean

all: run_diff_eq

run_diff_eq: $(OBJ)
	$(FC) $(LDFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJ) run_diff_eq *.mod
