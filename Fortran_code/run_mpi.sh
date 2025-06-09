#!/bin/bash
# ------------------------------------------------------------------
# Script to compile and run MPI version of diffusive equilibrium code
# For macOS with Homebrew-installed MPI
# ------------------------------------------------------------------

# Check if MPI is installed
if ! command -v mpif90 &> /dev/null; then
    echo "MPI not found. Installing via Homebrew..."
    # Check if Homebrew is installed
    if ! command -v brew &> /dev/null; then
        echo "Homebrew not found. Please install Homebrew first:"
        echo "/bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""
        exit 1
    fi
    # Install Open MPI
    brew install open-mpi
fi

# Check for required Fortran compiler
if ! command -v gfortran &> /dev/null; then
    echo "gfortran not found. Installing via Homebrew..."
    brew install gcc
fi

# Clean previous builds
echo "Cleaning previous builds..."
make clean

# Compile the MPI version
echo "Compiling MPI version..."
make -f Makefile

# Check if compilation was successful
if [ ! -f run_diff_eq_mpi ]; then
    echo "Compilation failed!"
    exit 1
fi

# Set number of processes (adjust based on your system)
# Get number of CPU cores on Mac
NCORES=$(sysctl -n hw.physicalcpu)
echo "Detected $NCORES physical CPU cores"

# Use all cores or specify a number
NPROCS=$NCORES
# Uncomment to use specific number:
# NPROCS=4

echo "Running with $NPROCS MPI processes..."

# Run the MPI program
mpirun -np $NPROCS ./run_diff_eq_mpi

# Check if run was successful
if [ $? -eq 0 ]; then
    echo "MPI run completed successfully!"
    echo "Output PNG files:"
    ls -la *.png
else
    echo "MPI run failed!"
    exit 1
fi
