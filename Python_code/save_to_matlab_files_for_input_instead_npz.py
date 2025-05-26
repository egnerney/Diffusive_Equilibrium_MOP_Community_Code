#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 26 01:19:11 2025

@author: edne8319
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NPZ → MAT-file converter (MATLAB v5 format)

Author : Edward (Eddie) G. Nerney
Updated: 2025-05-26
"""

import numpy as np
from scipy.io import savemat

# -------------------------------------------------------------
species_names = [
    'O+', 'O++', 'S+', 'S++', 'S+++',
    'H+', 'Na+', 'O+(hot)', 'eh-', 'e-'
]

def f_array(arr):
    """Return arr as Fortran-contiguous float64."""
    return np.asfortranarray(arr, dtype='float64')

# -------------------------------------------------------------
# 1. Field-line trace  →  MAT
# -------------------------------------------------------------
in_fl  = 'Io_aligned_field_line_trace_JRM33+con2020_dlat=0.1_degree.npz'
out_fl = 'Io_aligned_field_line_trace_JRM33+con2020_dlat=0.1_degree.mat'

with np.load(in_fl) as src:
    mdict = {key: f_array(src[key]) for key in
             ['s', 'x', 'y', 'z', 'rho', 'r', 'lat', 'wlong', 'B']}
    savemat(out_fl, mdict, oned_as='column')
print(f'✅  wrote {out_fl}')

# -------------------------------------------------------------
# 2. Centrifugal-equator reference model  →  MAT
# -------------------------------------------------------------
in_ref  = 'Nerney2025_1D_radial_ceq_reference_model_4-10RJ_601pts_dr=0.01RJ.npz'
out_ref = 'Nerney2025_reference_model_4-10RJ.mat'

with np.load(in_ref, allow_pickle=True) as src:
    def stack(d):
        return f_array(np.column_stack([d[name] for name in species_names]))

    matdict = {
        'n0'      : stack(src['n_0'].item()),
        'T0'      : stack(src['T_0'].item()),
        'kappa0'  : stack(src['kappa_0'].item()),
        'kappaT0' : stack(src['kappa_temps_0'].item()),
        'rho_ceq' : f_array(0.01*np.arange(601)+4.0),
        # MATLAB cell array of char vectors
        'species' : np.array([[s] for s in species_names], dtype=object)
    }
    savemat(out_ref, matdict, oned_as='column')
print(f'✅  wrote {out_ref}')

