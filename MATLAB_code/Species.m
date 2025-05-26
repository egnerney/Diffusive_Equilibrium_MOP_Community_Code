classdef Species
    %{
    MATLAB class mirroring your Python "Species" container.

    Attributes:
    -----------
    name : str
        Species name (e.g., 'O+', 'e-', 'S++', etc.).
    m : float
        Mass in AMU (or in kg if you convert it beforehand).
    q : float
        Charge in Coulombs (can be negative for electrons).
    T : float
        Characteristic temperature in eV or Joules (must be consistent with usage).
    A : float
        Anisotropy factor A0 = T_perp / T_par at reference location, if relevant.
    kappa : float
        Kappa parameter. For Maxwellian, kappa -> infinity in theory.
    lam : float
        λ = kappa_perp / kappa_par (only needed in product-kappa or exotic models).
    n : float
        Number density at reference location (in cm^-3 or consistent units).
    type : str
        Distribution type. Must match a key we call in the "plasma" class
        (e.g. 'Maxwellian', 'Aniso_kappa', 'Fried_Egg', etc.)
    %}
    properties
        name   (1,:) char
        m      double
        q      double
        T      double
        A      double
        kappa  double
        lam    double
        n      double
        type   (1,:) char
    end

    methods
        function obj = Species(name, m, q, T, A, kappa, lam, n, dist_type)
            % Constructor to replicate Python's __init__
            if nargin > 0
                obj.name  = name;
                obj.m     = m;
                obj.q     = q;
                obj.T     = T;
                obj.A     = A;
                obj.kappa = kappa;
                obj.lam   = lam;
                obj.n     = n;
                obj.type  = dist_type;
            end
        end
    end
end
