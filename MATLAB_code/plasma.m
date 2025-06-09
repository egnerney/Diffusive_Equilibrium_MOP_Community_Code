classdef plasma
    %{
    A single "toolbox" class containing static methods.

    Includes:
      - define_planet
      - find_local_maxima
      - calc_deltaU
      - diff_eq
      - calc_temps
      - calc_kappa_vals

    And private static methods that match each distribution function:
      maxwellian_density, aniso_maxwellian_density, aniso_kappa_density, ...
      plus "fried_egg_density" with all the if-else expansions,
      plus the product-kappa expansions,
      plus parallel T, perpendicular T, and kappa subroutines.

    We assume you have hyp2f1.m and hyperu.m in the same folder. The code calls
    them to handle hypergeometric integrals
    %}
    methods (Static)
        %==================================================================
        %                PUBLIC API  (equivalents of your Python)
        %==================================================================

        function [RP, Omega, GM] = define_planet(planet)
            %{
            Returns planet-specific parameters: radius (m), rotation rate (rad/s),
            and GM (m^3/s^2). 
            %}
            switch lower(planet)
                case 'earth'   % NASA Earth Fact Sheet
                    RP = 3.3781e6;
                    Omega = 7.29210e-5;    % rad/s
                    GM = 3.9860e14;       % m^3/s^2
                case 'jupiter' % NASA Jupiter Fact Sheet
                    RP = 7.1492e7;
                    Omega = 1.7585e-4;    % rad/s
                    GM = 1.26687e17;
                case 'saturn'
                    RP = 6.0268e7;
                    Omega = 1.6379e-4;
                    GM = 3.7931e16;
                case 'uranus'
                    RP = 2.5559e7;
                    Omega = 1.012e-4;
                    GM = 5.7940e15;
                case 'neptune'
                    RP = 2.4764e7;
                    Omega = 1.083e-4;
                    GM = 6.8351e15;
                otherwise
                    error('Planet %s is not supported in define_planet().', planet);
            end
        end


        function [max_indices, max_lats, max_deltaU] = find_local_maxima(deltaU, lat)
            %{
            Find local maxima in deltaU as a function of lat

            deltaU : 1D array
            lat    : 1D array (degrees)
            returns:
               max_indices -> indices of local maxima
               max_lats    -> latitudes of maxima
               max_deltaU  -> values of deltaU at those maxima
            %}
            if numel(deltaU) ~= numel(lat)
                error('deltaU and lat must have same length');
            end

            % Convert to double row vectors
            deltaU = double(deltaU(:));
            lat    = double(lat(:));

            % Sort by latitude
            [lat_sorted, sort_idx] = sort(lat);
            deltaU_sorted = deltaU(sort_idx);

            if length(deltaU_sorted) < 3
                max_indices = [];
                max_lats    = [];
                max_deltaU  = [];
                return
            end

            max_sorted_indices = [];
            for i = 2:length(deltaU_sorted)-1
                if deltaU_sorted(i) > deltaU_sorted(i-1) && ...
                   deltaU_sorted(i) > deltaU_sorted(i+1)
                    max_sorted_indices(end+1) = i; %#ok<AGROW>
                end
            end

            if isempty(max_sorted_indices)
                max_indices = [];
                max_lats    = [];
                max_deltaU  = [];
            else
                max_indices = sort_idx(max_sorted_indices);
                max_lats    = lat(max_indices);
                max_deltaU  = deltaU(max_indices);
            end
        end


        function dU = calc_deltaU(r0_in, lat0_in, r_in, lat_in, planet, fcor)
            %{
            Replicates your Python calc_deltaU. 
            r0_in, r_in : dimensionless radial coords (in planet radii).
            lat0_in, lat_in : latitudes in degrees.
            planet: string
            fcor: fraction of corotation (1 => full corotation)

            returns:
              dU : difference in gravitational + centrifugal potential
                   per unit mass, from (r0,lat0) to (r,lat).
            %}
            [RP, Omega, GM] = plasma.define_planet(planet);
            Omega = Omega * fcor;

            % Convert dimensionless radial coords to meters
            r0  = r0_in .* RP;     % if r0_in is array or scalar use '.*'
            r   = r_in  .* RP;
        
            % Convert lat from degrees to radians
            lat0 = deg2rad(lat0_in);
            lat  = deg2rad(lat_in);
        
            % Potential at reference location (scalar or array)
            U0 = -GM ./ r0 ...
                 - 0.5 .* (r0.^2) .* (cos(lat0).^2) .* (Omega.^2);
        
            % Potential at new location(s)
            U  = -GM ./ r ...
                 - 0.5 .* (r.^2)  .* (cos(lat).^2 ) .* (Omega.^2);

            dU = U - U0;
        end


        function [n, phi, dU_out] = diff_eq(Pos_in, Pos0_in, Bratio, species0, deltaU, planet, fcor)
            %{
            Solves for local electrostatic potential Phi (via charge neutrality)
            and returns densities n for each species. 

            1) Enforce neutrality at reference location by adjusting electron density.
            2) net_charge_density(Phi) -> bisection -> find root => Phi_root
            3) Return that Phi and densities n.

            We step outward from Phi=0 in +/- directions
            to bracket the root, then calling bisect.
            %}
            nspec = numel(species0);
            if nspec < 2
                error('Need at least 2 species (ions + electrons).');
            end

            %--- enforce neutrality at reference location s0
            total_ion_charge_density = 0.0;
            for i = 1:nspec-1
                total_ion_charge_density = total_ion_charge_density + ...
                                           species0(i).q * species0(i).n;
            end
            % last species is electrons => adjust n
            species0(nspec).n = -total_ion_charge_density / species0(nspec).q;

            % nested function for net charge density
            function nq = net_charge_density(Phi_test)
                n_local_temp = zeros(1, nspec);
                for j = 1:nspec
                    n_local_temp(j) = plasma.call_density_func(species0(j), deltaU, Phi_test, Bratio);
                end
                nq_temp = zeros(1, nspec);
                for j = 1:nspec
                    nq_temp(j) = species0(j).q * n_local_temp(j);
                end
                nq = sum(nq_temp);
            end

            %--- bracket the root
            Phi1 = 0.0;  Phi2 = 0.0;
            dPhi = 1.0;
            nq1 = net_charge_density(Phi1);
            nq2 = net_charge_density(Phi2);

            max_iterations = 10000;
            iter_count = 0;
            while (nq1 * nq2 > 0) && (iter_count < max_iterations)
                Phi1 = Phi1 - dPhi;
                Phi2 = Phi2 + dPhi;
                nq1 = net_charge_density(Phi1);
                nq2 = net_charge_density(Phi2);
                iter_count = iter_count + 1;
            end
            if iter_count >= max_iterations
                error('Failed to bracket root for net_charge_density in diff_eq.');
            end

            %--- use bisection to find the zero
            % We'll use fzero with bracket [Phi1, Phi2]
            f = @(P) net_charge_density(P);
            Phi_root = fzero(f, [Phi1, Phi2]);

            %--- final densities
            n = zeros(1, nspec);
            for i = 1:nspec
                n(i) = plasma.call_density_func(species0(i), deltaU, Phi_root, Bratio);
            end
            phi = Phi_root;
            dU_out = deltaU;
        end


        function [Tpar, Tperp] = calc_temps(species_list, deltaU, phi, Bratio)
            %{
            For each species, call the appropriate temperature function 
            to get Tpar, Tperp at location s given (deltaU, phi, Bratio).
            %}
            nspec = numel(species_list);
            Tpar  = zeros(1, nspec);
            Tperp = zeros(1, nspec);

            for i = 1:nspec
                [Tpar(i), Tperp(i)] = plasma.call_temperature_func(species_list(i), ...
                    deltaU, phi, Bratio);
            end
        end


        function [kappa_par, kappa_perp] = calc_kappa_vals(species_list, deltaU, phi, Bratio)
            %{
            For each species, find the local kappa_par and kappa_perp parameters 
            %}
            nspec = numel(species_list);
            kappa_par  = zeros(1, nspec);
            kappa_perp = zeros(1, nspec);

            for i = 1:nspec
                [kappa_par(i), kappa_perp(i)] = plasma.call_kappa_func(species_list(i), ...
                    deltaU, phi, Bratio);
            end
        end


        %==================================================================
        %              DISPATCHERS to the private kernel methods
        %==================================================================
        function val = call_density_func(sp, deltaU, Phi, Br)
            switch sp.type
                case 'Maxwellian'
                    val = plasma.maxwellian_density(sp, deltaU, Phi, Br);
                case 'Aniso_Maxwellian'
                    val = plasma.aniso_maxwellian_density(sp, deltaU, Phi, Br);
                case 'Aniso_Maxwellian_const_Tperp'
                    val = plasma.aniso_maxwellian_const_Tperp_density(sp, deltaU, Phi, Br);
                case 'Aniso_kappa'
                    val = plasma.aniso_kappa_density(sp, deltaU, Phi, Br);
                case 'Aniso_product_kappa'
                    val = plasma.aniso_product_kappa_density(sp, deltaU, Phi, Br);
                case 'Fried_Egg'
                    val = plasma.fried_egg_density(sp, deltaU, Phi, Br);
                otherwise
                    error('Unknown distribution type: %s', sp.type);
            end
        end

        function [Tp, Tperp] = call_temperature_func(sp, deltaU, Phi, Br)
            switch sp.type
                case 'Maxwellian'
                    [Tp, Tperp] = plasma.maxwellian_temperature(sp, deltaU, Phi, Br);
                case 'Aniso_Maxwellian'
                    [Tp, Tperp] = plasma.aniso_maxwellian_temperature(sp, deltaU, Phi, Br);
                case 'Aniso_Maxwellian_const_Tperp'
                    [Tp, Tperp] = plasma.aniso_maxwellian_const_Tperp_temperature(sp, deltaU, Phi, Br);
                case 'Aniso_kappa'
                    [Tp, Tperp] = plasma.aniso_kappa_temperature(sp, deltaU, Phi, Br);
                case 'Aniso_product_kappa'
                    [Tp, Tperp] = plasma.aniso_product_kappa_temperature(sp, deltaU, Phi, Br);
                case 'Fried_Egg'
                    [Tp, Tperp] = plasma.fried_egg_temperature(sp, deltaU, Phi, Br);
                otherwise
                    error('Unknown distribution type: %s', sp.type);
            end
        end

        function [kp_par, kp_perp] = call_kappa_func(sp, deltaU, Phi, Br)
            switch sp.type
                case 'Maxwellian'
                    [kp_par, kp_perp] = plasma.maxwellian_kappa(sp, deltaU, Phi, Br);
                case 'Aniso_Maxwellian'
                    [kp_par, kp_perp] = plasma.aniso_maxwellian_kappa(sp, deltaU, Phi, Br);
                case 'Aniso_Maxwellian_const_Tperp'
                    [kp_par, kp_perp] = plasma.aniso_maxwellian_const_Tperp_kappa(sp, deltaU, Phi, Br);
                case 'Aniso_kappa'
                    [kp_par, kp_perp] = plasma.aniso_kappa_kappa(sp, deltaU, Phi, Br);
                case 'Aniso_product_kappa'
                    [kp_par, kp_perp] = plasma.aniso_product_kappa_kappa(sp, deltaU, Phi, Br);
                case 'Fried_Egg'
                    [kp_par, kp_perp] = plasma.fried_egg_kappa(sp, deltaU, Phi, Br);
                otherwise
                    error('Unknown distribution type: %s', sp.type);
            end
        end

    end % public static methods


    methods (Static, Access=private)
        %==================================================================
        %                DENSITY FUNCTIONS (private)
        %==================================================================

        function n_local = maxwellian_density(species, deltaU, Phi, Bratio)
            % Isotropic Maxwellian
            exponent_term = -(species.m*deltaU + species.q*Phi) / species.T;
            n_local = species.n * exp(exponent_term);
        end


        function n_local = aniso_maxwellian_density(species, deltaU, Phi, Bratio)
            % Anisotropic Maxwellian distribution
            B0_over_B = 1./Bratio;
            denom = species.A + (1 - species.A)*B0_over_B;
            if denom <= 0, denom = eps; end

            exponent_term = -(species.m*deltaU + species.q*Phi) / species.T;
            factor_exp = exp(exponent_term);
            n_local = (species.n / denom) * factor_exp;
        end


        function n_local = aniso_maxwellian_const_Tperp_density(species, deltaU, Phi, Bratio)
            % Anisotropic Maxwellian with constant T_perp
            exponent_term = -(species.m*deltaU + species.q*Phi) / species.T;
            % n(s) ~ n0 * (B/B0)^(1 - A) * exp[-(m dU + q Phi)/T_par]
            power_factor = Bratio^(1 - species.A);
            n_local = species.n * power_factor .* exp(exponent_term);
        end


        function n_local = aniso_kappa_density(species, deltaU, Phi, Bratio)
            % Standard anisotropic kappa distribution
            B0_over_B = 1./Bratio;
            denom = species.A + (1 - species.A)*B0_over_B;
            if denom <= 0, denom = eps; end

            PEratio = (species.m * deltaU + species.q * Phi) / (species.kappa * species.T);
            if PEratio <= -1
                PEratio = -1 + eps;
            end

            exponent_factor = (1 + PEratio)^(0.5 - species.kappa);
            n_local = (species.n / denom) * exponent_factor;
        end


        function n_local = aniso_product_kappa_density(species, deltaU, Phi, Bratio)
            %{
            Product kappa distribution with possibly different kappa_par & kappa_perp.
            hypergeometric approach if (k0<17 & k_perp0<17).
            Otherwise large kappa fallback.
            %}
            k0 = species.kappa;
            k_perp0 = k0 * species.lam;
            k_tot0 = k0 + k_perp0;
            n0 = species.n;
            T_par0 = species.T;
            A0 = species.A;

            if Bratio < 1.0
                Bratio = 1.0 + eps;
            end
            beta = (species.m*deltaU + species.q*Phi)/(k0*T_par0);
            if beta <= -1
                beta = -1 + eps;
            end
            x = A0*(Bratio - 1.0);

            if Bratio ~= 1.0
                if (k0 < 17) && (k_perp0 < 17)
                    % Use hypergeometric approach 
                    delta_val = (species.lam * species.A) * ((Bratio - 1.) / (1. + beta));
                    % we define
                    % hf = hyp2f1(1., k_perp0+1.0, k_tot0+1.5, 1.-1./delta_val)
                    hf = hyp2f1(1.0, k_perp0+1.0, k_tot0+1.5, 1.0 - 1./delta_val);

                    factor = k0 / (k_tot0 + 0.5);
                    outer = (n0 / x) * ( (1.0 + beta)^(0.5 - k0) );
                    n_local = Bratio * outer * factor * hf;
                else
                    % fallback for large kappa
                    exponent_factor = Bratio * (1.0 + beta)^(0.5 - k0);
                    % n_local ~ n0 * exponent_factor / (x + 1 + beta)
                    hypfac = 1.0 / (x + 1.0 + beta);
                    n_local = n0 * exponent_factor * hypfac;
                end
            else
                n_local = n0 * (1.0 + beta)^(0.5 - k0);
            end
        end


        function n_local = fried_egg_density(species, deltaU, Phi, Bratio)
            %{
            "Fried-Egg" distribution: Maxwellian in parallel, kappa in perpendicular.
            %}
            if Bratio < 1.0
                fprintf('Warning: B/B0 < 1 for fried_egg_density => set Bratio=1+eps.\n');
                Bratio = 1.0 + eps;
            end

            k0 = species.kappa;    % kappa_perp
            A0 = species.A;        % Tperp0 / Tpar0
            x  = A0*(Bratio - 1.0);
            % exponent = -(m*deltaU + q*Phi)/Tpar
            exponent = -(species.m*deltaU + species.q*Phi)/species.T;
            % clamp exponent if huge
            if abs(exponent) > 700
                exponent = sign(exponent)*700;
            end
            exp_factor = exp(exponent);

            if Bratio ~= 1.0
                if x < 0.0001
                    % tiny x expansion
                    factor = 1.0;
                    n_local = species.n * Bratio * factor * exp_factor;
                elseif (k0 >= 15.0) && (x <= 1e8)
                    factor = (1.0/(1.0 + x)) - x/(k0*((1.0 + x)^3));
                    n_local = species.n * Bratio * factor * exp_factor;
                elseif (k0 < 15.0) && (x > 15.0) && (x < 30.0)
                    % 3rd order expansion in 1/k0
                    factor = 1.0/(1.0 + x) ...
                             - x/(k0*((1.0 + x)^3)) ...
                             + x*(2.0*x - 1.0)/( (k0^2)*((1.0 + x)^5) ) ...
                             - (6.0*(x^3) - 8.0*(x^2) + x)/( (k0^3)*((1.0 + x)^7) );
                    n_local = species.n * Bratio * factor * exp_factor;
                elseif (k0 < 15.0) && (x >= 30.0)
                    % 4th order expansion in 1/x
                    numer = -6.0 + k0*(-11.0 + 2.0*x - (6.0 + (-3.0 + x)*x)*k0 + (-1.0 + x)*(1.0 + x^2)*(k0^2));
                    denom = (x^4)*(k0^3);
                    factor = numer/denom;
                    n_local = species.n * Bratio * factor * exp_factor;
                else
                    % fallback => call hyperu approach 
                    z = k0*x;
                    func_value = k0 * hyperu(1.0, 1.0 - k0, z);
                    n_local = species.n * Bratio * func_value * exp_factor;
                end
            else
                % B/B0=1 => x=0 => n=n0 * exp_factor
                n_local = species.n * exp_factor;
            end
        end


        %==================================================================
        %         TEMPERATURE FUNCTIONS (private)
        %==================================================================

        function [T_par_local, T_perp_local] = maxwellian_temperature(species, deltaU, Phi, Bratio)
            % Isotropic Maxwellian => both Tpar, Tperp = T
            T_par_local  = species.T;
            T_perp_local = species.T;
        end


        function [T_par_local, T_perp_local] = aniso_maxwellian_temperature(species, deltaU, Phi, Bratio)
            % Anisotropic Maxwellian
            T_par_local = species.T;

            B0_over_B = 1./Bratio;
            denom = species.A + (1.0 - species.A)*B0_over_B;
            if denom <= 0, denom = eps; end

            T_perp_local = species.A * species.T / denom;
        end


        function [T_par_local, T_perp_local] = aniso_maxwellian_const_Tperp_temperature(species, deltaU, Phi, Bratio)
            % T_perp also constant => T_par also constant in your code
            T_par_local  = species.T;
            T_perp_local = species.A * species.T;
        end


        function [T_par_local, T_perp_local] = aniso_kappa_temperature(species, deltaU, Phi, Bratio)
            % standard anisotropic kappa
            PEratio = (species.m*deltaU + species.q*Phi) / (species.kappa * species.T);
            if PEratio <= -1
                PEratio = -1 + eps;
            end
            beta = PEratio;

            T_par_local = species.T * (1.0 + beta);

            B0_over_B = 1./Bratio;
            denom = species.A + (1.0 - species.A)*B0_over_B;
            if denom <= 0, denom = eps; end

            T_perp_local = species.A*species.T*(1.0 + beta)/denom;
        end


        function [T_par_local, T_perp_local] = aniso_product_kappa_temperature(species, deltaU, Phi, Bratio)
           
            if Bratio < 1.0
                fprintf('Warning: B/B0 < 1 in aniso_product_kappa_temperature => set Bratio=1+eps.\n');
                Bratio = 1.0 + eps;
            end
            beta = (species.m*deltaU + species.q*Phi)/(species.kappa*species.T);
            if beta <= -1
                beta = -1 + eps;
            end

            kappa_par = species.kappa;
            kappa_perp = kappa_par*species.lam;
            kappa_tot = kappa_par + kappa_perp;
            A0 = species.A;
            T_par0 = species.T;
            T_perp0 = A0*T_par0;
            x = A0*(Bratio - 1.0);
            g = (1.0 + beta);

            if Bratio ~= 1.0
                if (kappa_par < 17.0) && (kappa_perp < 17.0)
                    % Full "hypergeometric" expansions 
                    
                    delta_val = (species.lam * species.A)*((Bratio -1.)/(1.+beta));

                    % We define a set of hypergeometric calls. 
                   
                    % F1 = hyp2f1(1.0, kappa_perp+1.0, kappa_tot+1.5, 1.0 - 1./delta_val)
                    F1 = hyp2f1(1.0, kappa_perp+1.0, kappa_tot+1.5, 1.0 - 1./delta_val);
                    F2 = hyp2f1(2.0, kappa_perp+1.0, kappa_tot+1.5, 1.0 - 1./delta_val);
                    F3q= hyp2f1(0.75, kappa_perp+1.0, kappa_tot+1.5, 1.0 - 1./delta_val);
                    F7q= hyp2f1(1.75,kappa_perp+1.0, kappa_tot+1.5, 1.0 - 1./delta_val);

                    % H1h, H3q, H3h, H7q similarly
                    H1h= hyp2f1(1.0, kappa_perp+1.0, kappa_tot+0.5,  1.0 - 1./delta_val)/gamma(kappa_tot+0.5);
                    H3q= hyp2f1(1.0, kappa_perp+1.0, kappa_tot+0.75, 1.0 - 1./delta_val)/gamma(kappa_tot+0.75);
                    H3h= F1/gamma(kappa_tot+1.5);
                    H7q= hyp2f1(1.0, kappa_perp+1.0, kappa_tot+1.75, 1.0 - 1./delta_val)/gamma(kappa_tot+1.75);

                    C1 = 4.0*kappa_tot - 1.0;
                    C2 = 2.0*kappa_tot - 1.0;
                    C3 = 4.0*kappa_par - 1.0;
                    C4 = 2.0*kappa_par - 1.0;

                    % T_perp
                    T_perp_num = (beta+1.0)*kappa_par*T_perp0*(A0 + x)/(3.0*A0*x);
                    Dperp1 = (C1*F3q)/(3.0*F7q);
                    Dperp2 = (C2*F1)/(2.0*F2);
                    T_perp_denom = Dperp1 - Dperp2;
                    T_perp_local = T_perp_num / T_perp_denom;

                    % T_par
                    T_par_num = 8.0*(beta+1.0)*kappa_par*T_par0;
                    Dpar1 = C3*C1*H7q / H3q;
                    Dpar2 = 2.0*C4*C2*H3h / H1h;
                    T_par_denom = Dpar1 - Dpar2;
                    T_par_local = T_par_num / T_par_denom;

                elseif (kappa_par >= 17.0) && (kappa_perp >= 17.0)
                    % Large kappa => simpler approach
                    T_par_local  = T_par0*g;
                    T_perp_local = Bratio*T_perp0*g / (x + g);
                else
                    % mid-range fallback
                    T_par_local  = T_par0*g;
                    T_perp_local = Bratio*T_perp0*g/(x + g);
                end
            else
                % if Bratio=1
                T_perp_local = T_perp0;
                T_par_local  = T_par0;
            end
        end


        function [T_par_local, T_perp_local] = fried_egg_temperature(species, deltaU, Phi, Bratio)
           
            if Bratio < 1.0
                fprintf('Warning: B/B0 < 1 in fried_egg_temperature => set Bratio=1+eps.\n');
                Bratio = 1.0 + eps;
            end

            % T_parallel is constant for Fried-Egg
            T_par_local = species.T;
            T_perp0 = species.T * species.A;

            k0 = species.kappa; 
            x = species.A*(Bratio - 1.0);
            z = k0*x;

            if Bratio ~= 1.0
                if k0 > 100.0
                    % 0th order approximation => same as aniso Maxwellian
                    factor = 1.0/(1.0 + x);
                    T_perp_local = Bratio*T_perp0*factor;
                elseif k0 >= 15.0
                    % 1st order in 1/k0
                    factor = 1.0/(1.0 + x) - x/(k0*((1.0 + x)^3));
                    T_perp_local = Bratio*T_perp0*factor;
                elseif (x >= 7.0) && (k0 < 15.0)
                    % For x>7
                    C1 = (3.0 + k0)*(4.0 + k0*(7.0 + k0));
                    C2 = -x*((1.0 + k0)^2)*(4.0 + k0);
                    C3 = -(x^3)*(k0^2)*(1.0 + k0) + (x^2)*k0*(1.0 + k0)*(2.0 + k0);
                    C4 = (x^4)*(k0^3) + C3 + C2 + C1;
                    numerator = -4.0 + k0*C4;
                    denominator = (x^5)*(k0^4);
                    factor = numerator/denominator;
                    T_perp_local = Bratio*T_perp0*factor;
                else
                    % final "full approach" => hyperu
                    U0 = hyperu(k0 + 1.0, k0, z);
                    U1 = hyperu(k0 + 1.0, k0 + 1.0, z);
                    U1q= hyperu(k0 + 1.0, k0 + 0.25, z);
                    U5q= hyperu(k0 + 1.0, k0 + 1.25, z);

                    denominator = 4.* (U5q/U1q) - 3. *(U1/U0);
                    factor = (1./(x*denominator));
                    T_perp_local = Bratio*T_perp0*factor;
                end
            else
                T_perp_local = T_perp0;
            end
        end


        %==================================================================
        %                   KAPPA FUNCTIONS (private)
        %==================================================================
        function [kp_par, kp_perp] = maxwellian_kappa(species, deltaU, Phi, Bratio)
            % Maxwellian => effectively infinite kappa
            kp_par = Inf;
            kp_perp = Inf;
        end

        function [kp_par, kp_perp] = aniso_maxwellian_kappa(species, deltaU, Phi, Bratio)
            % same => infinite
            kp_par  = Inf;
            kp_perp = Inf;
        end

        function [kp_par, kp_perp] = aniso_maxwellian_const_Tperp_kappa(species, deltaU, Phi, Bratio)
            kp_par = Inf;
            kp_perp= Inf;
        end

        function [kp_par, kp_perp] = aniso_kappa_kappa(species, deltaU, Phi, Bratio)
            % standard anisotropic kappa => same parallel/perp k
            kp_par  = species.kappa;
            kp_perp = species.kappa;
        end

        function [kp_par, kp_perp] = aniso_product_kappa_kappa(species, deltaU, Phi, Bratio)
           
            if Bratio < 1.0
                fprintf('Warning: B/B0 < 1 in aniso_product_kappa_kappa => set Bratio=1+eps.\n');
                Bratio = 1.0 + eps;
            end
            beta = (species.m*deltaU + species.q*Phi)/(species.kappa*species.T);
            if beta <= -1
                beta = -1 + eps;
            end

            kappa_par = species.kappa;
            kappa_perp = kappa_par*species.lam;
            kappa_tot = kappa_par + kappa_perp;

            if Bratio ~= 1.0
                if (kappa_par >= 30.0) && (kappa_perp >= 30.0)
                    kp_par  = kappa_par;
                    kp_perp = kappa_perp;
                else
                    % Full hypergeom approach from your Python code
                    delta_val = (species.lam * species.A)*((Bratio -1.)/(1.+beta));
                    F1 = hyp2f1(1.0, kappa_perp+1.0, kappa_tot+1.5, 1.0 - 1./delta_val);
                    F2 = hyp2f1(2.0, kappa_perp+1.0, kappa_tot+1.5, 1.0 - 1./delta_val);
                    F3q= hyp2f1(0.75, kappa_perp+1.0, kappa_tot+1.5, 1.0 - 1./delta_val);
                    F7q= hyp2f1(1.75, kappa_perp+1.0, kappa_tot+1.5, 1.0 - 1./delta_val);

                    H1h= hyp2f1(1.0, kappa_perp+1.0, kappa_tot+0.5,  1.0 - 1./delta_val)/gamma(kappa_tot+0.5);
                    H3q= hyp2f1(1.0, kappa_perp+1.0, kappa_tot+0.75, 1.0 - 1./delta_val)/gamma(kappa_tot+0.75);
                    H3h= F1/gamma(kappa_tot+1.5);
                    H7q= hyp2f1(1.0, kappa_perp+1.0, kappa_tot+1.75, 1.0 - 1./delta_val)/gamma(kappa_tot+1.75);

                    C1 = 4.0*kappa_tot - 1.0;
                    C2 = 2.0*kappa_tot - 1.0;
                    C3 = 4.0*kappa_par - 1.0;
                    C4 = 2.0*kappa_par - 1.0;

                    
                    % num_kappa_perp = 2*C1*F3q*F2/(C2*F1*F7q) - 4
                    % kappa_perp_local = 1/num_kappa_perp +1
                    num_kappa_perp = 2.0*C1*F3q*F2/(C2*F1*F7q) - 4.0;
                    kp_perp = 1.0/num_kappa_perp + 1.0;

                    fac1 = (C3*C1)/(C4*C2);
                    fac2 = fac1*H1h*H7q/(H3q*H3h);
                    kp_par = (fac2 - 2.0) / (2.0*fac2 - 8.0);
                end
            else
                kp_par  = kappa_par;
                kp_perp = kappa_perp;
            end
        end

        function [kp_par, kp_perp] = fried_egg_kappa(species, deltaU, Phi, Bratio)
            % Fried-Egg => no parallel kappa (Maxwellian in parallel => infinite).
            if Bratio < 1.0
                fprintf('Warning: B/B0 < 1 in fried_egg_kappa => set Bratio=1+eps.\n');
                Bratio = 1.0 + eps;
            end
            kp_par = Inf;
            k0 = species.kappa;
            if Bratio ~= 1.0
                 
                kp_perp = k0;  % or expansions.
                % you had big expansions for large x, calling hyperu, etc.

                x = species.A*(Bratio - 1.0);
                z = k0*x;
                if k0 >= 50.0
                    kp_perp = k0;
                elseif x >= 15.0
                    % series expansion about x->infinity
                    term1 = (x^2)*(k0^2)/(1.0 + k0);
                    term2 = (289*(2.0 + k0)) / (16.0*x*k0*(1.0 + k0));
                    term3 = ( x*k0*(27.0 + 8.0*k0 ) ) / (4.0*(1.0 + k0));
                    kp_perp = term1 + term2 + term3;
                else
                    valU0 = hyperu(k0+1.0, k0, z);
                    valU1 = hyperu(k0+1.0, k0+1.0, z);
                    valU1q= hyperu(k0+1.0, k0+0.25, z);
                    valU5q= hyperu(k0+1.0, k0+1.25, z);

                    fac = (valU0*valU5q)/(valU1q*valU1);
                    kp_perp = (0.75 - fac)/(1.0 - fac);
                end
            else
                kp_perp = k0;
            end
        end

    end % private methods
end
