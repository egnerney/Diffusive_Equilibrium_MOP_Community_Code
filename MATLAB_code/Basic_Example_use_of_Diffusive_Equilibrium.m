% BASIC_EXAMPLE_USE_OF_DIFFUSIVE_EQUILIBRIUM_CODE.M
% ----------------------------------------------------------------------------
%
% Requirements:
%   - Io_aligned_field_line_trace_JRM33+con2020_dlat=0.1_degree.mat
%   - Nerney2025_reference_model_4-10RJ.mat
%   - Species.m, plasma.m
%
% This produces the same results and plots, but the species fields in the
% struct are now: 'Op', 'O2p', 'Sp', 'S2p', 'S3p', 'Hp', 'Nap', 'Oph', 'ehm', 'em'.
% ----------------------------------------------------------------------------

clear; clc; close all;

%% Fundamental constants
kg_per_u = 1.66053906892e-27;  % 1 u in kg
m_e = 5.485799090441e-4;       % electron mass in u

m_H  = 1.0078;  
m_O  = 15.999;
m_S  = 32.065;
m_Na = 22.990;

m_Hp  = m_H  - m_e;
m_Op  = m_O  - m_e;
m_O2p = m_O  - 2*m_e;
m_Sp  = m_S  - m_e;
m_S2p = m_S  - 2*m_e;
m_S3p = m_S  - 3*m_e;
m_Nap = m_Na - m_e;

ELEMENTARY_CHARGE = 1.602176634e-19;  
EV_TO_JOULE       = ELEMENTARY_CHARGE;
JOULE_TO_EV       = 1 / ELEMENTARY_CHARGE;


%% Planet, corotation fraction
planet = 'Jupiter';
fcor   = 1.0;

[RP, Omega, GM] = plasma.define_planet(planet);
Omega = Omega * fcor;

%% Load field-line trace data from MAT
traceData = load('Io_aligned_field_line_trace_JRM33+con2020_dlat=0.1_degree.mat');
s     = traceData.s;
x     = traceData.x;
y     = traceData.y;
z     = traceData.z;
rho   = traceData.rho;
r     = traceData.r;
lat   = traceData.lat;
wlong = traceData.wlong;
B     = traceData.B;

npointsfl = length(s);

%% Load reference model data from MAT
refData = load('Nerney2025_reference_model_4-10RJ.mat');
% refData.n0, refData.T0, refData.kappa0, refData.kappaT0 => size [601 x 10]
% refData.species => 10x1 cell array

iofl_idx = 191;
nrfls    = 601;  % matches shape of n0, T0, etc.

% species names to avoid +/- issues in names
species_names = {
    'Op',   ...  % was 'O+'
    'O2p',  ...  % was 'O++'
    'Sp',   ...  % was 'S+'
    'S2p',  ...  % was 'S++'
    'S3p',  ...  % was 'S+++'
    'Hp',   ...  % was 'H+'
    'Nap',  ...  % was 'Na+'
    'Oph',  ...  % was 'O+(hot)'
    'ehm',  ...  % was 'eh-'
    'em'    ...  % was 'e-'
};
nspec = numel(species_names);



% Prepare arrays
species_n    = zeros(1, nspec);
species_T    = zeros(1, nspec);
species_kval = zeros(1, nspec);
species_ktmp = zeros(1, nspec);

for i = 1:nspec
    species_n(i)    = refData.n0(iofl_idx, i);
    species_T(i)    = refData.T0(iofl_idx, i);
    species_kval(i) = refData.kappa0(iofl_idx, i);
    species_ktmp(i) = refData.kappaT0(iofl_idx, i);
end

%% Build mass & charge arrays in same order

species_m = [
    m_Op,  ...  % 'Op' 
    m_O2p, ...  % 'O2p'
    m_Sp,  ...  % 'Sp'
    m_S2p, ...  % 'S2p'
    m_S3p, ...  % 'S3p'
    m_Hp,  ...  % 'Hp'
    m_Nap, ...  % 'Nap'
    m_Op,  ...  % 'Oph' -> we can just reuse m_Op
    m_e,   ...  % 'ehm'
    m_e    ...  % 'em'
];
species_q = [
    1.0, 2.0, 1.0, 2.0, 3.0, ...
    1.0, 1.0, 1.0, -1.0, -1.0
];

% For all we set A=2, lam=1 (adjust as needed)
A0_array  = 2.0 * ones(1, nspec);
lam_array = 1.0 * ones(1, nspec);
%| **Swap distribution function** | each species can differ but in this test case they are all the same.
%Avaiable Distribution Types:
%dist_type = repmat({'Maxwellian'}, 1, nspec); % Bagenal & Sullivan (1981) etc. Basic default Assumed isotropic always, ignores anisotropy
dist_type = repmat({'Aniso_Maxwellian'}, 1, nspec); % Huang & Birmingham (1992) Anisotropic Maxwellian
%dist_type = repmat({'Aniso_Maxwellian_const_Tperp'}, 1, nspec); % Bagenal & Sullivan (1981), Bagenal (1994), Mei, Thorne, Bagenal (1995) Assumed Constant Anisotrpy along field line, inconsistent but good enough often 
%dist_type = repmat({'Aniso_kappa'}, 1, nspec); % Meyer-Vernet et al. (1995), Moncuquet et al. (2002) Anisotropic Standard Kappa distribution
%dist_type = repmat({'Aniso_product_kappa'}, 1, nspec); % Nerney (2025) Anisotropic Product Kappa Distribution 
%dist_type = repmat({'Fried_Egg'}, 1, nspec); % Nerney (2025) "Fried-Egg" Distribution Function (Maxwellian limit in parallel of aniso product kappa and still kappa in perp direction)


%% Create the 10 Species objects
species_0 = cell(1, nspec);
for i = 1:nspec
    sp = Species( ...
        species_names{i}, ...
        species_m(i)*kg_per_u, ...
        species_q(i)*ELEMENTARY_CHARGE, ...
        species_T(i)*ELEMENTARY_CHARGE, ...
        A0_array(i), ...
        species_kval(i), ...
        lam_array(i), ...
        species_n(i), ...
        dist_type{i} ...
    );
    species_0{i} = sp;
end

%% Prepare struct outputs
% Now the species_names are valid as struct fields, e.g. 'Op' not 'O+'
n          = struct();
T_par      = struct();
T_perp     = struct();
kappa_par  = struct();
kappa_perp = struct();

for i = 1:nspec
    n.(species_names{i})          = zeros(1, npointsfl);
    T_par.(species_names{i})      = zeros(1, npointsfl);
    T_perp.(species_names{i})     = zeros(1, npointsfl);
    kappa_par.(species_names{i})  = zeros(1, npointsfl);
    kappa_perp.(species_names{i}) = zeros(1, npointsfl);
end

%% Find "centrifugal equator" at max(rho)
[~, ceq_idx] = max(rho);

Bratio = B ./ B(ceq_idx);
Pos  = [r(:), lat(:), wlong(:)];
Pos0 = [r(ceq_idx), lat(ceq_idx), wlong(ceq_idx)];

deltaU = plasma.calc_deltaU(r(ceq_idx), lat(ceq_idx), r, lat, planet, fcor);
[max_idx, max_lats, max_vals] = plasma.find_local_maxima(deltaU, lat);

%% Plot deltaU vs lat
figure('Units','inches','Position',[1 1 10 6]);
plot(lat, deltaU, 'LineWidth',1.2); hold on;
for j = 1:length(max_lats)
    xline(max_lats(j), '--r', 'LineWidth',1);
end
xlabel('Latitude (^o)');
ylabel('\DeltaU (J/kg)');
title('Centrifugal + Gravitational Potential vs Latitude');
grid on; hold off;

%% Allocate arrays
n_out         = zeros(nspec, npointsfl);
Tpar_out      = zeros(nspec, npointsfl);
Tperp_out     = zeros(nspec, npointsfl);
kappapar_out  = zeros(nspec, npointsfl);
kappaperp_out = zeros(nspec, npointsfl);
phi           = zeros(1,npointsfl);

%% Loop over field line points
for ipt = 1:npointsfl
    if mod(ipt,100)==0
        disp(ipt);
    end
    
    [densVec, phiVal] = plasma.diff_eq(Pos(ipt,:), Pos0, Bratio(ipt), ...
                                       [species_0{:}], deltaU(ipt), planet, fcor);
    n_out(:, ipt) = densVec;
    phi(ipt)      = phiVal;

    [TparVec, TperpVec] = plasma.calc_temps([species_0{:}], deltaU(ipt), phiVal, Bratio(ipt));
    Tpar_out(:, ipt)  = TparVec;
    Tperp_out(:, ipt) = TperpVec;

    [kp_parVec, kp_perpVec] = plasma.calc_kappa_vals([species_0{:}], ...
                                       deltaU(ipt), phiVal, Bratio(ipt));
    kappapar_out(:, ipt)  = kp_parVec;
    kappaperp_out(:, ipt) = kp_perpVec;
end

%% Create pretty labels for plotting
species_labels = {
    'O^{+}', 'O^{++}', 'S^{+}', 'S^{++}', 'S^{+++}', ...
    'H^{+}', 'Na^{+}', 'O^{+}(hot)', 'eh^{-}', 'e^{-}'
};

%% 1) Plot densities vs s
figure; hold on;
for i = 1:nspec
    n.(species_names{i}) = n_out(i,:);
    plot(s, n.(species_names{i}), 'DisplayName', species_labels{i}, 'LineWidth',1.2);

    T_par.(species_names{i})  = Tpar_out(i,:) / EV_TO_JOULE;
    T_perp.(species_names{i}) = Tperp_out(i,:) / EV_TO_JOULE;

    kappa_par.(species_names{i})  = kappapar_out(i,:);
    kappa_perp.(species_names{i}) = kappaperp_out(i,:);
end
set(gca,'YScale','log');
title([strrep(dist_type{1}, '_', ' ') ' Diffusive Equilibrium, Aligned Io Field Line']);
xlabel('s (R_J)'); ylabel('n (cm^{-3})');
legend('Location','best','FontSize',8);
grid on; hold off;
saveas(gcf, [dist_type{1}, '_densities_vs_s.png']);

%% 2) Plot densities vs latitude
figure; hold on;
for i = 1:nspec
    plot(lat, n.(species_names{i}), 'DisplayName', species_labels{i}, 'LineWidth',1.2);
end
set(gca,'YScale','log');
title([strrep(dist_type{1}, '_', ' ') ' Diffusive Equilibrium, Aligned Io Field Line']);
xlabel('Latitude (^o)'); ylabel('n (cm^{-3})');
legend('Location','best','FontSize',8,'NumColumns',2);
grid on; hold off;
saveas(gcf, [dist_type{1}, '_densities_vs_lat.png']);

%% 3) Tparallel vs latitude
figure; hold on;
for i = 1:nspec
    plot(lat, T_par.(species_names{i}), 'DisplayName', species_labels{i}, 'LineWidth',1.2);
end
set(gca,'YScale','log');
title([strrep(dist_type{1}, '_', ' ') ' Diffusive Equilibrium, Aligned Io Field Line']);
xlabel('Latitude (^o)'); ylabel('T_{||} (eV)');
legend('Location','best','FontSize',8,'NumColumns',2);
grid on; hold off;
saveas(gcf, [dist_type{1}, '_Tpar_vs_lat.png']);

%% 4) Tperp vs latitude
figure; hold on;
for i = 1:nspec
    plot(lat, T_perp.(species_names{i}), 'DisplayName', species_labels{i}, 'LineWidth',1.2);
end
set(gca,'YScale','log');
title([strrep(dist_type{1}, '_', ' ') ' Diffusive Equilibrium, Aligned Io Field Line']);
xlabel('Latitude (^o)'); ylabel('T_{\perp} (eV)');
legend('Location','best','FontSize',8,'NumColumns',2);
grid on; hold off;
saveas(gcf, [dist_type{1}, '_Tperp_vs_lat.png']);

%% 5) kappa_par vs latitude
figure; hold on;
for i = 1:nspec
    plot(lat, kappa_par.(species_names{i}), 'DisplayName', species_labels{i}, 'LineWidth',1.2);
end
set(gca,'YScale','log');
title([strrep(dist_type{1}, '_', ' ') ' Diffusive Equilibrium, Aligned Io Field Line']);
xlabel('Latitude (^o)'); ylabel('\kappa_{||}');
legend('Location','best','FontSize',8,'NumColumns',2);
grid on; hold off;
saveas(gcf, [dist_type{1}, '_kappa_par_vs_lat.png']);

%% 6) kappa_perp vs latitude
figure; hold on;
for i = 1:nspec
    plot(lat, kappa_perp.(species_names{i}), 'DisplayName', species_labels{i}, 'LineWidth',1.2);
end
set(gca,'YScale','log');
title([strrep(dist_type{1}, '_', ' ') ' Diffusive Equilibrium, Aligned Io Field Line']);
xlabel('Latitude (^o)'); ylabel('\kappa_{\perp}');
legend('Location','best','FontSize',8,'NumColumns',2);
grid on; hold off;
saveas(gcf, [dist_type{1}, '_kappa_perp_vs_lat.png']);

disp('Done - all diffusive equilibrium plots generated.');


%% (Optional) Save final results to .mat
% outname = [dist_type{1}, '_DiffEq_IoAligned.mat'];
% save(outname, 'n','T_par','T_perp','kappa_par','kappa_perp','phi','deltaU','-v7');
