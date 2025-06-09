(* ::Package:: *)

(*
Highly Optimized Diffusive Equilibrium Package
Created on Thu Feb 20 17:10:10 2025
@author: Edward (Eddie) G. Nerney

Nerney (2025)

Diffusive Equilibrium Code for finding steady state plasma Densities
at extrapolated location s given plasma parameters at s along same 
magnetic field line trace. s is arc length along field line.
In assumed non inertial rotating frame with planet at \[CapitalOmega] in rad/s,
or with fraction of corotation fcor given
along fixed magnetic field lines given field line traces,
including gravitational potential energy (only important close to planet), 
Centrifugal potential energy, and electrostatic potential energy.

Assumes dynamic accessibility to 
all locations along field line in steady state,
plasma frozen to field line,
So No cross-field drift (beyond small gyration around B),
No significant scattering to transport particles across field lines,
equivalently for MHD Alfv\[EAcute]n's theorem frozen-in flux theorem or large magnetic reynolds number,
no collisions or sources or losses of particles (no loss cone), 
no drift or transport (though boundary condition can be informed by),
the dist. function same at reference location s0 as s,
energy conserved (no sources or losses of energy or dissipation),
magnetic moment conserved (adiabatic) or timescales all >>> gyroperiod,
Assumes neutrality spatial scales >>> debye length,
and assumes plasma frame is the same as the rotating frame.

T parallel or Tpar is given T parameter 
with T perpendicular or Tperp given by A*T.
or A = Tperp/Tpar (Anisotropy Factor)
For Kappa distributions Temps are characteristic temps or T_c not "real" temps. T
Where real temps are defined as as Tperp = pperp/n and Tpar = Tpar/n
For standard kappas the relation is
Tpar = Tpar_c/(1-3/(2*kappa)), Tperp = Tperp_c/(1-3/(2*kappa))  
For product kappas there can be different kappa_par and kappa_perp
Tpar = Tpar_c/(1-1/(2*kappa_par)), Tperp = Tperp_c/(1 - 1/(kappa_perp))  
kappa parameter given is assumed to be parallel kappa or kappa_par
kappa_par = kappa
then perp. kappa value is kappa_perp = lambda * kappa
lambda = kappa_perp/kappa_par
*)
(* ::Package:: *)

BeginPackage["DiffusiveEquilibriumOptimized`"]

(* Public symbols *)
Species::usage = "Species[name, m, q, T, A, kappa, lam, n, type] creates a species object.";
definePlanet::usage = "definePlanet[planet] returns planet-specific parameters.";
maxwellianDensity::usage = "Isotropic Maxwellian distribution density.";
anisoMaxwellianDensity::usage = "Anisotropic Maxwellian distribution density.";
anisoMaxwellianConstTperpDensity::usage = "Anisotropic Maxwellian with constant Tperp density.";
anisoKappaDensity::usage = "Standard anisotropic kappa distribution density.";
anisoProductKappaDensity::usage = "Product kappa distribution density.";
friedEggDensity::usage = "Fried-Egg distribution density.";
maxwellianTemperature::usage = "Temperature function for Maxwellian.";
anisoMaxwellianTemperature::usage = "Temperature function for anisotropic Maxwellian.";
anisoMaxwellianConstTperpTemperature::usage = "Temperature function for anisotropic Maxwellian with constant Tperp.";
anisoKappaTemperature::usage = "Temperature function for anisotropic kappa.";
anisoProductKappaTemperature::usage = "Temperature function for product kappa.";
friedEggTemperature::usage = "Temperature function for Fried-Egg.";
maxwellianKappa::usage = "Kappa function for Maxwellian.";
anisoMaxwellianKappa::usage = "Kappa function for anisotropic Maxwellian.";
anisoMaxwellianConstTperpKappa::usage = "Kappa function for anisotropic Maxwellian with constant Tperp.";
anisoKappaKappa::usage = "Kappa function for anisotropic kappa.";
anisoProductKappaKappa::usage = "Kappa function for product kappa.";
friedEggKappa::usage = "Kappa function for Fried-Egg.";
densityFunctions::usage = "Association of density functions.";
temperatureFunctions::usage = "Association of temperature functions.";
kappaFunctions::usage = "Association of kappa functions.";
diffEq::usage = "Compute densities and find phi from neutrality condition.";
findLocalMaxima::usage = "Find local maxima in data.";
calcDeltaU::usage = "Calculate deltaU.";
calcTemps::usage = "Calculate temperatures.";
calcKappaVals::usage = "Calculate kappa values.";
compiledDeltaU::usage = "Compiled version of deltaU calculation.";
compiledDensity::usage = "Compiled density calculations.";
compiledTemperature::usage = "Compiled temperature calculations.";
compiledKappa::usage = "Compiled kappa calculations.";
parallelDiffEq::usage = "Parallel computation of diffusive equilibrium.";

Begin["`Private`"]

(* CODATA RECOMMENDED VALUES OF THE FUNDAMENTAL PHYSICAL CONSTANTS: 2022, physics.nist.gov/constants *)
(* Unnecessary precision *)
kgPerU = 1.66053906892*10^-27; (* Atomic mass unit (u or daltons, defined as 1/12 carbon 12 mass) in kg/u *)
uPerKg = 1./kgPerU; (* u/kg *)
mE = 5.485799090441*10^-4; (* electron mass in u ~ 1.0/1822.9 u *)

(* Standard Periodic table values taking into account isotopic weighting and differences due to nuclear binding energy *)
mH = 1.0078; (* Hydrogen mass in u *)
mO = 15.999; (* Oxygen mass in u *)
mS = 32.065; (* Sulfur mass in u *)
mNa = 22.990; (* Sodium mass in u *)

(* Even more Unwarranted Unnecessary precision (ignoring electron binding energy decrease too) *)
mHp = mH - mE; (* Singly (Fully) Ionized Hydrogen or H^{+} or proton mass in u *)
mOp = mO - mE; (* Singly Ionized Oxygen or O^{+} mass in u *)
mO2p = mO - 2*mE; (* Doubly Ionized Oxygen or O^{++} mass in u *)
mSp = mS - mE; (* Singly Ionized Sulfur or S^{+} mass in u *)
mS2p = mS - 2*mE; (* Doubly Ionized Oxygen or S^{++} mass in u *)
mS3p = mS - 3*mE; (* Triply Ionized Oxygen or S^{+++} mass in u *)
mNap = mNa - mE; (* Singly Ionized Sodium or Na^{+} mass in u *)

ELEMENTARYCHARGE = 1.602176634*10^-19; (* Elementary charge in Coulombs *)
EVTOJOULE = ELEMENTARYCHARGE; (* Conversion from eV to Joules *)
JOULETOEV = 1./ELEMENTARYCHARGE;

(* Species data structure *)
(* Container for plasma species parameters *)
Species[name_, m_, q_, T_, A_, kappa_, lam_, n_, type_] := Association[
  "name" -> name,
  "m" -> m,
  "q" -> q,
  "T" -> T,
  "A" -> A,
  "kappa" -> kappa,
  "lam" -> lam,
  "n" -> n,
  "type" -> type
]

(* Planet definition function *)
(* Returns planet-specific parameters: radius (m), rotation rate (rad/s), and GM (m^3/s^2) *)
definePlanet[planet_String] := Module[{p = ToLowerCase[planet]},
  Switch[p,
    "earth", (* From NASA Earth Fact Sheet *)
    {3.3781*10^6, 7.29210*10^-5, 3.9860*10^14}, (* Equatorial Radius (1 bar level) in meters, Rotation rate in rad/s, G*M (m^3/s^2) *)
    "jupiter", (* From NASA Jupiter Fact Sheet *)
    {7.1492*10^7, 1.7585*10^-4, 1.26687*10^17}, (* Equatorial Radius (1 bar level) in meters, Rotation rate in rad/s, G*M (m^3/s^2) *)
    "saturn", (* From NASA Saturn Fact Sheet *)
    {6.0268*10^7, 1.6379*10^-4, 3.7931*10^16}, (* Equatorial Radius (1 bar level) in meters, Rotation rate in rad/s, G*M (m^3/s^2) *)
    "uranus", (* From NASA Uranus Fact Sheet *)
    {2.5559*10^7, 1.012*10^-4, 5.7940*10^15}, (* Equatorial Radius (1 bar level) in meters, Rotation rate in rad/s, G*M (m^3/s^2) *)
    "neptune", (* From NASA Neptune Fact Sheet *)
    {2.4764*10^7, 1.083*10^-4, 6.8351*10^15}, (* Equatorial Radius (1 bar level) in meters, Rotation rate in rad/s, G*M (m^3/s^2) *)
    _,
    (Print["Planet ", planet, " is not supported in definePlanet. Exiting."];
     Abort[])
  ]
]

(* Compiled versions of density functions for speed *)
(* These are used internally by the main density functions *)
compiledMaxwellianDensity = Compile[{{n0, _Real}, {m, _Real}, {deltaU, _Real}, {q, _Real}, {Phi, _Real}, {T, _Real}},
  n0 * Exp[-(m * deltaU + q * Phi)/T],
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed"
];

compiledAnisoMaxwellianDensity = Compile[{{n0, _Real}, {m, _Real}, {deltaU, _Real}, {q, _Real}, {Phi, _Real}, {T, _Real}, {A, _Real}, {Bratio, _Real}},
  Module[{B0OverB, denom, factorExp},
    B0OverB = 1.0/Bratio;
    denom = A + (1.0 - A) * B0OverB;
    If[denom <= 0, denom = 2.2251*10^-308]; (* machine epsilon *)
    factorExp = Exp[-(m * deltaU + q * Phi)/T];
    (n0/denom) * factorExp
  ],
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed"
];

compiledAnisoMaxwellianConstTperpDensity = Compile[{{n0, _Real}, {m, _Real}, {deltaU, _Real}, {q, _Real}, {Phi, _Real}, {T, _Real}, {A, _Real}, {Bratio, _Real}},
  Module[{factor, exponentTerm},
    factor = Bratio^(1.0 - A);
    exponentTerm = Exp[-(m * deltaU + q * Phi)/T];
    n0 * factor * exponentTerm
  ],
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed"
];

compiledAnisoKappaDensity = Compile[{{n0, _Real}, {m, _Real}, {deltaU, _Real}, {q, _Real}, {Phi, _Real}, {T, _Real}, {A, _Real}, {kappa, _Real}, {Bratio, _Real}},
  Module[{B0OverB, denom, PEratio, exponentFactor},
    B0OverB = 1.0/Bratio;
    denom = A + (1.0 - A) * B0OverB;
    If[denom <= 0, denom = 2.2251*10^-308];
    PEratio = (m * deltaU + q * Phi)/(kappa * T);
    If[PEratio <= -1., PEratio = -1.0 + 2.2251*10^-308];
    exponentFactor = (1.0 + PEratio)^(0.5 - kappa);
    (n0/denom) * exponentFactor
  ],
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed"
];

(* Density functions *)
(* Isotropic Maxwellian distribution: n(s) ~ exp[-(m*deltaU + q*Phi)/T ] *)
(* B-field ratio (Bratio) is not used in the isotropic Maxwellian formula, but is included for uniform interface *)
(* References: Bagenal & Sullivan (1981), eq. for isotropic Maxwellian *)
maxwellianDensity[species_, deltaU_, Phi_, Bratio_] := 
  compiledMaxwellianDensity[species["n"], species["m"], deltaU, species["q"], Phi, species["T"]]

(* Anisotropic Maxwellian distribution (standard form), where T_par is constant *)
(* along the field line, but T_perp(s) ~ 1 / [ A + (1 - A)*B0/B ] *)
(* n(s) = n0 / [ A0 + (1 - A0)*(B0/B) ] * exp[ - (m*deltaU + q*Phi)/T_par ] *)
(* References: Huang & Birmingham (1992) *)
anisoMaxwellianDensity[species_, deltaU_, Phi_, Bratio_] := 
  compiledAnisoMaxwellianDensity[species["n"], species["m"], deltaU, species["q"], Phi, species["T"], species["A"], Bratio]

(* Alternate anisotropic Maxwellian variant where T_perp is also assumed constant *)
(* (along with T_par). This leads to a different expression: *)
(* n(s) ~ (B(s)/B0)^(1 - A0) * exp[-(m*deltaU + q*Phi)/T_par ] *)
(* References: Bagenal & Sullivan (1981), Bagenal (1994), Mei, Thorne, Bagenal (1995) *)
anisoMaxwellianConstTperpDensity[species_, deltaU_, Phi_, Bratio_] := 
  compiledAnisoMaxwellianConstTperpDensity[species["n"], species["m"], deltaU, species["q"], Phi, species["T"], species["A"], Bratio]

(* Standard anisotropic kappa distribution with single kappa (same parallel & perp): *)
(* n(s) = [ n0 / ( A0 + (1-A0)(B0/B) ) ] * [ 1 + (m deltaU + q Phi) / (kappa T_par) ]^(0.5 - kappa) *)
(* References: Meyer-Vernet et al. (1995); Moncuquet et al. (2002) *)
anisoKappaDensity[species_, deltaU_, Phi_, Bratio_] := 
  compiledAnisoKappaDensity[species["n"], species["m"], deltaU, species["q"], Phi, species["T"], species["A"], species["kappa"], Bratio]

(* Product kappa distribution with possibly different kappa_par and kappa_perp: *)
(* f ~ [1 + (m v_par^2)/(2 kappa_par T_par)]^(-kappa_par - a) * *)
(*      [1 + (m v_perp^2)/(2 kappa_perp T_perp)]^(-kappa_perp - b) *)
(* The integral solution typically involves hypergeometric functions. *)
(* This function includes a rough check for B/B0 >= 1. If B < B0, we set Bratio = 1 + eps *)
(* as integral diverges for B/B0<1 *)
(* References: Nerney (2025) *)
anisoProductKappaDensity[species_, deltaU_, Phi_, Bratio_] := Module[
  {k0, kPerp0, kTot0, n0, TPar0, A0, l0, beta, x, delta, hf, factor, outer, nLocal, BratioUse},
  
  k0 = species["kappa"];
  kPerp0 = k0 * species["lam"];
  kTot0 = k0 + kPerp0;
  n0 = species["n"];
  TPar0 = species["T"];
  A0 = species["A"];
  l0 = species["lam"];
  
  (* For convergence of integral and positive density (physical) *)
  BratioUse = Bratio;
  If[BratioUse < 1.0,
    BratioUse = 1.0 + $MachineEpsilon
  ];
  
  (* For convergence of integral and positive density (physical) *)
  beta = (species["m"] * deltaU + species["q"] * Phi)/(k0 * TPar0);
  If[beta <= -1.,
    beta = -1.0 + $MachineEpsilon
  ];
  
  x = A0*(BratioUse - 1.0);
  
  If[BratioUse != 1.0,
    (* For smaller kappa (< ~17), attempt hypergeometric approach. Otherwise fallback. *)
    If[k0 < 17. && kPerp0 < 17.,
      delta = (l0 * A0) * ((BratioUse - 1.)/(1. + beta));
      (* Using 2F1 from scipy.special *)
      hf = Hypergeometric2F1[1., kPerp0 + 1.0, kTot0 + 1.5, 1. - 1./delta];
      
      factor = k0/(kTot0 + 0.5);
      outer = (n0/x) * ((1.0 + beta)^(0.5 - k0));
      nLocal = BratioUse * outer * factor * hf,
      
      (* Large kappa fallback: simpler expression ignoring the hypergeometric detail *)
      (* Use approximation for large kappa_perp +k_par *)
      (* k0*hf/(k0*(1+lambda0)+1/2) ~ 1 / (x + 1 + beta) *)
      (* and good to a couple % for kappa_perp and kappa_par >17 (better for larger kappa) *)
      nLocal = n0 * BratioUse * ((1. + beta)^(0.5 - k0)) / (x + 1. + beta)
    ],
    (* Else Bratio = 1 *)
    nLocal = n0*((1. + beta)^(0.5 - k0))
  ];
  
  nLocal
]

(* "Fried-Egg" distribution: Maxwellian in parallel (kappa_par -> inf), *)
(* kappa in the perpendicular direction. The velocity-space integral leads *)
(* to exponential integrals E_n(x). *)
(* n(s) ~ [some prefactor]* exp[-(m*deltaU + q*Phi)/T_par + ... ] * E_{kappa+1}(kappa * x) *)
(* with x = A0*(B/B0 - 1). *)
(* For physically consistent solutions in the derived formula, we require B >= B0. *)
(* References: Nerney (2025) *)
friedEggDensity[species_, deltaU_, Phi_, Bratio_] := Module[
  {k0, A0, x, z, exponent, expFactor, factor, funcValue, nLocal, BratioUse},
  
  BratioUse = Bratio;
  If[BratioUse < 1.0,
    Print["Warning: B/B0 < 1 for friedEggDensity => non-convergent integral. Setting Bratio=1."];
    BratioUse = 1.0 + $MachineEpsilon
  ];
  
  k0 = species["kappa"]; (* k0 = kappa_0 = kappa_perp0 as there is no kappa_par0 for Fried Egg *)
  A0 = species["A"]; (* A0 = Tperp0/Tpar0 *)
  x = A0 * (BratioUse - 1.0);
  z = k0*x;
  
  (* exponent = \Delta PE / Tpar0 *)
  exponent = -(species["m"] * deltaU + species["q"] * Phi)/species["T"];
  If[Abs[exponent] > 700,
    exponent = Sign[exponent]*700.
  ];
  expFactor = Exp[exponent];
  
  (* T assumed parallel T here *)
  (* Though we can use mpmath our approximations are highly accurate and faster *)
  (* In addition scipy only allows positive integer n for E_n(x) *)
  (* and scipy special incomplete upper gamma function only allows positive arguments so that respresntation won't work *)
  (* Scipy special allows for real valued a, b, and x for hyperu(a,b,x) and checks against mpmath at dps=100 show good aggremment for kappa0 and x<30 *)
  (* e^z * E_{k0+1}[z] = U[1, 1 - k0, z] *)
  
  If[BratioUse != 1.0,
    Which[
      (* Small x or Bratio expansion to second order in x *)
      (* good to 0.9%. Much better away from k0 = 1 *)
      (* but near k0 = 1 this expansion doesn't do well so keep to within x< 0.001 *)
      (* Prevents nans from small Bratio>1 and x>0 in confluent hypergeometric scipy special call *)
      x < 0.0001,
      factor = 1.;
      nLocal = species["n"] * BratioUse * factor * expFactor,
      
      (* if kappa=kappa_perp>30 and x>30. 0th order Approximation good to 0.1% that is same as Aniso-Maxwellian *)
      (* But easier and better approximation to just use first order approximation for k0>15. *)
      k0 >= 15.,
      (* First order in 1/k0 expansion *)
      (* k0 e^(k0*x)*E_{k0+1}[k0*x] ~1 / (1 + x) - x/(k0*(1+ x)^3) *)
      (* Good to 0.035% for k0>15. and all x values >0 => For Bratio>1 *)
      factor = 1./(1. + x) - x/(k0*((x + 1.)^3.));
      nLocal = species["n"] * BratioUse * factor * expFactor,
      
      (* Third order in 1/k0 expansion *)
      k0 < 15. && 15. < x < 30.,
      (* k0 e^(k0*x)*E_{k0+1}[k0*x] ~1 / (1 + x) - x/(k0*(1+ x)^3) + x(2x-1)/(k0^2*(1+ x)^5)- (6x^3 - 8x^2 +x)/(k0^3*(1+ x)^7) *)
      (* Didn't expand in x but 1/k0 but still works but for large x that x^7 gets huge so use below for x>30 *)
      (* if x> 15. then for any k0 (though refer to above approximation if k0>15 for any x) *)
      (* Good to 0.02% for k0=1.001, better for larger k0 and x *)
      factor = 1./(1. + x) - x/(k0*((x + 1.)^3.)) + 
               x*(2.*x - 1.)/((k0^2.)*((1. + x)^5.)) - 
               (6.*(x^3.) - 8.*(x^2.) + x)/((k0^3.)*((1. + x)^7.));
      nLocal = species["n"] * BratioUse * factor * expFactor,
      
      (* fourth order in 1/x expansion *)
      k0 < 15. && x >= 30.,
      (* For x>30 use this so that we don't have to calculate the x^7 above which gets large *)
      (* if x> 30. then for any k0 (though refer to above approximation if k0>15) *)
      (* Good to 0.01% *)
      Module[{numer, denom},
        numer = -6. + k0*(-11. + 2.*x - (6. + (-3. + x)*x)*k0 + 
                (-1. + x)*(1 + x^2.)*(k0^2.));
        denom = (x^4.)*(k0^3.);
        factor = numer/denom;
        nLocal = species["n"] * BratioUse * factor * expFactor
      ],
      
      (* General case using HypergeometricU *)
      True,
      funcValue = k0 * HypergeometricU[1., 1. - k0, z];
      nLocal = species["n"] * BratioUse * funcValue * expFactor
    ],
    (* If Bratio = 1, x= 0, and we are at reference location s0 *)
    (* but just in case \Delta PE does not equal 0 where Bratio = 1 *)
    (* then we use limiting form for non zero \delta PE same as Isotropic Maxwellian *)
    nLocal = species["n"]*expFactor
  ];
  
  nLocal
]

(* Temperature functions *)
(* Temperature function for isotropic Maxwellian distribution. *)
(* Temperature is constant along field line. *)
maxwellianTemperature[species_, deltaU_, Phi_, Bratio_] := {species["T"], species["T"]}

(* Temperature function for anisotropic Maxwellian distribution. *)
(* T_parallel is constant, T_perpendicular varies with B. *)
anisoMaxwellianTemperature[species_, deltaU_, Phi_, Bratio_] := Module[
  {TParLocal, B0OverB, denom, TPerpLocal},
  
  TParLocal = species["T"];
  B0OverB = 1.0/Bratio;
  denom = species["A"] + (1.0 - species["A"]) * B0OverB;
  
  If[denom <= 0,
    Print["Warning: Non-physical denom in anisoMaxwellianTemperature. Setting denom -> eps."];
    denom = $MachineEpsilon
  ];
  
  TPerpLocal = species["A"] * species["T"] / denom;
  {TParLocal, TPerpLocal}
]

(* Temperature function for anisotropic Maxwellian with constant T_perp. *)
(* Both T_parallel and T_perpendicular are constant. *)
anisoMaxwellianConstTperpTemperature[species_, deltaU_, Phi_, Bratio_] := 
  {species["T"], species["A"] * species["T"]}

(* Temperature function for standard anisotropic kappa distribution. *)
anisoKappaTemperature[species_, deltaU_, Phi_, Bratio_] := Module[
  {PEratio, beta, TParLocal, B0OverB, denom, TPerpLocal},
  
  PEratio = (species["m"] * deltaU + species["q"] * Phi)/(species["kappa"] * species["T"]);
  
  (* For physical consistency, clamp if PEratio <= -1 *)
  If[PEratio <= -1.,
    PEratio = -1.0 + $MachineEpsilon
  ];
  
  beta = PEratio;
  TParLocal = species["T"] * (1.0 + beta);
  
  B0OverB = 1.0/Bratio;
  denom = species["A"] + (1.0 - species["A"]) * B0OverB;
  
  If[denom <= 0,
    Print["Warning: Non-physical denom in anisoKappaTemperature. Setting denom -> eps."];
    denom = $MachineEpsilon
  ];
  
  TPerpLocal = species["A"] * species["T"] * (1.0 + beta) / denom;
  {TParLocal, TPerpLocal}
]

(* Temperature function for anisotropic product kappa distribution. *)
anisoProductKappaTemperature[species_, deltaU_, Phi_, Bratio_] := Module[
  {beta, kappaPar, kappaPerp, A0, TPar, TPerp, x, g, TParLocal, TPerpLocal,
   deltaVal, kappaTot, F1, F2, F3q, F7q, H1h, H3q, H3h, H7q, C1, C2, C3, C4,
   TPerpNum, Dperp1, Dperp2, TPerpDenom, TParNum, Dpar1, Dpar2, TParDenom, BratioUse},
  
  (* Check if B/B0 >= 1 *)
  BratioUse = Bratio;
  If[BratioUse < 1.0,
    Print["Warning: B/B0 < 1 in anisoProductKappaTemperature. Setting Bratio -> 1 + eps."];
    BratioUse = 1.0 + $MachineEpsilon
  ];
  
  beta = (species["m"] * deltaU + species["q"] * Phi)/(species["kappa"] * species["T"]);
  If[beta <= -1.,
    beta = -1.0 + $MachineEpsilon
  ];
  
  (* Reference Values at s0 *)
  kappaPar = species["kappa"];
  kappaPerp = kappaPar * species["lam"];
  A0 = species["A"];
  TPar = species["T"];
  TPerp = A0*TPar;
  x = A0*(BratioUse - 1.);
  g = (1.0 + beta);
  
  If[BratioUse != 1.0,
    (* For larger kappa values, use simpler approximations *)
    If[kappaPar >= 17. && kappaPerp >= 17.,
      (* Approximation for large kappa *)
      TParLocal = TPar * g;
      TPerpLocal = BratioUse * TPerp * g / (x + g),
      
      (* Full calculation with hypergeometric functions *)
      deltaVal = (species["lam"] * species["A"]) * ((BratioUse - 1.)/(1. + beta));
      kappaTot = kappaPar + kappaPerp;
      
      (* Calculate F_q values for Tperp *)
      (* F_q = 2_F_1(q,kappa_perp + 0.5, kappa_tot + 1.5, 1.0 - 1./delta_val) *)
      F1 = Hypergeometric2F1[1.0, kappaPerp + 1.0, kappaTot + 1.5, 1.0 - 1./deltaVal];
      F2 = Hypergeometric2F1[2.0, kappaPerp + 1.0, kappaTot + 1.5, 1.0 - 1./deltaVal];
      F3q = Hypergeometric2F1[0.75, kappaPerp + 1.0, kappaTot + 1.5, 1.0 - 1./deltaVal];
      F7q = Hypergeometric2F1[1.75, kappaPerp + 1.0, kappaTot + 1.5, 1.0 - 1./deltaVal];
      
      (* Calculate regularized 2_F_1 values for Tpar *)
      (* H_q = 2_F_1(1.0, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val) *)
      H1h = Hypergeometric2F1[1.0, kappaPerp + 1.0, kappaTot + 0.5, 1.0 - 1./deltaVal]/
            Gamma[kappaTot + 0.5];
      H3q = Hypergeometric2F1[1.0, kappaPerp + 1.0, kappaTot + 0.75, 1.0 - 1./deltaVal]/
            Gamma[kappaTot + 0.75];
      H3h = F1/Gamma[kappaTot + 1.5];
      H7q = Hypergeometric2F1[1.0, kappaPerp + 1.0, kappaTot + 1.75, 1.0 - 1./deltaVal]/
            Gamma[kappaTot + 1.75];
      
      (* Calculate constants *)
      C1 = 4.0 * kappaTot - 1.0;
      C2 = 2.0 * kappaTot - 1.0;
      C3 = 4.0 * kappaPar - 1.0;
      C4 = 2.0 * kappaPar - 1.0;
      
      (* Calculate temperatures *)
      TPerpNum = (beta + 1.0)*kappaPar*TPerp*(A0 + x)/(3. * A0 * x);
      Dperp1 = C1*F3q / (3.*F7q);
      Dperp2 = C2*F1 / (2.*F2);
      TPerpDenom = Dperp1 - Dperp2;
      
      (* Value of Perp Temperature at s *)
      TPerpLocal = TPerpNum / TPerpDenom;
      
      TParNum = 8.0*(beta + 1.0)*kappaPar*TPar;
      Dpar1 = C3*C1*H7q / H3q;
      Dpar2 = 2.*C4*C2*H3h / H1h;
      TParDenom = Dpar1 - Dpar2;
      
      (* Value of Parallel Temperature at s *)
      TParLocal = TParNum / TParDenom
    ],
    (* If B/B0 = 1 then we are at reference location s0 *)
    TPerpLocal = TPerp;
    TParLocal = TPar
  ];
  
  {TParLocal, TPerpLocal}
]

(* Temperature function for "Fried-Egg" distribution. *)
friedEggTemperature[species_, deltaU_, Phi_, Bratio_] := Module[
  {TParLocal, TPerp0, k0, x, z, factor, TPerpLocal, C1, C2, C3, C4,
   numerator, denominator, U5q, U1q, U1, U0, BratioUse},
  
  (* Check if B/B0 >= 1 *)
  BratioUse = Bratio;
  If[BratioUse < 1.0,
    Print["Warning: B/B0 < 1 in friedEggTemperature. Setting Bratio -> 1 + eps."];
    BratioUse = 1.0 + $MachineEpsilon
  ];
  
  (* T_parallel is constant for Fried-Egg *)
  TParLocal = species["T"];
  TPerp0 = species["T"] * species["A"];
  
  k0 = species["kappa"];
  (* For T_perpendicular, we use the derived formula with x and z *)
  x = species["A"] * (BratioUse - 1.0);
  z = k0 * x;
  
  If[BratioUse != 1.0,
    Which[
      (* 0th order approximation is just *)
      (* Tperp result is same as Aniso Maxwellian (Huang & Birmingham) *)
      (* factor = 1./(1. + x) as Bratio/(1. + x) = 1/(A0+(1-A0)B0/B) *)
      (* Good to 0.25 % for k0>100 for all x>0 *)
      (* Better for larger k0 and larger x *)
      k0 > 100.,
      factor = 1./(1. + x);
      TPerpLocal = BratioUse*TPerp0*factor,
      
      (* Good to 0.05% for all x>0 *)
      (* Better for larger kappa and larger x *)
      k0 >= 15.,
      factor = 1./(1. + x) - x/(((1. + x)^3.)*k0);
      TPerpLocal = BratioUse*TPerp0*factor,
      
      (* For k0> 1, Bratio>1, and x>7. this 5th order approximation *)
      (* Good to better than 0.22% for all all k0>1 *)
      (* better for larger x or k0 *)
      x >= 7. && k0 < 15.,
      C1 = (3. + k0)*(4. + k0*(7. + k0));
      C2 = -x*((1. + k0)^2.)*(4. + k0);
      C3 = -(x^3.)*(k0^2.)*(1. + k0) + (x^2.)*k0*(1. + k0)*(2. + k0);
      C4 = (x^4.)*(k0^3.) + C3 + C2 + C1;
      numerator = -4. + k0*C4;
      denominator = (x^5.)*(k0^4.);
      factor = numerator/denominator;
      TPerpLocal = BratioUse*TPerp0*factor,
      
      True,
      (* for x<7 and k0<15 *)
      U5q = HypergeometricU[k0 + 1., k0 + 5./4., z];
      U1q = HypergeometricU[k0 + 1., k0 + 1./4., z];
      U1 = HypergeometricU[k0 + 1., k0 + 1., z];
      U0 = HypergeometricU[k0 + 1., k0, z];
      denominator = 4.*(U5q/U1q) - 3.*(U1/U0);
      factor = (1./(x*denominator));
      TPerpLocal = BratioUse * TPerp0 * factor
    ],
    (* If B/B0 = 1 then we are at reference location s0 *)
    TPerpLocal = TPerp0
  ];
  
  {TParLocal, TPerpLocal}
]

(* Kappa functions *)
(* Kappa function for Maxwellian distribution. *)
(* Not applicable but included for interface consistency. *)
(* Though actually effectivly Infinite as Maxwellian is infinite kappa limit *)
maxwellianKappa[species_, deltaU_, Phi_, Bratio_] := {Null, Null}

(* Kappa function for anisotropic Maxwellian distribution. *)
(* Not applicable but included for interface consistency. *)
(* Though actually effectivly Infinite as Maxwellian is infinite kappa limit *)
anisoMaxwellianKappa[species_, deltaU_, Phi_, Bratio_] := {Null, Null}

(* Kappa function for anisotropic Maxwellian with constant T_perp. *)
(* Not applicable but included for interface consistency. *)
(* Though actually effectivly Infinite as Maxwellian is infinite kappa limit *)
anisoMaxwellianConstTperpKappa[species_, deltaU_, Phi_, Bratio_] := {Null, Null}

(* Kappa function for standard anisotropic kappa distribution. *)
(* Kappa parameter is constant along field line. *)
(* there is no par or perp in standard kappa only kappa *)
(* but include both for interface consistency. *)
anisoKappaKappa[species_, deltaU_, Phi_, Bratio_] := {species["kappa"], species["kappa"]}

(* Kappa function for anisotropic product kappa distribution. *)
anisoProductKappaKappa[species_, deltaU_, Phi_, Bratio_] := Module[
  {beta, kappaPar, kappaPerp, kappaParLocal, kappaPerpLocal, deltaVal, kappaTot,
   F1, F2, F3q, F7q, H1h, H3q, H3h, H7q, C1, C2, C3, C4, numKappaPerp, fac1, fac2, BratioUse},
  
  (* Check if B/B0 >= 1 *)
  BratioUse = Bratio;
  If[BratioUse < 1.0,
    Print["Warning: B/B0 < 1 in anisoProductKappaKappa. Setting Bratio -> 1 + eps."];
    BratioUse = 1.0 + $MachineEpsilon
  ];
  
  beta = (species["m"] * deltaU + species["q"] * Phi)/(species["kappa"] * species["T"]);
  If[beta <= -1.,
    beta = -1.0 + $MachineEpsilon
  ];
  
  (* kappa values at reference location s0 *)
  (* kappa_par0 *)
  (* kappa_perp0 *)
  kappaPar = species["kappa"];
  kappaPerp = kappaPar * species["lam"];
  
  If[BratioUse != 1.0,
    (* For larger kappa values, return constant kappa *)
    (* (approximation for large kappa matches standard kappa behavior) *)
    If[kappaPar >= 30. && kappaPerp >= 30.,
      kappaParLocal = kappaPar;
      kappaPerpLocal = kappaPerp,
      
      (* Full calculation with hypergeometric functions *)
      deltaVal = (species["lam"] * species["A"]) * ((BratioUse - 1.)/(1. + beta));
      kappaTot = kappaPar + kappaPerp;
      
      (* Calculate F_q values *)
      F1 = Hypergeometric2F1[1.0, kappaPerp + 1.0, kappaTot + 1.5, 1.0 - 1./deltaVal];
      F2 = Hypergeometric2F1[2.0, kappaPerp + 1.0, kappaTot + 1.5, 1.0 - 1./deltaVal];
      F3q = Hypergeometric2F1[0.75, kappaPerp + 1.0, kappaTot + 1.5, 1.0 - 1./deltaVal];
      F7q = Hypergeometric2F1[1.75, kappaPerp + 1.0, kappaTot + 1.5, 1.0 - 1./deltaVal];
      
      (* Calculate regularized values *)
      H1h = Hypergeometric2F1[1.0, kappaPerp + 1.0, kappaTot + 0.5, 1.0 - 1./deltaVal]/
            Gamma[kappaTot + 0.5];
      H3q = Hypergeometric2F1[1.0, kappaPerp + 1.0, kappaTot + 0.75, 1.0 - 1./deltaVal]/
            Gamma[kappaTot + 0.75];
      H3h = F1/Gamma[kappaTot + 1.5];
      H7q = Hypergeometric2F1[1.0, kappaPerp + 1.0, kappaTot + 1.75, 1.0 - 1./deltaVal]/
            Gamma[kappaTot + 1.75];
      
      (* Calculate constants *)
      C1 = 4.0 * kappaTot - 1.0;
      C2 = 2.0 * kappaTot - 1.0;
      C3 = 4.0 * kappaPar - 1.0;
      C4 = 2.0 * kappaPar - 1.0;
      
      numKappaPerp = 2.0 * C1 * F3q * F2 / (C2 * F1 * F7q) - 4.0;
      (* Calculate kappa_perp at s *)
      kappaPerpLocal = 1.0/numKappaPerp + 1.0;
      
      (* Calculate kappa_par at s *)
      fac1 = C3*C1/(C4*C2);
      fac2 = fac1*H1h*H7q/(H3q*H3h);
      kappaParLocal = (fac2 - 2.0)/(2*fac2 - 8.0)
    ],
    (* If Bratio = 1 then at reference location s0 *)
    kappaParLocal = kappaPar;
    kappaPerpLocal = kappaPerp
  ];
  
  {kappaParLocal, kappaPerpLocal}
]

(* Kappa function for "Fried-Egg" distribution. *)
(* No parallel kappa parameter (Maxwellian in parallel direction) *)
friedEggKappa[species_, deltaU_, Phi_, Bratio_] := Module[
  {kappaParLocal, k0, x, z, kappaPerpLocal, term1, term2, term3, C1, C2, C3, C4,
   numerator, denominator, U0, U1, U1q, U5q, fac, BratioUse},
  
  (* Check if B/B0 >= 1 *)
  BratioUse = Bratio;
  If[BratioUse < 1.0,
    Print["Warning: B/B0 < 1 in friedEggKappa. Setting Bratio -> 1 + eps."];
    BratioUse = 1.0 + $MachineEpsilon
  ];
  
  (* No parallel kappa parameter (Maxwellian in parallel direction so parallel kappa effectivly infinite) *)
  kappaParLocal = Null;
  
  (* k0 = kappa_perp0 *)
  k0 = species["kappa"];
  (* For kappa_perpendicular, use the derived formula with x and z *)
  x = species["A"] * (BratioUse - 1.0);
  z = k0 * x;
  
  If[BratioUse != 1.0,
    Which[
      k0 >= 50.0,
      kappaPerpLocal = k0,
      
      (* Series expansion about x = infinity to 2nd order *)
      (* Good to 0.4% for k0 >1 and x>15 *)
      x >= 15.0,
      term1 = ((x^2.)*(k0^2.))/(1. + k0);
      term2 = (289.*(2. + k0))/(16.*x*k0*(1. + k0));
      term3 = (x*k0*(27. + 8.*k0))/(4.*(1. + k0));
      kappaPerpLocal = term1 + term2 + term3,
      
      True,
      (* Via scipy special direct calculation using confluent hypergeometric function U(a,b,z) *)
      U0 = HypergeometricU[k0 + 1., k0, z];
      U1 = HypergeometricU[k0 + 1., k0 + 1., z];
      U1q = HypergeometricU[k0 + 1., k0 + 1./4., z];
      U5q = HypergeometricU[k0 + 1., k0 + 5./4., z];
      
      fac = U0*U5q/(U1q*U1);
      kappaPerpLocal = (0.75 - fac)/(1.0 - fac)
    ],
    (* If Bratio = 1 we are at reference location *)
    kappaPerpLocal = species["kappa"]
  ];
  
  {kappaParLocal, kappaPerpLocal}
]

(* Create associations for density, temperature, and kappa functions *)
densityFunctions = Association[
  "Maxwellian" -> maxwellianDensity,
  "Aniso_Maxwellian" -> anisoMaxwellianDensity,
  "Aniso_Maxwellian_const_Tperp" -> anisoMaxwellianConstTperpDensity,
  "Aniso_kappa" -> anisoKappaDensity,
  "Aniso_product_kappa" -> anisoProductKappaDensity,
  "Fried_Egg" -> friedEggDensity
];

temperatureFunctions = Association[
  "Maxwellian" -> maxwellianTemperature,
  "Aniso_Maxwellian" -> anisoMaxwellianTemperature,
  "Aniso_Maxwellian_const_Tperp" -> anisoMaxwellianConstTperpTemperature,
  "Aniso_kappa" -> anisoKappaTemperature,
  "Aniso_product_kappa" -> anisoProductKappaTemperature,
  "Fried_Egg" -> friedEggTemperature
];

kappaFunctions = Association[
  "Maxwellian" -> maxwellianKappa,
  "Aniso_Maxwellian" -> anisoMaxwellianKappa,
  "Aniso_Maxwellian_const_Tperp" -> anisoMaxwellianConstTperpKappa,
  "Aniso_kappa" -> anisoKappaKappa,
  "Aniso_product_kappa" -> anisoProductKappaKappa,
  "Fried_Egg" -> friedEggKappa
];

calcDeltaU[r0In_, lat0In_, rIn_, latIn_, planet_, fcor_] := Module[
  {RP, Omega, GM, r0, lat0, r, lat, U0, U, dU},
  
  (* Set planet-specific variables *)
  {RP, Omega, GM} = definePlanet[planet];
  Omega *= fcor;
  
  (* Calculate gravitational-centrifugal potentials per mass at s0 and s along field line *)
  (* Convert radius to meters *)
  r0 = r0In * RP;
  (* convert to radians *)
  lat0 = lat0In * Degree;
  (* Convert radius to meters *)
  r = rIn * RP;
  (* convert to radians *)
  lat = latIn * Degree;
  
  U0 = -GM/r0 - 0.5 * r0^2. * Cos[lat0]^2. * Omega^2.;
  U = -GM/r - 0.5 * r^2. * Cos[lat]^2. * Omega^2.;
  
  (* difference in gravitational + centrifdugal potential energy per mass *)
  (* along field line as a function of s *)
  dU = U - U0;
  dU
]

(* Compiled version for faster execution - properly vectorized *)
compiledDeltaU = Compile[{{r0, _Real}, {lat0, _Real}, {rVec, _Real, 1}, 
  {latVec, _Real, 1}, {RP, _Real}, {Omega, _Real}, {GM, _Real}},
  Module[{r0m, lat0rad, rm, latrad, U0, U},
    r0m = r0 * RP;
    lat0rad = lat0 * Degree;
    rm = rVec * RP;
    latrad = latVec * Degree;
    
    U0 = -GM/r0m - 0.5 * r0m^2. * Cos[lat0rad]^2. * Omega^2.;
    U = -GM/rm - 0.5 * rm^2. * Cos[latrad]^2. * Omega^2.;
    
    U - U0
  ],
  RuntimeAttributes -> {Listable},
  Parallelization -> True,
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed"
];

(* Find local maxima in deltaU as a function of lat *)
findLocalMaxima[deltaU_, lat_] := Module[
  {deltaUArray, latArray, sortIndices, deltaUSorted, maxSortedIndices, 
   maxIndices, maxLats, maxDeltaU},
  
  (* Ensure inputs are lists *)
  deltaUArray = Flatten[{deltaU}];
  latArray = Flatten[{lat}];
  
  (* Check if arrays have the same length *)
  If[Length[deltaUArray] != Length[latArray],
    Throw["deltaU and lat arrays must have the same length"]
  ];
  
  (* Sort by latitude (x-axis) *)
  sortIndices = Ordering[latArray];
  deltaUSorted = deltaUArray[[sortIndices]];
  
  (* Find local maxima *)
  maxSortedIndices = {};
  
  (* We need at least 3 points to find local maxima *)
  If[Length[deltaUSorted] < 3,
    Return[{{}, {}, {}}]
  ];
  
  (* Check each point (except first and last) *)
  Do[
    If[deltaUSorted[[i]] > deltaUSorted[[i-1]] && 
       deltaUSorted[[i]] > deltaUSorted[[i+1]],
      AppendTo[maxSortedIndices, i]
    ],
    {i, 2, Length[deltaUSorted] - 1}
  ];
  
  (* No maxima found *)
  If[Length[maxSortedIndices] == 0,
    Return[{{}, {}, {}}]
  ];
  
  (* Convert sorted indices to original indices *)
  maxIndices = sortIndices[[maxSortedIndices]];
  
  (* Get values at maxima *)
  maxLats = latArray[[maxIndices]];
  maxDeltaU = deltaUArray[[maxIndices]];
  
  {maxIndices, maxLats, maxDeltaU}
]

(* Compiled net charge density function for bisection - optimized for common cases *)
compiledNetChargeDensity = Compile[
  {{phi, _Real}, {deltaU, _Real}, {Bratio, _Real}, 
   {masses, _Real, 1}, {charges, _Real, 1}, {temps, _Real, 1}, 
   {ns, _Real, 1}, {As, _Real, 1}, {kappas, _Real, 1}},
  Module[{nLocal, totalCharge, i, nspec = Length[masses], expArg, PEratio, B0OverB, denom},
    totalCharge = 0.;
    B0OverB = 1.0/Bratio;
    
    (* Compute densities for each species - using simple compiled expressions *)
    For[i = 1, i <= nspec, i++,
      (* Use simplified expressions for speed in bisection *)
      (* This is for Maxwellian/Aniso_Maxwellian approximation *)
      denom = As[[i]] + (1.0 - As[[i]]) * B0OverB;
      If[denom <= 0, denom = 2.2251*10^-308];
      expArg = -(masses[[i]] * deltaU + charges[[i]] * phi)/temps[[i]];
      If[expArg < -700, expArg = -700];
      If[expArg > 700, expArg = 700];
      nLocal = (ns[[i]]/denom) * Exp[expArg];
      totalCharge += charges[[i]] * nLocal;
    ];
    totalCharge
  ],
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed"
];

(* Net charge density function for bisection *)
netChargeDensity[species0_, deltaU_, phi_, Bratio_] := Module[
  {nspec, nLocal, densityFunc, totalCharge, masses, charges, temps, ns, As, kappas, useCompiled},
  
  nspec = Length[species0];
  
  (* Check if we can use compiled version (for simple distributions) *)
  useCompiled = And @@ Table[
    MemberQ[{"Maxwellian", "Aniso_Maxwellian", "Aniso_Maxwellian_const_Tperp", "Aniso_kappa"}, species0[[i]]["type"]],
    {i, nspec}
  ];
  
  If[useCompiled,
    (* Use compiled version for speed *)
    masses = Table[species0[[i]]["m"], {i, nspec}];
    charges = Table[species0[[i]]["q"], {i, nspec}];
    temps = Table[species0[[i]]["T"], {i, nspec}];
    ns = Table[species0[[i]]["n"], {i, nspec}];
    As = Table[species0[[i]]["A"], {i, nspec}];
    kappas = Table[species0[[i]]["kappa"], {i, nspec}];
    
    totalCharge = compiledNetChargeDensity[phi, deltaU, Bratio, masses, charges, temps, ns, As, kappas],
    
    (* Fall back to general version for complex distributions *)
    nLocal = Table[0., {nspec}];
    Do[
      densityFunc = densityFunctions[species0[[i]]["type"]];
      nLocal[[i]] = densityFunc[species0[[i]], deltaU, phi, Bratio],
      {i, nspec}
    ];
    totalCharge = Sum[species0[[i]]["q"] * nLocal[[i]], {i, nspec}]
  ];
  
  totalCharge
]

(* Local relative error function for bisection *)
localRelativeError[species0_, deltaU_, thisPhi_, br_] := Module[
  {nQ, nspec, distFuncE, neVal, qe},
  
  nQ = netChargeDensity[species0, deltaU, thisPhi, br];
  nspec = Length[species0];
  
  (* cold electron density is species0[[nspec]] *)
  distFuncE = densityFunctions[species0[[nspec]]["type"]];
  neVal = distFuncE[species0[[nspec]], deltaU, thisPhi, br];
  
  (* e- charge *)
  qe = species0[[nspec]]["q"];
  
  (* If neVal ~0 => avoid divide by zero => return huge *)
  If[Abs[neVal] < 1.*10^-30,
    Return[1.*10^30]
  ];
  
  Abs[nQ/(qe * neVal)]
]

(* Main diffEq function using bisection method *)
(* Solve for local densities "n" and potential "phi" to ensure net charge=0 *)
diffEq[posIn_, pos0In_, bratio_, species0_, deltaU_, planet_, fcor_] := Module[
  {nspec, totalIonCharge, phi1, phi2, dphi, nq1, nq2, maxiters, iter,
   tolerance, it2, phiRoot, rel1, rel2, relMid, nqMid, n, phi, densityFunc,
   modSpecies, newElectronDensity},
  
  nspec = Length[species0];
  
  (* Force neutrality at reference => calculate cold-e density *)
  (* The last species must be cold e- => solve n_e = - (ion charge)/q_e *)
  totalIonCharge = Sum[species0[[i]]["q"] * species0[[i]]["n"], {i, nspec - 1}];
  newElectronDensity = -totalIonCharge / species0[[nspec]]["q"];
  
  (* Create modified species list with updated electron density *)
  modSpecies = Table[
    If[i == nspec,
      (* For electron, create new association with updated density *)
      Association[
        species0[[i]],
        "n" -> newElectronDensity
      ],
      (* For other species, keep as is *)
      species0[[i]]
    ],
    {i, nspec}
  ];
  
  (* Bracket the root in phi *)
  phi1 = 0.;
  phi2 = 0.;
  dphi = 1.;
  
  nq1 = netChargeDensity[modSpecies, deltaU, phi1, bratio];
  nq2 = netChargeDensity[modSpecies, deltaU, phi2, bratio];
  
  maxiters = 10000;
  iter = 0;
  
  (* Expand outward until we find nq1 & nq2 with opposite signs *)
  While[nq1*nq2 > 0 && iter < maxiters,
    phi1 -= dphi;
    phi2 += dphi;
    nq1 = netChargeDensity[modSpecies, deltaU, phi1, bratio];
    nq2 = netChargeDensity[modSpecies, deltaU, phi2, bratio];
    iter += 1
  ];
  
  If[iter >= maxiters,
    Print["Failed to bracket root in diffEq."];
    Return[Null]
  ];
  
  (* Bisection until net charge ~ 0 or maxiters *)
  tolerance = 1.*10^-8;
  it2 = 0;
  phiRoot = 0.5*(phi1 + phi2);
  
  rel1 = localRelativeError[modSpecies, deltaU, phi1, bratio];
  rel2 = localRelativeError[modSpecies, deltaU, phi2, bratio];
  relMid = localRelativeError[modSpecies, deltaU, phiRoot, bratio];
  
  While[relMid > tolerance && it2 < maxiters,
    (* Evaluate sign of netCharge at midpoint *)
    nqMid = netChargeDensity[modSpecies, deltaU, phiRoot, bratio];
    
    If[nq1 * nqMid < 0,
      phi2 = phiRoot;
      nq2 = nqMid,
      (* Else *)
      phi1 = phiRoot;
      nq1 = nqMid
    ];
    
    phiRoot = 0.5*(phi1 + phi2);
    relMid = localRelativeError[modSpecies, deltaU, phiRoot, bratio];
    it2 += 1
  ];
  
  (* Compute densities at phiRoot *)
  n = Table[0., {nspec}];
  Do[
    densityFunc = densityFunctions[modSpecies[[i]]["type"]];
    n[[i]] = densityFunc[modSpecies[[i]], deltaU, phiRoot, bratio],
    {i, nspec}
  ];
  
  phi = phiRoot;
  
  (* Return densities and phi *)
  {n, phi, deltaU}
]

calcTemps[species0_, deltaU_, phi_, Bratio_] := Module[
  {nspec, Tpar, Tperp, tempFunc, result},
  
  nspec = Length[species0];
  Tpar = Table[0., {nspec}];
  Tperp = Table[0., {nspec}];
  
  Do[
    tempFunc = temperatureFunctions[species0[[i]]["type"]];
    result = tempFunc[species0[[i]], deltaU, phi, Bratio];
    Tpar[[i]] = result[[1]];
    Tperp[[i]] = result[[2]],
    {i, nspec}
  ];
  
  {Tpar, Tperp}
]

calcKappaVals[species0_, deltaU_, phi_, Bratio_] := Module[
  {nspec, kappaPar, kappaPerp, kappaFunc, result},
  
  nspec = Length[species0];
  kappaPar = Table[0., {nspec}];
  kappaPerp = Table[0., {nspec}];
  
  Do[
    kappaFunc = kappaFunctions[species0[[i]]["type"]];
    result = kappaFunc[species0[[i]], deltaU, phi, Bratio];
    kappaPar[[i]] = result[[1]];
    kappaPerp[[i]] = result[[2]],
    {i, nspec}
  ];
  
  {kappaPar, kappaPerp}
]

(* Parallel computation function for field line calculations *)
parallelDiffEq[Pos_, Pos0_, Bratio_, species0_, deltaU_, planet_, fcor_, npointsfl_] := Module[
  {nspec = Length[species0], results, nOut, TparOut, TperpOut, kappaparOut, kappaperpOut, phi},
  
  (* Distribute work across kernels *)
  DistributeDefinitions[diffEq, calcTemps, calcKappaVals, densityFunctions, temperatureFunctions, kappaFunctions];
  
  (* Parallel computation with progress monitoring *)
  results = ParallelTable[
    If[Mod[i, 100] == 0, Print["Processing point ", i, " of ", npointsfl]];
    
    Module[{result, tempResult, kappaResult},
      result = diffEq[Pos[[i]], Pos0, Bratio[[i]], species0, deltaU[[i]], planet, fcor];
      
      If[result =!= Null,
        tempResult = calcTemps[species0, deltaU[[i]], result[[2]], Bratio[[i]]];
        kappaResult = calcKappaVals[species0, deltaU[[i]], result[[2]], Bratio[[i]]];
        {result[[1]], result[[2]], tempResult[[1]], tempResult[[2]], kappaResult[[1]], kappaResult[[2]]},
        {Table[0., {nspec}], Missing["NotAvailable"], Table[0., {nspec}], Table[0., {nspec}], Table[0., {nspec}], Table[0., {nspec}]}
      ]
    ],
    {i, 1, npointsfl},
    Method -> "FinestGrained"
  ];
  
  (* Unpack results *)
  nOut = Table[0., {nspec}, {npointsfl}];
  phi = Table[0., {npointsfl}];
  TparOut = Table[0., {nspec}, {npointsfl}];
  TperpOut = Table[0., {nspec}, {npointsfl}];
  kappaparOut = Table[0., {nspec}, {npointsfl}];
  kappaperpOut = Table[0., {nspec}, {npointsfl}];
  
  Do[
    nOut[[All, i]] = results[[i, 1]];
    phi[[i]] = results[[i, 2]];
    TparOut[[All, i]] = results[[i, 3]];
    TperpOut[[All, i]] = results[[i, 4]];
    kappaparOut[[All, i]] = results[[i, 5]];
    kappaperpOut[[All, i]] = results[[i, 6]],
    {i, npointsfl}
  ];
  
  {nOut, phi, TparOut, TperpOut, kappaparOut, kappaperpOut}
]

End[]

EndPackage[]
