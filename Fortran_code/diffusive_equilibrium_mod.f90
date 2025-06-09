!=======================================================================
!  diffusive_equilibrium_mod.f90
!Created May 2025
!
!@author: Edward (Eddie) G. Nerney
!
!Nerney (2025)
!
!Diffusive Equilibrium Code for finding steady state plasma Densities
!at extrapolated location s given plasma parameters at s along same 
!magnetic field line trace. s is arc length along field line.
!In assumed non inertial rotating frame with planet at \Omega in rad/s,
!or with fraction of corotation fcor given
!along fixed magnetic field lines given field line traces,
!including gravitational potential energy (only important close to planet), 
!Centrifugal potential energy, and electrostatic potenial energy.
!
!Assumes dynamic accesibility to 
!all locations along field line in steady state,
!plasma frozen to field line,
!So No cross‐field drift (beyond small gyration around B),
!No significant scattering to transport particles across field lines,
!equivalently for MHD Alfvén's theorem frozen-in flux theorem or large magnetic reynolds number,
!no collisions or sources or losses of particles (no loss cone), 
!no drift or transport (though boundary condition can be informed by),
!the dist. function same at reference location s0 as s,
!energy conserved (no sources or losses of energy or dissipation),
!magnetic moment conserved (adiabatic) or timescales all >>> gyroperiod,
!Assumes neutrality spatial scales >>> debye length,
!and assumes plasma frame is the same as the rotating frame.
!
!T parallel or Tpar is given T parameter 
!with T perpendicular or Tperp given by A*T.
!or A = Tperp/Tpar (Anistropy Factor)
!For Kappa distributions Temps are characterstic temps or T_c not "real" temps. T
!Where real temps are defined as as Tperp = pperp/n and Tpar = Tpar/n
!For standard kappas the relation is
!Tpar = Tpar_c/(1-3/(2*kappa)), Tperp = Tperp_c/(1-3/(2*kappa))  
!For product kappas there can be different kappa_par and kappa_perp
!Tpar = Tpar_c/(1-1/(2*kappa_par)), Tperp = Tperp_c/(1 - 1/(kappa_perp))  
!kappa parameter given is assumed to be parralel kappa or kappa_par
!kappa_par = kappa
!then perp. kappa value is kappa_perp = lambda * kappa
!lambda = kappa_perp/kappa_par
!=======================================================================
module diffusive_equilibrium_mod
   use, intrinsic :: iso_fortran_env, only : dp => real64               
   implicit none
   private

!---------------- public API -------------------------------------------
   public :: species_type, species_init, define_planet
   public :: maxwellian_density, aniso_maxwellian_density
   public :: aniso_maxwellian_const_tperp_density
   public :: aniso_kappa_density, aniso_product_kappa_density
   public :: fried_egg_density
   public :: maxwellian_temperature, aniso_maxwellian_temperature
   public :: aniso_maxwellian_const_tperp_temperature
   public :: aniso_kappa_temperature, aniso_product_kappa_temperature
   public :: fried_egg_temperature
   public :: maxwellian_kappa, aniso_maxwellian_kappa
   public :: aniso_maxwellian_const_tperp_kappa
   public :: aniso_kappa_kappa, aniso_product_kappa_kappa
   public :: fried_egg_kappa
   public :: calc_deltau, find_local_maxima
   public :: maxima_type
   public :: diff_eq, calc_temps, calc_kappa_vals
   public :: diffeq_result_type
   public :: get_temperature, get_kappa

!---------------- external special functions ---------------------------
   interface                             ! hypergeometric function used for 2F1(a,b,c,z)  is in hyper_2f1.f90
      !https://people.sc.fsu.edu/~jburkardt/f_src/hyper_2f1/hyper_2f1.html
      !https://people.sc.fsu.edu/~jburkardt/f_src/hyper_2f1/hyper_2f1.f90
      !N.J.~Higham, ``Accuracy and Stability of Numerical Algorithms'',
      ! SIAM, Philadelphia, 1996 for expm1 implementation.
      ! log1p follows instantly.
      !hyper_2f1, a Fortran90 code which evaluates the hypergeometric functions 2F1(a,b,c;x) for real or complex parameters a, b, c, and argument x, by N. Michel and M. Stoitsov.
      function hyper_2f1(a,b,c,z) result(res)
         import :: dp
         complex(dp)             :: res
         complex(dp), intent(in) :: a,b,c,z
      end function hyper_2f1
   end interface

   interface                             !Confluent Hyper. Geo. U Function chgu or U(a,b,x) used is in  (special_functions.f90 as chgu())
      !https://people.math.sc.edu/Burkardt/f_src/special_functions/special_functions.html
      !https://people.math.sc.edu/Burkardt/f_src/special_functions/special_functions.f90
      !Shanjie Zhang, Jianming Jin,
      !    Computation of Special Functions,
      !    Wiley, 1996,
      !    ISBN: 0-471-11963-6,
      !    LC: QA351.C45.
      ! 
      subroutine chgu(a,b,x,hu,md)
         import :: dp
         real(dp), intent(in)  :: a,b,x
         real(dp), intent(out) :: hu ! Output U(a,b,x)
         integer, intent(out)  :: md ! method used, integer MD, the method code. see chgu in special_functions.f90
      end subroutine chgu
   end interface
!-----------------------------------------------------------------------

   type :: species_type
   !Type for plasma species parameters.
   !Attributes:
   ! -----------
   ! name : str
   !     Species name (e.g., 'O+', 'e-', 'S++', etc.).
   ! m : float
   !     Mass in AMU (if prior to conversion) or in kg (if SI).
   ! q : float
   !     Charge in elementary charges (if prior to conversion) or in Coulombs (if SI).
   ! T : float
   !     Characteristic temperature in eV (if prior to conversion) or Joules (if SI).
   ! A : float
   !     Anisotropy factor A0 = T_perp / T_par at reference location, if relevant.
   ! kappa : float
   !     Kappa parameter. For Maxwellian, kappa -> infinity.
   ! lam : float
   !     λ = kappa_perp / kappa_par (only needed in product-kappa or more exotic models).
   ! n : float
   !     Number density at reference location (in cm^-3 or other consistent units).
   ! type : str
   !     Distribution type. Must match a key in density_functions dict (e.g., 'Maxwellian',
   !     'Aniso_Maxwellian', 'Aniso_kappa','Product_kappa', 'Fried_Egg').
      character(len=:), allocatable :: name
      real(dp)                      :: m, q, T, A, kappa, lam, n
      character(len=:), allocatable :: dist_type
   end type species_type

   type :: maxima_type
      integer,  allocatable :: idx(:)
      real(dp), allocatable :: lat(:)    ! latitude (deg)
      real(dp), allocatable :: val(:)    ! ΔU value
   end type maxima_type

   type :: diffeq_result_type
      real(dp), allocatable :: n(:)   ! densities for each species (m⁻³)
      real(dp)               :: phi   ! electrostatic potential (V)
      real(dp)               :: dU    ! ΔU supplied to solver (J kg⁻¹)
   end type diffeq_result_type
   

contains
!-----------------------------------------------------------------------
   pure function to_lower(str) result(out)
      character(len=*), intent(in) :: str
      character(len=len(str))      :: out
      integer :: i
      do i=1,len(str)
         select case(str(i:i))
         case('A':'Z'); out(i:i)=achar(iachar(str(i:i))+32)
         case default ; out(i:i)=str(i:i)
         end select
      end do
   end function to_lower
!-----------------------------------------------------------------------
   pure function species_init(name,m,q,T,A,kappa,lam,n,dist_type) result(sp)
   !Type for plasma species parameters.
   !Attributes:
   ! -----------
   ! name : str
   !     Species name (e.g., 'O+', 'e-', 'S++', etc.).
   ! m : float
   !     Mass in AMU (if prior to conversion) or in kg (if SI).
   ! q : float
   !     Charge in elementary charges (if prior to conversion) or in Coulombs (if SI).
   ! T : float
   !     Characteristic temperature in eV (if prior to conversion) or Joules (if SI).
   ! A : float
   !     Anisotropy factor A0 = T_perp / T_par at reference location, if relevant.
   ! kappa : float
   !     Kappa parameter. For Maxwellian, kappa -> infinity.
   ! lam : float
   !     λ = kappa_perp / kappa_par (only needed in product-kappa or more exotic models).
   ! n : float
   !     Number density at reference location (in cm^-3 or other consistent units).
   ! type : str
   !     Distribution type. Must match a key in density_functions dict (e.g., 'Maxwellian',
   !     'Aniso_Maxwellian', 'Aniso_kappa','Product_kappa', 'Fried_Egg').
     
      character(len=*), intent(in) :: name, dist_type
      real(dp), intent(in)         :: m,q,T,A,kappa,lam,n
      type(species_type)           :: sp
      sp%name=name; sp%m=m; sp%q=q; sp%T=T; sp%A=A
      sp%kappa=kappa; sp%lam=lam; sp%n=n; sp%dist_type=dist_type
   end function species_init
!-----------------------------------------------------------------------
   function define_planet(planet) result(p)
    !Returns planet-specific parameters: radius (m), rotation rate (rad/s),
    !and GM (m^3/s^2).
    !
    !Parameters
    !----------
    !planet : str
    !    Name of the planet
    !
    !Returns
    ! -------
    !RP : double precision real
    !    Planetary Equatorial radius (1 bar level) in meters.
    !Omega : double precision real
    !    Rotation rate in rad/s. 
    !GM : double precision real
    !    Gravitational parameter G*M_planet in m^3/s^2.
      character(len=*), intent(in) :: planet
      real(dp) :: p(3)
      select case(to_lower(trim(planet)))
      case('earth')
         p=[3.3781e6_dp,7.29210e-5_dp,3.9860e14_dp]
      case('jupiter')
         p=[7.1492e7_dp,1.7585e-4_dp,1.26687e17_dp]
      case('saturn')
         p=[6.0268e7_dp,1.6379e-4_dp,3.7931e16_dp]
      case('uranus')
         p=[2.5559e7_dp,1.0120e-4_dp,5.7940e15_dp]
      case('neptune')
         p=[2.4764e7_dp,1.0830e-4_dp,6.8351e15_dp]
      case default
         write(*,'(a)') 'Planet "'//trim(planet)//'" not supported in DEFINE_PLANET.'
         p=[-1.0_dp,-1.0_dp,-1.0_dp]
      end select
   end function define_planet

!====================  Stage-2 density helpers  ========================
!-----------------------------------------------------------------------
   pure real(dp) function safe_exp(x)
      real(dp), intent(in) :: x
      if     (x> 700.0_dp) then; safe_exp=exp(700.0_dp)
      elseif (x<-700.0_dp) then; safe_exp=exp(-700.0_dp)
      else                       ; safe_exp=exp(x)
      end if
   end function safe_exp
!-----------------------------------------------------------------------
   function maxwellian_density(sp,deltaU,phi,Bratio) result(n_local)
!!$     Isotropic Maxwellian distribution: n(s) ~ exp[-(m*deltaU + q*Phi)/T ].
!!$
!!$     B-field ratio (Bratio) is not used in the isotropic Maxwellian formula,
!!$     but is included for uniform interface.
!!$
!!$     References:
!!$     -----------
!!$     Bagenal & Sullivan (1981), eq. for isotropic Maxwellian
!!$
!!$     Parameters:
!!$     -----------
!!$     species : Species
!!$     Contains n, m, q, T (SI or consistent units).
!!$     deltaU : double precision real
!!$     Gravitational + centrifugal potential energy difference per unit mass.
!!$     Phi : double precision real
!!$     Electrostatic potential difference (voltage).
!!$     Bratio :  double precision real
!!$     B(s) / B0. (Unused here but included for interface consistency.)
!!$
!!$     Returns:
!!$     --------
!!$     n_local : double precision real
!!$     Number density at location s.
     type(species_type), intent(in) :: sp
      real(dp), intent(in) :: deltaU,phi,Bratio
      real(dp) :: n_local
      n_local = sp%n * safe_exp(-(sp%m*deltaU+sp%q*phi)/sp%T)
   end function maxwellian_density
!-----------------------------------------------------------------------
   function aniso_maxwellian_density(sp,deltaU,phi,Bratio) result(n_local)
!!$    Anisotropic Maxwellian distribution (standard form), where T_par is constant
!!$    along the field line, but T_perp(s) ~ 1 / [ A + (1 - A)*B0/B ].
!!$
!!$    n(s) = n0 / [ A0 + (1 - A0)*(B0/B) ] * exp[ - (m*deltaU + q*Phi)/T_par ]
!!$
!!$    References:
!!$    -----------
!!$    Huang & Birmingham (1992)
!!$
!!$    Parameters:
!!$    -----------
!!$    species : Species
!!$        Must have .A = T_perp0 / T_par0 at reference location.
!!$    deltaU : double precision real
!!$        Potential energy difference per unit mass (grav + cent).
!!$    Phi : double precision real
!!$        Electrostatic potential difference.
!!$    Bratio : double precision real
!!$        B(s)/B0.
!!$
!!$    Returns:
!!$    --------
!!$    n_local : double precision real
!!$        Computed number density at s.
      type(species_type), intent(in) :: sp
      real(dp), intent(in) :: deltaU,phi,Bratio
      real(dp) :: n_local,br,denom
      br=max(Bratio,1.0e-15_dp)
      denom = sp%A + (1.0_dp-sp%A)/br
      denom = max(denom,1.0e-15_dp)
      n_local = (sp%n/denom)*safe_exp(-(sp%m*deltaU+sp%q*phi)/sp%T)
   end function aniso_maxwellian_density
!-----------------------------------------------------------------------
   function aniso_maxwellian_const_tperp_density(sp,deltaU,phi,Bratio) result(n_local)
!!$    Alternate anisotropic Maxwellian variant where T_perp is also assumed constant
!!$    (along with T_par). This leads to a different expression:
!!$
!!$    n(s) ~ (B(s)/B0)^(1 - A0) * exp[-(m*deltaU + q*Phi)/T_par ].
!!$
!!$    References:
!!$    -----------
!!$    Bagenal & Sullivan (1981), Bagenal (1994), Mei, Thorne, Bagenal (1995)
!!$
!!$    Parameters:
!!$    -----------
!!$    species : Species
!!$        .A is T_perp/T_par (constant).
!!$    deltaU  : double precision real
!!$        Potential energy difference per unit mass.
!!$    Phi : double precision real
!!$        Electrostatic potential difference.
!!$    Bratio : double precision real
!!$        B(s)/B0.
!!$
!!$    Returns:
!!$    --------
!!$    n_local : double precision real
!!$        Number density at s.
      type(species_type), intent(in) :: sp
      real(dp), intent(in) :: deltaU,phi,Bratio
      real(dp) :: n_local
      n_local = sp%n * Bratio**(1.0_dp-sp%A) * safe_exp(-(sp%m*deltaU+sp%q*phi)/sp%T)
   end function aniso_maxwellian_const_tperp_density
!-----------------------------------------------------------------------
   function aniso_kappa_density(sp,deltaU,phi,Bratio) result(n_local)
!!$    Standard anisotropic kappa distribution with single kappa (same parallel & perp):
!!$      n(s) = [ n0 / ( A0 + (1-A0)(B0/B) ) ] * [ 1 + (m deltaU + q Phi) / (kappa T_par) ]^(0.5 - kappa)
!!$
!!$    References:
!!$    -----------
!!$    Meyer-Vernet et al. (1995); Moncuquet et al. (2002)
!!$
!!$    Parameters:
!!$    -----------
!!$    species : Species
!!$        Has .kappa, .A, .T, .m, .q
!!$    deltaU : double prec. real
!!$        Potential difference per unit mass (grav + cent).
!!$    Phi : double prec. real
!!$        Electrostatic potential difference.
!!$    Bratio : double prec. real
!!$        B(s)/B0.
!!$
!!$    Returns:
!!$    --------
!!$    n_local : double prec. real
!!$        Number density at s.
      type(species_type), intent(in) :: sp
      real(dp), intent(in) :: deltaU,phi,Bratio
      real(dp) :: n_local,br,denom,PE
      br = max(Bratio,1.0e-15_dp)
      denom = sp%A + (1.0_dp-sp%A)/br
      denom = max(denom,1.0e-15_dp)
      PE = (sp%m*deltaU+sp%q*phi)/(sp%kappa*sp%T)
      PE = max(PE,-1.0_dp+1.0e-15_dp)
      n_local = (sp%n/denom)*(1.0_dp+PE)**(0.5_dp-sp%kappa)
   end function aniso_kappa_density
!-----------------------------------------------------------------------
   function aniso_product_kappa_density(sp,deltaU,phi,Bratio) result(n_local)
!!$    Product kappa distribution with possibly different kappa_par and kappa_perp:
!!$       f ~ [1 + (m v_par^2)/(2 kappa_par T_par)]^(-kappa_par - a) *
!!$            [1 + (m v_perp^2)/(2 kappa_perp T_perp)]^(-kappa_perp - b)
!!$    The integral solution typically involves hypergeometric functions.
!!$
!!$    This function includes a rough check for B/B0 >= 1. If B < B0, we set Bratio = 1 + eps
!!$    as integral diverges for B/B0<1
!!$
!!$    References:
!!$    -----------
!!$    Nerney (2025)
!!$     
!!$      hypergeometric function used for 2F1(a,b,c,z)  is in hyper_2f1.f90
!!$      https://people.sc.fsu.edu/~jburkardt/f_src/hyper_2f1/hyper_2f1.html
!!$      https://people.sc.fsu.edu/~jburkardt/f_src/hyper_2f1/hyper_2f1.f90
!!$      N.J.~Higham, ``Accuracy and Stability of Numerical Algorithms'',
!!$      SIAM, Philadelphia, 1996 for expm1 implementation.
!!$      log1p follows instantly.
!!$      hyper_2f1, a Fortran90 code which evaluates the hypergeometric functions 2F1(a,b,c;x)
!!$      for real or complex parameters a, b, c, and argument x, by N. Michel and M. Stoitsov.
!!$
!!$    Parameters:
!!$    -----------
!!$    species : double prec. real
!!$        Has .kappa (interpreted as kappa_par), .lam ( = kappa_perp / kappa_par ), etc.
!!$    deltaU : double prec. real
!!$    Phi : double prec. real
!!$    Bratio : double prec. real
!!$        B(s)/B0.
!!$
!!$    Returns:
!!$    --------
!!$    n_local : double prec. real
!!$        Number density at s (via hypergeometric for small kappa; a simpler expression for large).
      type(species_type), intent(in) :: sp
      real(dp), intent(in) :: deltaU,phi,Bratio
      real(dp) :: n_local
      real(dp) :: br,k0,kperp0,ktot,n0,Tpar0,A0,x,beta_,delta_val,hf,factor,outer,expf,hypfac
      br=max(Bratio,1.0_dp)
      k0=sp%kappa; kperp0=k0*sp%lam; ktot=k0+kperp0
      n0=sp%n; Tpar0=sp%T; A0=sp%A
      beta_=(sp%m*deltaU+sp%q*phi)/(k0*Tpar0); beta_=max(beta_,-1.0_dp+1.0e-15_dp)
      x=A0*(br-1.0_dp)

      if(br==1.0_dp) then
         n_local=n0*(1.0_dp+beta_)**(0.5_dp-k0)
      else
         if(k0<17.0_dp .and. kperp0<17.0_dp) then
            delta_val=(sp%lam*A0)*(br-1.0_dp)/(1.0_dp+beta_)
            hf=real(hyper_2f1( cmplx(1.0_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                               cmplx(ktot+1.5_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) ))
            factor=k0/(ktot+0.5_dp)
            outer =(n0/x)*(1.0_dp+beta_)**(0.5_dp-k0)
            n_local=br*outer*factor*hf
         else
            expf=br*(1.0_dp+beta_)**(0.5_dp-k0)
            hypfac=1.0_dp/(x+1.0_dp+beta_)
            n_local=n0*expf*hypfac
         end if
      end if
   end function aniso_product_kappa_density
!-----------------------------------------------------------------------
   function fried_egg_density(sp,deltaU,phi,Bratio) result(n_local)
!!$    "Fried-Egg" distribution: Maxwellian in parallel (kappa_par -> inf),
!!$    kappa in the perpendicular direction. The velocity-space integral leads
!!$    to exponential integrals E_n(x).
!!$
!!$    n(s) ~ [some prefactor]* exp[-(m*deltaU + q*Phi)/T_par + ... ] * E_{kappa+1}(kappa * x)
!!$    with x = A0*(B/B0 - 1).
!!$
!!$    For physically consistent solutions in the derived formula, we
!!$    require B >= B0.
!!$
!!$    References:
!!$    -----------
!!$    Nerney (2025)
!!$
!!$    Tricommi Confluent Hyper. Function U(a,b,x) used is in  (special_functions.f90 as chgu())
!!$       https://people.math.sc.edu/Burkardt/f_src/special_functions/special_functions.html
!!$       https://people.math.sc.edu/Burkardt/f_src/special_functions/special_functions.f90
!!$       Shanjie Zhang, Jianming Jin,
!!$           Computation of Special Functions,
!!$           Wiley, 1996,
!!$           ISBN: 0-471-11963-6,
!!$           LC: QA351.C45. 
!!$
!!$    Parameters:
!!$    -----------
!!$    species : Species
!!$        Has .kappa (the perp kappa), .A = T_perp0/T_par0, .T = T_par in eV or J, etc.
!!$    deltaU : double prec. real
!!$        Potential difference per unit mass
!!$    Phi : double prec. real
!!$        Electrostatic potential difference
!!$    Bratio : double prec. real
!!$        B(s)/B0
!!$
!!$    Returns:
!!$    --------
!!$    n_local : double prec. real
!!$        Number density at s, using the "fried-egg" limit.
      type(species_type), intent(in) :: sp
      real(dp), intent(in) :: deltaU,phi,Bratio
      real(dp) :: n_local
      real(dp) :: br,k0,A0,x,z,expf,factor,fu
      integer  :: md

      br=max(Bratio,1.0_dp)
      k0=sp%kappa; A0=sp%A
      x=A0*(br-1.0_dp); z=k0*x
      expf=safe_exp(-(sp%m*deltaU+sp%q*phi)/sp%T)

      if(br==1.0_dp) then
         n_local=sp%n*expf; return
      end if

      if(x<1.0e-4_dp) then
         factor=1.0_dp
      elseif(k0>=15.0_dp) then
         factor=1.0_dp/(1.0_dp+x)-x/(k0*(1.0_dp+x)**3)
      elseif(x>=30.0_dp) then
         factor=(-6.0_dp+k0*(-11.0_dp+2.0_dp*x-(6.0_dp+(-3.0_dp+x)*x)*k0+(-1.0_dp+x)*(1.0_dp+x**2)*k0**2)) &
                 /(x**4*k0**3)
      else
         call chgu(1.0_dp,1.0_dp-k0,z,fu,md)
         factor=k0*fu
      end if
      n_local=sp%n*br*factor*expf
   end function fried_egg_density

!#######################################################################
   !                            TEMPERATURES
!#######################################################################

!-----------------------------------------------------------------------
!  4.1  isotropic Maxwellian  (T∥ = T⊥ = const)
!-----------------------------------------------------------------------
   pure function maxwellian_temperature(sp,deltaU,phi,Bratio) result(Tout)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: deltaU,phi,Bratio
      real(dp)                       :: Tout(2)
      Tout = [ sp%T, sp%T ]
   end function maxwellian_temperature

!-----------------------------------------------------------------------
!  4.2  anisotropic Maxwellian  (T∥ const, T⊥ varies with B)
!-----------------------------------------------------------------------
   pure function aniso_maxwellian_temperature(sp,deltaU,phi,Bratio) result(Tout)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: deltaU,phi,Bratio
      real(dp)                       :: Tout(2), br, denom
      br    = max(Bratio,1.0e-15_dp)
      denom = sp%A + (1.0_dp-sp%A)/br
      denom = max(denom,1.0e-15_dp)
      Tout  = [ sp%T, sp%A*sp%T/denom ]
   end function aniso_maxwellian_temperature

!-----------------------------------------------------------------------
!  4.3  anisotropic Maxwellian with constant T⊥
!-----------------------------------------------------------------------
   pure function aniso_maxwellian_const_tperp_temperature(sp,deltaU,phi,Bratio) result(Tout)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: deltaU,phi,Bratio
      real(dp)                       :: Tout(2)
      Tout = [ sp%T , sp%A*sp%T ]
   end function aniso_maxwellian_const_tperp_temperature

!-----------------------------------------------------------------------
!  4.4  standard anisotropic κ
!-----------------------------------------------------------------------
   pure function aniso_kappa_temperature(sp,deltaU,phi,Bratio) result(Tout)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: deltaU,phi,Bratio
      real(dp)                       :: Tout(2), br, beta_, denom
      beta_ = (sp%m*deltaU+sp%q*phi)/(sp%kappa*sp%T)
      beta_ = max(beta_,-1.0_dp+1.0e-15_dp)

      br = max(Bratio,1.0e-15_dp)
      denom = sp%A + (1.0_dp-sp%A)/br
      denom = max(denom,1.0e-15_dp)

      Tout(1) = sp%T*(1.0_dp+beta_)                         ! T∥
      Tout(2) = sp%A*sp%T*(1.0_dp+beta_)/denom              ! T⊥
   end function aniso_kappa_temperature

!-----------------------------------------------------------------------
   !  4.5  product-κ  (impure: needs ₂F₁) Nerney (2025)
!!$     
!!$      hypergeometric function used for 2F1(a,b,c,z)  is in hyper_2f1.f90
!!$      https://people.sc.fsu.edu/~jburkardt/f_src/hyper_2f1/hyper_2f1.html
!!$      https://people.sc.fsu.edu/~jburkardt/f_src/hyper_2f1/hyper_2f1.f90
!!$      N.J.~Higham, ``Accuracy and Stability of Numerical Algorithms'',
!!$      SIAM, Philadelphia, 1996 for expm1 implementation.
!!$      log1p follows instantly.
!!$      hyper_2f1, a Fortran90 code which evaluates the hypergeometric functions 2F1(a,b,c;x)
!!$      for real or complex parameters a, b, c, and argument x, by N. Michel and M. Stoitsov.
!!$
!-----------------------------------------------------------------------
   function aniso_product_kappa_temperature(sp,deltaU,phi,Bratio) result(Tout)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: deltaU,phi,Bratio
      real(dp)                       :: Tout(2)
      ! ---- locals mirroring IDL variables ----------------------------
      real(dp) :: br,beta_,kpar0,kperp0,ktot,A0,x,Tpar0,Tperp0,delta_val
      real(dp) :: F1,F2,F3q,F7q,H1h,H3q,H3h,H7q,C1,C2,C3,C4
      real(dp) :: numP,denP,numT,denT
      ! ----------------------------------------------------------------
      br = max(Bratio,1.0_dp)
      if (br == 1.0_dp) then
         Tout = [ sp%T, sp%A*sp%T ]
         return
      end if

      beta_  = (sp%m*deltaU+sp%q*phi)/(sp%kappa*sp%T)
      beta_  = max(beta_,-1.0_dp+1.0e-15_dp)
      kpar0  = sp%kappa
      kperp0 = kpar0*sp%lam
      ktot   = kpar0+kperp0
      A0     = sp%A
      x      = A0*(br-1.0_dp)
      Tpar0  = sp%T
      Tperp0 = A0*Tpar0

      ! --- large-κ shortcut ------------------------------------------
      if (kpar0>=17.0_dp .and. kperp0>=17.0_dp) then
         ! Approximation for large kappa
         Tout(1)= Tpar0*(1.0_dp + beta_)
         Tout(2)= br*Tperp0*(1.0_dp + beta_)/(x + 1.0_dp + beta_)
         return
      end if

      ! --- full calculation ------------------------------------------
      delta_val = (sp%lam*A0)*(br-1.0_dp)/(1.0_dp+beta_)
      F1   = real(hyper_2f1( cmplx(1.0_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                             cmplx(ktot+1.5_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) ))
      F2   = real(hyper_2f1( cmplx(2.0_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                             cmplx(ktot+1.5_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) ))
      F3q  = real(hyper_2f1( cmplx(0.75_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                             cmplx(ktot+1.5_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) ))
      F7q  = real(hyper_2f1( cmplx(1.75_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                             cmplx(ktot+1.5_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) ))

      H1h = real(hyper_2f1( cmplx(1.0_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                            cmplx(ktot+0.5_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) )) &
            / gamma(ktot+0.5_dp)
      H3q = real(hyper_2f1( cmplx(1.0_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                            cmplx(ktot+0.75_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) )) &
            / gamma(ktot+0.75_dp)
      H3h = F1 / gamma(ktot+1.5_dp)
      H7q = real(hyper_2f1( cmplx(1.0_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                            cmplx(ktot+1.75_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) )) &
            / gamma(ktot+1.75_dp)

      C1 = 4.0_dp*ktot - 1.0_dp
      C2 = 2.0_dp*ktot - 1.0_dp
      C3 = 4.0_dp*kpar0 - 1.0_dp
      C4 = 2.0_dp*kpar0 - 1.0_dp

      ! ---- T⊥ --------------------------------------------------------
      numP = (beta_+1.0_dp)*kpar0*Tperp0*(A0+x)/(3.0_dp*A0*x)
      denP = C1*F3q/(3.0_dp*F7q) - C2*F1/(2.0_dp*F2)
      Tout(2) = numP/denP

      ! ---- T∥ --------------------------------------------------------
      numT = 8.0_dp*(beta_+1.0_dp)*kpar0*Tpar0
      denT = C3*C1*H7q/H3q - 2.0_dp*C4*C2*H3h/H1h
      Tout(1) = numT/denT
   end function aniso_product_kappa_temperature

!-----------------------------------------------------------------------
   !  4.6  “Fried-Egg” distribution  (impure: needs U) (Nerney 2025)
!!$
!!$    Tricommi Confluent Hyper. Function U(a,b,x) used is in  (special_functions.f90 as chgu())
!!$       https://people.math.sc.edu/Burkardt/f_src/special_functions/special_functions.html
!!$       https://people.math.sc.edu/Burkardt/f_src/special_functions/special_functions.f90
!!$       Shanjie Zhang, Jianming Jin,
!!$           Computation of Special Functions,
!!$           Wiley, 1996,
!!$           ISBN: 0-471-11963-6,
!!$           LC: QA351.C45. 
!!$
!-----------------------------------------------------------------------
   function fried_egg_temperature(sp,deltaU,phi,Bratio) result(Tout)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: deltaU,phi,Bratio
      real(dp)                       :: Tout(2)
      real(dp):: br,k0,A0,x,z,Tperp0,factor,U0,U1,U1q,U5q,C1,C2,C3,C4,numer,denom
      integer :: md

      br = max(Bratio,1.0_dp)
      k0 = sp%kappa; A0=sp%A
      x  = A0*(br-1.0_dp); z=k0*x
      Tperp0 = A0*sp%T
      Tout(1)= sp%T                              ! T∥ constant in Fried-Egg

      if(br==1.0_dp) then
         ! If B/B0 = 1 then we are at reference location s0
         Tout(2)=Tperp0
         return
      end if
      

      if(x<1.0e-4_dp) then
         factor=1.0_dp
      elseif(k0>=100.0_dp) then
         ! 0th order approximation is just
         ! Tperp result is same as Aniso Maxwellian (Huang & Birmingham)
         ! factor = 1./(1. + x) as Bratio/(1. + x) = 1/(A0+(1-A0)B0/B)
         ! Good to 0.25 % for k0>100 for all x>0
         ! Better for larger k0 and larger x 
         factor = 1.0_dp/(1.0_dp + x)
      elseif(k0>=15.0_dp) then
         ! Good to 0.05% for all x>0 
         ! Better for larger kappa and larger x 
         factor = 1.0_dp/(1.0_dp + x) - x/(k0*(1.0_dp + x)**3)
      elseif(x>=7.0_dp) then
         ! For k0> 1, Bratio>1, and x>7. this 5th order approximation
         ! Good to better than 0.22% for all all k0>1 
         ! better for larger x or k0
         C1 = (3.0_dp + k0)*(4.0_dp + k0 * (7.0_dp + k0))
         C2 = -x*((1.0_dp + k0) ** 2)*(4.0_dp + k0)
         C3 = -( x**3) * (k0**2) *(1.0_dp + k0) + ( x**2) * k0 *(1.0_dp + k0) *(2.0_dp + k0)
         C4 = ( x**4.) * (k0**3.) + C3 + C2 + C1
         numer =  -4.0_dp + k0 *C4
         denom = ( x**5) * (k0**4)
         factor = numer/denom
      else
         !for x<7 and k0 <15
         call chgu(k0 + 1.0_dp,k0,z,U0,md)
         call chgu(k0 + 1.0_dp,k0 + 1.0_dp,z,U1,md)
         call chgu(k0 + 1.0_dp,k0 + 0.25_dp,z,U1q,md)
         call chgu(k0 + 1.0_dp,k0 + 1.25_dp,z,U5q,md)

         !denominator = 4.0_dp* (U5q/U1q) - 3.0_dp *(U1/U0)
         factor = (1.0_dp/(x*( 4.0_dp* (U5q/U1q) - 3.0_dp *(U1/U0))))
      end if
      Tout(2)= br*Tperp0*factor
   end function fried_egg_temperature

!#######################################################################
!                           KAPPA 
!#######################################################################

!-----------------------------------------------------------------------
!  5.1  Maxwellian  →  no κ (return NaNs)
!-----------------------------------------------------------------------
   pure function maxwellian_kappa(sp,deltaU,phi,Bratio) result(kout)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: deltaU,phi,Bratio
      real(dp)                       :: kout(2)
      kout = [ real( huge(1.0_dp), dp ), real( huge(1.0_dp), dp ) ]   ! ∞ → represent by huge()
   end function maxwellian_kappa

!-----------------------------------------------------------------------
!  5.2  anisotropic Maxwellian  (same)
!-----------------------------------------------------------------------
   pure function aniso_maxwellian_kappa(sp,deltaU,phi,Bratio) result(kout)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: deltaU,phi,Bratio
      real(dp)                       :: kout(2)
      kout = [ real( huge(1.0_dp), dp ), real( huge(1.0_dp), dp ) ]
   end function aniso_maxwellian_kappa

!-----------------------------------------------------------------------
!  5.3  anisotropic Maxwellian const-T⊥  (same)
!-----------------------------------------------------------------------
   pure function aniso_maxwellian_const_tperp_kappa(sp,deltaU,phi,Bratio) result(kout)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: deltaU,phi,Bratio
      real(dp)                       :: kout(2)
      kout = [ real( huge(1.0_dp), dp ), real( huge(1.0_dp), dp ) ]
   end function aniso_maxwellian_const_tperp_kappa

!-----------------------------------------------------------------------
!  5.4  standard κ  (κ constant)
!-----------------------------------------------------------------------
   pure function aniso_kappa_kappa(sp,deltaU,phi,Bratio) result(kout)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: deltaU,phi,Bratio
      real(dp)                       :: kout(2)
      kout = [ sp%kappa, sp%kappa ]
   end function aniso_kappa_kappa

!-----------------------------------------------------------------------
   !  5.5  product-κ   (impure – uses ₂F₁) Nerney (2025)
!!$     
!!$      hypergeometric function used for 2F1(a,b,c,z)  is in hyper_2f1.f90
!!$      https://people.sc.fsu.edu/~jburkardt/f_src/hyper_2f1/hyper_2f1.html
!!$      https://people.sc.fsu.edu/~jburkardt/f_src/hyper_2f1/hyper_2f1.f90
!!$      N.J.~Higham, ``Accuracy and Stability of Numerical Algorithms'',
!!$      SIAM, Philadelphia, 1996 for expm1 implementation.
!!$      log1p follows instantly.
!!$      hyper_2f1, a Fortran90 code which evaluates the hypergeometric functions 2F1(a,b,c;x)
!!$      for real or complex parameters a, b, c, and argument x, by N. Michel and M. Stoitsov.
!!$
!-----------------------------------------------------------------------
   function aniso_product_kappa_kappa(sp,deltaU,phi,Bratio) result(kout)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: deltaU,phi,Bratio
      real(dp)                       :: kout(2)
      real(dp):: br,kpar0,kperp0,k_tot,A0,x,beta_,delta_val
      real(dp):: F1,F2,F3q,F7q,H1h,H3q,H3h,H7q,C1,C2,C3,C4,fac1,fac2
      br = max(Bratio,1.0_dp)
      kpar0 = sp%kappa
      kperp0 = kpar0*sp%lam
      if (br==1.0_dp) then
         ! If Bratio = 1 then at reference location s0
         kout=[kpar0,kperp0]; return
      end if

      beta_ = (sp%m*deltaU+sp%q*phi)/(sp%kappa*sp%T)
      beta_ = max(beta_,-1.0_dp+1.0e-15_dp)
      if (kpar0>=30.0_dp .and. kperp0>=30.0_dp) then
         ! For larger kappa values, return constant kappa 
         ! (approximation for large kappa matches standard kappa behavior)
         kout=[kpar0,kperp0]; return
      end if
      ! Full calculation with hypergeometric functions

      A0=sp%A; x=A0*(br-1.0_dp); k_tot=kpar0+kperp0
      delta_val=(sp%lam*A0)*(br-1.0_dp)/(1.0_dp+beta_)

      F1  = real(hyper_2f1( cmplx(1.0_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                            cmplx(k_tot+1.5_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) ))
      F2  = real(hyper_2f1( cmplx(2.0_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                            cmplx(k_tot+1.5_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) ))
      F3q = real(hyper_2f1( cmplx(0.75_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                            cmplx(k_tot+1.5_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) ))
      F7q = real(hyper_2f1( cmplx(1.75_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                            cmplx(k_tot+1.5_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) ))

      H1h = real(hyper_2f1( cmplx(1.0_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                            cmplx(k_tot+0.5_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) )) &
            / gamma(k_tot+0.5_dp)
      H3q = real(hyper_2f1( cmplx(1.0_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                            cmplx(k_tot+0.75_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) )) &
            / gamma(k_tot+0.75_dp)
      H3h = F1 / gamma(k_tot+1.5_dp)
      H7q = real(hyper_2f1( cmplx(1.0_dp,0.0_dp,dp), cmplx(kperp0+1.0_dp,0.0_dp,dp), &
                            cmplx(k_tot+1.75_dp,0.0_dp,dp), cmplx(1.0_dp-1.0_dp/delta_val,0.0_dp,dp) )) &
            / gamma(k_tot+1.75_dp)

      C1=4.0_dp*k_tot-1.0_dp ;  C2=2.0_dp*k_tot-1.0_dp
      C3=4.0_dp*kpar0-1.0_dp ;  C4=2.0_dp*kpar0-1.0_dp

      kout(2) = 1.0_dp/( 2.0_dp*C1*F3q*F2/(C2*F1*F7q) - 4.0_dp ) + 1.0_dp  ! κ⊥
      fac1    = C3*C1/(C4*C2)
      fac2    = fac1*H1h*H7q/(H3q*H3h)
      kout(1) = (fac2-2.0_dp)/(2.0_dp*fac2-8.0_dp)                          ! κ∥
   end function aniso_product_kappa_kappa

!-----------------------------------------------------------------------
   !  5.6  Fried-Egg distribution (impure – uses U) (Nerney 2025)
!!$
!!$    Tricommi Confluent Hyper. Function U(a,b,x) used is in  (special_functions.f90 as chgu())
!!$       https://people.math.sc.edu/Burkardt/f_src/special_functions/special_functions.html
!!$       https://people.math.sc.edu/Burkardt/f_src/special_functions/special_functions.f90
!!$       Shanjie Zhang, Jianming Jin,
!!$           Computation of Special Functions,
!!$           Wiley, 1996,
!!$           ISBN: 0-471-11963-6,
!!$           LC: QA351.C45. 
!!$
!-----------------------------------------------------------------------
   function fried_egg_kappa(sp,deltaU,phi,Bratio) result(kout)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: deltaU,phi,Bratio
      real(dp)                       :: kout(2)
      real(dp):: br,k0,A0,x,z,U0,U1,U1q,U5q,fac
      integer :: md
      br = max(Bratio,1.0_dp)
      k0 = sp%kappa
      if (br==1.0_dp) then
         ! If Bratio = 1 we are at reference location
         kout=[ huge(1.0_dp), k0 ]   ! κ∥ → ∞, κ⊥ constant
         return
      end if
      A0=sp%A; x=A0*(br-1.0_dp); z=k0*x

      if(x<1.0e-4_dp) then
         kout=[ huge(1.0_dp), k0 ]
      elseif (k0>=50.0_dp) then
         kout=[ huge(1.0_dp), k0 ]
      elseif (x>=15.0_dp) then
         !x>=15.
         ! Series expansion about x = infinity to 2nd order
         ! Good to 0.4% for k0 >1 and x>15
         kout=[ huge(1.0_dp), &
                ( (x**2*k0**2)/(1.0_dp+k0) + (x*k0*(27.0_dp+8.0_dp*k0))/(4.0_dp*(1.0_dp+k0)) + &
                  (289.0_dp*(2.0_dp+k0))/(16.0_dp*x*k0*(1.0_dp+k0)) ) ]
      else
         call chgu(k0+1.0_dp,k0,z,U0,md)
         call chgu(k0+1.0_dp,k0 + 1.0_dp,z,U1,md)
         call chgu(k0+1.0_dp,k0 + 0.25_dp,z,U1q,md)
         call chgu(k0+1.0_dp,k0 + 1.25_dp,z,U5q,md)
         fac=U0*U5q/(U1q*U1)
         kout=[ huge(1.0_dp), (0.75_dp-fac)/(1.0_dp-fac) ]
      end if
   end function fried_egg_kappa

!#######################################################################
!                            Other
!#######################################################################

!-----------------------------------------------------------------------
! 6)  ΔU  (gravitational + centrifugal)  — pure
!-----------------------------------------------------------------------
   function calc_deltau(r0_in,lat0_in,r_in,lat_in,planet,fcor) result(dU)
      real(dp), intent(in) :: r0_in, lat0_in, r_in, lat_in  ! planet-radii & deg
      character(len=*), intent(in) :: planet
      real(dp), intent(in) :: fcor                           ! corotation fraction
      real(dp) :: dU                                         ! J kg⁻¹

      real(dp) :: p(3), rp, omega0, omega, gm
      real(dp) :: r0, r, lat0, lat
      real(dp), parameter :: deg2rad = acos(-1.0_dp)/180.0_dp

      p = define_planet(planet)
      if (any(p < 0.0_dp)) then
         dU = 0.0_dp
         return
      end if

      rp      = p(1);  omega0 = p(2);  gm = p(3)
      omega   = omega0 * fcor

      r0   =  r0_in * rp
      r    =     r_in * rp
      lat0 = lat0_in * deg2rad
      lat  =  lat_in * deg2rad

      dU = ( -gm/r  - 0.5_dp*r **2 *cos(lat )**2 *omega**2 ) - &
           ( -gm/r0 - 0.5_dp*r0**2*cos(lat0)**2*omega**2 )
    end function calc_deltau

!-----------------------------------------------------------------------
!  sort_perm  – returns permutation that sorts array x(:) ascending
!               (simple insertion sort; O(n²) but fine for n≲10⁴)
!-----------------------------------------------------------------------
pure subroutine sort_perm(x, perm)
   use, intrinsic :: iso_fortran_env, only : dp => real64
   real(dp), intent(in)  :: x(:)
   integer , intent(out) :: perm(size(x))

   integer :: i, j, tmp

   perm = [(i, i=1, size(x))]        ! identity permutation

   do i = 2, size(x)
      tmp = perm(i)
      j   = i - 1
      do while (j >= 1 .and. x(tmp) < x(perm(j)))
         perm(j+1) = perm(j)
         j = j - 1
      end do
      perm(j+1) = tmp
   end do
end subroutine sort_perm


!=======================================================================
!  find_local_maxima – identify all local peaks in ΔU versus latitude
!                      (output already sorted by latitude ascending)
!=======================================================================
   function find_local_maxima(du, lat) result(peaks)
   use, intrinsic :: iso_fortran_env, only : dp => real64
   real(dp), intent(in) :: du(:), lat(:)
   type(maxima_type)    :: peaks                ! result structure

   integer                :: n, i
   integer, allocatable   :: perm(:), peakPos(:)
   real(dp), allocatable  :: duS(:), latS(:)

   !---------------- sanity checks -------------------------------------
   n = size(du)
   if (n /= size(lat) .or. n < 3) then
      allocate(peaks%idx(0), peaks%lat(0), peaks%val(0))
      return
   end if

   !---------------- sort by latitude ----------------------------------
   allocate(perm(n));  call sort_perm(lat, perm)
   allocate(duS(n), latS(n))
   duS  = du (perm)
   latS = lat(perm)

   !---------------- locate peaks --------------------------------------
   allocate(peakPos(0))
   do i = 2, n-1
      if (duS(i) > duS(i-1) .and. duS(i) > duS(i+1)) peakPos = [peakPos, i]
   end do

   if (size(peakPos) == 0) then             ! no maxima
      allocate(peaks%idx(0), peaks%lat(0), peaks%val(0))
      return
   end if

   !---------------- fill result (already latitude-sorted) -------------
   allocate(peaks%idx(size(peakPos)), peaks%lat(size(peakPos)), peaks%val(size(peakPos)))
   peaks%idx = perm(peakPos)     ! original indices
   peaks%lat = latS(peakPos)
   peaks%val = duS (peakPos)
 end function find_local_maxima


!#######################################################################
!                           Solver
!#######################################################################

!--------------------- generic dispatch helpers ------------------------
   function get_density(sp, dU, phi, br) result(nloc)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: dU, phi, br
      real(dp)                       :: nloc
      select case (to_lower(sp%dist_type))
      case ('maxwellian')
         nloc = maxwellian_density(sp,dU,phi,br)
      case ('aniso_maxwellian')
         nloc = aniso_maxwellian_density(sp,dU,phi,br)
      case ('aniso_maxwellian_const_tperp')
         nloc = aniso_maxwellian_const_tperp_density(sp,dU,phi,br)
      case ('aniso_kappa')
         nloc = aniso_kappa_density(sp,dU,phi,br)
      case ('aniso_product_kappa')
         nloc = aniso_product_kappa_density(sp,dU,phi,br)
      case ('fried_egg')
         nloc = fried_egg_density(sp,dU,phi,br)
      case default
         nloc = 0.0_dp
      end select
   end function get_density

   function get_temperature(sp, dU, phi, br) result(Tpair)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: dU, phi, br
      real(dp)                       :: Tpair(2)
      select case (to_lower(sp%dist_type))
      case ('maxwellian')
         Tpair = maxwellian_temperature(sp,dU,phi,br)
      case ('aniso_maxwellian')
         Tpair = aniso_maxwellian_temperature(sp,dU,phi,br)
      case ('aniso_maxwellian_const_tperp')
         Tpair = aniso_maxwellian_const_tperp_temperature(sp,dU,phi,br)
      case ('aniso_kappa')
         Tpair = aniso_kappa_temperature(sp,dU,phi,br)
      case ('aniso_product_kappa')
         Tpair = aniso_product_kappa_temperature(sp,dU,phi,br)
      case ('fried_egg')
         Tpair = fried_egg_temperature(sp,dU,phi,br)
      end select
   end function get_temperature

   function get_kappa(sp, dU, phi, br) result(kpair)
      type(species_type), intent(in) :: sp
      real(dp), intent(in)           :: dU, phi, br
      real(dp)                       :: kpair(2)
      select case (to_lower(sp%dist_type))
      case ('maxwellian')
         kpair = maxwellian_kappa(sp,dU,phi,br)
      case ('aniso_maxwellian')
         kpair = aniso_maxwellian_kappa(sp,dU,phi,br)
      case ('aniso_maxwellian_const_tperp')
         kpair = aniso_maxwellian_const_tperp_kappa(sp,dU,phi,br)
      case ('aniso_kappa')
         kpair = aniso_kappa_kappa(sp,dU,phi,br)
      case ('aniso_product_kappa')
         kpair = aniso_product_kappa_kappa(sp,dU,phi,br)
      case ('fried_egg')
         kpair = fried_egg_kappa(sp,dU,phi,br)
      end select
   end function get_kappa
!-----------------------------------------------------------------------

!-------------------- charge density helpers ---------------------------
   function net_charge_density(species, dU, phi, br) result(rho)
      type(species_type), intent(in) :: species(:)
      real(dp), intent(in)           :: dU, phi, br
      real(dp)                       :: rho
      integer :: i
      rho = 0.0_dp
      do i = 1, size(species)
         rho = rho + species(i)%q * get_density(species(i),dU,phi,br)
      end do
   end function net_charge_density

   function local_relative_error(species, dU, phi, br, eIdx) result(err)
      type(species_type), intent(in) :: species(:)
      real(dp), intent(in)           :: dU, phi, br
      integer , intent(in)           :: eIdx           ! index of cold electrons
      real(dp)                       :: err, ne
      ne = get_density(species(eIdx), dU, phi, br)
      if (abs(ne) < 1.0e-30_dp) then
         err = huge(1.0_dp)
      else
         err = abs(net_charge_density(species,dU,phi,br) / (species(eIdx)%q*ne))
      end if
   end function local_relative_error
!-----------------------------------------------------------------------

!-------------------- DIFF_EQ root solver ------------------------------
   function diff_eq(species, dU, br, eIdx) result(res)
      !  species(:)  : array of species (electron must be at eIdx)
      !  dU          : ΔU for this location (J kg⁻¹)
      !  br          : B(s)/B₀
      !  eIdx        : index of cold electron species in the input array
      type(species_type), intent(inout) :: species(:)
      real(dp), intent(in)              :: dU, br
      integer , intent(in)              :: eIdx
      type(diffeq_result_type)          :: res

      integer :: nSpec,i, maxIter, iter
      real(dp) :: phi1,phi2,dphi,nq1,nq2,phiMid,relMid, tol
      real(dp) :: ionCharge   
      real(dp), allocatable :: nvals(:)

      nSpec = size(species)

      !---------- Step-1 : adjust cold-electron reference density -------
      ionCharge = 0.0_dp
      do i = 1, nSpec
         if (i /= eIdx) ionCharge = ionCharge + species(i)%q * species(i)%n
      end do
      species(eIdx)%n = - ionCharge / species(eIdx)%q   ! enforce neutrality at s₀

      !---------- Step-2 : bracket the root ----------------------------
      phi1 = 0.0_dp;  phi2 = 0.0_dp;  dphi = 1.0_dp
      nq1 = net_charge_density(species,dU,phi1,br)
      nq2 = net_charge_density(species,dU,phi2,br)
      maxIter = 10000; iter = 0
      do while (nq1*nq2 > 0.0_dp .and. iter < maxIter)
         phi1 = phi1 - dphi
         phi2 = phi2 + dphi
         nq1 = net_charge_density(species,dU,phi1,br)
         nq2 = net_charge_density(species,dU,phi2,br)
         iter = iter + 1
      end do
      if (iter == maxIter) then
         write(*,*) 'DIFF_EQ: failed to bracket root.'
         allocate(res%n(0)); res%phi = 0.0_dp; res%dU = dU
         return
      end if

      !---------- Step-3 : bisection to required tolerance -------------
      tol = 1.0e-8_dp ! Relative neutrality to electron density tolerance, change if you want tighter tolerance though then code is slower
      phiMid = 0.5_dp*(phi1+phi2)
      relMid = local_relative_error(species,dU,phiMid,br,eIdx)
      iter = 0
      do while (relMid > tol .and. iter < maxIter)
         if ( net_charge_density(species,dU,phiMid,br) * nq1 < 0.0_dp ) then
            phi2 = phiMid ; nq2 = net_charge_density(species,dU,phi2,br)
         else
            phi1 = phiMid ; nq1 = net_charge_density(species,dU,phi1,br)
         end if
         phiMid = 0.5_dp*(phi1+phi2)
         relMid = local_relative_error(species,dU,phiMid,br,eIdx)
         iter = iter + 1
      end do

      !---------- Step-4 : final densities -----------------------------
      allocate(nvals(nSpec))
      do i = 1, nSpec
         nvals(i) = get_density(species(i), dU, phiMid, br)
      end do

      res%n   = nvals
      res%phi = phiMid
      res%dU  = dU
   end function diff_eq
!-----------------------------------------------------------------------

!---------------- convenience vector wrappers -------------------------
   subroutine calc_temps(species, dU, phi, br, Tpar, Tperp)
     type(species_type), intent(in)  :: species(:)
     real(dp), intent(in)            :: dU, phi, br
     real(dp), intent(out)           :: Tpar(:), Tperp(:)

     integer :: i
     real(dp) :: pair(2)                       

     do i = 1, size(species)
        pair     = get_temperature(species(i), dU, phi, br)
        Tpar(i)  = pair(1)
        Tperp(i) = pair(2)
     end do
   end subroutine calc_temps


 subroutine calc_kappa_vals(species, dU, phi, br, kpar, kperp)
   type(species_type), intent(in)  :: species(:)
   real(dp), intent(in)            :: dU, phi, br
   real(dp), intent(out)           :: kpar(:), kperp(:)

   integer :: i
   real(dp) :: pair(2)                       

   do i = 1, size(species)
      pair     = get_kappa(species(i), dU, phi, br)
      kpar(i)  = pair(1)
      kperp(i) = pair(2)
   end do
 end subroutine calc_kappa_vals

end module diffusive_equilibrium_mod
