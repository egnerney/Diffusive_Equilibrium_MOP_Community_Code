!=======================================================================
!  basic_example_use_of_diffusive_equilibrium_code.f90
!Created May 2025
!
!@author: Edward (Eddie) G. Nerney
!
!Nerney (2025)
!
!Diffusive Equilibrium Code for basic example use along Io aligned JRM33+Con2020 field line
!finding steady state plasma Densities
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
program basic_example_use_of_diffusive_equilibrium_code
   !--------------------------------------------------------------------
   !  Modules
   !--------------------------------------------------------------------
   use, intrinsic :: iso_fortran_env , only : dp => real64
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use npy_reader
   use diffusive_equilibrium_mod
   implicit none
   !--------------------------------------------------------------------
   !  0.  Physical constants
   !--------------------------------------------------------------------
   real(dp), parameter :: kg_per_u    = 1.66053906892d-27
   real(dp), parameter :: e_charge    = 1.602176634d-19
   real(dp), parameter :: EV_TO_JOULE = e_charge
   real(dp), parameter :: JOULE_TO_EV = 1.0_dp / e_charge
   !--------------------------------------------------------------------
   !  1.  Field-line arrays
   !--------------------------------------------------------------------
   real(dp), allocatable :: s(:),x(:),y(:),z(:),rho(:),r(:),lat(:),wlon(:),B(:)
   integer :: npoints
   !--------------------------------------------------------------------
   !  2.  Reference-model matrices (601×10)
   !--------------------------------------------------------------------
   real(dp), allocatable :: n0_mat(:,:),T0_mat(:,:),k0_mat(:,:),kappaT0(:,:),tmp2(:,:)
   !--------------------------------------------------------------------
   !  3.  Species meta-data
   !--------------------------------------------------------------------
   integer, parameter         :: nspec      = 10
   integer, parameter         :: iofl_idx   = 191   ! Io field-line in ρ grid
   !--------------------------------------------------------------------
   !  Species display labels (uniform length = 7)
   !--------------------------------------------------------------------
   character(len=7), dimension(nspec), parameter :: spname =  &
        [ character(len=7) ::                              &
        'O+',      &  ! 2 chars, padded with blanks
        'O++',     &  ! 3
        'S+',      &  ! 2
        'S++',     &  ! 3
        'S+++',    &  ! 4
        'H+',      &  ! 2
        'Na+',     &  ! 3
        'O+(hot)', &  ! 7 (longest) so we set len to 7 adjust accordingly in len above 
        'eh-',     &  ! 3
        'e-' ]        ! 2

   real(dp), parameter :: m_e  = 5.485799090441d-4
   real(dp), parameter :: m_H  =  1.0078_dp
   real(dp), parameter :: m_O  = 15.999_dp
   real(dp), parameter :: m_S  = 32.065_dp
   real(dp), parameter :: m_Na = 22.990_dp
   real(dp)                       :: m_u(nspec), q_unit(nspec)
   real(dp)                       :: A_arr(nspec), lam_arr(nspec)
   !Avaiable Distribution Types:
   !character(len=32)              :: dist_tag = 'Maxwellian' ! Bagenal & Sullivan (1981) etc. Basic default Assumed isotropic always, ignores anisotropy
   !character(len=32)              :: dist_tag = 'Aniso_Maxwellian' ! Huang & Birmingham (1992) Anisotropic Maxwellian
   !character(len=32)              :: dist_tag = 'Aniso_Maxwellian_Const_Tperp' ! Bagenal & Sullivan (1981), Bagenal (1994), Mei, Thorne, Bagenal (1995) Assumed Constant Anisotrpy along field line, inconsistent but good enough often 
   !character(len=32)              :: dist_tag = 'Aniso_Kappa' ! Meyer-Vernet et al. (1995), Moncuquet et al. (2002) Anisotropic Standard Kappa distribution
   !character(len=32)              :: dist_tag = 'Aniso_Product_Kappa' ! Nerney (2025) Anisotropic Product Kappa Distribution 
   !character(len=32)              :: dist_tag = 'Fried_Egg' ! Nerney (2025) "Fried-Egg" Distribution Function (Maxwellian limit in parallel of aniso product kappa and still kappa in perp direction) 
  
   character(len=32)              :: dist_tag = 'Aniso_Maxwellian' ! Huang & Birmingham (1992) Anisotropic Maxwellian
   
   ! User-friendly title   ("Aniso_Maxwellian" → "Aniso Maxwellian") or whatever dist_tag is set to above raplces underscores with space for titles
   character(len=:), allocatable :: nice_tag
   integer :: pos
   

   ! When building titles / filenames use nice_tag or trimmed tag
   type(species_type)             :: sp(nspec)
   integer, parameter             :: cold_e_idx = 10
   integer                         :: j
   !--------------------------------------------------------------------
   !  4.  Planet geometry
   !--------------------------------------------------------------------
   character(len=*), parameter :: planet_ = 'Jupiter'
   real(dp), parameter         :: fcor    = 1.0_dp
   real(dp)                    :: pvals(3), B0
   integer                     :: ceq_idx
   real(dp), allocatable       :: Bratio(:), dU(:)
   !--------------------------------------------------------------------
   !  5.  Output arrays
   !--------------------------------------------------------------------
   real(dp), allocatable :: n_out(:,:), Tpar_out(:,:), Tperp_out(:,:)
   real(dp), allocatable :: kpar_out(:,:), kperp_out(:,:), phi_out(:)
   !--------------------------------------------------------------------
   !  Working scalars
   !--------------------------------------------------------------------
   integer :: i
   real(dp) :: Ttmp(2), Ktmp(2)
   type(diffeq_result_type) :: sol
   real(dp), parameter :: BIG = 1.0e100_dp
   real(dp) :: nan_dp
   integer :: npeak, ii
   real(dp), allocatable :: lat_peak(:), dU_row(:)

   nan_dp = ieee_value(0.0_dp, ieee_quiet_nan)

   dist_tag = trim(adjustl(dist_tag))
   nice_tag = dist_tag
   do pos = 1, len_trim(nice_tag)
      if (nice_tag(pos:pos) == '_') nice_tag(pos:pos) = ' '
   end do
   !====================================================================
   !  Load field-line vectors
   !====================================================================
   call read_npy_1d(s,   's.npy');      npoints = size(s)
   call read_npy_1d(x,   'x.npy');      call read_npy_1d(y,   'y.npy')
   call read_npy_1d(z,   'z.npy');      call read_npy_1d(rho, 'rho.npy')
   call read_npy_1d(r,   'r.npy');      call read_npy_1d(lat, 'lat.npy')
   call read_npy_1d(wlon,'wlong.npy');  call read_npy_1d(B,   'B.npy')
   !====================================================================
   !  Load reference model matrices
   !====================================================================
   call read_npy_2d(tmp2,'n0.npy');            n0_mat  = tmp2
   call read_npy_2d(tmp2,'T0.npy');            T0_mat  = tmp2
   call read_npy_2d(tmp2,'kappa0.npy');        k0_mat  = tmp2
   call read_npy_2d(tmp2,'kappa_temps_0.npy'); kappaT0 = tmp2
   !====================================================================
   !  Build species array
   !====================================================================
   m_u = [ (m_O-m_e), (m_O-2*m_e), (m_S-m_e), (m_S-2*m_e), (m_S-3*m_e), &
           (m_H-m_e), (m_Na-m_e),  (m_O-m_e), m_e, m_e ]
   q_unit = [ 1,2,1,2,3, 1,1,1,-1,-1 ]
   A_arr  = 2.0_dp
   lam_arr= 1.0_dp

   do j = 1, nspec
      sp(j) = species_init( trim(spname(j)),                      &
                            m_u(j)*kg_per_u,                      &
                            q_unit(j)*e_charge,                   &
                            T0_mat(iofl_idx,j)*EV_TO_JOULE,       &
                            A_arr(j), k0_mat(iofl_idx,j),         &
                            lam_arr(j), n0_mat(iofl_idx,j),       &
                            dist_tag )
   end do
   !====================================================================
   !  Planet constants & ΔU
   !====================================================================
   pvals  = define_planet(planet_)
   ceq_idx = maxloc(rho, dim=1)
   B0   = B(ceq_idx)
   allocate(Bratio(npoints), dU(npoints))
   Bratio = B / B0
   do i = 1, npoints
      dU(i) = calc_deltau( r(ceq_idx), lat(ceq_idx), r(i), lat(i), planet_, fcor )
   end do
   !====================================================================
   !  Quick ΔU vs latitude plot
   !====================================================================

   ! --- compute local maxima (simple 3-point criterion) ---------------
   allocate(lat_peak(npoints))
   allocate(dU_row(npoints)); dU_row = dU 
   npeak = 0
   do ii = 2, npoints-1
      if ( dU(ii) > dU(ii-1) .and. dU(ii) > dU(ii+1) ) then
         npeak = npeak + 1
         lat_peak(npeak) = lat(ii)
      end if
   end do
   if (npeak>0) lat_peak = lat_peak(:npeak)          ! trim

   call make_plot_gnuplot(lat, reshape(dU_row,[1,npoints]),         &
        ['ΔU'], 'Latitude (deg)', 'ΔU (J/kg)',                      &
        'Io FL: centrifugal+grav ΔU', 'deltaU_vs_lat.png',          &
        vlines = lat_peak)          
   !====================================================================
   !  Allocate outputs
   !====================================================================
   allocate( n_out(nspec,npoints), Tpar_out(nspec,npoints), Tperp_out(nspec,npoints))
   allocate( kpar_out(nspec,npoints),kperp_out(nspec,npoints), phi_out(npoints))
   !====================================================================
   !  Main field-line loop, could easily be vectorised or parallelized as each operation only depends on reference values
   !====================================================================
   do i = 1, npoints
      if (mod(i-1,100)==0) write(*,'(a,i0,a,i0)') 'point ',i-1,' / ',npoints-1
      sol = diff_eq(sp, dU(i), Bratio(i), cold_e_idx)
      n_out(:,i) = sol%n
      phi_out(i) = sol%phi
      do j = 1, nspec
         Ttmp = get_temperature(sp(j), dU(i), phi_out(i), Bratio(i))
         Tpar_out(j,i)  = Ttmp(1)
         Tperp_out(j,i) = Ttmp(2)
         Ktmp = get_kappa(sp(j), dU(i), phi_out(i), Bratio(i))
         kpar_out(j,i)  = Ktmp(1)
         kperp_out(j,i) = Ktmp(2)
      end do
   end do
   write(*,*) 'All points solved.'
   !====================================================================
   !  Produce six PNG plots
   !====================================================================
   call make_plot_gnuplot(s, n_out, spname,                      &
        's (R_J)', 'n (m^{-3})', nice_tag//' n vs s',            &
        trim(dist_tag)//'_n_vs_s.png', logy=.true.)
   call make_plot_gnuplot(lat, n_out, spname,                    &
        'Latitude (deg)', 'n (m^{-3})', nice_tag//' n vs lat',    &
        trim(dist_tag)//'_n_vs_lat.png', logy=.true.)
   call make_plot_gnuplot(lat, Tpar_out*JOULE_TO_EV, spname,     &
        'Latitude (deg)', 'T‖ (eV)', nice_tag//' T‖ vs lat',      &
        trim(dist_tag)//'_Tpar_vs_lat.png', logy=.true.)
   call make_plot_gnuplot(lat, Tperp_out*JOULE_TO_EV, spname,    &
        'Latitude (deg)', 'T⊥ (eV)', nice_tag//' T⊥ vs lat',      &
        trim(dist_tag)//'_Tperp_vs_lat.png', logy=.true.)

   ! should have been done in diffusive_equilbrium_mod.f90 kappa value funcs for Maxwellians but for now just done here, though in reality the Maxwellian limit is
   !  actually the infinite kappa limit for consistency accross languages we set it to NaNs
   where (kpar_out  > BIG) kpar_out  = nan_dp   ! NaN mask
   where (kperp_out > BIG) kperp_out = nan_dp
   
   call make_plot_gnuplot(lat, kpar_out, spname,                 &
        'Latitude (deg)', 'κ‖', nice_tag//' κ‖ vs lat',           &
        trim(dist_tag)//'_kappa_par_vs_lat.png', logy=.true.)
   call make_plot_gnuplot(lat, kperp_out, spname,                &
        'Latitude (deg)', 'κ⊥', nice_tag//' κ⊥ vs lat',           &
        trim(dist_tag)//'_kappa_perp_vs_lat.png', logy=.true.)
   write(*,*) 'PNG files written – done.'


 contains
   ! Attempt to get first 10 matplotlib default like colors for plots in gnuplot
   pure function tab10_hex(idx) result(col)
     !! index 0…9 → Matplotlib “tab10” hex string
     integer, intent(in) :: idx
     character(len=9)    :: col
     character(len=9), dimension(10), parameter :: tab10 = [ &
          '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', &
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf' ]
     col = tab10(mod(idx,10) + 1)
   end function tab10_hex

   pure function to_string(r) result(txt)
      real(dp), intent(in) :: r
      character(len=32) :: txt
      write(txt,'(es16.8)') r
   end function to_string

!-----------------------------------------------------------------------
!  gnuplot helper  –  skips all-NaN curves and can draw vertical lines
!-----------------------------------------------------------------------
subroutine make_plot_gnuplot(x, y, lbl, xt, yt, ttl, fname, logy, vlines)
   use, intrinsic :: iso_fortran_env ,  only : dp=>real64
   use, intrinsic :: ieee_arithmetic , only : ieee_is_nan
   implicit none
   real(dp),          intent(in)            :: x(:)
   real(dp),          intent(in)            :: y(:,:)
   character(*),      intent(in)            :: lbl(:), xt, yt, ttl, fname
   logical, optional, intent(in)            :: logy
   real(dp), optional,intent(in)            :: vlines(:)

   integer            :: ns, npts, i, j, k, udat, ugp, first
   character(len=9)   :: c
   logical, allocatable :: good(:)
   real(dp) :: ymin, ymax

   ns   = size(lbl);  npts = size(x)
   allocate(good(ns))

   ! mark curves that contain at least one finite value
   do j = 1, ns
      good(j) = any(.not. ieee_is_nan(y(j,:)))
   end do

   ! safety range for empty plots
   ymin = minval( merge(y, huge(1.0_dp), ieee_is_nan(y) ) )
   ymax = maxval( merge(y,-huge(1.0_dp), ieee_is_nan(y) ) )
   if (ymin >= ymax) then
      ymin =  1.0_dp
      ymax = 10.0_dp
   end if
   !---------------- data file -----------------------------------------
   open (newunit = udat, file = '__tmp_data.dat', status = 'replace')

   do i = 1, npts
      ! We keep ALL rows and let gnuplot treat the ‘NaN’ strings as gaps.
      write(udat,'(*(es20.10,1x))') x(i), ( y(j,i), j = 1, ns )
   end do
   close(udat)
   !---------------- gnuplot script ------------------------------------
   open(newunit=ugp,file='__tmp_plot.gp',status='replace')
   write(ugp,*) "set term pngcairo size 2500,1875"
   write(ugp,*) "set output '"//trim(fname)//"'"
   write(ugp,*) "set encoding utf8"
   write(ugp,*) "set xlabel '"//trim(xt)//"'"
   write(ugp,*) "set ylabel '"//trim(yt)//"'"
   write(ugp,*) "set title  '"//trim(ttl)//"'"
   write(ugp,*) "set datafile missing 'NaN'"

   if (present(logy) .and. logy .and. any(good)) then
      write(ugp,*) "set logscale y"
      write(ugp,*) "set yrange [1e-3:*]"
   end if                !! ← no manual yrange for linear plots
   ! optional vertical dashed lines
   if (present(vlines)) then
      do i = 1, size(vlines)
         write(ugp,'("set arrow from ",es14.6,", graph 0 to ",es14.6,", graph 1 nohead dt 2 lw 1")') &
              vlines(i), vlines(i)
      end do
   end if

   ! style lines
   do j=1,ns
      c = tab10_hex(j-1)
      write(ugp,'("set style line ",i0," lt 1 lw 2 lc rgb ''",a,"''")') j, trim(c)
   end do
   !---------------- plot command ---------------------------------------
   if (.not. any(good)) then                       ! every curve is NaN
      ! Establish finite axis limits
      write(ugp,*) "set xrange ["//trim(to_string(minval(x)))//":"// &
           trim(to_string(maxval(x)))//"]"
      write(ugp,*) "set yrange [0.0:1.0]"

      ! Draw one invisible line to force the frame
      write(ugp,*) "plot '-' using 1:2 with lines lc rgb '#ffffff' notitle"
      write(ugp,*) x(1), 0.0_dp
      write(ugp,*) x(npts), 0.0_dp
      write(ugp,*) "e"
   else
      do first = 1, ns; if (good(first)) exit; end do
         write(ugp,'("plot ''__tmp_data.dat'' using 1:",i0," with lines ls ",i0,&
              " title ''",a,"''",$)') first+1, first, trim(lbl(first))
         do k = first + 1, ns
            if (good(k)) then
               write(ugp,'(", ''__tmp_data.dat'' using 1:",i0," with lines ls ",i0,&
                    " title ''",a,"''",$)') k+1, k, trim(lbl(k))
            end if
         end do
         write(ugp,*)
   end if
   close(ugp)
   !---------------- run & tidy -----------------------------------------
   call execute_command_line("gnuplot __tmp_plot.gp")
   call execute_command_line("rm -f __tmp_data.dat __tmp_plot.gp")
end subroutine make_plot_gnuplot


end program basic_example_use_of_diffusive_equilibrium_code
