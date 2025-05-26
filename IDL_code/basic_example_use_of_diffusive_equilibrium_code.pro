
;------------------------------------------------------------------
; make_plot  –  Function-Graphics log-y plot with legend
;------------------------------------------------------------------
PRO make_plot, x, yarr, lbl, xtitle, ytitle, title, file

  ;-----------------------------------------------------------------
  ;  Matplotlib’s first 10 “tab10” colors (hex strings are fine in
  ;  IDL 8.4+ function graphics).
  ;-----------------------------------------------------------------
  mp_colors = [  $
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', $
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf' ]

  ns       = N_ELEMENTS(lbl)
  handles  = OBJARR(ns)

  ;-----------------------------------------------------------------
  ;  First curve ───────────────────────────────────────────────────
  ;-----------------------------------------------------------------
  handles[0] = PLOT( x, yarr[0,*], $
                     YLOG      = 1, $
                     THICK     = 2, $
                     COLOR     = mp_colors[0], $
                     NAME      = lbl[0], $
                     XTITLE    = xtitle, $
                     YTITLE    = ytitle, $
                     TITLE     = title )

  ;-----------------------------------------------------------------
  ;  Additional curves (over-plot) ─────────────────────────────────
  ;-----------------------------------------------------------------
  FOR i = 1, ns-1 DO BEGIN
    handles[i] = PLOT( x, yarr[i,*], /OVERPLOT, $
                       YLOG  = 1, $
                       THICK = 2, $
                       COLOR = mp_colors[i], $
                       NAME  = lbl[i] )
  ENDFOR
  
  leg = LEGEND( TARGET               = handles, $
                /RELATIVE, $
                POSITION             = [0.98, 0.3], $
                HORIZONTAL_ALIGNMENT = 'center', $
                VERTICAL_ALIGNMENT   = 'center', $
                ORIENTATION          = 0, $   
                /AUTO_TEXT_COLOR )

  ;-----------------------------------------------------------------
  ;  Save to PNG at 500 dpi
  ;-----------------------------------------------------------------
  handles[0].Save, file, /PNG, RESOLUTION=500
  PRINT, 'Saved ', file

END

@diffusive_equilibrium

PRO BASIC_EXAMPLE_USE_OF_DIFFUSIVE_EQUILIBRIUM_CODE
  ;*************************************************************
  ;  BASIC_EXAMPLE_USE_OF_DIFFUSIVE_EQUILIBRIUM_CODE.pro
  ;*************************************************************
  ;COMPILE_OPT idl2

  ;-----------------------------------------------------------
  ; 0.  Constants
  ;-----------------------------------------------------------
  kg_per_u         = 1.66053906892d-27
  ELEMENTARY_CHARGE = 1.602176634d-19
  EV_TO_JOULE       = ELEMENTARY_CHARGE
  JOULE_TO_EV       = 1d / ELEMENTARY_CHARGE

  ;-----------------------------------------------------------
  ; 1.  READ FIELD-LINE TRACE (Fortran-order arrays)
  ;-----------------------------------------------------------
  f1  = 'Io_aligned_field_line_trace_JRM33+con2020_dlat=0.1_degree.h5'
  fid = H5F_OPEN(f1)
  s      = H5D_READ(H5D_OPEN(fid,'s'))
  x      = H5D_READ(H5D_OPEN(fid,'x'))
  y      = H5D_READ(H5D_OPEN(fid,'y'))
  z      = H5D_READ(H5D_OPEN(fid,'z'))
  rho    = H5D_READ(H5D_OPEN(fid,'rho'))
  r      = H5D_READ(H5D_OPEN(fid,'r'))
  lat    = H5D_READ(H5D_OPEN(fid,'lat'))
  wlong  = H5D_READ(H5D_OPEN(fid,'wlong'))
  Bmag   = H5D_READ(H5D_OPEN(fid,'B'))
  H5F_CLOSE, fid
  npoints = N_ELEMENTS(s)

  ;-----------------------------------------------------------
  ; 2.  READ CENTRIFUGAL-EQUATOR REFERENCE MODEL
  ;-----------------------------------------------------------
  f2  = 'Nerney2025_reference_model_4-10RJ.h5'
  fid = H5F_OPEN(f2)
  rho_ceq  = H5D_READ(H5D_OPEN(fid,'rho_ceq'))
  
  ; transpose each 2-D matrix so dim0 = 601 (ρ-grid), dim1 = 10 (species)
  n0_mat  = TRANSPOSE( H5D_READ(H5D_OPEN(fid,'n0')) )
  T0_mat  = TRANSPOSE( H5D_READ(H5D_OPEN(fid,'T0')) )
  k0_mat  = TRANSPOSE( H5D_READ(H5D_OPEN(fid,'kappa0')) )
  kT0_mat = TRANSPOSE( H5D_READ(H5D_OPEN(fid,'kappaT0')) )
  
  sps_byte = H5D_READ(H5D_OPEN(fid,'species'))
  

  H5F_CLOSE, fid
  species_names = STRTRIM(STRING(sps_byte),2)
  nspec   = N_ELEMENTS(species_names)

  ; 601 field lines between 4.00 and 10.00 at 0.01 RJ resolution
  ; so index that corresponds to Io field line (rho_ceq = 5.91 RJ) is 191
  iofl_idx = 191 

  ;-----------------------------------------------------------
  ; 3.  DEFINE SPECIES PARAMETERS  (same ordering as Python)
  ;-----------------------------------------------------------
  m_e = 5.485799090441d-4
  m_H  = 1.0078d & m_O  = 15.999d & m_S = 32.065d & m_Na = 22.990d
  m_Hp  = m_H  - m_e & m_Op = m_O - m_e & m_O2p = m_O - 2*m_e
  m_Sp  = m_S  - m_e & m_S2p= m_S - 2*m_e & m_S3p= m_S - 3*m_e
  m_Nap = m_Na - m_e
  species_m = [m_Op,m_O2p,m_Sp,m_S2p,m_S3p,m_Hp,m_Nap,m_Op,m_e,m_e]
  species_q = [ 1,2,1,2,3, 1,1,1,-1,-1 ]

  species_A   = REPLICATE(2.0d, nspec)
  species_lam = REPLICATE(1.0d, nspec)
  ;Avaiable Distribution Types:
  ;type_ = 'Maxwellian' ; Bagenal & Sullivan (1981) etc. Basic default Assumed isotropic always, ignores anisotropy
  type_ = 'Aniso_Maxwellian' ; Huang & Birmingham (1992) Anisotropic Maxwellian
  ;type_ = 'Aniso_Maxwellian_Const_Tperp' ; Bagenal & Sullivan (1981), Bagenal (1994), Mei, Thorne, Bagenal (1995) Assumed Constant Anisotrpy along field line, inconsistent but good enough often 
  ;type_ = 'Aniso_Kappa' ; Meyer-Vernet et al. (1995); Moncuquet et al. (2002) Anisotropic Standard Kappa distribution
  ;type_ = 'Aniso_Product_Kappa' ; Nerney (2025) Anisotropic Product Kappa Distribution 
  ;type_ = 'Fried_Egg' ; Nerney (2025) "Fried-Egg" Distribution Function (Maxwellian limit in parallel of aniso product kappa and still kappa in perp direction) 
  species_type = REPLICATE(type_, nspec)
  ; -- build an empty template that has the right tag list ----
  template = SPECIES_INIT('', 0d, 0d, 0d, 0d, 0d, 0d, 0d, '')

  ; -- replicate that template nspec times --------------------
  species_arr = REPLICATE(template, nspec)
  
  FOR i=0, nspec-1 DO BEGIN
    n0   = n0_mat[iofl_idx,i]
    T0ev = T0_mat[iofl_idx,i]
    k0   = k0_mat[iofl_idx,i]

    species_arr[i] = SPECIES_INIT( $
      species_names[i], $
      species_m[i]*kg_per_u, $
      species_q[i]*ELEMENTARY_CHARGE, $
      T0ev*EV_TO_JOULE, $
      species_A[i], $
      k0, $
      species_lam[i], $
      n0, $
      species_type[i] )
  ENDFOR

  ;-----------------------------------------------------------
  ; 4.  PLANET DATA & GEOMETRY ARRAYS
  ;-----------------------------------------------------------
  ; Set planet to Jupiter and and perfect corotation with fcor =1.0
  planet_ = 'Jupiter' & fcor = 1.0d
  pvals  = DEFINE_PLANET(planet_)      ; [RP, Ω, GM]
  ceq_idx = WHERE(rho EQ MAX(rho), /NULL) ; Find Centrifugal equator along field line
  B0      = Bmag[ceq_idx[0]]            ; make it a true scalar
  B_ratio  = Bmag / B0                   ; now Bratio is length = npoints
  deltaU  = CALC_DELTAU(r[ceq_idx[0]], lat[ceq_idx[0]], r, lat, planet_ , fcor)

  ;-----------------------------------------------------------
  ; 5.  QUICK ΔU ( Centrifugal + Gravitational PE/mass) vs LAT PLOT
  ;-----------------------------------------------------------
  p0 = PLOT( lat, deltaU, TITLE='Io Field Line, Centrifugal + Gravitational $\Delta$PE/m vs Latitude', $
    XTITLE='Latitude ($\deg$)', YTITLE='$\Delta$U (J/kg)',THICK=2)
    
  ; draw vertical lines at local maxima
  maxS = FIND_LOCAL_MAXIMA(deltaU, lat)
  IF SIZE(maxS, /TYPE) EQ 8 THEN BEGIN   ; 8 → structure
    FOR j=0,N_ELEMENTS(maxS.lat)-1 DO $
      p1 = Plot( [maxS.lat[j],maxS.lat[j]], [MIN(deltaU),MAX(deltaU)], LINESTYLE=2,/overplot)
  ENDIF
  
  ;p1.save, 'deltaU_vs_lat_Io_aligned_FL.png', RESOLUTION=500



  ;-----------------------------------------------------------
  ; 6.  OUTPUT ARRAYS
  ;-----------------------------------------------------------
  n_out       = DBLARR(nspec, npoints)
  Tpar_out    = DBLARR(nspec, npoints)
  Tperp_out   = DBLARR(nspec, npoints)
  kpar_out    = DBLARR(nspec, npoints)
  kperp_out   = DBLARR(nspec, npoints)
  phi_out     = DBLARR(npoints)

  ;-----------------------------------------------------------
  ; 7.  MAIN LOOP ALONG FIELD LINE
  ;-----------------------------------------------------------
  FOR i=0, npoints-1 DO BEGIN
    IF (i MOD 100 EQ 0) THEN PRINT, 'point ',i,' / ',npoints-1

    res = DIFF_EQ([r[i],lat[i],wlong[i]], [r[ceq_idx[0]],lat[ceq_idx[0]],wlong[ceq_idx[0]]], $
      B_ratio[i], species_arr, deltaU[i], planet_ , fcor)

    n_out[*,i] = res.n
    phi_out[i] = res.PHI

    CALC_TEMPS,  species_arr, deltaU[i], phi_out[i], B_ratio[i], Tpl, Tpp
    CALC_KAPPA_VALS, species_arr, deltaU[i], phi_out[i], B_ratio[i], kpl, kpp
    Tpar_out[*,i]  = Tpl &  Tperp_out[*,i] = Tpp
    kpar_out[*,i]  = kpl &  kperp_out[*,i] = kpp
  ENDFOR

  ;-----------------------------------------------------------
  ; 8.  PLOTS  (densities vs s)
  ;-----------------------------------------------------------
  species_lbl = ['O!U+!N','O!U++!N','S!U+!N','S!U++!N','S!U+++!N', $
                 'H!U+!N','Na!U+!N','O!U+!N (hot)','e!U-!N (hot)','e!U-!N (cold)']

  
  ;-------------------------------------------------------------
  ; 1.  n(s)
  ;-------------------------------------------------------------
  make_plot, s, n_out, species_lbl, 's (R!DJ!N)', 'n (cm!U-3!N)', $
    STRJOIN( STRSPLIT( type_, '_', /EXTRACT ), ' ' ) + ' n vs s, Aligned Io FL', type_ + '_n_vs_s_Io_aligned_FL.png'

  ;-------------------------------------------------------------
  ; 2.  n(lat)
  ;-------------------------------------------------------------
  make_plot, lat, n_out, species_lbl, 'Latitude!DIII!N (deg)', 'n (cm!U-3!N)', $
    STRJOIN( STRSPLIT( type_, '_', /EXTRACT ), ' ' ) + ' n vs Latitude, Aligned Io FL', type_ + '_n_vs_lat_Io_aligned_FL.png'

  ;-------------------------------------------------------------
  ; 3.  T‖(lat)  [convert to eV]
  ;-------------------------------------------------------------
  make_plot, lat, Tpar_out*JOULE_TO_EV, species_lbl, 'Latitude!DIII!N (deg)', 'T!D‖!N (eV)', $
    STRJOIN( STRSPLIT( type_, '_', /EXTRACT ), ' ' ) + ' T!D‖!N vs Latitude, Aligned Io FL', type_ + '_Tpar_vs_lat_Io_aligned_FL.png'

  ;-------------------------------------------------------------
  ; 4.  T⊥(lat)
  ;-------------------------------------------------------------
  make_plot, lat, Tperp_out*JOULE_TO_EV, species_lbl, 'Latitude!DIII!N (deg)', 'T!D⊥!N (eV)', $
    STRJOIN( STRSPLIT( type_, '_', /EXTRACT ), ' ' ) + ' T!D⊥!N vs Latitude, Aligned Io FL', type_ + '_Tperp_vs_lat_Io_aligned_FL.png'

  ;-------------------------------------------------------------
  ; 5.  κ‖(lat)
  ;-------------------------------------------------------------
  make_plot, lat, kpar_out, species_lbl, 'Latitude!DIII!N (deg)', 'κ!D‖!N', $
    STRJOIN( STRSPLIT( type_, '_', /EXTRACT ), ' ' ) + ' κ!D‖!N vs Latitude, Aligned Io FL', type_ + '_kappa_par_vs_lat_Io_aligned_FL.png'

  ;-------------------------------------------------------------
  ; 6.  κ⊥(lat)
  ;-------------------------------------------------------------
  make_plot, lat, kperp_out, species_lbl, 'Latitude!DIII!N (deg)', 'κ!D⊥!N', $
    STRJOIN( STRSPLIT( type_, '_', /EXTRACT ), ' ' ) + ' κ!D⊥!N vs Latitude, Aligned Io FL', type_ + '_kappa_perp_vs_lat_Io_aligned_FL.png'
  
  PRINT, 'All computations and PNG files finished.'
END
