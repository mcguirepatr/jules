! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Routine to calculate the spectral albedo of the land surface and the profile
! of PAR through the canopy used by some canopy radiation models using
! the two stream approach of Sellers, 1995. It is called both to calculate the
! albedo if l_spec_albedo is true and to calculate the absorption of sunlight
! for photosynthesis.
! *********************************************************************
MODULE albpft_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALBPFT_MOD'
CONTAINS
SUBROUTINE albpft(                                                             &
  !INTENT(IN)
    l_getprofile, land_pts, ilayers, albpft_call,                              &
    surft_pts, surft_index,                                                    &
    cosz_gb, lai_physical, albudir, albudif,                                   &
  !INTENT(INOUT)
    alb_type,                                                                  &
  !INTENT(OUT)
    fapar_dir, fapar_dif, fapar_dir2dif,                                       &
    fapar_dif2dif, fapar_dir2dir, fsun,                                        &
    !New arguments replacing USE statements
    !jules_mod (IN OUT)
    albobs_scaling_surft)

USE jules_surface_types_mod, ONLY: npft, ntype
USE pftparm, ONLY: omega,  omnir,  alpar,  alnir,                              &
                   omegal, omnirl, alparl, alnirl,                             &
                   omegau, omniru, alparu, alniru,                             &
                   can_struct_a, orient

USE stochastic_physics_run_mod, ONLY: l_rp2, i_rp_scheme, i_rp2b,              &
                                      alnir_rp, alpar_rp, omega_rp, omnir_rp

USE jules_radiation_mod, ONLY: l_spec_albedo, l_albedo_obs, l_niso_direct

USE jules_vegetation_mod, ONLY: can_rad_mod

USE ancil_info, ONLY: rad_nband

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE ereport_mod, ONLY: ereport

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER, INTENT(IN) ::                                                         &
 land_pts                 &! Number of land points.
,ilayers                    &! Number of layers over which
                             ! the PAR absorption profile is to
                             ! be calculated.
,albpft_call                 ! 0 - No scaling to obs
                             ! 1 - receive an input scaling to obs
                             ! and apply a correction to it for the
                             ! non-linearity (and store that scaling)
                             ! 2 - use the stored scaling to obs

LOGICAL, INTENT(IN) :: l_getprofile
                             ! TRUE means profiles through the
!                                    canopy are calculated
!                                  ! FALSE means no profiles calculated

!   Array arguments with intent(in):
INTEGER, INTENT(IN) ::                                                         &
 surft_pts(ntype)                                                              &
                             ! Number of land points which
!                                  ! include the surface type.
,surft_index(land_pts,ntype)
                             ! Indices of land points which
!                                  ! include the surface type.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 albudir(land_pts,2,npft)                                                      &
                             ! Direct albedo of the underlying
!                            ! surface.
,albudif(land_pts,2,npft)                                                      &
                             ! Diffuse albedo of the underlying
!                            ! surface.
,cosz_gb(land_pts)                                                             &
                             ! Cosine of the zenith angle.
,lai_physical(land_pts,npft) ! Leaf area index.

!   Array arguments with intent(out):
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
fapar_dir(land_pts,npft,ilayers)                                               &
!                                  ! Profile of absorbed PAR
!                                  !     - Direct (fraction/LAI)
,fapar_dif(land_pts,npft,ilayers)                                              &
!                                  ! Profile of absorbed PAR
!                                  !     - Diffuse (fraction/LAI)
,fapar_dir2dif(land_pts,npft,ilayers)                                          &
!                                  ! Profile of scattered PAR
!                                  !   -direct to diffuse
!                                  !      (fraction/LAI)
,fapar_dif2dif(land_pts,npft,ilayers)                                          &
!                                  ! Profile of absorbed PAR
!                                  !     - Diffuse (fraction/LAI)
,fapar_dir2dir(land_pts,npft,ilayers)                                          &
!                                  ! Profile of absorbed only direct PAR
!                                  !     -  (fraction/LAI)
,fsun(land_pts,npft,ilayers)
                             ! fraction of sunlit leaves

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 alb_type(land_pts,ntype,4)
!                                  ! Albedos for surface types.
!                                  !   (*,*,1) - Direct beam visible
!                                  !   (*,*,2) - Diffuse visible
!                                  !   (*,*,3) - Direct beam near-IR
!                                  !   (*,*,4) - Diffuse near-IR

!New arguments replacing USE statements
!jules_mod (IN OUT)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: albobs_scaling_surft(land_pts,ntype,rad_nband)

! Local arrays:
REAL(KIND=real_jlslsm) ::                                                      &
 alpl                                                                          &
                            ! Leaf reflection coeff
,om                                                                            &
                            ! Leaf scattering coeff
,rnet_dir(0:ilayers)                                                           &
!                                 ! Net downward PAR at layer
!                                 ! boundaries - incident direct
,rnet_dif(0:ilayers)                                                           &
!                                 ! Net downward PAR at layer
!                                 ! boundaries - incident diffuse
,taul                                                                          &
                                  ! Leaf transmission coefficient.
,lai(land_pts,npft)
                                  ! effective leaf area index
! Local scalars:
REAL(KIND=real_jlslsm) ::                                                      &
 betadir                                                                       &
                            ! Upscatter parameter for direct beam.
,betadif                                                                       &
                            ! Upscatter parameter for diffuse beam
,coszm                                                                         &
                            ! Mean value of COSZ.
,k                                                                             &
                            ! Optical depth per unit leaf area.
,g                                                                             &
                            ! Relative projected leaf area in
                            ! direction cosz.
,salb                                                                          &
                            ! Single scattering albedo.
,sqcost                                                                        &
                            ! Cosine squared of the mean leaf
                            ! angle to the horizontal.
,albdif_tmp                                                                    &
                            ! Temporary diffuse albedo
,albscl_tmp                                                                    &
                            ! Temporary scaling value
,albscl_new                                                                    &
                            ! Temporary scaling value
,b,c,ca,d,f,h,u1                                                               &
                            ! Miscellaneous variables from
,p1,p2,p3,p4,d1                                                                &
                            ! Sellers (1985).
,h1,h2,h3,h7,h8                                                                &
                            !
,s1,s2,sig                  !

!-----------------------------------------------------------------------
! Additional work variables to calculate profiles of absorbed PAR
!-----------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
 dlai                                                                          &
                                ! LAI increment.
,la                                                                            &
                                ! Cumulative LAI through canopy.
,drdird_dlai,drdiru_dlai                                                       &
,drdifd_dlai,drdifu_dlai                                                       &
,u2,u3,d2,h4,h5,h6,h9,h10                                                      &
                                ! Rate of change of fluxes with LAI
                                ! (W/m2/LAI).
,rdiru,rdird,rdifu,rdifd                                                       &
                                ! work variable
,dabeer_dla(0:ilayers)          ! Attenuation of direct solar
                                ! beam per unit LAI

!Optimisation variables to reduce the number of EXP evaluations
REAL(KIND=real_jlslsm) :: exp_hla, exp_minus_hla, exp_minus_kla,               &
        exp_k_dlai_min_la

!Optimisation variables to reduce the number of IF evaluations
REAL(KIND=real_jlslsm) :: omega_arr(2),  alpa_arr(2),                          &
        omegal_arr(2), alpal_arr(2),                                           &
        omegau_arr(2), alpau_arr(2)

!variables for layeredCanopyTwoStream
REAL :: lower_boundary_r, Ru_last, Rc_last, Rd_last                           &
,mu, muBar, Gfunc, Zfunc, Kfunc, w, leaf_r, leaf_t                            &
,ssa, B0, scalefunc, g1, g2, g3, g4                                           &
,tau, tauD, a1, a2, Rd, Td, Rc, Ru, Tu                                        &
,Tc, U, Z, RuNew, TuNew, RcNew, TcNew                                         &
,Zdif, rdirup(ilayers), rdirdn(ilayers)                                       &
,rdirup_u(ilayers), rdirup_c(ilayers), rdirdn_u(ilayers)                      &
,rdirdn_c(ilayers), rdifup(ilayers), rdifdn(ilayers)                          &
,F11,G11,H11,H112 ! change from H to H11, etc., to avoid duplicate variables

!,d, D, B, k, G, F, H, Tc
REAL :: Tdif, S, rdirdn_above, rdifdn_above, rdirup_below                     &
,rdifup_below, rdn, rup, t2

INTEGER ::                                                                    &
 band,i,j,l,n               ! Loop counters.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALBPFT'

INTEGER ::              errcode            ! Error code
CHARACTER(LEN=80) :: ERRMSG             ! Error message text
CHARACTER(LEN=50000) :: lineBuffer
CHARACTER(LEN=50000) :: lineBuffer

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check the consistency of the options.
IF ( albpft_call < 0 .OR. albpft_call > 2) THEN
  errcode = 1
  WRITE(ERRMSG,'(A,I0)') 'albpft_call must be 0, 1 or 2. '                     &
      // 'Called with albpft_call = ', albpft_call
  CALL ereport (routinename,errcode,ERRMSG)
END IF
IF ( albpft_call > 0 ) THEN
  IF ( .NOT. (l_albedo_obs .AND. l_spec_albedo) ) THEN
    errcode = 1
    ERRMSG  = 'invalid call to albpft_jls, needs l_albedo_obs and l_spec_albedo'
    CALL ereport ('albfpt',errcode,ERRMSG)
  END IF
END IF

IF ( l_rp2 .AND. i_rp_scheme == i_rp2b ) THEN
  DO n = 1,npft
    alnir(n) = alnir_rp(n)
    alpar(n) = alpar_rp(n)
    omega(n) = omega_rp(n)
    omnir(n) = omnir_rp(n)
  END DO
END IF

!$OMP PARALLEL DEFAULT(NONE)                                                   &
!$OMP PRIVATE(n, l, i, j, albdif_tmp, albscl_new, albscl_tmp, b, betadif,      &
!$OMP         betadir, c, ca, coszm, d, d1, d2, dabeer_dla, dlai,              &
!$OMP         drdird_dlai, drdiru_dlai, drdifu_dlai, drdifd_dlai, f, g, h,     &
!$OMP         h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, k, la, p1, p2,          &
!$OMP         p3, p4, rdifd, rdifu, rdird, rdiru, rnet_dif, rnet_dir, s1,      &
!$OMP         s2, salb, sig, sqcost, taul, u1, u2, u3, exp_hla,                &
!$OMP         exp_minus_hla, exp_minus_kla, om, alpl, band, omega_arr,         &
!$OMP         alpa_arr, omegal_arr, alpal_arr, omegau_arr,alpau_arr,           &
!$OMP         exp_k_dlai_min_la)                                               &
!$OMP SHARED(npft,omega,omnir,alpar,alnir,omegal,omnirl,alparl,alnirl,         &
!$OMP        omegau,omniru,alparu,alniru,lai_physical,surft_pts,               &
!$OMP        surft_index, albpft_call, land_pts,                               &
!$OMP        albobs_scaling_surft,  orient, l_getprofile,                      &
!$OMP        alb_type, cosz_gb, lai, albudif, ilayers, fapar_dir,              &
!$OMP        fapar_dif, fapar_dir2dir, fapar_dir2dif, fapar_dif2dif,           &
!$OMP        albudir, can_rad_mod, can_struct_a, fsun, l_niso_direct,          &
!$OMP        l_albedo_obs)

! Initialisations
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
DO n = 1,npft
  DO l = 1,land_pts
    lai(l,n) = 0.0
  END DO
END DO
!$OMP END DO NOWAIT

! Initialise to zero at top.
!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
DO i = 1,ilayers
  DO n = 1,npft
    DO l = 1, land_pts
      fapar_dir2dir(l,n,i) = 0.0
      fapar_dir2dif(l,n,i) = 0.0
      fapar_dif2dif(l,n,i) = 0.0
      fsun(l,n,i)          = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

! We need an explicit OMP barrier to ensure that the arrays are fully
! initialised before use.
!$OMP BARRIER

DO n = 1,npft

  !Package up variables to fit in with looping over bands. They come as pairs:
  !omega, omnir:   Leaf scattering coefficient for PAR.
  omega_arr(1)  = omega(n)
  omega_arr(2)  = omnir(n)
  !alpar, alnir Leaf reflection coefficient for PAR.
  alpa_arr(1)   = alpar(n)
  alpa_arr(2)   = alnir(n)
  ! IF l_albedo_obs = F then the following will be set to missing data
  IF ( l_albedo_obs ) THEN
    !omegal, omnirl: lower limit on omega, when scaled to albedo obs
    omegal_arr(1) = omegal(n)
    omegal_arr(2) = omnirl(n)
    !alparl, alnirl: lower limit on alpar, when scaled to albedo obs
    alpal_arr(1)  = alparl(n)
    alpal_arr(2)  = alnirl(n)
    !omegau, omniru: upper limit on omega, when scaled to albedo obs
    omegau_arr(1) = omegau(n)
    omegau_arr(2) = omniru(n)
    !alparu, alniru: upper limit on alpar, when scaled to albedo obs
    alpau_arr(1)  = alparu(n)
    alpau_arr(2)  = alniru(n)
  END IF

!$OMP DO SCHEDULE(STATIC)
  DO j = 1,surft_pts(n)
    l = surft_index(j,n)
    lai(l,n) = lai_physical(l,n) * can_struct_a(n)
  END DO
!$OMP END DO NOWAIT

  DO band = 1,2  ! Visible and near-IR bands

    ! If albpft_call is 0, then either there's no scaling to albedo obs,
    ! or this is the first call in that process to make a first guess:
    IF (albpft_call == 0 ) THEN
      om   = omega_arr(band)
      alpl = alpa_arr(band)
    END IF

!$OMP DO SCHEDULE(STATIC)
    DO j = 1,surft_pts(n)
      l = surft_index(j,n)

      IF ( albpft_call == 1 ) THEN

        !Firstly apply the scaling for this point to the inputs:
        om   = omega_arr(band) * albobs_scaling_surft(l,n,band)
        alpl = alpa_arr(band)  * albobs_scaling_surft(l,n,band)

        ! Test to see if the scaling puts us in totally unphysical part of
        !phase-space and set to the upper/lower values if it is:
        IF (om <= 0.0 .OR. alpl <= 0.0 ) THEN
          om   = omegal_arr(band)
          alpl = alpal_arr(band)

          ! Don't do anything with the scaling, leave it to make unphysical values
          ! so the later call will just use the upper/lower range values
        ELSE IF (om >= 1.0 .OR. alpl >= 1.0 ) THEN
          om   = omegau_arr(band)
          alpl = alpau_arr(band)
        ELSE
          ! Now work out a temporary diffuse albedo for the band
          !(a predictor):
          taul = om - alpl
          IF (orient(n) == 0) THEN
            sqcost = 1.0 / 3.0
            coszm  = 1.0
          ELSE IF (orient(n) == 1) THEN
            sqcost = 1.0
            coszm  = 1.0
          END IF

          c           = 0.5 * ( alpl + taul + (alpl - taul) * sqcost )
          betadif     = c / om
          b           = 1.0 - (1.0 - betadif) * om
          h           = SQRT(b * b - c * c) / coszm
          u1          = b - c / albudif(l,band,n)
          s1          = EXP(-h * lai(l,n))
          p1          = b + coszm * h
          p2          = b - coszm * h
          d1          = p1 * (u1 - coszm * h) / s1 - p2 * (u1 + coszm * h) * s1
          h7          =  (c / d1) * (u1 - coszm * h) / s1
          h8          = -(c / d1) * (u1 + coszm * h) * s1
          albdif_tmp  = h7 + h8

          ! compare it to the albedo from the previous first guess call:
          albscl_tmp = albdif_tmp / alb_type(l,n,2 * band)
          ! and apply a corrector to the scaling factor to account for
          !non-linearity:
          albscl_new = albobs_scaling_surft(l,n,band)                          &
                       * (albobs_scaling_surft(l,n,band) / albscl_tmp)
          albobs_scaling_surft(l,n,band) = albscl_new

          om   = omega_arr(band) * albscl_new
          alpl = alpa_arr(band)  * albscl_new

        END IF ! ends test of totally unphysical 'predictor' scaled constants.

      ELSE IF ( albpft_call == 2 ) THEN
        ! Using pre-corrected scalings, so just apply scaling:
        om   = omega_arr(band) * albobs_scaling_surft(l,n,band)
        alpl = alpa_arr(band)  * albobs_scaling_surft(l,n,band)
      END IF

      ! If some scaling has been applied, check to see if any of the values are
      ! unphyical, then check the values are in the lower to upper value range:
      ! (consistent with the unphsyical values checks above)
      IF (albpft_call > 0) THEN
        IF (om <= 0.0 .OR. alpl <= 0.0 ) THEN
          om   = omegal_arr(band)
          alpl = alpal_arr(band)
        ELSE IF (om >= 1.0 .OR. alpl >= 1.0 ) THEN
          om   = omegau_arr(band)
          alpl = alpau_arr(band)
        ELSE
          om   = MIN( MAX(om, omegal_arr(band)), omegau_arr(band) )
          alpl = MIN( MAX(alpl, alpal_arr(band)), alpau_arr(band) )
        END IF
      END IF

      ! Now continue with rest of the subroutine with the om and alpl:
      taul = om - alpl

      IF (orient(n) == 0) THEN
        sqcost = 1.0 / 3.0
        g      = 0.5
        coszm  = 1.0
        k      = g / 0.01
        salb   = 0.5 * om

        IF (cosz_gb(l) > 0.01) THEN
          k = g / cosz_gb(l)

          IF (l_niso_direct) THEN
            !Use the full non-isotropic expression for the albedo Asymmetry
            salb = (4.0 / 9.0) * (om-2.0 * alpl) / om
            !Albedo
            salb = 0.5 * om                                                    &
                   * ((1.0+3.0 * salb *  cosz_gb(l)**2)                        &
                      * (1.0 - cosz_gb(l)                                      &
                         * LOG((cosz_gb(l) + 1.0) / cosz_gb(l)))               &
                      - 1.5 * salb * cosz_gb(l) )
          ELSE
            !Original form directly from Table 2 of Sellers (1985);
            !this is valid only for isotropic scattering.
            salb = 0.5 * om *                                                  &
                   (1.0 - cosz_gb(l) * LOG((cosz_gb(l) + 1.0) / cosz_gb(l)))
          END IF
        END IF

      ELSE IF (orient(n) == 1) THEN
        sqcost = 1.0
        coszm  = 1.0
        !K = G / COSZ(I) = 1.0 because G = COSZ(I) for horizontal leaves
        k = 1.0
        IF (l_niso_direct) THEN
          !Full non-isotropic form.
          salb = alpl * 0.5
        ELSE
          !Original form directly from Table 2 of Sellers (1985);
          !this is valid only for isotropic scattering.
          salb = om * 0.25
        END IF
      END IF

      betadir = (1.0 + coszm * k) / (om * coszm * k) * salb
      c       = 0.5 * ( alpl + taul + (alpl - taul) * sqcost )
      betadif = c / om
      b       = 1.0 - (1.0 - betadif) * om
      d       = om * coszm * k * betadir
      f       = om * coszm * k * (1.0 - betadir)
      h       = SQRT(b * b - c * c) / coszm
      sig     = (coszm * k)**2 + c * c - b * b

      IF ( ABS(sig) < 1.0e-4 ) THEN
        sig = SIGN( 1.0e-4,sig )
      END IF

      u1      = b - c / albudif(l,band,n)
      ca      = c * albudir(l,band,n) / albudif(l,band,n)
      s1      = EXP(-h * lai(l,n))
      s2      = EXP(-k * lai(l,n))
      p1      = b + coszm * h
      p2      = b - coszm * h
      p3      = b + coszm * k
      p4      = b - coszm * k
      d1      = p1 * (u1 - coszm * h) / s1 - p2 * (u1 + coszm * h) * s1
      h1      = -d * p4 - c * f
      h2      =  ((d - p3 * h1 / sig) * (u1 - coszm *  h ) / s1                &
                  - (d - ca - (u1 + coszm * k) * h1 / sig) * p2 * s2) / d1
      h3      = -((d - p3 * h1 / sig) * (u1 + coszm * h) * s1                  &
                  - (d - ca - (u1 + coszm * k) * h1 / sig) * p1 * s2) / d1
      h7      =  (c / d1) * (u1 - coszm * h) / s1
      h8      = -(c / d1) * (u1 + coszm * h) * s1

      alb_type(l,n,2 * band-1) = h1 / sig + h2 + h3   ! Direct beam
      alb_type(l,n,2 * band)   = h7 + h8              ! Diffuse

      !-----------------------------------------------------------------------
      ! If required calculate the profile of absorbed PAR through the canopy
      ! of direct and diffuse beams (BEWARE: assumes PAR is band 1)
      !-----------------------------------------------------------------------
      IF ( l_getprofile .AND. band == 1 ) THEN

        u2  =  b - c * albudif(l,band,n)
        u3  =  f + c * albudif(l,band,n)
        d2  =  (u2 + coszm * h) / s1 - (u2 - coszm * h) * s1
        h4  = -f * p3 - c * d
        h5  = -1.0 / d2 * (h4 / sig * (u2 + coszm * h)                         &
               / s1 + (u3 - h4 / sig * (u2 - coszm * k)) * s2)
        h6  =  1.0 / d2 * (h4 / sig * (u2 - coszm * h)                         &
               * s1 + (u3 - h4 / sig * (u2 - coszm * k)) * s2)
        h9  =  1.0 / d2 * (u2 + coszm * h) / s1
        h10 = -1.0 / d2 * (u2 - coszm * h) * s1
        !-----------------------------------------------------------------------
        ! Two-stream equations for direct and diffuse upward and downward beams:
        !             RDIRU(I)=H1/SIG*EXP(-K*LA)+H2*EXP(-H*LA)+H3*EXP(H*LA)
        !             RDIRD(I)=(H4/SIG+1)*EXP(-K*LA)+H5*EXP(-H*LA)+H6*EXP(H*LA)
        !             RDIFU(I)=H7*EXP(-H*LA)+H8*EXP(H*LA)
        !             RDIFD(I)=H9*EXP(-H*LA)+H10*EXP(H*LA)
        ! Note: The original equations in Sellers (1985) include just the
        ! diffuse radiation scattered from the direct beam. The total downward
        ! flux is represented by RDIRD, with the extra 1 added to h4/sig to
        ! represent direct solar radiation.
        !-----------------------------------------------------------------------
        dlai = lai(l,n) / REAL(ilayers)

        IF ( can_rad_mod /= 5 .AND. can_rad_mod /= 6 .AND.                    &
             can_rad_mod /= 7 .AND. can_rad_mod /= 8 ) THEN
          !-----------------------------------------------------------------------
          ! Differentiate these equations to calculate PAR absorption per unit
          ! LAI down through the canopy. Centre derivatives in the centre of each
          ! LAI layer.
          !
          ! With these options fapar_dir represents the total PAR at each level
          ! normalized by the incident PAR in the case when the incident
          ! flux is entirely direct. fapar_dif is the normalized PAR at
          ! each level when the incident flux is entirely diffuse.
          !-----------------------------------------------------------------------
          la = 0.5 * dlai
          DO i = 1,ilayers

            !This is a more efficient version but will cause a change in KGO
            exp_hla       = EXP(h * la)
            exp_minus_hla = EXP(-h * la)  !We could 1/exp_hla too
            exp_minus_kla = EXP(-k * la)

            drdiru_dlai   = -k *  h1 / sig    * exp_minus_kla                  &
                             - h * h2 * exp_minus_hla + h * h3 * exp_hla
            drdird_dlai   = -k * (h4 / sig+1) * exp_minus_kla                  &
                             - h * h5 * exp_minus_hla + h * h6 * exp_hla
            drdifu_dlai   = -h * h7 * exp_minus_hla + h * h8  * exp_hla
            drdifd_dlai   = -h * h9 * exp_minus_hla + h * h10 * exp_hla

            fapar_dir(l,n,i) = (-drdird_dlai + drdiru_dlai) * can_struct_a(n)
            fapar_dif(l,n,i) = (-drdifd_dlai + drdifu_dlai) * can_struct_a(n)
            !WRITE(lineBuffer,*)'ZZZR2 ',l,n,i,fapar_dir(l,n,i),fapar_dif(l,n,i) 
            !CALL jules_print('albpft',lineBuffer)

            la = la + dlai
          END DO  !layers

        ELSE IF ( can_rad_mod == 5 .OR. can_rad_mod == 6 ) THEN

          !-----------------------------------------------------------------------
          !               can_rad_mod = 5 or 6.
          !               Don't use derivatives.
          !-----------------------------------------------------------------------

          la = dlai
          rnet_dir(0) = (h4 - h1) / sig + (h5 - h2) + (h6 - h3)
          rnet_dif(0) = (h9 - h7) + (h10 - h8)

          DO i = 1,ilayers

            exp_hla       = EXP(h * la)
            exp_minus_hla = EXP(-h * la)  !We could 1/exp_hla too
            exp_minus_kla = EXP(-k * la)
            exp_k_dlai_min_la = EXP(k * (dlai - la))

            rdiru = h1 / sig * exp_minus_kla + h2 * exp_minus_hla + h3 * exp_hla
            rdird = h4 / sig * exp_minus_kla + h5 * exp_minus_hla + h6 * exp_hla
            rdifu = h7 * exp_minus_hla + h8 * exp_hla
            rdifd = h9 * exp_minus_hla + h10 * exp_hla

            rnet_dir(i) = rdird - rdiru
            rnet_dif(i) = rdifd - rdifu
            WRITE(lineBuffer,*)'ZZZS1 ',i,rdird,rdiru,rnet_dir(i-1),rnet_dir(i) 
            CALL jules_print('albpft',lineBuffer)

            fapar_dir(l,n,i) = (rnet_dir(i-1) - rnet_dir(i)) / dlai
            fapar_dif(l,n,i) = (rnet_dif(i-1) - rnet_dif(i)) / dlai
            WRITE(lineBuffer,*)'ZZZS2 ',l,n,i,fapar_dir(l,n,i),fapar_dif(l,n,i) 
            CALL jules_print('albpft',lineBuffer)

            dabeer_dla(i) = (exp_k_dlai_min_la - exp_minus_kla) / dlai

            fapar_dir2dir(l,n,i) = (1.0 - om) * dabeer_dla(i)
            fapar_dir2dif(l,n,i) =  om * dabeer_dla(i) + fapar_dir(l,n,i)
            fapar_dif2dif(l,n,i) = fapar_dif(l,n,i)

            fsun(l,n,i) = (exp_k_dlai_min_la - exp_minus_kla) / (k * dlai)
            WRITE(lineBuffer,*)'ZZZS3 ',k,la,exp_minus_kla,i,dabeer_dla(i)
            CALL jules_print('albpft',lineBuffer)
            WRITE(lineBuffer,*)'ZZZS4 ',om,fapar_dir(l,n,i)
            CALL jules_print('albpft',lineBuffer)
            WRITE(lineBuffer,*)'ZZZS5 ',fapar_dir2dir(l,n,i),fapar_dir2dif(l,n,i)
            CALL jules_print('albpft',lineBuffer)
            WRITE(lineBuffer,*)'ZZZS6 ',fapar_dif(l,n,i),fsun(l,n,i)
            CALL jules_print('albpft',lineBuffer)
            !WRITE(lineBuffer,*)'ZZZS2 ',l,n,i,fapar_dir(l,n,i),fapar_dif(l,n,i) 
            !CALL jules_print('albpft',lineBuffer)
            !WRITE(lineBuffer,*)'ZZZS3 ',k,la,exp_minus_kla,i,dabeer_dla(i)
            !CALL jules_print('albpft',lineBuffer)
            !WRITE(lineBuffer,*)'ZZZS4 ',om,fapar_dir(l,n,i)
            !CALL jules_print('albpft',lineBuffer)
            !WRITE(lineBuffer,*)'ZZZS5 ',fapar_dir2dir(l,n,i),fapar_dir2dif(l,n,i)
            !CALL jules_print('albpft',lineBuffer)
            !WRITE(lineBuffer,*)'ZZZS6 ',fapar_dif(l,n,i)
            !CALL jules_print('albpft',lineBuffer)

            fapar_dir(l,n,i)     = fapar_dir(l,n,i) * can_struct_a(n)
            fapar_dif(l,n,i)     = fapar_dif(l,n,i) * can_struct_a(n)
            fapar_dir2dir(l,n,i) = fapar_dir2dir(l,n,i) * can_struct_a(n)
            fapar_dir2dif(l,n,i) = fapar_dir2dif(l,n,i) * can_struct_a(n)
            fapar_dif2dif(l,n,i) = fapar_dif2dif(l,n,i) * can_struct_a(n)

            la = la + dlai
          END DO  !layers
        ELSE IF ( can_rad_mod == 7 .OR. can_rad_mod == 8 ) THEN

          !Calculate the fluxes between each layer for the 
          !whole canopy under both DIRECT & DIFFUSE illumination
    
          !n.b., for DIRECT illumination: 
    
          !Ru = hemispheric reflectance arising from uncollided beam 
          !Rc = hemispheric reflectance arising from collided 
          !etc... 
    
          !work from lowest layer up using adding to
          !calculate values of reflectance and transmission        
          la = lai(l,n)
          lower_boundary_r = albudif(l,band,n)
          Ru_last = lower_boundary_r
          Rc_last = lower_boundary_r
          Rd_last = lower_boundary_r
          mu   = cosz_gb(l)

          !WRITE(lineBuffer,*)'BA0 = ',l,mu,lower_boundary_r,la
          !CALL jules_print('albpft',lineBuffer)

          !for layer in xrange(self.nLayers-1,-1,-1):
          !need to avoid division by 0; no else case (mu < 0.001) currently
          IF ( mu .GE. 0.001 ) THEN !setting this limit to 0.01 causes overflows with EXP(tau/mu)
           DO i = ilayers,1,-1

            !self.JULES_lad='uniform' !this is the setting currently used

            !Average inverse diffuse optical depth per unit leaf area
            !as calculated in the JULES implementation 
            !n.b. is unity in both cases
    
            !if self.JULES_lad=='uniform':
            muBar  = 1.0
            !elif self.JULES_lad=='horizontal':
            !mubar  = 1.0

            Gfunc  = 0.5 ! Gfunc=G(mu)   
            Zfunc  = 1.0 ! Zfunc=Z(mu)   
            Kfunc  = Gfunc/mu
            leaf_r = 0.1
            leaf_t = 0.1
            w      = leaf_r + leaf_t
            !def volssa_JULES(self):
            !The volume single scattering albedo as defined in JULES.     
            !if self.JULES_lad=='uniform':
            ssa=0.5*w*(1.0-mu*LOG((mu+1.0)/mu)) !this was also log in Python,
                                   ! which is natural log, same as in FORTRAN
            !elif self.JULES_lad=='horizontal':
            !ssa=w/4.

            !self.gamma_scaling_diffuse="delta"   !default setting for JULES  
            !self.gamma_scaling_direct="delta"    
            !  def getGamma(self,method="delta"):
            !Get the g1,g2,g3,g4 coefficients for the Meador and Weaver
            !two stream model that are consistent with the 
            !Sellers model.
    
            !B=self.B_diffuse()
            !def B_diffuse_JULES(self):
            !The Diffuse upscatter as calculated in JULES
            d     = leaf_r - leaf_t
            !if self.JULES_lad=='uniform':
            sqcost=1.0/3.0
            !elif self.JULES_lad=='horizontal':
            !sqcost=1.0
            B = 0.5*(w+d*sqcost)/w

            !B0=self.B_direct()
            !def B_direct_Dickinson(self,ssa):
            !Direct upscatter parameter as defined by Dickinson
            !both JULES and CLM use this formulation but differ
            !in the calculation of the single scattering albedo.

            B0 = (1.0/w)*ssa*(1.0+muBar*Kfunc)/(muBar*Kfunc)

            !WRITE(lineBuffer,*)'BA1 = ',i,leaf_r,leaf_t,ssa,d,B,B0
            !CALL jules_print('albpft',lineBuffer)
            !if method=="delta":   
            scalefunc = 1.0/(Gfunc*muBar)
            !elif method=="quad": 
            !   !modified quadrature
            !  scalefunc = SQRT(3.0)
            !else:
            !  raise TypeError("Unknown gamma scaling type:%s",method)
      
    
            g1 = (1.0 - (1.0 - B) * w) * scalefunc
            g2 = w * B * scalefunc
            g3 = B0
            g4 = 1.0 - g3
    

            !Ru,Tu=self.layers[layer].rtDirect()
            !Get reflectance and transmittance for
            !a single canopy layer under direct illumination
    
            !g1,g2,g3,g4=self.getGamma(method=self.gamma_scaling_direct)

            tau  = la*Gfunc*Zfunc
            tauD = tau
    
            a1 = g1*g4+g2*g3
            a2 = g1*g3+g2*g4
            k  = SQRT(g1*g1-g2*g2)
  
            D  = (1.0-k*k*mu*mu)*((k+g1)*EXP(k*tauD)+(k-g1)*EXP(-k*tauD))

            !reflectance:
            F11 = (1.0-k*mu)*(a2+k*g3)*EXP(k*tauD)
            G11 = (1.0+k*mu)*(a2-k*g3)*EXP(-k*tauD)
            H11 = 2.0*k*(g3-a2*mu)*EXP(-tau/mu)
      
            Ru = w/D*(F11-G11-H11)
            !{MPI Task 1} [INFO] albpft:  BA1 =            8  0.1000000      0.1000000 3.9199978E-02  0.0000000E+00  0.5000000      0.4538153
            !{MPI Task 1} [INFO] albpft:  BA1b =    1.788854      0.6576928 3.9999997E-03  3.9999997E-03  7.9999994E-03  0.5000000       1.000000 1.073895      0.5461847
    
            !WRITE(lineBuffer,*)'BA1b = ',k,mu,tau,tauD,la,Gfunc,Zfunc,        &
            !      a1,g4
            !CALL jules_print('albpft',lineBuffer)
            !transmittance
            F11 = (1.0+k*mu)*(a1+k*g4)*EXP(k*tauD)
            G11 = (1.0-k*mu)*(a1-k*g4)*EXP(-k*tauD)
            !WRITE(lineBuffer,*)'BA1c = ',F11,G11
            !CALL jules_print('albpft',lineBuffer)
            !WRITE(lineBuffer,*)'BA1c2 = ',tau/mu
            !CALL jules_print('albpft',lineBuffer)
            !WRITE(lineBuffer,*)'BA1c3 = ',EXP(1.1)
            !CALL jules_print('albpft',lineBuffer)
            !WRITE(lineBuffer,*)'BA1d = ',EXP(tau/mu)
            !CALL jules_print('albpft',lineBuffer)
            !H11 = 2.0*k*(g4+a1*mu)*EXP(tau/mu)
            H112 = 2.0*k*(g4+a1*mu)
            !WRITE(lineBuffer,*)'BA1e = ',H
            !CALL jules_print('albpft',lineBuffer)

            !Tu = EXP(-tau/mu)*(1.0-w/D*(F11-G11-H11)) 
            Tu = EXP(-tau/mu)*(1.0-w/D*(F11-G11))+w*H112/D

            !WRITE(lineBuffer,*)'BA2 = ',g1,g2,g3,g4,a1,a2,Ru,Tu
            !CALL jules_print('albpft',lineBuffer)
    
  
            !Get reflectance and transmittance for
            !   a single canopy layer under diffuse illumination
            !Rc,Tc=self.layers[layer].rtDiffuse()
            !with the default standard settings in the Python code, these
            !g1,g2,g3,g4 values are the same as previously calculated (above)
            !g1,g2,g3,g4=self.getGamma(method=self.gamma_scaling_diffuse)
            tau = la*Gfunc*Zfunc    
            !tau= la/muBar_uniformG_pintyZ_analytic()
            k   = SQRT(g1*g1-g2*g2)
            D   = k+g1+(k-g1)*EXP(-2.0*k*tau)
            Rc  = (g2*(1.0-EXP(-2.0*k*tau)))/D
            Tc  = (2.0*k*EXP(-k*tau))/D
    
            !U=self.layers[layer].tUncollidedDirect()
            !!U= EXP(-G(mu)*Z(mu)*lai/mu)
            !!self.JULES_lad='uniform'
            !!if self.JULES_lad=='uniform':
            !!return self.G_uniform(mu)
            !!return 0.5
            !!  self.Z=self.zeta_noStruct
            !! def zeta_noStruct(self,mu):
            !! return 1.0
            !!Given these settings, we hardwire G(mu)=0.5 and Z(mu)=1.0 for now:
            U           = EXP(-Gfunc*Zfunc*la/mu)

            !!This looks identical to Collided radiation
            !!R,T=self.layers[layer].rtDiffuse()
            !!g1,g2,g3,g4 are the same as previously calculated
            !!g1,g2,g3,g4=self.getGamma(method=self.gamma_scaling_diffuse)
            tau = la*Gfunc*Zfunc    
            !tau= la/muBar_uniformG_pintyZ_analytic()
            k   = SQRT(g1*g1-g2*g2)
            D   = k+g1+(k-g1)*EXP(-2.0*k*tau)
            Rd  = (g2*(1.0-EXP(-2.0*k*tau)))/D
            Td  = (2.0*k*EXP(-k*tau))/D
    
      
            Z = Rc_last * Rc  
            !WRITE(lineBuffer,*)'BA3 = ',tau,Rc,Tc,U,Rd,Td,Rc_last,Z
            !CALL jules_print('albpft',lineBuffer)
    
            !Collimated beam
            RuNew = Ru+U*Ru_last*Tc*(1.0+Z/(1.0-Z)) +                           &
                    (Tu-U)*Rc_last*Tc*(1.0+Z/(1.0-Z))
            TuNew = U+U*Ru_last*Rc*(1.+Z/(1.-Z)) + (Tu-U)*(1.+Z/(1.-Z))
      
            !Collided radiation
            RcNew = Rc+Tc*Rc_last*Tc*(1.0+Z/(1.0-Z))
            TcNew = Tc*(1.0+Z/(1.0-Z))
      
            Ru_last=RuNew
            Rc_last=RcNew

            Zdif=Rd_last*Rd  
            !Diffuse radiation 
            Rd=Rd+Td*Td*Rd_last*(1.+Zdif/(1.-Zdif))
            Td=Td*(1.+Zdif/(1.-Zdif))
            Rd_last=Rd

            !self.layers[layer].Iup_dir_u=RuNew
            !self.layers[layer].Iup_dir_c=RcNew
            !self.layers[layer].Idn_dir_u=TuNew
            !self.layers[layer].Idn_dir_c=TcNew
            !self.layers[layer].Iup_dif=Rd_last
            !self.layers[layer].Idn_dif=Td
            rdirup_u(i) = RuNew
            rdirup_c(i) = RcNew
            rdirdn_u(i) = TuNew
            rdirdn_c(i) = TcNew
            rdifup(i)   = Rd_last
            rdifdn(i)   = Td
            !WRITE(lineBuffer,*)'BB = ',i,RuNew,RcNew,TuNew,TcNew,Rd_last,Td
            !CALL jules_print('albpft',lineBuffer)

            la = la - dlai ! go down one layer
            !WRITE(lineBuffer,*)'ZZ la = ',la
            !CALL jules_print('albpft',lineBuffer)
           END DO

           !these alb_type albedos might need to be returned, they might need the
           !top-of-canopy values
           !alb_type(l,n,2 * band-1) =    ! Direct beam
           !alb_type(l,n,2 * band)   =    ! Diffuse
      


           !work back down through the layers to get
           !fluxes by normalising by total transmission
           la = dlai
           !WRITE(lineBuffer,*)'AA1 la = ',la
           !CALL jules_print('albpft',lineBuffer)
      
           Tdif = 1.0
           Tu   = 1.0
           Tc   = 0.0
           DO i = 1,ilayers
            Gfunc = 0.5 ! Gfunc=G(mu)   
            Zfunc = 1.0 ! Zfunc=Z(mu)   

            !self.layers[i].Iup_dir_u*=Tu
            !self.layers[i].Idn_dir_u*=Tu
            !self.layers[i].Iup_dir_c*=Tc
            !self.layers[i].Idn_dir_c*=Tc
            !self.layers[i].Iup_dir=self.layers[i].Iup_dir_c+self.layers[i].Iup_dir_u
            !self.layers[i].Idn_dir=self.layers[i].Idn_dir_c+self.layers[i].Idn_dir_u
            rdirup_u(i) = rdirup_u(i) * Tu
            rdirdn_u(i) = rdirdn_u(i) * Tu
            rdirup_c(i) = rdirup_c(i) * Tc
            rdirdn_c(i) = rdirdn_c(i) * Tc
            rdirup(i)   = rdirup_c(i) + rdirup_u(i) 
            rdirdn(i)   = rdirdn_c(i) + rdirdn_u(i) 
            !t=self.layers[layer].Idn_dif
            !self.layers[layer].Iup_dif*=T
            !self.layers[layer].Idn_dif*=T
            !T*=t
            !! FORTRAN is not case sensitive, so we redefine T as Tdif to
            !! distinguish T from t; we also change t to t2 for a longer
            !! variable name
            t2          = rdifdn(i)
            rdifup(i)   = rdifup(i) * Tdif
            rdifdn(i)   = rdifdn(i) * Tdif
            Tdif        = Tdif * t2
      
            !Keep track of how much uncollided
            !radiation there is incident on next
            !layer down:
            !U=self.layers[i].tUncollidedDirect()

            !! spelling out function: 
            !!U=np.exp(-self.G(self.mu)*self.Z(self.mu)*self.lai/self.mu)
            Gfunc = 0.5 ! Gfunc=G(mu)   
            Zfunc = 1.0 ! Zfunc=Z(mu)   
            U           = EXP(-Gfunc*Zfunc*la/mu)
            !Uncollided transmission of collimated radiation
            !Tu*=U
            Tu          = Tu * U

            !downward source term (note that Idn_dir_u 
            !has already been normalised by the incoming
            !collimated beam).
            !S=self.layers[i].Idn_dir_u-Tu           
            S           = rdirdn_u(i) - Tu           
            !Total collided radiation (again, note
            !that Idn_dir_c is normalised already):
            !Tc=self.layers[i].Idn_dir_c+S
            Tc          = rdirdn_c(i) + S

            la = la + dlai
            !WRITE(lineBuffer,*)'ZZ1 la = ',la
            !CALL jules_print('albpft',lineBuffer)
           END DO  !layers

           la = dlai
           !calculate absorption by solving energy balance
           DO i = 1,ilayers
            !if layer==0:
            !  Idn_above=1
            !else:
            !  Idn_above=self.layers[layer-1].Idn_dir
            !if layer==0:
            !  Idn_above=1
            !else:
            !  Idn_above=self.layers[layer-1].Idn_dif
            IF ( i==1 ) THEN
              rdirdn_above = 1.0
              rdifdn_above = 1.0
            ELSE
              rdirdn_above = rdirdn(i-1)
              rdifdn_above = rdifdn(i-1)
            END IF
            !if layer==(self.nLayers-1):
            !  Iup_below=self.layers[layer].Idn_dir*self.lower_boundary_r
            !else:
            !  Iup_below=self.layers[layer+1].Iup_dir
            !if layer==(self.nLayers-1):
            !  Iup_below=self.layers[layer].Idn_dif*self.lower_boundary_r
            !else:
            !  Iup_below=self.layers[layer+1].Iup_dif
            IF ( i==ilayers ) THEN
              rdirup_below = rdirdn(i)*lower_boundary_r
              rdifup_below = rdifdn(i)*lower_boundary_r
            ELSE
              rdirup_below = rdirup(i+1)
              rdifup_below = rdifup(i+1)
            END IF
            !Idn=self.layers[layer].Idn_dir
            !Iup=self.layers[layer].Iup_dir
            !self.layers[layer].Iab_dir=np.abs(Iup_below-Iup+Idn_above-Idn)
            !Idn=self.layers[layer].Idn_dif
            !Iup=self.layers[layer].Iup_dif
            !self.layers[layer].Iab_dif=np.abs(Iup_below-Iup+Idn_above-Idn)
            rdn = rdirdn(i)
            rup = rdirup(i)
            fapar_dir(l,n,i) = ABS(rdirup_below-rup+rdirdn_above-rdn)
            rnet_dir(i) = rdn - rup
            !WRITE(lineBuffer,*)'ZZZZ1 ',i,rdn,rup,rnet_dir(i-1),rnet_dir(i) 
            !CALL jules_print('albpft',lineBuffer)
            rdn = rdifdn(i)
            rup = rdifup(i)
            fapar_dif(l,n,i) = ABS(rdifup_below-rup+rdifdn_above-rdn)
            rnet_dif(i) = rdn - rup
            !WRITE(lineBuffer,*)'ZZZZ2 ',l,n,i,fapar_dir(l,n,i),fapar_dif(l,n,i) 
            !CALL jules_print('albpft',lineBuffer)

            !the following lines of code (6 equations) are from the can_rad_mod == 5 code above.
            !This assumes that can_struct_a == 1.0
            exp_minus_kla = EXP(-k * la)

            dabeer_dla(i) = (EXP(-k * (la - dlai)) - exp_minus_kla) / dlai
            !WRITE(lineBuffer,*)'ZZZZ3 ',k,la,exp_minus_kla,i,dabeer_dla(i)
            !CALL jules_print('albpft',lineBuffer)

            fapar_dir2dir(l,n,i) = (1.0 - om) * dabeer_dla(i)
            fapar_dir2dif(l,n,i) = om * dabeer_dla(i) + fapar_dir(l,n,i)
            fapar_dif2dif(l,n,i) = fapar_dif(l,n,i)
            !WRITE(lineBuffer,*)'ZZZZ4 ',om,fapar_dir(l,n,i)
            !CALL jules_print('albpft',lineBuffer)

            fsun(l,n,i) = exp_minus_kla * (EXP(k * dlai) - 1.0) / (k * dlai)
            !WRITE(lineBuffer,*)'ZZZZ5 ',fapar_dir2dir(l,n,i),fapar_dir2dif(l,n,i)
            !CALL jules_print('albpft',lineBuffer)
            !WRITE(lineBuffer,*)'ZZZZ6 ',fapar_dif(l,n,i)
            !CALL jules_print('albpft',lineBuffer)

            la = la + dlai

           END DO  !layers
          END IF
      
        END IF  !  can_rad_mod

      END IF  !  abs. par (band=1+getProfile)

    END DO  !  points
!$OMP END DO NOWAIT

  END DO  !  bands

END DO  !  pfts

!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE albpft
END MODULE albpft_mod
