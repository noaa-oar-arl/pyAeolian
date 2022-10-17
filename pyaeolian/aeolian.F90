subroutine fengsha(rho_phy,smois,ssm,xland,ust,clay,sand,rdrag,u_ts0,emis_dust)
  IMPLICIT NONE

  REAL, INTENT(OUT) :: emis_dust
  REAL, INTENT(IN) :: smois         ! volumetric soil moisture m3/m3
  REAL, INTENT(IN) :: ssm           ! sediment supply map
  REAL, INTENT(IN) :: xland,      & ! land=1 or water=0
       ust,        & ! friction velocity (m/s)
       clay,       & ! clay fractions
       sand,       & ! sand fraction
       rdrag,      & ! drag partition (1/m)
       u_ts0,       & ! dry threshold friction velocity (m/s)
       rho_phy       ! air density [kg/m3]

  REAL, PARAMETER :: cmb=1.0
  ! Local variables

  integer :: ilwi
  real, parameter :: g0 = 9.81 * 100 ! gravity
  real :: kvh
  real :: rhoa
  real :: u_ts
  real :: H
  real :: Q

  ! ! threshold values
  ! conver=1.e-9
  ! converi=1.e9

  ! Don't do dust over water!!!
  ilwi = 0
  if (xland.eq.1) then
     ilwi = 1
  end if
  ! dont do where ssm says not possible
  if (ssm > 0) then
     ilwi = 1
  end if

  IF (ilwi .eq. 0) THEN
     emis_dust = 0.
  ELSE
     rhoa = rho_phy * 1.0E-3
     ! soil moisture correction
     call fecan_moisture_correction(smois,sand,clay, H)

     ! vertical to horizontal mass flux
     call MB95_kvh(clay,kvh)

     ! modified threshold velocity
     call  modified_threshold(u_ts0, H, rdrag, u_ts)

     ! horizontal mass flux
     call fengsha_hflux(ust,u_ts, Q)

     emis_dust = cmb * ssm * rho_phy / g0 * kvh * Q

  end if

  return

end subroutine fengsha

subroutine GinouxDustEmission(radius, smois, xland, w10, rho_phy, emis_dust)

  IMPLICIT NONE

  ! output
  REAL, dimension(:), INTENT(OUT) :: emis_dust

  ! input
  REAL, dimension(:), INTENT(IN) :: radius
  REAL, INTENT(IN) :: smois         ! volumetric soil moisture m3/m3
  REAL, INTENT(IN) :: xland         ! land=1 or water=0
  REAL, INTENT(IN) :: w10           ! friction velocity (m/s)
  REAL, INTENT(IN) :: rho_phy       ! air density [kg/m3] (GOCART2G uses rho_phy=1.25)

  ! Local variables
  integer :: nbins, n
  real :: diameter, u_thresh0, u_thresh
  real, parameter :: g0 = 9.81 ! gravity
  real, parameter :: soil_density = 2650.  ! km m-3

  nbins = size(radius)

  do n = 1, nbins
    diameter = 2. * radius(n)
    u_thresh0 = 0.13 * sqrt(soil_density*g0*diameter/rho_phy) &
                       * sqrt(1.+6.e-7/(soil_density*g0*diameter**2.5)) &
              / sqrt(1.928*(1331.*(100.*diameter)**1.56+0.38)**0.092 - 1.)

    if ( xland /= 1 ) cycle ! done emit over water

    if(smois .lt. 0.5) then
        u_thresh = amax1(0.,u_thresh0* (1.2+0.2*alog10(max(1.e-3,smois))))

        if(w10 .gt. u_thresh) then
  !       Emission of dust [kg m-2 s-1]
          emis_dust(n) = w10**2. * (w10-u_thresh)
        endif
    endif


  end do
  return

end subroutine GinouxDustEmission

subroutine mackinnon_drag(z0,R)

  IMPLICIT NONE

  real, intent(in) :: z0
  real, intent(out) :: R
  real, parameter :: z0s = 1.0e-04 !Surface roughness for ideal bare surface [m]
  ! ------------------------------------------------------------------------
  ! Function: Calculates the MacKinnon et al. 2004 Drag Partition Correction
  !
  !   R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)
  !
  !--------------------------------------------------------------------------
  ! Drag partition correction. See MacKinnon et al. (2004),
  !     doi:10.1016/j.geomorph.2004.03.009
  R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)

  return

end subroutine mackinnon_drag


subroutine mb95_drag(z0,R)

  IMPLICIT NONE

  real, intent(in) :: z0
  real, intent(out) :: R
  real, parameter :: z0s = 1.0e-04 !Surface roughness for ideal bare surface [m]

  ! Drag partition correction. See MacKinnon et al. (2004),
  !     doi:10.1016/j.geomorph.2004.03.009
  !  R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)

  ! Drag partition correction. See Marticorena et al. (1997),
  !     doi:10.1029/96JD02964

  R = 1.0 - log(z0 / z0s) / log( 0.7 * (10./z0s) ** 0.8)

  return

end subroutine mb95_drag


subroutine fengsha_hflux(ust,utst, Q)
  !---------------------------------------------------------------------
  ! Function: Calculates the Horizontal Saltation Flux, Q, and then
  !           calculates the vertical flux.
  !
  ! formula of Draxler & Gillette (2001) Atmos. Environ.
  ! F   =  K A (r/g) U* ( U*^2 - Ut*^2 )
  !
  ! where:
  !     F   = vertical emission flux  [g/m**2-s]
  !     K   = constant 2.0E-04                      [1/m]
  !     A   = 0~3.5  mean = 2.8  (fudge factor)
  !     U*  = friction velocity                     [m/s]
  !     Ut* = threshold friction velocity           [m/s]
  !
  !--------------------------------------------------------------------
  real, intent(in) :: ust ! friction velocity
  real, intent(in) ::utst ! threshold friction velocity

  real, intent(out) :: Q
  Q = max(0.,ust * (ust * ust - utst * utst))

  return

end subroutine fengsha_hflux


subroutine MB95_kvh(clay,kvh)
  !---------------------------------------------------------------------
  ! Function: Calculates the vertical to horizontal mass flux ratio.
  !
  ! formula of MB95
  ! kvh = 10.0**(13.4*clay-6.0)
  !
  ! where:
  !     kvh   = vertical to hoizontal mass flux ratio [-]
  !     clay  = fractional clay content [-]
  !
  !--------------------------------------------------------------------
  real, intent(in) :: clay ! fractional clay content [-]

  real, intent(out) :: kvh
  if (clay <= 0.2) then
     kvh=10.0**(13.4*clay-6.0)
  else
     kvh = 2.E-4
  endif

  return

end subroutine MB95_kvh

subroutine fecan_moisture_correction(vol_soil_moisture,sand,clay, H)
  !---------------------------------------------------------------------
  ! Function: calculates the fecan soil moisture
  ! drylimit = 14.0*clay*clay+17.0*clay
  ! H = sqrt(1.0 + 1.21*(gravsm-drylimit)**0.68)
  !---------------------------------------------------------------------
  real, intent(in) :: vol_soil_moisture ! fractional clay content [-]
  real, intent(in) :: sand ! fractional sand content [-]
  real, intent(in) :: clay ! fractional clay content [-]
  real, parameter :: soil_dens = 2650. ! soil density [kg/m3]
  real, parameter :: water_dens = 1000. ! water density [kg/m3]
  real :: GRAVSM
  real :: drylimit
  real, intent(out) :: H ! fecan soil moisture adjustment

  H = 0.

  call volumetric_to_gravimetric(vol_soil_moisture,sand,gravsm)

  ! fecan dry limit
  drylimit=14.0*clay*clay+17.0*clay

  ! fecan soil moisture correction
  IF (gravsm > drylimit) THEN
     H = sqrt(1.0 + 1.21*(gravsm-drylimit)**0.68)
  ELSE
     H = 1.0
  END IF

  return

end subroutine fecan_moisture_correction

subroutine shao_1996_soil_moisture(w, H)

  Implicit None

  ! inputs
  real, intent(in) :: w ! volumetric soil moisture [m3/m3]

  !outputs
  real, intent(out) :: H

  H = 0.

  H = exp(22.7 * w)

  return

end subroutine shao_1996_soil_moisture


subroutine shao_2004_soil_moisture(w, H)

  Implicit None

  ! inputs
  real, intent(in) :: w ! volumetric soil moisture [m3/m3]

  !outputs
  real, intent(out) :: H

  H = 0.
  if (w <= 0.03) then
     H = exp(22.7 * w)
  else
     H = exp(95.3 * w - 2.029)
  end if

  return

end subroutine shao_2004_soil_moisture

subroutine fecan_dry_limit(clay,drylimit)

  IMPLICIT NONE

  !Inputs
  real, intent(in) :: clay ! fractional clay content
  !Outputs
  real, intent(out) :: drylimit ! fecan dry limit [kg/kg]

  drylimit = 0.
  if (clay <= 0) then
     drylimit = 14. * 1e-4 * 1e-4 + 17. * 1e-4
  else
     drylimit = 14 * clay * clay + 17 * clay
  end if

  return
end subroutine fecan_dry_limit

subroutine volumetric_to_gravimetric(vsoil, sandfrac, grav_soil)

  IMPLICIT NONE

  ! Inputs
  real, intent(in) :: vsoil ! Volumetric Soil Moisture [m3/m3]
  real, intent(in) :: sandfrac ! fractional Sand content
  ! outputs
  real, intent(out) :: grav_soil ! gravimetric soil moisture [kg/kg]
  ! Local
  real, parameter :: soil_dens = 2650.
  real, parameter :: water_dens = 1000.
  real :: vsat
  grav_soil = 0.
  vsat = 0.

  ! Saturated volumetric water content (sand-dependent) ! [m3 m-3]
  vsat = 0.489 - 0.00126 * ( sandfrac * 100 )
  ! gravimetric soil content
  grav_soil = vsoil * water_dens / (soil_dens * (1. - vsat))

  return
end subroutine volumetric_to_gravimetric

subroutine modified_threshold(u_ts0, H, drag, u_ts)

  IMPLICIT NONE

  real, intent(in) :: u_ts0 ! dry threshold velocity
  real, intent(in) :: H ! fecan soil moisture correction
  real, intent(in) :: drag ! drag partition
  real, intent(out) :: u_ts ! modified threshold velocity
  u_ts = 0.
  u_ts = u_ts0 * H / drag

  return
end subroutine modified_threshold

subroutine GinouxDryThreshold(radius, u_thresh0)

    Implicit NONE

    ! Input parameter
    real, dimension(:), intent(in)  :: radius      ! Dry particle bin effective radius [um]

    ! output parameters
    real, dimension(:), intent(out) :: u_thresh0

    ! local variables
    integer :: nbins, n
    real, parameter ::  air_dens = 1.25  ! Air density = 1.25 kg m-3
    real, parameter ::  soil_density  = 2650.  ! km m-3
    real            ::  diameter         ! dust effective diameter [m]
    real, parameter ::  grav = 9.81

    nbins = size(radius)

    do n = 1, nbins
      diameter = 2. * radius(n)

      u_thresh0(n) = 0.13 * sqrt(soil_density*grav*diameter/air_dens) &
                       * sqrt(1.+6.e-7/(soil_density*grav*diameter**2.5)) &
              / sqrt(1.928*(1331.*(100.*diameter)**1.56+0.38)**0.092 - 1.)
    end do

    return

  end subroutine GinouxDryThreshold




subroutine DustAerosolDistributionKok ( radius, rLow, rUp, distribution )

  ! !USES:
  implicit NONE

  ! !INPUT PARAMETERS:
  real, dimension(:), intent(in)  :: radius      ! Dry particle bin effective radius [um]
  real, dimension(:), intent(in)  :: rLow, rUp   ! Dry particle bin edge radii [um]

  ! !OUTPUT PARAMETERS:
  real, dimension(:), intent(out) :: distribution    ! Normalized dust aerosol distribution [1]

  ! !DESCRIPTION: Computes lognormal aerosol size distribution for dust bins according to
  !               J.F.Kok, PNAS, Jan 2011, 108 (3) 1016-1021; doi:10.1073/pnas.1014798108
  !
  ! !REVISION HISTORY:
  !
  ! 22Feb2020 B.Baker/NOAA    - Original implementation
  ! 01Apr2021 R.Montuoro/NOAA - Refactored for GOCART process library
  !

  ! !Local Variables
  integer :: n, nbins
  real    :: diameter, dlam, dvol

  ! !CONSTANTS
  real, parameter    :: mmd    = 3.4          ! median mass diameter [um]
  real, parameter    :: stddev = 3.0          ! geometric standard deviation [1]
  real, parameter    :: lambda = 12.0         ! crack propagation length [um]
  real, parameter    :: factor = 1.e0 / (sqrt(2.e0) * log(stddev))  ! auxiliary constant


  distribution = 0.

  nbins = size(radius)

  dvol = 0.
  do n = 1, nbins
    diameter = 2 * radius(n)
    dlam = diameter/lambda
    distribution(n) = diameter * (1. + erf(factor * log(diameter/mmd))) &
                    * exp(-dlam * dlam * dlam) * log(rUp(n)/rLow(n))
    dvol = dvol + distribution(n)
  end do

  !  Normalize distribution
  do n = 1, nbins
    distribution(n) = distribution(n) / dvol
  end do

  return
  end subroutine DustAerosolDistributionKok
