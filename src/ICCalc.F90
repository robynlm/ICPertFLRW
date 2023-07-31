#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine ICPertFLRW_ICCalc (CCTK_ARGUMENTS)
  use DefineRc
  use TensorCalc
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: i,j,k
  integer, parameter :: dp = 8
  real(dp), parameter :: pi = 4._dp*atan(1._dp)
  CCTK_REAL :: H0, a0, t0_EdS, Omega_lambda0
  CCTK_REAL :: t, aa, a2, ta2, mta2, kappa, Hprop
  CCTK_REAL :: Omega_matter, f_Loc, F, iFH, iFH2, m2iFH2
  CCTK_REAL :: rhoflrw

  CCTK_REAL :: Rc, dxRc, dyRc, dzRc
  CCTK_REAL :: dxdxRc, dxdyRc, dxdzRc, dydyRc, dydzRc, dzdzRc
  CCTK_REAL :: dxdxdxRc, dxdxdyRc, dxdxdzRc, dxdydyRc, dxdydzRc, dxdzdzRc
  CCTK_REAL :: dydydyRc, dydydzRc, dydzdzRc, dzdzdzRc
  CCTK_REAL :: dxdxdxdxRc, dxdxdxdyRc, dxdxdxdzRc, dxdxdydyRc, dxdxdydzRc
  CCTK_REAL :: dxdxdzdzRc, dxdydydyRc, dxdydydzRc, dxdydzdzRc, dxdzdzdzRc
  CCTK_REAL :: dydydydyRc, dydydydzRc, dydydzdzRc, dydzdzdzRc, dzdzdzdzRc

  CCTK_REAL :: dxgxx, dxgxy, dxgxz, dxgyy, dxgyz, dxgzz
  CCTK_REAL :: dygxx, dygxy, dygxz, dygyy, dygyz, dygzz
  CCTK_REAL :: dzgxx, dzgxy, dzgxz, dzgyy, dzgyz, dzgzz
  CCTK_REAL :: dxdxgxx, dxdxgxy, dxdxgxz, dxdxgyy, dxdxgyz, dxdxgzz
  CCTK_REAL :: dxdygxx, dxdygxy, dxdygxz, dxdygyy, dxdygyz, dxdygzz
  CCTK_REAL :: dxdzgxx, dxdzgxy, dxdzgxz, dxdzgyy, dxdzgyz, dxdzgzz
  CCTK_REAL :: dydygxx, dydygxy, dydygxz, dydygyy, dydygyz, dydygzz
  CCTK_REAL :: dydzgxx, dydzgxy, dydzgxz, dydzgyy, dydzgyz, dydzgzz
  CCTK_REAL :: dzdzgxx, dzdzgxy, dzdzgxz, dzdzgyy, dzdzgyz, dzdzgzz

  CCTK_REAL :: gdet, dxgdet, dygdet, dzgdet
  CCTK_REAL :: gxxu, gxyu, gxzu, gyyu, gyzu, gzzu
  CCTK_REAL :: dxgxxu, dxgxyu, dxgxzu, dxgyyu, dxgyzu, dxgzzu
  CCTK_REAL :: dygxxu, dygxyu, dygxzu, dygyyu, dygyzu, dygzzu
  CCTK_REAL :: dzgxxu, dzgxyu, dzgxzu, dzgyyu, dzgyzu, dzgzzu

  CCTK_REAL :: Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz
  CCTK_REAL :: Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz
  CCTK_REAL :: Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz
  CCTK_REAL :: dxGxxx, dxGxxy, dxGxxz, dxGxyy, dxGxyz, dxGxzz
  CCTK_REAL :: dxGyxx, dxGyxy, dxGyxz, dxGyyy, dxGyyz, dxGyzz
  CCTK_REAL :: dxGzxx, dxGzxy, dxGzxz, dxGzyy, dxGzyz, dxGzzz
  CCTK_REAL :: dyGxxx, dyGxxy, dyGxxz, dyGxyy, dyGxyz, dyGxzz
  CCTK_REAL :: dyGyxx, dyGyxy, dyGyxz, dyGyyy, dyGyyz, dyGyzz
  CCTK_REAL :: dyGzxx, dyGzxy, dyGzxz, dyGzyy, dyGzyz, dyGzzz
  CCTK_REAL :: dzGxxx, dzGxxy, dzGxxz, dzGxyy, dzGxyz, dzGxzz
  CCTK_REAL :: dzGyxx, dzGyxy, dzGyxz, dzGyyy, dzGyyz, dzGyzz
  CCTK_REAL :: dzGzxx, dzGzxy, dzGzxz, dzGzyy, dzGzyz, dzGzzz
  CCTK_REAL :: Gxxxu, Gxxyu, Gxxzu, Gxyyu, Gxyzu, Gxzzu
  CCTK_REAL :: Gyxxu, Gyxyu, Gyxzu, Gyyyu, Gyyzu, Gyzzu
  CCTK_REAL :: Gzxxu, Gzxyu, Gzxzu, Gzyyu, Gzyzu, Gzzzu
  CCTK_REAL :: dnGnxx, dnGnxy, dnGnxz, dnGnyy, dnGnyz, dnGnzz
  CCTK_REAL :: dxGnxn, dxGnyn, dxGnzn
  CCTK_REAL :: dyGnxn, dyGnyn, dyGnzn
  CCTK_REAL :: dzGnxn, dzGnyn, dzGnzn
  CCTK_REAL :: Gnxn, Gnyn, Gnzn, Gxnn, Gynn, Gznn

  CCTK_REAL :: RicciScalar, R1, R2, R3, R4

  CCTK_REAL :: kxxu, kxyu, kxzu, kyyu, kyzu, kzzu
  CCTK_REAL :: K_Loc, KijKji
  
  logical :: want_pflrw, want_Lambda

  want_pflrw = CCTK_EQUALS (my_initial_data, "ICPertFLRW")
  if (want_pflrw) then
    want_Lambda = CCTK_EQUALS (ICPertFLRW_Lambda, "yes")
    H0 = ICPertFLRW_h * ICPertFLRW_c / 2997.9_dp ! Units are Mpc
    a0 = 1._dp + ICPertFLRW_z_comoving_ref
    t0_EdS = 2._dp / ( 3._dp * H0 )
    Omega_lambda0 = 1._dp - ICPertFLRW_Omega_matter0

    t = cctk_time
    if (want_Lambda) then
      aa = a0 * ( ICPertFLRW_Omega_matter0 / Omega_lambda0 )**(1._dp/3._dp) &
           * sinh( sqrt(Omega_lambda0) * t / t0_EdS )**(2._dp/3._dp)
      Hprop = H0 * sqrt( ICPertFLRW_Omega_matter0 * ( aa / a0 )**(-3._dp) + Omega_lambda0 )
      Omega_matter = ICPertFLRW_Omega_matter0 / ( ICPertFLRW_Omega_matter0 &
                                                 + Omega_lambda0 * ( aa / a0 )**3._dp )
      Lambda = 3._dp * Omega_lambda0 * H0**2._dp / ICPertFLRW_c**2._dp
    else
      aa = a0 * ( t / t0_EdS )**(2._dp/3._dp)
      Hprop = 2._dp / ( 3._dp * t )
      Omega_matter = 1.0
      Lambda = 0._dp
    endif

    a2 = aa**2._dp
    ta2 = 2._dp * a2
    mta2 = - ta2

    f_Loc = Omega_matter**(6._dp/11._dp)
    F = f_Loc + (3._dp/2._dp) * Omega_matter
    iFH = 1._dp / ( F * Hprop )
    iFH2 = 1._dp / ( F * Hprop**2._dp )
    m2iFH2 = - 2._dp * iFH2

    kappa = 8._dp * pi * ICPertFLRW_G 
    rhoflrw = 3._dp * Hprop**2._dp * Omega_matter / kappa

    

    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

          call ICPertFLRW_DefineRc ( x(i,j,k), y(i,j,k), z(i,j,k), &
                                     Rc, dxRc, dyRc, dzRc, &
                                     dxdxRc, dxdyRc, dxdzRc, dydyRc, dydzRc, dzdzRc, &
                                     dxdxdxRc, dxdxdyRc, dxdxdzRc, dxdydyRc, dxdydzRc, dxdzdzRc, &
                                     dydydyRc, dydydzRc, dydzdzRc, dzdzdzRc, &
                                     dxdxdxdxRc, dxdxdxdyRc, dxdxdxdzRc, dxdxdydyRc, dxdxdydzRc, &
                                     dxdxdzdzRc, dxdydydyRc, dxdydydzRc, dxdydzdzRc, dxdzdzdzRc, &
                                     dydydydyRc, dydydydzRc, dydydzdzRc, dydzdzdzRc, dzdzdzdzRc )

          !============================================================================================
          ! Metric
          !============================================================================================

          gxx(i,j,k) = a2 * ( 1._dp - 2._dp * Rc ) + m2iFH2 * dxdxRc
          gxy(i,j,k) = - 2._dp * dxdyRc * iFH2
          gxz(i,j,k) = - 2._dp * dxdzRc * iFH2
          gyy(i,j,k) = a2 * ( 1._dp - 2._dp * Rc ) + m2iFH2 * dydyRc
          gyz(i,j,k) = - 2._dp * dydzRc * iFH2
          gzz(i,j,k) = a2 * ( 1._dp - 2._dp * Rc ) + m2iFH2 * dzdzRc

          ! Derivatives
          dxgxx = mta2 * dxRc + m2iFH2 * dxdxdxRc
          dygxx = mta2 * dyRc + m2iFH2 * dxdxdyRc
          dzgxx = mta2 * dzRc + m2iFH2 * dxdxdzRc
          dxgxy = m2iFH2 * dxdxdyRc
          dygxy = m2iFH2 * dxdydyRc
          dzgxy = m2iFH2 * dxdydzRc
          dxgxz = m2iFH2 * dxdxdzRc
          dygxz = m2iFH2 * dxdydzRc
          dzgxz = m2iFH2 * dxdzdzRc
          dxgyy = mta2 * dxRc + m2iFH2 * dxdydyRc
          dygyy = mta2 * dyRc + m2iFH2 * dydydyRc
          dzgyy = mta2 * dzRc + m2iFH2 * dydydzRc
          dxgyz = m2iFH2 * dxdydzRc
          dygyz = m2iFH2 * dydydzRc
          dzgyz = m2iFH2 * dydzdzRc
          dxgzz = mta2 * dxRc + m2iFH2 * dxdzdzRc
          dygzz = mta2 * dyRc + m2iFH2 * dydzdzRc
          dzgzz = mta2 * dzRc + m2iFH2 * dzdzdzRc

          dxdxgxx = mta2 * dxdxRc + m2iFH2 * dxdxdxdxRc
          dxdygxx = mta2 * dxdyRc + m2iFH2 * dxdxdxdyRc
          dxdzgxx = mta2 * dxdzRc + m2iFH2 * dxdxdxdzRc
          dxdxgxy = m2iFH2 * dxdxdxdyRc
          dxdygxy = m2iFH2 * dxdxdydyRc
          dxdzgxy = m2iFH2 * dxdxdydzRc
          dxdxgxz = m2iFH2 * dxdxdxdzRc
          dxdygxz = m2iFH2 * dxdxdydzRc
          dxdzgxz = m2iFH2 * dxdxdzdzRc
          dxdxgyy = mta2 * dxdxRc + m2iFH2 * dxdxdydyRc
          dxdygyy = mta2 * dxdyRc + m2iFH2 * dxdydydyRc
          dxdzgyy = mta2 * dxdzRc + m2iFH2 * dxdydydzRc
          dxdxgyz = m2iFH2 * dxdxdydzRc
          dxdygyz = m2iFH2 * dxdydydzRc
          dxdzgyz = m2iFH2 * dxdydzdzRc
          dxdxgzz = mta2 * dxdxRc + m2iFH2 * dxdxdzdzRc
          dxdygzz = mta2 * dxdyRc + m2iFH2 * dxdydzdzRc
          dxdzgzz = mta2 * dxdzRc + m2iFH2 * dxdzdzdzRc

          dydygxx = mta2 * dydyRc + m2iFH2 * dxdxdydyRc
          dydzgxx = mta2 * dydzRc + m2iFH2 * dxdxdydzRc
          dydygxy = m2iFH2 * dxdydydyRc
          dydzgxy = m2iFH2 * dxdydydzRc
          dydygxz = m2iFH2 * dxdydydzRc
          dydzgxz = m2iFH2 * dxdydzdzRc
          dydygyy = mta2 * dydyRc + m2iFH2 * dydydydyRc
          dydzgyy = mta2 * dydzRc + m2iFH2 * dydydydzRc
          dydygyz = m2iFH2 * dydydydzRc
          dydzgyz = m2iFH2 * dydydzdzRc
          dydygzz = mta2 * dydyRc + m2iFH2 * dydydzdzRc
          dydzgzz = mta2 * dydzRc + m2iFH2 * dydzdzdzRc

          dzdzgxx = mta2 * dzdzRc + m2iFH2 * dxdxdzdzRc
          dzdzgxy = m2iFH2 * dxdydzdzRc
          dzdzgxz = m2iFH2 * dxdzdzdzRc
          dzdzgyy = mta2 * dzdzRc + m2iFH2 * dydydzdzRc
          dzdzgyz = m2iFH2 * dydzdzdzRc
          dzdzgzz = mta2 * dzdzRc + m2iFH2 * dzdzdzdzRc

          ! Determinant
          gdet = gxx(i,j,k) * ( gyy(i,j,k)*gzz(i,j,k) - gyz(i,j,k)*gyz(i,j,k) ) &
               - gxy(i,j,k) * ( gxy(i,j,k)*gzz(i,j,k) - gxz(i,j,k)*gyz(i,j,k) ) &
               + gxz(i,j,k) * ( gxy(i,j,k)*gyz(i,j,k) - gxz(i,j,k)*gyy(i,j,k) )

          dxgdet = dxgxx * ( gyy(i,j,k)*gzz(i,j,k) - gyz(i,j,k)*gyz(i,j,k) ) &
                 - dxgxy * ( gxy(i,j,k)*gzz(i,j,k) - gxz(i,j,k)*gyz(i,j,k) ) &
                 + dxgxz * ( gxy(i,j,k)*gyz(i,j,k) - gxz(i,j,k)*gyy(i,j,k) ) &
                 + gxx(i,j,k) * ( dxgyy*gzz(i,j,k) - dxgyz*gyz(i,j,k) ) &
                 - gxy(i,j,k) * ( dxgxy*gzz(i,j,k) - dxgxz*gyz(i,j,k) ) &
                 + gxz(i,j,k) * ( dxgxy*gyz(i,j,k) - dxgxz*gyy(i,j,k) ) &
                 + gxx(i,j,k) * ( gyy(i,j,k)*dxgzz - gyz(i,j,k)*dxgyz ) &
                 - gxy(i,j,k) * ( gxy(i,j,k)*dxgzz - gxz(i,j,k)*dxgyz ) &
                 + gxz(i,j,k) * ( gxy(i,j,k)*dxgyz - gxz(i,j,k)*dxgyy )

          dygdet = dygxx * ( gyy(i,j,k)*gzz(i,j,k) - gyz(i,j,k)*gyz(i,j,k) ) &
                 - dygxy * ( gxy(i,j,k)*gzz(i,j,k) - gxz(i,j,k)*gyz(i,j,k) ) &
                 + dygxz * ( gxy(i,j,k)*gyz(i,j,k) - gxz(i,j,k)*gyy(i,j,k) ) &
                 + gxx(i,j,k) * ( dygyy*gzz(i,j,k) - dygyz*gyz(i,j,k) ) &
                 - gxy(i,j,k) * ( dygxy*gzz(i,j,k) - dygxz*gyz(i,j,k) ) &
                 + gxz(i,j,k) * ( dygxy*gyz(i,j,k) - dygxz*gyy(i,j,k) ) &
                 + gxx(i,j,k) * ( gyy(i,j,k)*dygzz - gyz(i,j,k)*dygyz ) &
                 - gxy(i,j,k) * ( gxy(i,j,k)*dygzz - gxz(i,j,k)*dygyz ) &
                 + gxz(i,j,k) * ( gxy(i,j,k)*dygyz - gxz(i,j,k)*dygyy )

          dzgdet = dzgxx * ( gyy(i,j,k)*gzz(i,j,k) - gyz(i,j,k)*gyz(i,j,k) ) &
                 - dzgxy * ( gxy(i,j,k)*gzz(i,j,k) - gxz(i,j,k)*gyz(i,j,k) ) &
                 + dzgxz * ( gxy(i,j,k)*gyz(i,j,k) - gxz(i,j,k)*gyy(i,j,k) ) &
                 + gxx(i,j,k) * ( dzgyy*gzz(i,j,k) - dzgyz*gyz(i,j,k) ) &
                 - gxy(i,j,k) * ( dzgxy*gzz(i,j,k) - dzgxz*gyz(i,j,k) ) &
                 + gxz(i,j,k) * ( dzgxy*gyz(i,j,k) - dzgxz*gyy(i,j,k) ) &
                 + gxx(i,j,k) * ( gyy(i,j,k)*dzgzz - gyz(i,j,k)*dzgyz ) &
                 - gxy(i,j,k) * ( gxy(i,j,k)*dzgzz - gxz(i,j,k)*dzgyz ) &
                 + gxz(i,j,k) * ( gxy(i,j,k)*dzgyz - gxz(i,j,k)*dzgyy )

          ! indices up
          gxxu = ( gyy(i,j,k)*gzz(i,j,k) - gyz(i,j,k)*gyz(i,j,k) ) / gdet
          gxyu = ( gyz(i,j,k)*gxz(i,j,k) - gxy(i,j,k)*gzz(i,j,k) ) / gdet
          gxzu = ( gxy(i,j,k)*gyz(i,j,k) - gyy(i,j,k)*gxz(i,j,k) ) / gdet
          gyyu = ( gxx(i,j,k)*gzz(i,j,k) - gxz(i,j,k)*gxz(i,j,k) ) / gdet
          gyzu = ( gxy(i,j,k)*gxz(i,j,k) - gxx(i,j,k)*gyz(i,j,k) ) / gdet
          gzzu = ( gxx(i,j,k)*gyy(i,j,k) - gxy(i,j,k)*gxy(i,j,k) ) / gdet

          dxgxxu = (( dxgyy*gzz(i,j,k) - dxgyz*gyz(i,j,k) &
                 + gyy(i,j,k)*dxgzz - gyz(i,j,k)*dxgyz )*gdet &
                 - ( gyy(i,j,k)*gzz(i,j,k) - gyz(i,j,k)*gyz(i,j,k) )*dxgdet) / (gdet*gdet)
          dxgxyu = (( dxgyz*gxz(i,j,k) - dxgxy*gzz(i,j,k) &
                 + gyz(i,j,k)*dxgxz - gxy(i,j,k)*dxgzz )*gdet &
                 - ( gyz(i,j,k)*gxz(i,j,k) - gxy(i,j,k)*gzz(i,j,k) )*dxgdet) / (gdet*gdet)
          dxgxzu = (( dxgxy*gyz(i,j,k) - dxgyy*gxz(i,j,k) &
                 + gxy(i,j,k)*dxgyz - gyy(i,j,k)*dxgxz )*gdet &
                 - ( gxy(i,j,k)*gyz(i,j,k) - gyy(i,j,k)*gxz(i,j,k) )*dxgdet) / (gdet*gdet)
          dxgyyu = (( dxgxx*gzz(i,j,k) - dxgxz*gxz(i,j,k) &
                 + gxx(i,j,k)*dxgzz - gxz(i,j,k)*dxgxz )*gdet &
                 - ( gxx(i,j,k)*gzz(i,j,k) - gxz(i,j,k)*gxz(i,j,k) )*dxgdet) / (gdet*gdet)
          dxgyzu = (( dxgxy*gxz(i,j,k) - dxgxx*gyz(i,j,k) &
                 + gxy(i,j,k)*dxgxz - gxx(i,j,k)*dxgyz )*gdet &
                 - ( gxy(i,j,k)*gxz(i,j,k) - gxx(i,j,k)*gyz(i,j,k) )*dxgdet) / (gdet*gdet)
          dxgzzu = (( dxgxx*gyy(i,j,k) - dxgxy*gxy(i,j,k) &
                 + gxx(i,j,k)*dxgyy - gxy(i,j,k)*dxgxy )*gdet &
                 - ( gxx(i,j,k)*gyy(i,j,k) - gxy(i,j,k)*gxy(i,j,k) )*dxgdet) / (gdet*gdet)

          dygxxu = (( dygyy*gzz(i,j,k) - dygyz*gyz(i,j,k) &
                 + gyy(i,j,k)*dygzz - gyz(i,j,k)*dygyz )*gdet &
                 - ( gyy(i,j,k)*gzz(i,j,k) - gyz(i,j,k)*gyz(i,j,k) )*dygdet) / (gdet*gdet)
          dygxyu = (( dygyz*gxz(i,j,k) - dygxy*gzz(i,j,k) &
                 + gyz(i,j,k)*dygxz - gxy(i,j,k)*dygzz )*gdet &
                 - ( gyz(i,j,k)*gxz(i,j,k) - gxy(i,j,k)*gzz(i,j,k) )*dygdet) / (gdet*gdet)
          dygxzu = (( dygxy*gyz(i,j,k) - dygyy*gxz(i,j,k) &
                 + gxy(i,j,k)*dygyz - gyy(i,j,k)*dygxz )*gdet &
                 - ( gxy(i,j,k)*gyz(i,j,k) - gyy(i,j,k)*gxz(i,j,k) )*dygdet) / (gdet*gdet)
          dygyyu = (( dygxx*gzz(i,j,k) - dygxz*gxz(i,j,k) &
                 + gxx(i,j,k)*dygzz - gxz(i,j,k)*dygxz )*gdet &
                 - ( gxx(i,j,k)*gzz(i,j,k) - gxz(i,j,k)*gxz(i,j,k) )*dygdet) / (gdet*gdet)
          dygyzu = (( dygxy*gxz(i,j,k) - dygxx*gyz(i,j,k) &
                 + gxy(i,j,k)*dygxz - gxx(i,j,k)*dygyz )*gdet &
                 - ( gxy(i,j,k)*gxz(i,j,k) - gxx(i,j,k)*gyz(i,j,k) )*dygdet) / (gdet*gdet)
          dygzzu = (( dygxx*gyy(i,j,k) - dygxy*gxy(i,j,k) &
                 + gxx(i,j,k)*dygyy - gxy(i,j,k)*dygxy )*gdet &
                 - ( gxx(i,j,k)*gyy(i,j,k) - gxy(i,j,k)*gxy(i,j,k) )*dygdet) / (gdet*gdet)

          dzgxxu = (( dzgyy*gzz(i,j,k) - dzgyz*gyz(i,j,k) &
                 + gyy(i,j,k)*dzgzz - gyz(i,j,k)*dzgyz )*gdet &
                 - ( gyy(i,j,k)*gzz(i,j,k) - gyz(i,j,k)*gyz(i,j,k) )*dzgdet) / (gdet*gdet)
          dzgxyu = (( dzgyz*gxz(i,j,k) - dzgxy*gzz(i,j,k) &
                 + gyz(i,j,k)*dzgxz - gxy(i,j,k)*dzgzz )*gdet &
                 - ( gyz(i,j,k)*gxz(i,j,k) - gxy(i,j,k)*gzz(i,j,k) )*dzgdet) / (gdet*gdet)
          dzgxzu = (( dzgxy*gyz(i,j,k) - dzgyy*gxz(i,j,k) &
                 + gxy(i,j,k)*dzgyz - gyy(i,j,k)*dzgxz )*gdet &
                 - ( gxy(i,j,k)*gyz(i,j,k) - gyy(i,j,k)*gxz(i,j,k) )*dzgdet) / (gdet*gdet)
          dzgyyu = (( dzgxx*gzz(i,j,k) - dzgxz*gxz(i,j,k) &
                 + gxx(i,j,k)*dzgzz - gxz(i,j,k)*dzgxz )*gdet &
                 - ( gxx(i,j,k)*gzz(i,j,k) - gxz(i,j,k)*gxz(i,j,k) )*dzgdet) / (gdet*gdet)
          dzgyzu = (( dzgxy*gxz(i,j,k) - dzgxx*gyz(i,j,k) &
                 + gxy(i,j,k)*dzgxz - gxx(i,j,k)*dzgyz )*gdet &
                 - ( gxy(i,j,k)*gxz(i,j,k) - gxx(i,j,k)*gyz(i,j,k) )*dzgdet) / (gdet*gdet)
          dzgzzu = (( dzgxx*gyy(i,j,k) - dzgxy*gxy(i,j,k) &
                 + gxx(i,j,k)*dzgyy - gxy(i,j,k)*dzgxy )*gdet &
                 - ( gxx(i,j,k)*gyy(i,j,k) - gxy(i,j,k)*gxy(i,j,k) )*dzgdet) / (gdet*gdet)

          !============================================================================================
          ! Christoffel symbols
          !============================================================================================

          Gxxx = (dxgxx + dxgxx - dxgxx) / 2._dp
          Gxxy = (dygxx + dxgxy - dxgxy) / 2._dp
          Gxxz = (dzgxx + dxgxz - dxgxz) / 2._dp
          Gxyy = (dygxy + dygxy - dxgyy) / 2._dp
          Gxyz = (dzgxy + dygxz - dxgyz) / 2._dp
          Gxzz = (dzgxz + dzgxz - dxgzz) / 2._dp

          Gyxx = (dxgxy + dxgxy - dygxx) / 2._dp
          Gyxy = (dygxy + dxgyy - dygxy) / 2._dp
          Gyxz = (dzgxy + dxgyz - dygxz) / 2._dp
          Gyyy = (dygyy + dygyy - dygyy) / 2._dp
          Gyyz = (dzgyy + dygyz - dygyz) / 2._dp
          Gyzz = (dzgyz + dzgyz - dygzz) / 2._dp

          Gzxx = (dxgxz + dxgxz - dzgxx) / 2._dp
          Gzxy = (dygxz + dxgyz - dzgxy) / 2._dp
          Gzxz = (dzgxz + dxgzz - dzgxz) / 2._dp
          Gzyy = (dygyz + dygyz - dzgyy) / 2._dp
          Gzyz = (dzgyz + dygzz - dzgyz) / 2._dp
          Gzzz = (dzgzz + dzgzz - dzgzz) / 2._dp

          ! Derivatives
          dxGxxx = (dxdxgxx + dxdxgxx - dxdxgxx) / 2._dp
          dxGxxy = (dxdygxx + dxdxgxy - dxdxgxy) / 2._dp
          dxGxxz = (dxdzgxx + dxdxgxz - dxdxgxz) / 2._dp
          dxGxyy = (dxdygxy + dxdygxy - dxdxgyy) / 2._dp
          dxGxyz = (dxdzgxy + dxdygxz - dxdxgyz) / 2._dp
          dxGxzz = (dxdzgxz + dxdzgxz - dxdxgzz) / 2._dp
          dxGyxx = (dxdxgxy + dxdxgxy - dxdygxx) / 2._dp
          dxGyxy = (dxdygxy + dxdxgyy - dxdygxy) / 2._dp
          dxGyxz = (dxdzgxy + dxdxgyz - dxdygxz) / 2._dp
          dxGyyy = (dxdygyy + dxdygyy - dxdygyy) / 2._dp
          dxGyyz = (dxdzgyy + dxdygyz - dxdygyz) / 2._dp
          dxGyzz = (dxdzgyz + dxdzgyz - dxdygzz) / 2._dp
          dxGzxx = (dxdxgxz + dxdxgxz - dxdzgxx) / 2._dp
          dxGzxy = (dxdygxz + dxdxgyz - dxdzgxy) / 2._dp
          dxGzxz = (dxdzgxz + dxdxgzz - dxdzgxz) / 2._dp
          dxGzyy = (dxdygyz + dxdygyz - dxdzgyy) / 2._dp
          dxGzyz = (dxdzgyz + dxdygzz - dxdzgyz) / 2._dp
          dxGzzz = (dxdzgzz + dxdzgzz - dxdzgzz) / 2._dp

          dyGxxx = (dxdygxx + dxdygxx - dxdygxx) / 2._dp
          dyGxxy = (dydygxx + dxdygxy - dxdygxy) / 2._dp
          dyGxxz = (dydzgxx + dxdygxz - dxdygxz) / 2._dp
          dyGxyy = (dydygxy + dydygxy - dxdygyy) / 2._dp
          dyGxyz = (dydzgxy + dydygxz - dxdygyz) / 2._dp
          dyGxzz = (dydzgxz + dydzgxz - dxdygzz) / 2._dp
          dyGyxx = (dxdygxy + dxdygxy - dydygxx) / 2._dp
          dyGyxy = (dydygxy + dxdygyy - dydygxy) / 2._dp
          dyGyxz = (dydzgxy + dxdygyz - dydygxz) / 2._dp
          dyGyyy = (dydygyy + dydygyy - dydygyy) / 2._dp
          dyGyyz = (dydzgyy + dydygyz - dydygyz) / 2._dp
          dyGyzz = (dydzgyz + dydzgyz - dydygzz) / 2._dp
          dyGzxx = (dxdygxz + dxdygxz - dydzgxx) / 2._dp
          dyGzxy = (dydygxz + dxdygyz - dydzgxy) / 2._dp
          dyGzxz = (dydzgxz + dxdygzz - dydzgxz) / 2._dp
          dyGzyy = (dydygyz + dydygyz - dydzgyy) / 2._dp
          dyGzyz = (dydzgyz + dydygzz - dydzgyz) / 2._dp
          dyGzzz = (dydzgzz + dydzgzz - dydzgzz) / 2._dp

          dzGxxx = (dxdzgxx + dxdzgxx - dxdzgxx) / 2._dp
          dzGxxy = (dydzgxx + dxdzgxy - dxdzgxy) / 2._dp
          dzGxxz = (dzdzgxx + dxdzgxz - dxdzgxz) / 2._dp
          dzGxyy = (dydzgxy + dydzgxy - dxdzgyy) / 2._dp
          dzGxyz = (dzdzgxy + dydzgxz - dxdzgyz) / 2._dp
          dzGxzz = (dzdzgxz + dzdzgxz - dxdzgzz) / 2._dp
          dzGyxx = (dxdzgxy + dxdzgxy - dydzgxx) / 2._dp
          dzGyxy = (dydzgxy + dxdzgyy - dydzgxy) / 2._dp
          dzGyxz = (dzdzgxy + dxdzgyz - dydzgxz) / 2._dp
          dzGyyy = (dydzgyy + dydzgyy - dydzgyy) / 2._dp
          dzGyyz = (dzdzgyy + dydzgyz - dydzgyz) / 2._dp
          dzGyzz = (dzdzgyz + dzdzgyz - dydzgzz) / 2._dp
          dzGzxx = (dxdzgxz + dxdzgxz - dzdzgxx) / 2._dp
          dzGzxy = (dydzgxz + dxdzgyz - dzdzgxy) / 2._dp
          dzGzxz = (dzdzgxz + dxdzgzz - dzdzgxz) / 2._dp
          dzGzyy = (dydzgyz + dydzgyz - dzdzgyy) / 2._dp
          dzGzyz = (dzdzgyz + dydzgzz - dzdzgyz) / 2._dp
          dzGzzz = (dzdzgzz + dzdzgzz - dzdzgzz) / 2._dp

          ! \partial_{l}(g^{ij} G_{ikj})
          dxGnxn = dxgxxu*Gxxx + gxxu*dxGxxx + dxgxyu*Gyxx + gxyu*dxGyxx + dxgxzu*Gzxx + gxzu*dxGzxx &
                 + dxgxyu*Gxxy + gxyu*dxGxxy + dxgyyu*Gyxy + gyyu*dxGyxy + dxgyzu*Gzxy + gyzu*dxGzxy &
                 + dxgxzu*Gxxz + gxzu*dxGxxz + dxgyzu*Gyxz + gyzu*dxGyxz + dxgzzu*Gzxz + gzzu*dxGzxz
          dxGnyn = dxgxxu*Gxxy + gxxu*dxGxxy + dxgxyu*Gyxy + gxyu*dxGyxy + dxgxzu*Gzxy + gxzu*dxGzxy &
                 + dxgxyu*Gxyy + gxyu*dxGxyy + dxgyyu*Gyyy + gyyu*dxGyyy + dxgyzu*Gzyy + gyzu*dxGzyy &
                 + dxgxzu*Gxyz + gxzu*dxGxyz + dxgyzu*Gyyz + gyzu*dxGyyz + dxgzzu*Gzyz + gzzu*dxGzyz
          dxGnzn = dxgxxu*Gxxz + gxxu*dxGxxz + dxgxyu*Gyxz + gxyu*dxGyxz + dxgxzu*Gzxz + gxzu*dxGzxz &
                 + dxgxyu*Gxyz + gxyu*dxGxyz + dxgyyu*Gyyz + gyyu*dxGyyz + dxgyzu*Gzyz + gyzu*dxGzyz &
                 + dxgxzu*Gxzz + gxzu*dxGxzz + dxgyzu*Gyzz + gyzu*dxGyzz + dxgzzu*Gzzz + gzzu*dxGzzz
          dyGnxn = dygxxu*Gxxx + gxxu*dyGxxx + dygxyu*Gyxx + gxyu*dyGyxx + dygxzu*Gzxx + gxzu*dyGzxx &
                 + dygxyu*Gxxy + gxyu*dyGxxy + dygyyu*Gyxy + gyyu*dyGyxy + dygyzu*Gzxy + gyzu*dyGzxy &
                 + dygxzu*Gxxz + gxzu*dyGxxz + dygyzu*Gyxz + gyzu*dyGyxz + dygzzu*Gzxz + gzzu*dyGzxz
          dyGnyn = dygxxu*Gxxy + gxxu*dyGxxy + dygxyu*Gyxy + gxyu*dyGyxy + dygxzu*Gzxy + gxzu*dyGzxy &
                 + dygxyu*Gxyy + gxyu*dyGxyy + dygyyu*Gyyy + gyyu*dyGyyy + dygyzu*Gzyy + gyzu*dyGzyy &
                 + dygxzu*Gxyz + gxzu*dyGxyz + dygyzu*Gyyz + gyzu*dyGyyz + dygzzu*Gzyz + gzzu*dyGzyz
          dyGnzn = dygxxu*Gxxz + gxxu*dyGxxz + dygxyu*Gyxz + gxyu*dyGyxz + dygxzu*Gzxz + gxzu*dyGzxz &
                 + dygxyu*Gxyz + gxyu*dyGxyz + dygyyu*Gyyz + gyyu*dyGyyz + dygyzu*Gzyz + gyzu*dyGzyz &
                 + dygxzu*Gxzz + gxzu*dyGxzz + dygyzu*Gyzz + gyzu*dyGyzz + dygzzu*Gzzz + gzzu*dyGzzz
          dzGnxn = dzgxxu*Gxxx + gxxu*dzGxxx + dzgxyu*Gyxx + gxyu*dzGyxx + dzgxzu*Gzxx + gxzu*dzGzxx &
                 + dzgxyu*Gxxy + gxyu*dzGxxy + dzgyyu*Gyxy + gyyu*dzGyxy + dzgyzu*Gzxy + gyzu*dzGzxy &
                 + dzgxzu*Gxxz + gxzu*dzGxxz + dzgyzu*Gyxz + gyzu*dzGyxz + dzgzzu*Gzxz + gzzu*dzGzxz
          dzGnyn = dzgxxu*Gxxy + gxxu*dzGxxy + dzgxyu*Gyxy + gxyu*dzGyxy + dzgxzu*Gzxy + gxzu*dzGzxy &
                 + dzgxyu*Gxyy + gxyu*dzGxyy + dzgyyu*Gyyy + gyyu*dzGyyy + dzgyzu*Gzyy + gyzu*dzGzyy &
                 + dzgxzu*Gxyz + gxzu*dzGxyz + dzgyzu*Gyyz + gyzu*dzGyyz + dzgzzu*Gzyz + gzzu*dzGzyz
          dzGnzn = dzgxxu*Gxxz + gxxu*dzGxxz + dzgxyu*Gyxz + gxyu*dzGyxz + dzgxzu*Gzxz + gxzu*dzGzxz &
                 + dzgxyu*Gxyz + gxyu*dzGxyz + dzgyyu*Gyyz + gyyu*dzGyyz + dzgyzu*Gzyz + gyzu*dzGzyz &
                 + dzgxzu*Gxzz + gxzu*dzGxzz + dzgyzu*Gyzz + gyzu*dzGyzz + dzgzzu*Gzzz + gzzu*dzGzzz

          !============================================================================================
          ! Ricci Scalar
          !============================================================================================

          ! R1
          call ICPertFLRW_derivative_of_Gudd ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                               dxgxxu, dxgxyu, dxgxzu, &
                                               dygxyu, dygyyu, dygyzu, dzgxzu, dzgyzu, dzgzzu, &
                                               Gxxx, Gyxx, Gzxx, dxGxxx, dxGyxx, dxGzxx, &
                                               dyGxxx, dyGyxx, dyGzxx, dzGxxx, dzGyxx, dzGzxx, &
                                               dnGnxx )
          call ICPertFLRW_derivative_of_Gudd ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                               dxgxxu, dxgxyu, dxgxzu, &
                                               dygxyu, dygyyu, dygyzu, dzgxzu, dzgyzu, dzgzzu, &
                                               Gxxy, Gyxy, Gzxy, dxGxxy, dxGyxy, dxGzxy, &
                                               dyGxxy, dyGyxy, dyGzxy, dzGxxy, dzGyxy, dzGzxy, &
                                               dnGnxy )
          call ICPertFLRW_derivative_of_Gudd ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                               dxgxxu, dxgxyu, dxgxzu, &
                                               dygxyu, dygyyu, dygyzu, dzgxzu, dzgyzu, dzgzzu, &
                                               Gxxz, Gyxz, Gzxz, dxGxxz, dxGyxz, dxGzxz, &
                                               dyGxxz, dyGyxz, dyGzxz, dzGxxz, dzGyxz, dzGzxz, &
                                               dnGnxz )
          call ICPertFLRW_derivative_of_Gudd ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                               dxgxxu, dxgxyu, dxgxzu, &
                                               dygxyu, dygyyu, dygyzu, dzgxzu, dzgyzu, dzgzzu, &
                                               Gxyy, Gyyy, Gzyy, dxGxyy, dxGyyy, dxGzyy, &
                                               dyGxyy, dyGyyy, dyGzyy, dzGxyy, dzGyyy, dzGzyy, &
                                               dnGnyy )
          call ICPertFLRW_derivative_of_Gudd ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                               dxgxxu, dxgxyu, dxgxzu, &
                                               dygxyu, dygyyu, dygyzu, dzgxzu, dzgyzu, dzgzzu, &
                                               Gxyz, Gyyz, Gzyz, dxGxyz, dxGyyz, dxGzyz, &
                                               dyGxyz, dyGyyz, dyGzyz, dzGxyz, dzGyyz, dzGzyz, &
                                               dnGnyz )
          call ICPertFLRW_derivative_of_Gudd ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                               dxgxxu, dxgxyu, dxgxzu, &
                                               dygxyu, dygyyu, dygyzu, dzgxzu, dzgyzu, dzgzzu, &
                                               Gxzz, Gyzz, Gzzz, dxGxzz, dxGyzz, dxGzzz, &
                                               dyGxzz, dyGyzz, dyGzzz, dzGxzz, dzGyzz, dzGzzz, &
                                               dnGnzz )
          call ICPertFLRW_take_trace ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                       dnGnxx, dnGnxy, dnGnxz, dnGnyy, dnGnyz, dnGnzz, R1 )

          ! R2
          call ICPertFLRW_take_trace_notsym ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                              dxGnxn, dxGnyn, dxGnzn, &
                                              dyGnxn, dyGnyn, dyGnzn, &
                                              dzGnxn, dzGnyn, dzGnzn, R2 )

          ! R3
          call ICPertFLRW_Christoffel_contraction_one ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                                        Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz, &
                                                        Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz, &
                                                        Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz, &
                                                        Gnxn, Gnyn, Gnzn )
          call ICPertFLRW_Christoffel_contraction_two ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                                        Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz, &
                                                        Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz, &
                                                        Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz, &
                                                        Gxnn, Gynn, Gznn )
          call ICPertFLRW_product_of_two_vectors ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                                   Gnxn, Gnyn, Gnzn, &
                                                   Gxnn, Gynn, Gznn, &
                                                   R3 )

          ! R4
          call ICPertFLRW_raise_three_indices ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                                Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz, &
                                                Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz, &
                                                Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz, &
                                                Gxxxu, Gxxyu, Gxxzu, Gxyyu, Gxyzu, Gxzzu, &
                                                Gyxxu, Gyxyu, Gyxzu, Gyyyu, Gyyzu, Gyzzu, &
                                                Gzxxu, Gzxyu, Gzxzu, Gzyyu, Gzyzu, Gzzzu )
          call ICPertFLRW_Christoffel_contraction_three ( Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz, &
                                                          Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz, &
                                                          Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz, &
                                                          Gxxxu, Gxxyu, Gxxzu, Gxyyu, Gxyzu, Gxzzu, &
                                                          Gyxxu, Gyxyu, Gyxzu, Gyyyu, Gyyzu, Gyzzu, &
                                                          Gzxxu, Gzxyu, Gzxzu, Gzyyu, Gzyzu, Gzzzu, &
                                                          R4 )

          RicciScalar = R1 - R2 + R3 - R4

          !============================================================================================
          ! Extrinsic Curvature
          !============================================================================================

          kxx(i,j,k) = - a2 * Hprop * ( 1._dp - 2._dp * Rc ) + ( 2._dp + f_Loc ) * dxdxRc * iFH
          kxy(i,j,k) = ( 2._dp + f_Loc ) * dxdyRc * iFH
          kxz(i,j,k) = ( 2._dp + f_Loc ) * dxdzRc * iFH
          kyy(i,j,k) = - a2 * Hprop * ( 1._dp - 2._dp * Rc ) + ( 2._dp + f_Loc ) * dydyRc * iFH
          kyz(i,j,k) = ( 2._dp + f_Loc ) * dydzRc * iFH
          kzz(i,j,k) = - a2 * Hprop * ( 1._dp - 2._dp * Rc ) + ( 2._dp + f_Loc ) * dzdzRc * iFH

          call ICPertFLRW_take_trace ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                       kxx(i,j,k), kxy(i,j,k), kxz(i,j,k), &
                                       kyy(i,j,k), kyz(i,j,k), kzz(i,j,k), K_Loc )

          call ICPertFLRW_raise_two_indices ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                              kxx(i,j,k), kxy(i,j,k), kxz(i,j,k), &
                                              kyy(i,j,k), kyz(i,j,k), kzz(i,j,k), &
                                              kxxu, kxyu, kxzu, kyyu, kyzu, kzzu )

          call ICPertFLRW_contract_over_two_indices ( kxx(i,j,k), kxy(i,j,k), kxz(i,j,k), &
                                                      kyy(i,j,k), kyz(i,j,k), kzz(i,j,k), &
                                                      kxxu, kxyu, kxzu, kyyu, kyzu, kzzu, &
                                                      KijKji )

          !============================================================================================
          ! Density
          !============================================================================================

          rho(i,j,k) = (RicciScalar + K_Loc**2._dp - KijKji - 2._dp*Lambda) / (2._dp * kappa)

          rhodp(i,j,k) = rho(i,j,k) * gdet ** (1._dp/2._dp)
          prs(i,j,k) = 0._dp
          eps(i,j,k) = 0._dp
          u1(i,j,k) = 0._dp
          u2(i,j,k) = 0._dp
          u3(i,j,k) = 0._dp
           
        enddo  !i
      enddo  !j
    enddo  !k
  endif  !if you want pflrw
end subroutine ICPertFLRW_ICCalc




