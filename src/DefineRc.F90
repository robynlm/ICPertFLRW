#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module DefineRc

  implicit none
  contains

    subroutine ICPertFLRW_DefineRc ( x, y, z, &
                                     Rc, dxRc, dyRc, dzRc, &
                                     dxdxRc, dxdyRc, dxdzRc, dydyRc, dydzRc, dzdzRc, &
                                     dxdxdxRc, dxdxdyRc, dxdxdzRc, dxdydyRc, dxdydzRc, dxdzdzRc, &
                                     dydydyRc, dydydzRc, dydzdzRc, dzdzdzRc, &
                                     dxdxdxdxRc, dxdxdxdyRc, dxdxdxdzRc, dxdxdydyRc, dxdxdydzRc, &
                                     dxdxdzdzRc, dxdydydyRc, dxdydydzRc, dxdydzdzRc, dxdzdzdzRc, &
                                     dydydydyRc, dydydydzRc, dydydzdzRc, dydzdzdzRc, dzdzdzdzRc )
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS
      ! Input
      CCTK_REAL :: x, y, z

      ! Output
      CCTK_REAL :: Rc, dxRc, dyRc, dzRc
      CCTK_REAL :: dxdxRc, dxdyRc, dxdzRc, dydyRc, dydzRc, dzdzRc
      CCTK_REAL :: dxdxdxRc, dxdxdyRc, dxdxdzRc, dxdydyRc, dxdydzRc, dxdzdzRc
      CCTK_REAL :: dydydyRc, dydydzRc, dydzdzRc, dzdzdzRc
      CCTK_REAL :: dxdxdxdxRc, dxdxdxdyRc, dxdxdxdzRc, dxdxdydyRc, dxdxdydzRc
      CCTK_REAL :: dxdxdzdzRc, dxdydydyRc, dxdydydzRc, dxdydzdzRc, dxdzdzdzRc
      CCTK_REAL :: dydydydyRc, dydydydzRc, dydydzdzRc, dydzdzdzRc, dzdzdzdzRc

      ! Local variables
      integer, parameter :: dp = 8
      integer :: m
      real(dp), parameter :: pi = 4._dp*atan(1._dp)
      CCTK_REAL :: twopi, kx, ky, kz
      CCTK_REAL :: sinx, siny, sinz
      CCTK_REAL :: cosx, cosy, cosz
      twopi = 2._dp * pi

      ! Initialise them to zero
      Rc = 0._dp
      ! Single derivatives
      dxRc = 0._dp
      dyRc = 0._dp
      dzRc = 0._dp
      ! Double derivatives
      dxdxRc = 0._dp
      dxdyRc = 0._dp
      dxdzRc = 0._dp
      dydyRc = 0._dp
      dydzRc = 0._dp
      dzdzRc = 0._dp
      ! Triple derivatives
      dxdxdxRc = 0._dp
      dxdxdyRc = 0._dp
      dxdxdzRc = 0._dp
      dxdydyRc = 0._dp
      dxdydzRc = 0._dp
      dxdzdzRc = 0._dp
      dydydyRc = 0._dp
      dydydzRc = 0._dp
      dydzdzRc = 0._dp
      dzdzdzRc = 0._dp
      ! Quadruple derivatives
      dxdxdxdxRc = 0._dp
      dxdxdxdyRc = 0._dp
      dxdxdxdzRc = 0._dp
      dxdxdydyRc = 0._dp
      dxdxdydzRc = 0._dp
      dxdxdzdzRc = 0._dp
      dxdydydyRc = 0._dp
      dxdydydzRc = 0._dp
      dxdydzdzRc = 0._dp
      dxdzdzdzRc = 0._dp
      dydydydyRc = 0._dp
      dydydydzRc = 0._dp
      dydydzdzRc = 0._dp
      dydzdzdzRc = 0._dp
      dzdzdzdzRc = 0._dp

      do m = 1, 20        
        kx = twopi / ICPertFLRW_lambda_x(m)
        ky = twopi / ICPertFLRW_lambda_y(m)
        kz = twopi / ICPertFLRW_lambda_z(m)

        sinx = sin(x * kx + ICPertFLRW_phi_x(m))
        siny = sin(y * ky + ICPertFLRW_phi_y(m))
        sinz = sin(z * kz + ICPertFLRW_phi_z(m))
        cosx = cos(x * kx + ICPertFLRW_phi_x(m))
        cosy = cos(y * ky + ICPertFLRW_phi_y(m))
        cosz = cos(z * kz + ICPertFLRW_phi_z(m))

        Rc = Rc + ICPertFLRW_Amp_x(m) * sinx &
                + ICPertFLRW_Amp_y(m) * siny &
                + ICPertFLRW_Amp_z(m) * sinz

        ! Single derivatives
        dxRc = dxRc + ICPertFLRW_Amp_x(m) * kx * cosx
        dyRc = dyRc + ICPertFLRW_Amp_y(m) * ky * cosy
        dzRc = dzRc + ICPertFLRW_Amp_z(m) * kz * cosz

        ! Double derivatives
        dxdxRc = dxdxRc - ICPertFLRW_Amp_x(m) * (kx**2._dp) * sinx
        dydyRc = dydyRc - ICPertFLRW_Amp_y(m) * (ky**2._dp) * siny
        dzdzRc = dzdzRc - ICPertFLRW_Amp_z(m) * (kz**2._dp) * sinz

        ! Triple derivatives
        dxdxdxRc = dxdxdxRc - ICPertFLRW_Amp_x(m) * (kx**3._dp) * cosx
        dydydyRc = dydydyRc - ICPertFLRW_Amp_y(m) * (ky**3._dp) * cosy
        dzdzdzRc = dzdzdzRc - ICPertFLRW_Amp_z(m) * (kz**3._dp) * cosz

        ! Quadruple derivatives
        dxdxdxdxRc = dxdxdxdxRc + ICPertFLRW_Amp_x(m) * (kx**4._dp) * sinx
        dydydydyRc = dydydydyRc + ICPertFLRW_Amp_y(m) * (ky**4._dp) * siny
        dzdzdzdzRc = dzdzdzdzRc + ICPertFLRW_Amp_z(m) * (kz**4._dp) * sinz
      enddo
  end subroutine ICPertFLRW_DefineRc
end module DefineRc
