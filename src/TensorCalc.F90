#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module TensorCalc

  implicit none
  contains

    ! The following subroutine computes A^{i}_{j} = g^{ik} A_{kj}
    subroutine ICPertFLRW_raise_one_indice ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                             Axx, Axy, Axz, Ayy, Ayz, Azz, &
                                             Axxm, Axym, Axzm, Ayym, Ayzm, Azzm )
      ! Input
      CCTK_REAL :: gxxu, gxyu, gxzu, gyyu, gyzu, gzzu
      CCTK_REAL :: Axx, Axy, Axz, Ayy, Ayz, Azz

      ! Output
      CCTK_REAL :: Axxm, Axym, Axzm, Ayym, Ayzm, Azzm

      Axxm = gxxu*Axx + gxyu*Axy + gxzu*Axz
      Axym = gxxu*Axy + gxyu*Ayy + gxzu*Ayz
      Axzm = gxxu*Axz + gxyu*Ayz + gxzu*Azz
      Ayym = gxyu*Axy + gyyu*Ayy + gyzu*Ayz
      Ayzm = gxyu*Axz + gyyu*Ayz + gyzu*Azz
      Azzm = gxzu*Axz + gyzu*Ayz + gzzu*Azz
      
    end subroutine ICPertFLRW_raise_one_indice

    ! The following subroutine computes A^{ij} = g^{ik} g^{jl} A_{kl}
    subroutine ICPertFLRW_raise_two_indices ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                              Axx, Axy, Axz, Ayy, Ayz, Azz, &
                                              Axxu, Axyu, Axzu, Ayyu, Ayzu, Azzu )
      ! Input
      CCTK_REAL :: gxxu, gxyu, gxzu, gyyu, gyzu, gzzu
      CCTK_REAL :: Axx, Axy, Axz, Ayy, Ayz, Azz

      ! Output
      CCTK_REAL :: Axxu, Axyu, Axzu, Ayyu, Ayzu, Azzu

      Axxu = gxxu*gxxu*Axx + gxxu*gxyu*Axy + gxxu*gxzu*Axz &
           + gxyu*gxxu*Axy + gxyu*gxyu*Ayy + gxyu*gxzu*Ayz &
           + gxzu*gxxu*Axz + gxzu*gxyu*Ayz + gxzu*gxzu*Azz
      Axyu = gxxu*gxyu*Axx + gxxu*gyyu*Axy + gxxu*gyzu*Axz &
           + gxyu*gxyu*Axy + gxyu*gyyu*Ayy + gxyu*gyzu*Ayz &
           + gxzu*gxyu*Axz + gxzu*gyyu*Ayz + gxzu*gyzu*Azz
      Axzu = gxxu*gxzu*Axx + gxxu*gyzu*Axy + gxxu*gzzu*Axz &
           + gxyu*gxzu*Axy + gxyu*gyzu*Ayy + gxyu*gzzu*Ayz &
           + gxzu*gxzu*Axz + gxzu*gyzu*Ayz + gxzu*gzzu*Azz

      Ayyu = gxyu*gxyu*Axx + gxyu*gyyu*Axy + gxyu*gyzu*Axz &
           + gyyu*gxyu*Axy + gyyu*gyyu*Ayy + gyyu*gyzu*Ayz &
           + gyzu*gxyu*Axz + gyzu*gyyu*Ayz + gyzu*gyzu*Azz
      Ayzu = gxyu*gxzu*Axx + gxyu*gyzu*Axy + gxyu*gzzu*Axz &
           + gyyu*gxzu*Axy + gyyu*gyzu*Ayy + gyyu*gzzu*Ayz &
           + gyzu*gxzu*Axz + gyzu*gyzu*Ayz + gyzu*gzzu*Azz

      Azzu = gxzu*gxzu*Axx + gxzu*gyzu*Axy + gxzu*gzzu*Axz &
           + gyzu*gxzu*Axy + gyzu*gyzu*Ayy + gyzu*gzzu*Ayz &
           + gzzu*gxzu*Axz + gzzu*gyzu*Ayz + gzzu*gzzu*Azz

    end subroutine ICPertFLRW_raise_two_indices

    ! The following subroutine computes G^{ijk} = g^{il} g^{jm} g^{kn} G_{lmn}
    ! With G having the Christoffel symbols symmetries
    subroutine ICPertFLRW_raise_three_indices ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                                Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz, &
                                                Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz, &
                                                Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz, &
                                                Gxxxu, Gxxyu, Gxxzu, Gxyyu, Gxyzu, Gxzzu, &
                                                Gyxxu, Gyxyu, Gyxzu, Gyyyu, Gyyzu, Gyzzu, &
                                                Gzxxu, Gzxyu, Gzxzu, Gzyyu, Gzyzu, Gzzzu )
      ! Input
      CCTK_REAL :: gxxu, gxyu, gxzu, gyyu, gyzu, gzzu
      CCTK_REAL :: Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz
      CCTK_REAL :: Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz
      CCTK_REAL :: Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz

      ! Output
      CCTK_REAL :: Gxxxu, Gxxyu, Gxxzu, Gxyyu, Gxyzu, Gxzzu
      CCTK_REAL :: Gyxxu, Gyxyu, Gyxzu, Gyyyu, Gyyzu, Gyzzu
      CCTK_REAL :: Gzxxu, Gzxyu, Gzxzu, Gzyyu, Gzyzu, Gzzzu

      ! Gxnnu
      Gxxxu = gxxu*gxxu*gxxu*Gxxx + gxxu*gxxu*gxyu*Gxxy + gxxu*gxxu*gxzu*Gxxz + gxxu*gxyu*gxxu*Gxxy &
            + gxxu*gxyu*gxyu*Gxyy + gxxu*gxyu*gxzu*Gxyz + gxxu*gxzu*gxxu*Gxxz + gxxu*gxzu*gxyu*Gxyz &
            + gxxu*gxzu*gxzu*Gxzz + gxyu*gxxu*gxxu*Gyxx + gxyu*gxxu*gxyu*Gyxy + gxyu*gxxu*gxzu*Gyxz &
            + gxyu*gxyu*gxxu*Gyxy + gxyu*gxyu*gxyu*Gyyy + gxyu*gxyu*gxzu*Gyyz + gxyu*gxzu*gxxu*Gyxz &
            + gxyu*gxzu*gxyu*Gyyz + gxyu*gxzu*gxzu*Gyzz + gxzu*gxxu*gxxu*Gzxx + gxzu*gxxu*gxyu*Gzxy &
            + gxzu*gxxu*gxzu*Gzxz + gxzu*gxyu*gxxu*Gzxy + gxzu*gxyu*gxyu*Gzyy + gxzu*gxyu*gxzu*Gzyz &
            + gxzu*gxzu*gxxu*Gzxz + gxzu*gxzu*gxyu*Gzyz + gxzu*gxzu*gxzu*Gzzz 

      Gxxyu = gxxu*gxxu*gxyu*Gxxx + gxxu*gxxu*gyyu*Gxxy + gxxu*gxxu*gyzu*Gxxz + gxxu*gxyu*gxyu*Gxxy &
            + gxxu*gxyu*gyyu*Gxyy + gxxu*gxyu*gyzu*Gxyz + gxxu*gxzu*gxyu*Gxxz + gxxu*gxzu*gyyu*Gxyz &
            + gxxu*gxzu*gyzu*Gxzz + gxyu*gxxu*gxyu*Gyxx + gxyu*gxxu*gyyu*Gyxy + gxyu*gxxu*gyzu*Gyxz &
            + gxyu*gxyu*gxyu*Gyxy + gxyu*gxyu*gyyu*Gyyy + gxyu*gxyu*gyzu*Gyyz + gxyu*gxzu*gxyu*Gyxz &
            + gxyu*gxzu*gyyu*Gyyz + gxyu*gxzu*gyzu*Gyzz + gxzu*gxxu*gxyu*Gzxx + gxzu*gxxu*gyyu*Gzxy &
            + gxzu*gxxu*gyzu*Gzxz + gxzu*gxyu*gxyu*Gzxy + gxzu*gxyu*gyyu*Gzyy + gxzu*gxyu*gyzu*Gzyz &
            + gxzu*gxzu*gxyu*Gzxz + gxzu*gxzu*gyyu*Gzyz + gxzu*gxzu*gyzu*Gzzz 

      Gxxzu = gxxu*gxxu*gxzu*Gxxx + gxxu*gxxu*gyzu*Gxxy + gxxu*gxxu*gzzu*Gxxz + gxxu*gxyu*gxzu*Gxxy &
            + gxxu*gxyu*gyzu*Gxyy + gxxu*gxyu*gzzu*Gxyz + gxxu*gxzu*gxzu*Gxxz + gxxu*gxzu*gyzu*Gxyz &
            + gxxu*gxzu*gzzu*Gxzz + gxyu*gxxu*gxzu*Gyxx + gxyu*gxxu*gyzu*Gyxy + gxyu*gxxu*gzzu*Gyxz &
            + gxyu*gxyu*gxzu*Gyxy + gxyu*gxyu*gyzu*Gyyy + gxyu*gxyu*gzzu*Gyyz + gxyu*gxzu*gxzu*Gyxz &
            + gxyu*gxzu*gyzu*Gyyz + gxyu*gxzu*gzzu*Gyzz + gxzu*gxxu*gxzu*Gzxx + gxzu*gxxu*gyzu*Gzxy &
            + gxzu*gxxu*gzzu*Gzxz + gxzu*gxyu*gxzu*Gzxy + gxzu*gxyu*gyzu*Gzyy + gxzu*gxyu*gzzu*Gzyz &
            + gxzu*gxzu*gxzu*Gzxz + gxzu*gxzu*gyzu*Gzyz + gxzu*gxzu*gzzu*Gzzz 

      Gxyyu = gxxu*gxyu*gxyu*Gxxx + gxxu*gxyu*gyyu*Gxxy + gxxu*gxyu*gyzu*Gxxz + gxxu*gyyu*gxyu*Gxxy &
            + gxxu*gyyu*gyyu*Gxyy + gxxu*gyyu*gyzu*Gxyz + gxxu*gyzu*gxyu*Gxxz + gxxu*gyzu*gyyu*Gxyz &
            + gxxu*gyzu*gyzu*Gxzz + gxyu*gxyu*gxyu*Gyxx + gxyu*gxyu*gyyu*Gyxy + gxyu*gxyu*gyzu*Gyxz &
            + gxyu*gyyu*gxyu*Gyxy + gxyu*gyyu*gyyu*Gyyy + gxyu*gyyu*gyzu*Gyyz + gxyu*gyzu*gxyu*Gyxz &
            + gxyu*gyzu*gyyu*Gyyz + gxyu*gyzu*gyzu*Gyzz + gxzu*gxyu*gxyu*Gzxx + gxzu*gxyu*gyyu*Gzxy &
            + gxzu*gxyu*gyzu*Gzxz + gxzu*gyyu*gxyu*Gzxy + gxzu*gyyu*gyyu*Gzyy + gxzu*gyyu*gyzu*Gzyz &
            + gxzu*gyzu*gxyu*Gzxz + gxzu*gyzu*gyyu*Gzyz + gxzu*gyzu*gyzu*Gzzz 

      Gxyzu = gxxu*gxyu*gxzu*Gxxx + gxxu*gxyu*gyzu*Gxxy + gxxu*gxyu*gzzu*Gxxz + gxxu*gyyu*gxzu*Gxxy &
            + gxxu*gyyu*gyzu*Gxyy + gxxu*gyyu*gzzu*Gxyz + gxxu*gyzu*gxzu*Gxxz + gxxu*gyzu*gyzu*Gxyz &
            + gxxu*gyzu*gzzu*Gxzz + gxyu*gxyu*gxzu*Gyxx + gxyu*gxyu*gyzu*Gyxy + gxyu*gxyu*gzzu*Gyxz &
            + gxyu*gyyu*gxzu*Gyxy + gxyu*gyyu*gyzu*Gyyy + gxyu*gyyu*gzzu*Gyyz + gxyu*gyzu*gxzu*Gyxz &
            + gxyu*gyzu*gyzu*Gyyz + gxyu*gyzu*gzzu*Gyzz + gxzu*gxyu*gxzu*Gzxx + gxzu*gxyu*gyzu*Gzxy &
            + gxzu*gxyu*gzzu*Gzxz + gxzu*gyyu*gxzu*Gzxy + gxzu*gyyu*gyzu*Gzyy + gxzu*gyyu*gzzu*Gzyz &
            + gxzu*gyzu*gxzu*Gzxz + gxzu*gyzu*gyzu*Gzyz + gxzu*gyzu*gzzu*Gzzz 

      Gxzzu = gxxu*gxzu*gxzu*Gxxx + gxxu*gxzu*gyzu*Gxxy + gxxu*gxzu*gzzu*Gxxz + gxxu*gyzu*gxzu*Gxxy &
            + gxxu*gyzu*gyzu*Gxyy + gxxu*gyzu*gzzu*Gxyz + gxxu*gzzu*gxzu*Gxxz + gxxu*gzzu*gyzu*Gxyz &
            + gxxu*gzzu*gzzu*Gxzz + gxyu*gxzu*gxzu*Gyxx + gxyu*gxzu*gyzu*Gyxy + gxyu*gxzu*gzzu*Gyxz &
            + gxyu*gyzu*gxzu*Gyxy + gxyu*gyzu*gyzu*Gyyy + gxyu*gyzu*gzzu*Gyyz + gxyu*gzzu*gxzu*Gyxz &
            + gxyu*gzzu*gyzu*Gyyz + gxyu*gzzu*gzzu*Gyzz + gxzu*gxzu*gxzu*Gzxx + gxzu*gxzu*gyzu*Gzxy &
            + gxzu*gxzu*gzzu*Gzxz + gxzu*gyzu*gxzu*Gzxy + gxzu*gyzu*gyzu*Gzyy + gxzu*gyzu*gzzu*Gzyz &
            + gxzu*gzzu*gxzu*Gzxz + gxzu*gzzu*gyzu*Gzyz + gxzu*gzzu*gzzu*Gzzz 

      ! Gynnu
      Gyxxu = gxyu*gxxu*gxxu*Gxxx + gxyu*gxxu*gxyu*Gxxy + gxyu*gxxu*gxzu*Gxxz + gxyu*gxyu*gxxu*Gxxy &
            + gxyu*gxyu*gxyu*Gxyy + gxyu*gxyu*gxzu*Gxyz + gxyu*gxzu*gxxu*Gxxz + gxyu*gxzu*gxyu*Gxyz &
            + gxyu*gxzu*gxzu*Gxzz + gyyu*gxxu*gxxu*Gyxx + gyyu*gxxu*gxyu*Gyxy + gyyu*gxxu*gxzu*Gyxz &
            + gyyu*gxyu*gxxu*Gyxy + gyyu*gxyu*gxyu*Gyyy + gyyu*gxyu*gxzu*Gyyz + gyyu*gxzu*gxxu*Gyxz &
            + gyyu*gxzu*gxyu*Gyyz + gyyu*gxzu*gxzu*Gyzz + gyzu*gxxu*gxxu*Gzxx + gyzu*gxxu*gxyu*Gzxy &
            + gyzu*gxxu*gxzu*Gzxz + gyzu*gxyu*gxxu*Gzxy + gyzu*gxyu*gxyu*Gzyy + gyzu*gxyu*gxzu*Gzyz &
            + gyzu*gxzu*gxxu*Gzxz + gyzu*gxzu*gxyu*Gzyz + gyzu*gxzu*gxzu*Gzzz 

      Gyxyu = gxyu*gxxu*gxyu*Gxxx + gxyu*gxxu*gyyu*Gxxy + gxyu*gxxu*gyzu*Gxxz + gxyu*gxyu*gxyu*Gxxy &
            + gxyu*gxyu*gyyu*Gxyy + gxyu*gxyu*gyzu*Gxyz + gxyu*gxzu*gxyu*Gxxz + gxyu*gxzu*gyyu*Gxyz &
            + gxyu*gxzu*gyzu*Gxzz + gyyu*gxxu*gxyu*Gyxx + gyyu*gxxu*gyyu*Gyxy + gyyu*gxxu*gyzu*Gyxz &
            + gyyu*gxyu*gxyu*Gyxy + gyyu*gxyu*gyyu*Gyyy + gyyu*gxyu*gyzu*Gyyz + gyyu*gxzu*gxyu*Gyxz &
            + gyyu*gxzu*gyyu*Gyyz + gyyu*gxzu*gyzu*Gyzz + gyzu*gxxu*gxyu*Gzxx + gyzu*gxxu*gyyu*Gzxy &
            + gyzu*gxxu*gyzu*Gzxz + gyzu*gxyu*gxyu*Gzxy + gyzu*gxyu*gyyu*Gzyy + gyzu*gxyu*gyzu*Gzyz &
            + gyzu*gxzu*gxyu*Gzxz + gyzu*gxzu*gyyu*Gzyz + gyzu*gxzu*gyzu*Gzzz 

      Gyxzu = gxyu*gxxu*gxzu*Gxxx + gxyu*gxxu*gyzu*Gxxy + gxyu*gxxu*gzzu*Gxxz + gxyu*gxyu*gxzu*Gxxy &
            + gxyu*gxyu*gyzu*Gxyy + gxyu*gxyu*gzzu*Gxyz + gxyu*gxzu*gxzu*Gxxz + gxyu*gxzu*gyzu*Gxyz &
            + gxyu*gxzu*gzzu*Gxzz + gyyu*gxxu*gxzu*Gyxx + gyyu*gxxu*gyzu*Gyxy + gyyu*gxxu*gzzu*Gyxz &
            + gyyu*gxyu*gxzu*Gyxy + gyyu*gxyu*gyzu*Gyyy + gyyu*gxyu*gzzu*Gyyz + gyyu*gxzu*gxzu*Gyxz &
            + gyyu*gxzu*gyzu*Gyyz + gyyu*gxzu*gzzu*Gyzz + gyzu*gxxu*gxzu*Gzxx + gyzu*gxxu*gyzu*Gzxy &
            + gyzu*gxxu*gzzu*Gzxz + gyzu*gxyu*gxzu*Gzxy + gyzu*gxyu*gyzu*Gzyy + gyzu*gxyu*gzzu*Gzyz &
            + gyzu*gxzu*gxzu*Gzxz + gyzu*gxzu*gyzu*Gzyz + gyzu*gxzu*gzzu*Gzzz 

      Gyyyu = gxyu*gxyu*gxyu*Gxxx + gxyu*gxyu*gyyu*Gxxy + gxyu*gxyu*gyzu*Gxxz + gxyu*gyyu*gxyu*Gxxy &
            + gxyu*gyyu*gyyu*Gxyy + gxyu*gyyu*gyzu*Gxyz + gxyu*gyzu*gxyu*Gxxz + gxyu*gyzu*gyyu*Gxyz &
            + gxyu*gyzu*gyzu*Gxzz + gyyu*gxyu*gxyu*Gyxx + gyyu*gxyu*gyyu*Gyxy + gyyu*gxyu*gyzu*Gyxz &
            + gyyu*gyyu*gxyu*Gyxy + gyyu*gyyu*gyyu*Gyyy + gyyu*gyyu*gyzu*Gyyz + gyyu*gyzu*gxyu*Gyxz &
            + gyyu*gyzu*gyyu*Gyyz + gyyu*gyzu*gyzu*Gyzz + gyzu*gxyu*gxyu*Gzxx + gyzu*gxyu*gyyu*Gzxy &
            + gyzu*gxyu*gyzu*Gzxz + gyzu*gyyu*gxyu*Gzxy + gyzu*gyyu*gyyu*Gzyy + gyzu*gyyu*gyzu*Gzyz &
            + gyzu*gyzu*gxyu*Gzxz + gyzu*gyzu*gyyu*Gzyz + gyzu*gyzu*gyzu*Gzzz 

      Gyyzu = gxyu*gxyu*gxzu*Gxxx + gxyu*gxyu*gyzu*Gxxy + gxyu*gxyu*gzzu*Gxxz + gxyu*gyyu*gxzu*Gxxy &
            + gxyu*gyyu*gyzu*Gxyy + gxyu*gyyu*gzzu*Gxyz + gxyu*gyzu*gxzu*Gxxz + gxyu*gyzu*gyzu*Gxyz &
            + gxyu*gyzu*gzzu*Gxzz + gyyu*gxyu*gxzu*Gyxx + gyyu*gxyu*gyzu*Gyxy + gyyu*gxyu*gzzu*Gyxz &
            + gyyu*gyyu*gxzu*Gyxy + gyyu*gyyu*gyzu*Gyyy + gyyu*gyyu*gzzu*Gyyz + gyyu*gyzu*gxzu*Gyxz &
            + gyyu*gyzu*gyzu*Gyyz + gyyu*gyzu*gzzu*Gyzz + gyzu*gxyu*gxzu*Gzxx + gyzu*gxyu*gyzu*Gzxy &
            + gyzu*gxyu*gzzu*Gzxz + gyzu*gyyu*gxzu*Gzxy + gyzu*gyyu*gyzu*Gzyy + gyzu*gyyu*gzzu*Gzyz &
            + gyzu*gyzu*gxzu*Gzxz + gyzu*gyzu*gyzu*Gzyz + gyzu*gyzu*gzzu*Gzzz 

      Gyzzu = gxyu*gxzu*gxzu*Gxxx + gxyu*gxzu*gyzu*Gxxy + gxyu*gxzu*gzzu*Gxxz + gxyu*gyzu*gxzu*Gxxy &
            + gxyu*gyzu*gyzu*Gxyy + gxyu*gyzu*gzzu*Gxyz + gxyu*gzzu*gxzu*Gxxz + gxyu*gzzu*gyzu*Gxyz &
            + gxyu*gzzu*gzzu*Gxzz + gyyu*gxzu*gxzu*Gyxx + gyyu*gxzu*gyzu*Gyxy + gyyu*gxzu*gzzu*Gyxz &
            + gyyu*gyzu*gxzu*Gyxy + gyyu*gyzu*gyzu*Gyyy + gyyu*gyzu*gzzu*Gyyz + gyyu*gzzu*gxzu*Gyxz &
            + gyyu*gzzu*gyzu*Gyyz + gyyu*gzzu*gzzu*Gyzz + gyzu*gxzu*gxzu*Gzxx + gyzu*gxzu*gyzu*Gzxy &
            + gyzu*gxzu*gzzu*Gzxz + gyzu*gyzu*gxzu*Gzxy + gyzu*gyzu*gyzu*Gzyy + gyzu*gyzu*gzzu*Gzyz &
            + gyzu*gzzu*gxzu*Gzxz + gyzu*gzzu*gyzu*Gzyz + gyzu*gzzu*gzzu*Gzzz

      ! Gznnu
      Gzxxu = gxzu*gxxu*gxxu*Gxxx + gxzu*gxxu*gxyu*Gxxy + gxzu*gxxu*gxzu*Gxxz + gxzu*gxyu*gxxu*Gxxy &
            + gxzu*gxyu*gxyu*Gxyy + gxzu*gxyu*gxzu*Gxyz + gxzu*gxzu*gxxu*Gxxz + gxzu*gxzu*gxyu*Gxyz &
            + gxzu*gxzu*gxzu*Gxzz + gyzu*gxxu*gxxu*Gyxx + gyzu*gxxu*gxyu*Gyxy + gyzu*gxxu*gxzu*Gyxz &
            + gyzu*gxyu*gxxu*Gyxy + gyzu*gxyu*gxyu*Gyyy + gyzu*gxyu*gxzu*Gyyz + gyzu*gxzu*gxxu*Gyxz &
            + gyzu*gxzu*gxyu*Gyyz + gyzu*gxzu*gxzu*Gyzz + gzzu*gxxu*gxxu*Gzxx + gzzu*gxxu*gxyu*Gzxy &
            + gzzu*gxxu*gxzu*Gzxz + gzzu*gxyu*gxxu*Gzxy + gzzu*gxyu*gxyu*Gzyy + gzzu*gxyu*gxzu*Gzyz &
            + gzzu*gxzu*gxxu*Gzxz + gzzu*gxzu*gxyu*Gzyz + gzzu*gxzu*gxzu*Gzzz 

      Gzxyu = gxzu*gxxu*gxyu*Gxxx + gxzu*gxxu*gyyu*Gxxy + gxzu*gxxu*gyzu*Gxxz + gxzu*gxyu*gxyu*Gxxy &
            + gxzu*gxyu*gyyu*Gxyy + gxzu*gxyu*gyzu*Gxyz + gxzu*gxzu*gxyu*Gxxz + gxzu*gxzu*gyyu*Gxyz &
            + gxzu*gxzu*gyzu*Gxzz + gyzu*gxxu*gxyu*Gyxx + gyzu*gxxu*gyyu*Gyxy + gyzu*gxxu*gyzu*Gyxz &
            + gyzu*gxyu*gxyu*Gyxy + gyzu*gxyu*gyyu*Gyyy + gyzu*gxyu*gyzu*Gyyz + gyzu*gxzu*gxyu*Gyxz &
            + gyzu*gxzu*gyyu*Gyyz + gyzu*gxzu*gyzu*Gyzz + gzzu*gxxu*gxyu*Gzxx + gzzu*gxxu*gyyu*Gzxy &
            + gzzu*gxxu*gyzu*Gzxz + gzzu*gxyu*gxyu*Gzxy + gzzu*gxyu*gyyu*Gzyy + gzzu*gxyu*gyzu*Gzyz &
            + gzzu*gxzu*gxyu*Gzxz + gzzu*gxzu*gyyu*Gzyz + gzzu*gxzu*gyzu*Gzzz 

      Gzxzu = gxzu*gxxu*gxzu*Gxxx + gxzu*gxxu*gyzu*Gxxy + gxzu*gxxu*gzzu*Gxxz + gxzu*gxyu*gxzu*Gxxy &
            + gxzu*gxyu*gyzu*Gxyy + gxzu*gxyu*gzzu*Gxyz + gxzu*gxzu*gxzu*Gxxz + gxzu*gxzu*gyzu*Gxyz &
            + gxzu*gxzu*gzzu*Gxzz + gyzu*gxxu*gxzu*Gyxx + gyzu*gxxu*gyzu*Gyxy + gyzu*gxxu*gzzu*Gyxz &
            + gyzu*gxyu*gxzu*Gyxy + gyzu*gxyu*gyzu*Gyyy + gyzu*gxyu*gzzu*Gyyz + gyzu*gxzu*gxzu*Gyxz &
            + gyzu*gxzu*gyzu*Gyyz + gyzu*gxzu*gzzu*Gyzz + gzzu*gxxu*gxzu*Gzxx + gzzu*gxxu*gyzu*Gzxy &
            + gzzu*gxxu*gzzu*Gzxz + gzzu*gxyu*gxzu*Gzxy + gzzu*gxyu*gyzu*Gzyy + gzzu*gxyu*gzzu*Gzyz &
            + gzzu*gxzu*gxzu*Gzxz + gzzu*gxzu*gyzu*Gzyz + gzzu*gxzu*gzzu*Gzzz 

      Gzyyu = gxzu*gxyu*gxyu*Gxxx + gxzu*gxyu*gyyu*Gxxy + gxzu*gxyu*gyzu*Gxxz + gxzu*gyyu*gxyu*Gxxy &
            + gxzu*gyyu*gyyu*Gxyy + gxzu*gyyu*gyzu*Gxyz + gxzu*gyzu*gxyu*Gxxz + gxzu*gyzu*gyyu*Gxyz &
            + gxzu*gyzu*gyzu*Gxzz + gyzu*gxyu*gxyu*Gyxx + gyzu*gxyu*gyyu*Gyxy + gyzu*gxyu*gyzu*Gyxz &
            + gyzu*gyyu*gxyu*Gyxy + gyzu*gyyu*gyyu*Gyyy + gyzu*gyyu*gyzu*Gyyz + gyzu*gyzu*gxyu*Gyxz &
            + gyzu*gyzu*gyyu*Gyyz + gyzu*gyzu*gyzu*Gyzz + gzzu*gxyu*gxyu*Gzxx + gzzu*gxyu*gyyu*Gzxy &
            + gzzu*gxyu*gyzu*Gzxz + gzzu*gyyu*gxyu*Gzxy + gzzu*gyyu*gyyu*Gzyy + gzzu*gyyu*gyzu*Gzyz &
            + gzzu*gyzu*gxyu*Gzxz + gzzu*gyzu*gyyu*Gzyz + gzzu*gyzu*gyzu*Gzzz 

      Gzyzu = gxzu*gxyu*gxzu*Gxxx + gxzu*gxyu*gyzu*Gxxy + gxzu*gxyu*gzzu*Gxxz + gxzu*gyyu*gxzu*Gxxy &
            + gxzu*gyyu*gyzu*Gxyy + gxzu*gyyu*gzzu*Gxyz + gxzu*gyzu*gxzu*Gxxz + gxzu*gyzu*gyzu*Gxyz &
            + gxzu*gyzu*gzzu*Gxzz + gyzu*gxyu*gxzu*Gyxx + gyzu*gxyu*gyzu*Gyxy + gyzu*gxyu*gzzu*Gyxz &
            + gyzu*gyyu*gxzu*Gyxy + gyzu*gyyu*gyzu*Gyyy + gyzu*gyyu*gzzu*Gyyz + gyzu*gyzu*gxzu*Gyxz &
            + gyzu*gyzu*gyzu*Gyyz + gyzu*gyzu*gzzu*Gyzz + gzzu*gxyu*gxzu*Gzxx + gzzu*gxyu*gyzu*Gzxy &
            + gzzu*gxyu*gzzu*Gzxz + gzzu*gyyu*gxzu*Gzxy + gzzu*gyyu*gyzu*Gzyy + gzzu*gyyu*gzzu*Gzyz &
            + gzzu*gyzu*gxzu*Gzxz + gzzu*gyzu*gyzu*Gzyz + gzzu*gyzu*gzzu*Gzzz 

      Gzzzu = gxzu*gxzu*gxzu*Gxxx + gxzu*gxzu*gyzu*Gxxy + gxzu*gxzu*gzzu*Gxxz + gxzu*gyzu*gxzu*Gxxy &
            + gxzu*gyzu*gyzu*Gxyy + gxzu*gyzu*gzzu*Gxyz + gxzu*gzzu*gxzu*Gxxz + gxzu*gzzu*gyzu*Gxyz &
            + gxzu*gzzu*gzzu*Gxzz + gyzu*gxzu*gxzu*Gyxx + gyzu*gxzu*gyzu*Gyxy + gyzu*gxzu*gzzu*Gyxz &
            + gyzu*gyzu*gxzu*Gyxy + gyzu*gyzu*gyzu*Gyyy + gyzu*gyzu*gzzu*Gyyz + gyzu*gzzu*gxzu*Gyxz &
            + gyzu*gzzu*gyzu*Gyyz + gyzu*gzzu*gzzu*Gyzz + gzzu*gxzu*gxzu*Gzxx + gzzu*gxzu*gyzu*Gzxy &
            + gzzu*gxzu*gzzu*Gzxz + gzzu*gyzu*gxzu*Gzxy + gzzu*gyzu*gyzu*Gzyy + gzzu*gyzu*gzzu*Gzyz &
            + gzzu*gzzu*gxzu*Gzxz + gzzu*gzzu*gyzu*Gzyz + gzzu*gzzu*gzzu*Gzzz 

    end subroutine ICPertFLRW_raise_three_indices

    ! The following subroutine computes A = g^{ij} A_{ij}
    subroutine ICPertFLRW_take_trace ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                       Axx, Axy, Axz, Ayy, Ayz, Azz, A)
      ! Input
      CCTK_REAL :: gxxu, gxyu, gxzu, gyyu, gyzu, gzzu
      CCTK_REAL :: Axx, Axy, Axz, Ayy, Ayz, Azz

      ! Output
      CCTK_REAL :: A

      A = gxxu*Axx + gxyu*Axy + gxzu*Axz &
        + gxyu*Axy + gyyu*Ayy + gyzu*Ayz &
        + gxzu*Axz + gyzu*Ayz + gzzu*Azz

    end subroutine ICPertFLRW_take_trace

    ! The following subroutine computes A = g^{ij} A_{ij}, and A_{ij} isn't symmetric
    subroutine ICPertFLRW_take_trace_notsym ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                              Axx, Axy, Axz, Ayx, Ayy, Ayz, Azx, Azy, Azz, A)
      ! Input
      CCTK_REAL :: gxxu, gxyu, gxzu, gyyu, gyzu, gzzu
      CCTK_REAL :: Axx, Axy, Axz
      CCTK_REAL :: Ayx, Ayy, Ayz
      CCTK_REAL :: Azx, Azy, Azz

      ! Output
      CCTK_REAL :: A

      A = gxxu*Axx + gxyu*Axy + gxzu*Axz &
        + gxyu*Ayx + gyyu*Ayy + gyzu*Ayz &
        + gxzu*Azx + gyzu*Azy + gzzu*Azz

    end subroutine ICPertFLRW_take_trace_notsym

    ! The following subroutine computes A^{ij}A_{ij}
    subroutine ICPertFLRW_contract_over_two_indices ( Axx, Axy, Axz, Ayy, Ayz, Azz, &
                                                      Axxu, Axyu, Axzu, Ayyu, Ayzu, Azzu, &
                                                      AijuAijd )
      ! Input
      CCTK_REAL :: Axx, Axy, Axz, Ayy, Ayz, Azz
      CCTK_REAL :: Axxu, Axyu, Axzu, Ayyu, Ayzu, Azzu

      ! Output
      CCTK_REAL :: AijuAijd

      AijuAijd = Axxu*Axx + Axyu*Axy + Axzu*Axz &
               + Axyu*Axy + Ayyu*Ayy + Ayzu*Ayz &
               + Axzu*Axz + Ayzu*Ayz + Azzu*Azz

    end subroutine ICPertFLRW_contract_over_two_indices

    ! The following subroutine computes g^{ij}G_{imj}
    subroutine ICPertFLRW_Christoffel_contraction_one ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                                        Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz, &
                                                        Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz, &
                                                        Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz, &
                                                        Gnxn, Gnyn, Gnzn )
      ! Input
      CCTK_REAL :: gxxu, gxyu, gxzu, gyyu, gyzu, gzzu
      CCTK_REAL :: Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz
      CCTK_REAL :: Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz
      CCTK_REAL :: Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz 
    
      ! Output
      CCTK_REAL :: Gnxn, Gnyn, Gnzn
    
      Gnxn = gxxu*Gxxx + gxyu*Gxxy + gxzu*Gxxz + &
             gxyu*Gyxx + gyyu*Gyxy + gyzu*Gyxz + &
             gxzu*Gzxx + gyzu*Gzxy + gzzu*Gzxz
    
      Gnyn = gxxu*Gxxy + gxyu*Gxyy + gxzu*Gxyz + &
             gxyu*Gyxy + gyyu*Gyyy + gyzu*Gyyz + &
             gxzu*Gzxy + gyzu*Gzyy + gzzu*Gzyz
    
      Gnzn = gxxu*Gxxz + gxyu*Gxyz + gxzu*Gxzz + &
             gxyu*Gyxz + gyyu*Gyyz + gyzu*Gyzz + &
             gxzu*Gzxz + gyzu*Gzyz + gzzu*Gzzz
    
    end subroutine ICPertFLRW_Christoffel_contraction_one

    ! The following subroutine computes g^{ij}G_{mij}
    subroutine ICPertFLRW_Christoffel_contraction_two ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                                        Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz, &
                                                        Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz, &
                                                        Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz, &
                                                        Gxnn, Gynn, Gznn )
      ! Input
      CCTK_REAL :: gxxu, gxyu, gxzu, gyyu, gyzu, gzzu
      CCTK_REAL :: Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz
      CCTK_REAL :: Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz
      CCTK_REAL :: Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz 
    
      ! Output
      CCTK_REAL :: Gxnn, Gynn, Gznn
    
      call ICPertFLRW_take_trace ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                   Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz, Gxnn)
      call ICPertFLRW_take_trace ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                   Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz, Gynn)
      call ICPertFLRW_take_trace ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                   Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz, Gznn)
    
    end subroutine ICPertFLRW_Christoffel_contraction_two

    ! The following subroutine computes g^{ij}G_{mij}
    subroutine ICPertFLRW_Christoffel_contraction_three ( Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz, &
                                                          Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz, &
                                                          Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz, &
                                                          Gxxxu, Gxxyu, Gxxzu, Gxyyu, Gxyzu, Gxzzu, &
                                                          Gyxxu, Gyxyu, Gyxzu, Gyyyu, Gyyzu, Gyzzu, &
                                                          Gzxxu, Gzxyu, Gzxzu, Gzyyu, Gzyzu, Gzzzu, &
                                                          G )
      ! Input
      CCTK_REAL :: Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz
      CCTK_REAL :: Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz
      CCTK_REAL :: Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz
      CCTK_REAL :: Gxxxu, Gxxyu, Gxxzu, Gxyyu, Gxyzu, Gxzzu
      CCTK_REAL :: Gyxxu, Gyxyu, Gyxzu, Gyyyu, Gyyzu, Gyzzu
      CCTK_REAL :: Gzxxu, Gzxyu, Gzxzu, Gzyyu, Gzyzu, Gzzzu
    
      ! Output
      CCTK_REAL :: G
    
      G = Gxxxu*Gxxx + Gxxyu*Gxxy + Gxxzu*Gxxz &
        + Gxxyu*Gyxx + Gxyyu*Gyxy + Gxyzu*Gyxz &
        + Gxxzu*Gzxx + Gxyzu*Gzxy + Gxzzu*Gzxz &
        + Gyxxu*Gxxy + Gyxyu*Gxyy + Gyxzu*Gxyz &
        + Gyxyu*Gyxy + Gyyyu*Gyyy + Gyyzu*Gyyz &
        + Gyxzu*Gzxy + Gyyzu*Gzyy + Gyzzu*Gzyz &
        + Gzxxu*Gxxz + Gzxyu*Gxyz + Gzxzu*Gxzz &
        + Gzxyu*Gyxz + Gzyyu*Gyyz + Gzyzu*Gyzz &
        + Gzxzu*Gzxz + Gzyzu*Gzyz + Gzzzu*Gzzz
    
    end subroutine ICPertFLRW_Christoffel_contraction_three

    ! The following subroutine computes g^{ij}A_{i}B_{j}
    subroutine ICPertFLRW_product_of_two_vectors ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                                   Ax, Ay, Az, &
                                                   Bx, By, Bz, &
                                                   AB )
      ! Input
      CCTK_REAL :: gxxu, gxyu, gxzu, gyyu, gyzu, gzzu
      CCTK_REAL :: Ax, Ay, Az
      CCTK_REAL :: Bx, By, Bz

      ! Output
      CCTK_REAL :: AB

      AB = gxxu*Ax*Bx + gxyu*Ax*By + gxzu*Ax*Bz &
         + gxyu*Ay*Bx + gyyu*Ay*By + gyzu*Ay*Bz &
         + gxzu*Az*Bx + gyzu*Az*By + gzzu*Az*Bz

    end subroutine ICPertFLRW_product_of_two_vectors

    ! The following subroutine computes derivative of Gudd: \partial_{i}(g^{im} G_{mkl})
    subroutine ICPertFLRW_derivative_of_Gudd ( gxxu, gxyu, gxzu, gyyu, gyzu, gzzu, &
                                               dxgxxu, dxgxyu, dxgxzu, &
                                               dygxyu, dygyyu, dygyzu, dzgxzu, dzgyzu, dzgzzu, &
                                               Gxkl, Gykl, Gzkl, dxGxkl, dxGykl, dxGzkl, &
                                               dyGxkl, dyGykl, dyGzkl, dzGxkl, dzGykl, dzGzkl, dGkl )

      ! Input
      CCTK_REAL :: gxxu, gxyu, gxzu, gyyu, gyzu, gzzu
      CCTK_REAL :: dxgxxu, dxgxyu, dxgxzu, dygxyu, dygyyu, dygyzu, dzgxzu, dzgyzu, dzgzzu
      CCTK_REAL :: Gxkl, Gykl, Gzkl
      CCTK_REAL :: dxGxkl, dxGykl, dxGzkl, dyGxkl, dyGykl, dyGzkl, dzGxkl, dzGykl, dzGzkl

      ! Output
      CCTK_REAL :: dGkl

      dGkl = dxgxxu*Gxkl + gxxu*dxGxkl + dxgxyu*Gykl + gxyu*dxGykl + dxgxzu*Gzkl + gxzu*dxGzkl &
           + dygxyu*Gxkl + gxyu*dyGxkl + dygyyu*Gykl + gyyu*dyGykl + dygyzu*Gzkl + gyzu*dyGzkl &
           + dzgxzu*Gxkl + gxzu*dzGxkl + dzgyzu*Gykl + gyzu*dzGykl + dzgzzu*Gzkl + gzzu*dzGzkl

    end subroutine ICPertFLRW_derivative_of_Gudd

end module TensorCalc



















