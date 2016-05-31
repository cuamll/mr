module common
  implicit none
  real, public :: q, lambda, volume, ebar_x, ebar_y, ebar_z, temp, beta
  real, public :: rot_delt
  integer, public :: L,seed,accepth,acceptr,acceptg,iterations
  integer, dimension(:), allocatable, public :: pos,neg
  integer, dimension(:,:,:), allocatable, public :: v
  real*8, dimension(:,:,:), allocatable, public :: mnphi_x, mnphi_y, mnphi_z
  real*8, dimension(:,:,:), allocatable, public :: e_x, e_y, e_z
  real*8, dimension(:,:,:), allocatable, public :: phi_lapack, e_x_lapack
  real*8, dimension(:,:,:), allocatable, public :: e_y_lapack, e_z_lapack
  real*8, dimension(:,:,:), allocatable, public :: e_rot_x, e_rot_y, e_rot_z
  real*8, dimension(:,:,:,:,:,:), allocatable, public :: lgf
  real*8, dimension(:), allocatable, public :: energy, sq_energy, energy_run
  ! probably more things need to go here
  character(len=11), public :: lattfile

  integer, public :: have_lgf = 0
  real*8, public :: rot_ratio, g_ratio, hop_ratio
  real, parameter, public :: eps_0=1.0
  real, parameter, public :: pi=3.141592653589793
  real, parameter, public :: twopi=6.283185307179586
  real, parameter, public :: e=2.718281828459045
  save

  contains

    subroutine PBCs
      integer :: i
      do i=1,L
        pos(i)=mod(i,L)+1
        neg(i)=mod(i+L-2,L)+1
      end do
      return
    end subroutine PBCs

    function one_to_three(x) result(coord)
      integer*8, intent(in) :: x
      integer*8 :: coord(3)
      coord(1) = (x - 1)/L**2 + 1
      coord(2) = modulo((x - 1)/L,L) + 1
      coord(3) = modulo(x - 1,L) + 1
    end function one_to_three

    function three_to_one(coord) result(x)
      integer*8, intent(in) :: coord(3)
      integer*8 :: x
      x = (coord(1) - 1) * L**2 + (coord(2) - 1) * L + coord(3)
    end function three_to_one



!----------------------------------------------------------------------C
!                                                                      C
!  Lagged Fibonacci random number generator RANMAR.                    C
!  Must be initialized with randinit() before use.                     C
!                                                                      C
!  See F. James, Comp. Phys. Comm. 60, 329 (1990), or                  C
!  G. Marsaglia et al., Stat. Prob. Lett. 9, 35 (1990).                C
!                                                                      C
!----------------------------------------------------------------------C


!----------------------------------------------------------------------C
!                                                                      C
! This is the initialization routine RMARIN for the random number      C
!     generator RANMAR                                                 C
!                                                                      C
! NOTE: The seed variables can have values between:  0 <= IJ <= 31328  C
!                                                    0 <= KL <= 30081  C
!----------------------------------------------------------------------C

      SUBROUTINE randinit(seed)
      IMPLICIT NONE
      INTEGER seed
      INTEGER ij,kl, i,j,k,l, ii,jj, m
      REAL*8 s,t
      INTEGER Maxseed
      PARAMETER (Maxseed = 900000000)
      REAL u(97), c, cd, cm
      INTEGER i97, j97, ivec
      COMMON /raset1/ u, c, cd, cm, i97, j97, ivec

      seed = mod(seed,Maxseed)
      ij = seed / 30082
      kl = seed - (30082 * ij)
      i = mod(ij/177, 177) + 2
      j = mod(ij    , 177) + 2
      k = mod(kl/169, 178) + 1
      l = mod(kl,     169)
      DO 2 ii = 1, 97
        s = 0.0
        t = 0.5
        DO 3 jj = 1, 24
          m = mod(mod(i*j, 179)*k, 179)
          i = j
          j = k
          k = m
          l = mod(53*l+1, 169)
          IF (mod(l*m, 64) .ge. 32) then
            s = s + t
          ENDIF
          t = 0.5 * t
    3   CONTINUE
        u(ii) = s
    2 CONTINUE
      c = 362436.0 / 16777216.0
      cd = 7654321.0 / 16777216.0
      cm = 16777213.0 /16777216.0
      i97 = 97
      j97 = 33
      RETURN
      END



!----------------------------------------------------------------------C
!                                                                      C
!  Lagged Fibonacci random number generator RANMAR().                  C
!                                                                      C
!----------------------------------------------------------------------C

      FUNCTION rand()
      IMPLICIT NONE
      REAL u(97), c, cd, cm, uni, rand
      INTEGER i97, j97, ivec
      COMMON /raset1/ u, c, cd, cm, i97, j97, ivec

      uni = u(i97) - u(j97)
      IF( uni .LT. 0.0 ) uni = uni + 1.0
      u(i97) = uni
      i97 = i97 - 1
      IF(i97 .EQ. 0) i97 = 97
      j97 = j97 - 1
      IF(j97 .EQ. 0) j97 = 97
      c = c - cd
      IF( c .LT. 0.0 ) c = c + cm
      uni = uni - c
      IF( uni .LT. 0.0 ) uni = uni + 1.0
      rand = uni
      RETURN
      END

end module common
