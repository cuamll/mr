module common
  use fftw
  implicit none

  integer, public, parameter :: iprec = 15
  integer, public, parameter :: prec = 15
  integer, public, parameter :: expo = 300
  integer, public, parameter :: rk = selected_real_kind(prec, expo)
  integer, public, parameter :: ik = selected_int_kind(iprec)
  integer, public, parameter :: bz=2

  real(kind=8), parameter, public :: pi=3.141592653589793
  real(kind=8), parameter, public :: twopi=6.283185307179586
  real(kind=8), parameter, public :: e=2.718281828459045

  logical, public :: do_corr = .false., verbose = .false., canon = .false.

  integer, public :: have_lgf = 0
  integer(kind=ik), public :: mu_tot = 0
  integer, public :: glob = 0
  integer, public :: MPI_NEW_INT, MPI_NEW_REAL, MPI_NEW_COMPLEX
  integer(kind=4), public :: seed, no_samples, no_threads, L, add_charges,&
  no_measurements, therm_sweeps, measurement_sweeps, sample_interval
  integer(kind=ik), dimension(6), public :: attempts, accepts
  integer(kind=ik), dimension(:), allocatable, public :: bin_count
  integer(kind=4), dimension(:), allocatable, public :: pos,neg
  integer(kind=4), dimension(:,:), allocatable, public :: v

  real(kind=8), public :: q, lambda, volume, temp, beta, u_tot, eps_0,&
  bin_size, rot_delt, g0, rot_ratio, g_ratio, hop_ratio, g_thr, e_c
  real(kind=8), dimension(:,:,:), allocatable, public :: e_field, mnphi
  real(kind=8), dimension(:,:), allocatable, public :: lgf
  real(kind=rk), public :: ener_tot_sum, ener_tot_sq_sum, rho_avg, div, divsq
  real(kind=rk), dimension(2), public :: ebar, ebar_dip, ebar_wind, ebar_sum,&
  ebar_sq_sum, ebar_dip_sum, ebar_dip_sq_sum, ebar_wind_sum, ebar_wind_sq_sum,&
  avg_field_total, avg_field_sq_total
  real(kind=rk), dimension(:), allocatable, public :: dist_r
  real(kind=rk), dimension(:,:), allocatable, public :: v_avg, dir_struc,&
  windings, windings_sq
  real(kind=rk), dimension(:,:,:), allocatable, public :: e_tot_avg

  complex(kind=rk), dimension(:,:,:), allocatable, public :: s_ab,&
  s_ab_rot, s_ab_irrot

  complex(kind=rk), dimension(:,:,:,:), allocatable, public :: s_ab_large,&
  s_ab_rot_large, s_ab_irrot_large

  complex(kind=rk), dimension(:,:), allocatable, public :: ch_ch,&
  rho_k_m,rho_k_p
  complex(kind=rk), public :: runtot

  character(:), allocatable :: lattfile, arg, charge_st_file,&
  dir_st_file, dir_dist_file, sphe_sus_file, s_ab_file,&
  e_field_file, avg_field_file, chi_ab_file,&
  windings_file, windings_sq_file, lgf_path

  character(6) :: charge_gen

  type(C_PTR) :: plan_x, plan_y, plan_ch
  real(C_DOUBLE), dimension(:,:), allocatable :: ch_in, e_in
  complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: chk
  complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: exk, eyk

  save

  contains

    subroutine PBCs
      integer(kind=4) :: i
      do i=1,L
        pos(i)=mod(i,L)+1
        neg(i)=mod(i+L-2,L)+1
      end do
      return
    end subroutine PBCs


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
      REAL(KIND=8) s,t
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
      END SUBROUTINE randinit



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
      END FUNCTION rand

end module common
