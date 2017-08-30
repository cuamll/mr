module common
  implicit none
  logical, public :: do_corr
  integer, public :: seed,no_samples,no_threads
  integer(kind=8), public :: L,add_charges,no_measurements
  integer(kind=8), public :: accepth,acceptr,acceptg
  integer(kind=8), public :: therm_sweeps,measurement_sweeps,sample_interval
  real(kind=8), public :: q, lambda, volume, temp, beta
  real(kind=8), public :: eps_0, bin_size, rot_delt, g0
  real(kind=16), public :: ener_tot_sum, ener_rot_sum, ener_irrot_sum
  real(kind=16), public :: ener_tot_sq_sum, ener_rot_sq_sum, ener_irrot_sq_sum
  real(kind=16), dimension(2), public :: ebar_sum, ebar_sq_sum
  integer(kind=8), dimension(:), allocatable, public :: bin_count
  integer, dimension(:), allocatable, public :: pos,neg
  integer, dimension(:,:), allocatable, public :: v
  real(kind=8), dimension(:), allocatable, public :: ebar, energy, sq_energy
  !real(kind=16), dimension(:), allocatable, public :: ebar_sum, ebar_sq_sum
  real(kind=16), dimension(:), allocatable, public :: dist_r
  real(kind=16), dimension(:,:), allocatable, public :: v_avg
  real(kind=8), dimension(:,:,:), allocatable, public :: e_field, mnphi
  real(kind=16), dimension(:,:,:), allocatable, public :: e_tot_avg, e_rot_avg, e_irrot_avg
  real(kind=16), dimension(:,:), allocatable, public :: dir_struc
  real(kind=8), dimension(:,:,:,:), allocatable, public :: lgf
  real(kind=16), dimension(:,:,:,:), allocatable, public :: s_ab, s_ab_rot, s_ab_irrot
  complex(kind=16), dimension(:,:), allocatable, public :: ch_ch
  complex(kind=16), dimension(:,:), allocatable, public :: rho_k_m,rho_k_p

    ! probably more things need to go here
  character(len=200), public :: lattfile_long, en_long, sq_en_long
  character(len=200), public :: e_field_long, arg_long, ch_st_l, fi_st_l
  character(len=200), public :: s_ab_l, s_p_l, dir_st_l, dir_d_s_l, fe_ch_l
  character(len=200), public :: ir_fe_l,ir_sab_l,ir_sp_l,r_fe_l,r_sab_l,r_sp_l
  character(len=200), public :: spa_l, r_spa_l, ir_spa_l, sp_su_l, av_fe_l
  character(:), allocatable :: lattfile, arg, charge_st_file, field_st_file
  character(:), allocatable :: dir_st_file, dir_dist_file, sphe_sus_file
  character(:), allocatable :: s_ab_file, s_perp_file, field_charge_file
  character(:), allocatable :: energy_file, sq_energy_file, e_field_file
  character(:), allocatable :: irrot_field_file, irrot_sab_file, irrot_sperp_file
  character(:), allocatable :: rot_field_file, rot_sab_file, rot_sperp_file
  character(:), allocatable :: spar_file, rot_spar_file, irrot_spar_file
  character(:), allocatable :: avg_field_file

  !type parameters
    !integer(kind=8), public :: L,add_charges,no_measurements
    !integer(kind=8), public :: therm_sweeps,measurement_sweeps,sample_interval,no_samples
    !real(kind=8), public :: q, lambda, volume, temp, beta
    !real(kind=8), public :: eps_0, bin_size, rot_delt, g0
    !real(kind=8) :: temp, lambda, q, rot_delt, hop_ratio, rot_ratio, g_ratio
    !real(kind=8), dimension(:,:,:,:), allocatable, public :: lgf

    !character(:), allocatable :: lattfile, arg, charge_st_file, field_st_file
    !character(:), allocatable :: dir_st_file, dir_dist_file, sphe_sus_file
    !character(:), allocatable :: s_ab_file, s_perp_file, field_charge_file
    !character(:), allocatable :: energy_file, sq_energy_file, e_field_file
    !character(:), allocatable :: irrot_field_file, irrot_sab_file, irrot_sperp_file
    !character(:), allocatable :: rot_field_file, rot_sab_file, rot_sperp_file
    !character(:), allocatable :: spar_file, rot_spar_file, irrot_spar_file
    !character(:), allocatable :: avg_field_file
  !end type parameters

  !type output

  !  real(kind=16), public :: ener_tot_sum, ener_rot_sum, ener_irrot_sum
  !  real(kind=16), public :: ener_tot_sq_sum, ener_rot_sq_sum, ener_irrot_sq_sum
  !  real(kind=16), dimension(2), public :: ebar_sum, ebar_sq_sum
  !  integer(kind=8), dimension(:), allocatable, public :: bin_count
  !  real(kind=16), dimension(:,:), allocatable, public :: v_avg
  !  real(kind=16), dimension(:,:,:), allocatable, public :: e_tot_avg, e_rot_avg, e_irrot_avg
  !  real(kind=16), dimension(:,:), allocatable, public :: dir_struc
  !  real(kind=16), dimension(:,:,:,:), allocatable, public :: s_ab, s_ab_rot, s_ab_irrot
  !  complex(kind=16), dimension(:,:), allocatable, public :: ch_ch
  !  complex(kind=16), dimension(:,:), allocatable, public :: rho_k_m,rho_k_p
  !  integer(kind=8) :: output_size

  !end type output

  integer, public :: have_lgf = 0
  integer, public, parameter :: bz=2
  real*8, public :: rot_ratio, g_ratio, hop_ratio
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
