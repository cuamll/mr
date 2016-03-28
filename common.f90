module common
  implicit none
  real, public :: q, lambda, volume, ebar_x, ebar_y, ebar_z
  integer, public :: L
  integer, dimension(:), allocatable, public :: pos,neg
  integer, dimension(:,:,:), allocatable, public :: v
  real*8, dimension(:,:,:), allocatable, public :: mnphi_x, mnphi_y, mnphi_z
  ! probably more things need to go here

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

end module common
