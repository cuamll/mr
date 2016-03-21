program config_linsolve

use common
use linear_solver
implicit none
integer :: i, j, k, row, col
character(len=11) :: filename
integer, dimension(:,:), allocatable :: v_temp

! assignment statements from common: these should be read in from a file, really
q = 1.0 ! charge magnitude
lambda = 1.0 ! lattice spacing
volume = lambda**3 ! duh

! note: getarg(pos, value) does argc/argv

! these parameters will be put in an input file eventually
filename='lattice.dat'

! read in lattice size (for now) then allocate arrays accordingly
write(*,*) "Lattice dimension L: "
read(*,*) L
allocate(v_temp(L**2,L))
allocate(v(L,L,L))

! irrotational part of E-field
allocate(mnphi_x(L,L,L))
allocate(mnphi_y(L,L,L))
allocate(mnphi_z(L,L,L))

open(unit=1, file=filename)
read(1,*)((v_temp(row,col),col=1,L),row=1,L**2)

!do i=1,L**2
!  write(*,*) (v_temp(i,j), j=1,L)
!enddo

v = reshape(v_temp, (/L,L,L/))
deallocate(v_temp) ! we don't need it anymore

do j=1,L
  do i=1,L
    write(*,*) (v(i,j,k), k=1,L)
  enddo
enddo

call linsol

end program config_linsolve
