program mr

  use common
  use linear_solver
  use io
  implicit none
  integer :: i, j, k, row, col
  integer, dimension(:,:), allocatable :: v_temp

  ebar_x = 0.0
  ebar_y = 0.0
  ebar_z = 0.0

  call read_input

  write(*,*) 'L = ',L
  write(*,*) 'lambda = ',lambda
  write(*,*) 'q = ',q
  write(*,*) 'lattice config file = ',lattfile

  allocate(v_temp(L**2,L))
  allocate(v(L,L,L))
  allocate(pos(L))
  allocate(neg(L))

  ! irrotational part of E-field
  allocate(mnphi_x(L,L,L))
  allocate(mnphi_y(L,L,L))
  allocate(mnphi_z(L,L,L))
  allocate(lgf(L,L,L,L,L,L))

  open(unit=2, file=lattfile)
  read(2,*)((v_temp(row,col),col=1,L),row=1,L**2)

  v = reshape(v_temp, (/L,L,L/), ORDER = (/2,3,1/))
  deallocate(v_temp) ! we don't need it anymore
  close(2)

  call linsol

  stop

end program mr
