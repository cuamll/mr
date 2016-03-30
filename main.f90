program config_linsolve

  use common
  use linear_solver
  implicit none
  integer :: i, j, k, row, col
  character(len=11) :: filename
  integer, dimension(:,:), allocatable :: v_temp

  ! assignment statements from common:
  ! these should be read in from a file, really
  q = 1.0 ! charge magnitude
  lambda = 1.0 ! lattice spacing
  volume = lambda**3 ! duh

  ebar_x = 0.0
  ebar_y = 0.0
  ebar_z = 0.0

  ! note: getarg(pos, value) does argc/argv

  ! these parameters will be put in an input file eventually
  filename='lattice.dat'

  ! read in lattice size (for now) then allocate arrays accordingly
  write(*,*) "Lattice dimension L: "
  read(*,*) L
  allocate(v_temp(L**2,L))
  allocate(v(L,L,L))
  allocate(pos(L))
  allocate(neg(L))

  ! irrotational part of E-field
  allocate(mnphi_x(L,L,L))
  allocate(mnphi_y(L,L,L))
  allocate(mnphi_z(L,L,L))
  allocate(lgf(L,L,L,L,L,L))

  open(unit=1, file=filename)
  read(1,*)((v_temp(row,col),col=1,L),row=1,L**2)

  v = reshape(v_temp, (/L,L,L/), ORDER = (/2,3,1/))
  deallocate(v_temp) ! we don't need it anymore

  call linsol

end program config_linsolve
