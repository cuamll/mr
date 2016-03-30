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

  call PBCs

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

  call randinit(seed)
  write(*,*) rand(seed) 

  call linsol

  stop

end program mr

subroutine update(iter)
  integer :: i,x,y,z
  real :: chooser
  real*8 :: old_e,new_e,delta_e
  logical :: x_choice, y_choice1, y_choice2, z_choice

  ! charge hop updates
  do i=1,L**3
    ! one sweep is L**3 attempts
    x=int(rand()*L)+1
    y=int(rand()*L)+1
    z=int(rand()*L)+1

    chooser=rand()
    if (chooser.lt.(1.0/3.0)) then
      ! x component
      ! pick x component of field here
      chooser=rand()
      if (chooser.lt.(0.5)) then
        ! decrease bond
      else
        ! increase bond
      end if
    else if (chooser.ge.(2.0/3.0)) then
      ! y component
    else
      ! z component
    end if
  end do


  ! plaquette rot. update



  ! e bar update probably not needed, vanishes



end subroutine update
