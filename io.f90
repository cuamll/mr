module io
  use common
  implicit none
  logical, private :: start_file_there
  character(len=8) :: in_file

  contains

  subroutine read_input
    in_file='start.in'
    inquire(file=in_file,exist=start_file_there)
    if (start_file_there) then
      open(unit=1, file=in_file)
      read(1,*) L
      read(1,*) lambda
      read(1,*) q
      read(1,*) lattfile
      read(1,*) seed
      ! more things can go here as needed

      volume=lambda**3

      allocate(v(L,L,L))
      allocate(pos(L))
      allocate(neg(L))
      allocate(mnphi_x(L,L,L))
      allocate(mnphi_y(L,L,L))
      allocate(mnphi_z(L,L,L))
      allocate(e_rot_x(L,L,L))
      allocate(e_rot_y(L,L,L))
      allocate(e_rot_z(L,L,L))
      allocate(e_x(L,L,L))
      allocate(e_y(L,L,L))
      allocate(e_z(L,L,L))
      allocate(lgf(L,L,L,L,L,L))

    else
      write (*,*) "can't find input file"
      stop
    end if

    close(1)

  end subroutine read_input

end module io
