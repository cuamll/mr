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
      read(1,*) iterations
      read(1,*) temp
      read(1,*) lambda
      read(1,*) q
      read(1,*) rot_delt
      read(1,*) lattfile
      read(1,*) seed
      read(1,*) rot_ratio
      read(1,*) g_ratio
      ! more things can go here as needed

      volume = lambda**3
      beta = 1.0 / temp

      write (*,*) 'L = ',L
      write (*,*) 'T = ',temp
      write (*,*) 'rot. update ratio = ',rot_ratio

    else
      write (*,*) "can't find input file"
      stop
    end if

    close(1)

  end subroutine read_input

end module io
