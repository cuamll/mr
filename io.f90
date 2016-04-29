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

      ! set up other variables, allocations etc.
      volume = lambda**3
      beta = 1.0 / temp
      allocate(energy(iterations + 1))
      allocate(energy_run(iterations + 1))
      allocate(sq_energy(iterations + 1))

      write (*,*) 'L = ',L
      write (*,*) 'T = ',temp
      write (*,*) 'iter = ',iterations
      write (*,*) 'rot. ratio = ',rot_ratio

    else
      write (*,*) "can't find input file"
      stop
    end if

    close(1)

  end subroutine read_input

  subroutine write_output
    use common
    implicit none
    character(10) :: energy_filename
    character(13) :: sq_energy_filename
    integer :: i
    real*8 :: avg_e, avg_e2, prefac, sp_he

    ! guess what this one does
    energy_filename = "energy.out"
    sq_energy_filename = "sq_energy.out"

    open(unit=2, file=energy_filename)
    open(unit=3, file=sq_energy_filename)

    ! write out the parameters in the energy file
    write (2,*) "L",L
    write (2,*) "T",temp
    write (2,*) "iter",iterations
    write (2,*) "rot. ratio",rot_ratio
    write (2,*) "ebar ratio",g_ratio

    do i = 1,iterations + 1
      write(2,*) i - 1, energy(i)
      write(3,*) i - 1, sq_energy(i)
      avg_e = avg_e + energy(i)
      avg_e2 = avg_e2 + sq_energy(i)
    end do

    avg_e = avg_e / (iterations * L**3)
    avg_e2 = avg_e2 / (iterations * L**3 * L**3)

    ! this prefactor is wrong somehow - gas constant???
    prefac = (L**3) / (temp**2)
    sp_he = prefac * (avg_e2 - (avg_e * avg_e))

    write (*,*) "specific heat = ",sp_he

    close(2)
    close(3)

  end subroutine write_output

end module io
