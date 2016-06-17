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
      read(1,*) hop_ratio
      read(1,*) rot_ratio
      read(1,*) g_ratio

      ! Check L is even
      if (modulo(L,2)==1) then
        write(*,*) "L is odd - try again"
        STOP
      end if

      ! set up other variables, allocations etc.
      volume = lambda**3
      beta = 1.0 / temp
      allocate(energy(iterations + 1))
      allocate(energy_run(iterations + 1))
      allocate(sq_energy(iterations + 1))

      write(*,*)
      write(*,*) " --- input parameters ---"
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

    ! guess what this one does
    use common
    implicit none
    character(10) :: energy_filename
    character(11) :: efield_filename
    character(13) :: sq_energy_filename
    integer*8 :: i, j, k, N
    real*8 :: avg_e, avg_e2, prefac, sp_he

    avg_e = 0.0
    avg_e2 = 0.0
    sp_he = 0.0
    N = L**3

    energy_filename = "energy.out"
    sq_energy_filename = "sq_energy.out"
    efield_filename = "e_field.out"

    open(unit=2, file=energy_filename)
    open(unit=3, file=sq_energy_filename)
    open(unit=4, file=efield_filename)

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

    do i = 1,L
      do j = 1,L
        do k = 1,L

        write (4,*) i, j, k, e_x(i,j,k), e_y(i,j,k), e_z(i,j,k)

        end do
      end do
    end do

    !write (*,*) "<U> = ", avg_e / iterations
    !write (*,*) "sum U^2 / iter = ",avg_e2 / iterations
    !write (*,*) "N^2 = ",N**2
    !write (*,*) "N T = ",N * temp

    avg_e = avg_e / (iterations * N)
    avg_e2 = avg_e2 / (iterations * N**2)

    ! prefactor troubles
    prefac = 1.0 * N / (temp**2)

    write(*,*)
    write(*,*) " --- averages and specific heat ---"
    write(*,*) "<U> norm. = ",avg_e
    !write(*,*) "<U>^2 norm. = ",avg_e**2
    !write(*,*) "<U^2> norm. = ",avg_e2
    !write(*,*) "prefactor = ",prefac
    write(*,*)

    sp_he = prefac * ((avg_e2) - ((avg_e)**2))

    write (*,*) "sp. heat (C) = ",sp_he
    !write (*,*) "C / (N) = ",sp_he / (L**3)
    !write (*,*) "C / sqrt(L) = ",sp_he / sqrt(float(L))
    !write(*,*) "N(<U^2> - <U>^2) = ",N * (avg_e2 - ((avg_e)**2))
    write(*,*)

    close(2)
    close(3)

    call correlations

  end subroutine write_output

  subroutine correlations

    ! calculate charge-charge correlations, I guess
    ! Fourier transform the field and then time average
    ! Project out transverse part, it should be flat
    ! Longitudinal part should display pinch points
    ! Width of pinch points --> charge density
    use common
    implicit none
    real*8, dimension(:,:,:), allocatable :: e_kx, e_ky, e_kz
    !real*8, dimension(:,:,:,:), allocatable :: real_vec, k_vec
    integer :: i, j, k, m, n, p
    complex*16 :: imag, kdotx

    imag = (0.0,1.0)

    ! I think L + 1, to allow for components with k_mu = 0?
    ! Need to be careful of K = (0,0,0) though
    ! can i allocate here? not actually sure
    allocate(e_kx(L+1,L+1,L+1))
    allocate(e_ky(L+1,L+1,L+1))
    allocate(e_kz(L+1,L+1,L+1))

    ! Fourier transform
    do i = -L/2, L/2
      do j = -L/2, L/2
        do k = -L/2, L/2

          e_kx(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) = 0.0
          e_ky(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) = 0.0
          e_kz(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) = 0.0

          ! vec(q) = vec(0) term:
          if (i.eq.0.and.j.eq.0.and.k.eq.0) then
            write (*,*) "q=0 term; need to work this out"
            write (*,*) "e_kx(",i,j,k,") = ",e_kx(i,j,k)
            CYCLE
          end if

          ! m,n,p are the real space coordinates
          do m = 1,L
            do n = 1,L
              do p = 1,L

                kdotx = ((-1)*imag*((2*pi*(m-1)*i/(L*lambda)) + &
                        (2*pi*(n-1)*j/(L*lambda)) + &
                        (2*pi*(p-1)*k/(L*lambda))))

                e_kx(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) = e_kx(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) + e**(kdotx)*e_x(m,n,p)
                e_ky(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) = e_ky(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) + e**(kdotx)*e_y(m,n,p)
                e_kz(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) = e_kz(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) + e**(kdotx)*e_z(m,n,p)

              end do
            end do
          end do

          write (*,*) "e_kx(",(2*pi*i/L*lambda),(2*pi*j/L*lambda),(2*pi*k/L*lambda),") = ",e_kx(i,j,k),e_ky(i,j,k),e_kz(i,j,k)

        end do
      end do
    end do

    deallocate(e_kx)
    deallocate(e_ky)
    deallocate(e_kz)

  end subroutine correlations

end module io
