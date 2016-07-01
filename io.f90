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
      read(1,*) add_charges
      read(1,*) seed
      read(1,*) hop_ratio
      read(1,*) rot_ratio
      read(1,*) g_ratio

      ! Check L is even
      if (modulo(L,2)==1) then
        write(*,*) "L is odd - try again"
        STOP
      end if
      if (modulo(add_charges,2)==1) then
        write(*,*) "can't add an odd number of charges - try again"
        STOP
      end if

      ! set up other variables, allocations etc.
      volume = lambda**3
      ! beta has dimensions of [e_0 * length]
      ! need to figure out how that works exactly
      eps_0 = 1.0
      beta = 1.0 * L / temp
      allocate(energy(iterations + 1))
      allocate(energy_run(iterations + 1))
      allocate(sq_energy(iterations + 1))

      write (*,*)
      write (*,*) " --- input parameters ---"
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
    !complex*16, dimension(:,:,:), allocatable :: e_kx, e_ky, e_kz
    real*8 :: dot, dot_avg
    integer :: i, j, k, m, n, p
    complex*16 :: imag, kdotx
    character(7) :: struc_filename
    character(74) :: struc_format_string

    !format_string = "(ES18.9, ES18.9, ES18.9, ES18.9,&
    !                & ES18.9, ES18.9, ES18.9, ES18.9, ES18.9)"
    struc_format_string = "(ES18.9, ES18.9, ES18.9, ES18.9)"

    struc_filename = "s_q.out"

    do i = 1,L + 1
      do j = 1,L + 1
        do k = 1,L + 1
          do n = 1,iterations

            struc(i,j,k) = struc(i,j,k) + ch_ch(i,j,k,n)

          end do
        end do
      end do
    end do
    
    struc = struc / iterations

    open(unit=5, file=struc_filename)

    do i = 1,L + 1
      do j = 1,L + 1
        do k = 1,L + 1

          ! output is kx, ky, kz, S(\vec{k})
          write (5, struc_format_string)&
          2*pi*(i - 1 - L/2)/L*lambda,&
          2*pi*(j - 1 - L/2)/L*lambda,&
          2*pi*(k - 1 - L/2)/L*lambda,&
          struc(i,j,k)

        end do
      end do
    end do

    close(5)

    !imag = (0.0,1.0)

    ! I think L + 1, to allow for components with k_mu = 0?
    ! Need to be careful of K = (0,0,0) though

    ! Fourier transform
    !do i = -L/2, L/2
    !  do j = -L/2, L/2
    !    do k = -L/2, L/2

    !      e_kx(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) = 0.0
    !      e_ky(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) = 0.0
    !      e_kz(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) = 0.0

    !      ! vec(q) = vec(0) term:
    !      if (i.eq.0.and.j.eq.0.and.k.eq.0) then
    !        write (*,*) "q=0 term; need to work this out"
    !        write (*,*) "e_kx(",i,j,k,") = ",e_kx(i,j,k)
    !        CYCLE
    !      end if

    !      ! m,n,p are the real space coordinates
    !      do m = 1,L
    !        do n = 1,L
    !          do p = 1,L

    !            kdotx = ((-1)*imag*2*pi*(((m-1)*i/(L*lambda)) + &
    !                    ((n-1)*j/(L*lambda)) + &
    !                    ((p-1)*k/(L*lambda))))

    !            e_kx(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) = &
    !                        e_kx(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) &
    !                        + e**(kdotx)*e_x(m,n,p)
    !            e_ky(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) = &
    !                        e_ky(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) &
    !                        + e**(kdotx)*e_y(m,n,p)
    !            e_kz(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) = &
    !                        e_kz(i + 1 + L/2,j + 1 + L/2,k + 1 + L/2) &
    !                        + e**(kdotx)*e_z(m,n,p)

    !          end do
    !        end do
    !      end do

    !      ! now we need to do a (?) thermal (?) average to get correlations
    !      ! this is the dot of the fourier space field with itself at i,j,k
    !      dot = (e_kx(i + 1 + L/2, j + 1 + L/2, k + 1 + L/2)&
    !            * CONJG(e_kx(i + 1 + L/2, j + 1 + L/2, k + 1 + L/2))) +&
    !            (e_ky(i + 1 + L/2, j + 1 + L/2, k + 1 + L/2) *&
    !            CONJG(e_ky(i + 1 + L/2, j + 1 + L/2, k + 1 + L/2))) +&
    !            (e_kz(i + 1 + L/2, j + 1 + L/2, k + 1 + L/2) *&
    !            CONJG(e_kz(i + 1 + L/2, j + 1 + L/2, k + 1 + L/2)))

    !      !write (*,*) "E \cdot E (",2*pi*i/(L*lambda),2*pi*j/(L*lambda),2*pi*k/(L*lambda),") = ",dot

    !      dot_avg = dot_avg + dot

    !    end do
    !  end do
    !end do

    ! need to check this is being done properly
    ! this isn't what we actually we want to do
    !dot_avg = dot_avg / (L + 1)**3
    !write (*,*) "avg. E dot E over k = ",dot_avg

    !open(unit=5, file=e_k_filename)

    !do i = 1,L + 1
    !  do j = 1,L + 1
    !    do k = 1,L + 1

    !      ! output is (kx, ky, kz), real(e_kx), imag(e_kx), etc.
    !      write (5, format_string)&
    !      2*pi*(i - 1 - L/2)/L*lambda,&
    !      2*pi*(j - 1 - L/2)/L*lambda,&
    !      2*pi*(k - 1 - L/2)/L*lambda,&
    !      REAL(e_kx(i,j,k)), AIMAG(e_kx(i,j,k)), REAL(e_ky(i,j,k)),&
    !      AIMAG(e_ky(i,j,k)), REAL(e_kz(i,j,k)), AIMAG(e_kz(i,j,k))

    !    end do
    !  end do
    !end do

    !close(5)
    !deallocate(e_kx)
    !deallocate(e_ky)
    !deallocate(e_kz)

  end subroutine correlations

end module io
