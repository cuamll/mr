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
      eps_0 = 1.0 / L
      beta = 1.0 / temp
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

    write(*,*) "sp. heat (C) = ",sp_he
    write(*,*) "C / (N) = ",sp_he / N
    write(*,*) "C / sqrt(L) = ",sp_he / sqrt(float(L))
    write(*,*) "N(<U^2> - <U>^2) = ",N * (avg_e2 - ((avg_e)**2))
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
    real*8 :: dot, dot_avg, knorm
    integer :: i, j, k, m, n, p, kx, ky, kz
    complex*16 :: imag, kdotx
    character(16) :: charge_struc_filename
    character(15) :: field_struc_filename
    character(10) :: s_perp_filename
    character(74) :: struc_format_string
    character(97) :: field_format_string

    field_format_string = "(ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9)"
    struc_format_string = "(ES18.9, ES18.9, ES18.9, ES18.9)"

    charge_struc_filename = "charge_struc.out"
    field_struc_filename = "field_struc.out"
    s_perp_filename = "s_perp.out"

    do i = 1,bz*(L+1)
      do j = 1,bz*(L+1)
        do k = 1,bz*(L+1)
          do n = 1,iterations

            !charge_struc(i,j,k) = charge_struc(i,j,k) +&
            !                      0.5 * (ch_ch(i,j,k,n) - 2 * ch_ch_pp(i,j,k,n))
            field_struc(i,j,k) = field_struc(i,j,k) + fe_fe(i,j,k,n)
            charge_struc(i,j,k) = charge_struc(i,j,k) +&
              abs(abs(ch_ch(i,j,k,n)) - abs(rho_k_p_t(i,j,k,n))&
              *abs(rho_k_m_t(i,j,k,n)))

            ! s_ab averaging
            do m = 1,3
              do p = 1,3
                s_ab(m,p,i,j,k) = s_ab(m,p,i,j,k) + s_ab_n(m,p,i,j,k,n)
              end do
            end do

            !s_ab(1,1,i,j,k) = s_ab(1,1,i,j,k) + s_ab_n(1,1,i,j,k,n)

            if (n.eq.iterations) then
              ! s_ab projection for perpendicular component
              ! for n = iterations we've just finished the sum
              ! keep this in the loop bc we still need dummy variables

              ! normalise s_ab
              s_ab = s_ab / iterations

              ! find kx, ky, kz from i,j,k and their norm
              kx = i - 1 - bz*(L/2)
              ky = j - 1 - bz*(L/2)
              kz = k - 1 - bz*(L/2)
              if (kx.eq.0.and.ky.eq.0.and.kz.eq.0) then
                knorm = 0
              else
                knorm = 1.0/(kx*kx + ky*ky + kz*kz)
              end if

              ! sum over alpha, beta
              s_perp(i,j,k) = (1 - kx*kx*knorm)*s_ab(1,1,i,j,k) +&
                              ((-1)*kx*ky*knorm)*s_ab(1,2,i,j,k)+&
                              ((-1)*kx*kz*knorm)*s_ab(1,3,i,j,k)+&
                              ((-1)*ky*kx*knorm)*s_ab(2,1,i,j,k)+&
                              (1 - ky*ky*knorm)*s_ab(2,2,i,j,k) +&
                              ((-1)*ky*kz*knorm)*s_ab(2,3,i,j,k)+&
                              ((-1)*kz*kx*knorm)*s_ab(3,1,i,j,k)+&
                              ((-1)*kz*ky*knorm)*s_ab(3,2,i,j,k)+&
                              (1 - kz*kz*knorm)*s_ab(3,3,i,j,k)

            end if ! if (n.eq.iterations) for s_ab and s_perp
          end do ! n iteration loop
        end do ! k
      end do ! j
    end do ! i

    charge_struc = charge_struc / iterations
    field_struc = field_struc / iterations

    open(unit=5, file=charge_struc_filename)
    open(unit=7, file=field_struc_filename)
    open(unit=9, file=s_perp_filename)

    do k = 1,bz*(L) + 1
      do i = 1,bz*(L) + 1
        do j = 1,bz*(L) + 1

          ! output is kx, ky, kz, S(\vec{k})
          write (5, struc_format_string)&
          2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
          2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
          2*pi*(k - 1 - bz*(L/2))/(L*lambda),&
          charge_struc(i,j,k)

          write (7, field_format_string)&
          2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
          2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
          2*pi*(k - 1 - bz*(L/2))/(L*lambda),&
          field_struc(i,j,k)

          write (9, field_format_string)&
          2*pi*(i - 1 - bz*(l/2))/(l*lambda),&
          2*pi*(j - 1 - bz*(l/2))/(l*lambda),&
          2*pi*(k - 1 - bz*(l/2))/(l*lambda),&
          s_perp(i,j,k)

        end do
      end do
    end do

    close(5)
    close(7)
    close(9)

  end subroutine correlations

end module io
