module io
  use common
  implicit none
  logical, private :: start_file_there
  character(len = 100) :: label, buffer
  integer :: posit
  integer :: ios = 0
  integer :: line = 0

  contains

  subroutine read_input

    call getarg(1, arg_long)
    arg = trim(arg_long)

    inquire(file=arg,exist=start_file_there)
    if (start_file_there) then
      write (*,*)
      write (*,*) "--- INPUT PARAMETERS: ---"
      write (*,*) "Input file: ",arg

      open(unit=10, file=arg)

      do while (ios == 0)
        read(10, '(A)', iostat=ios) buffer
        if (ios == 0) then
          line = line + 1

          ! Split at whitespace
          posit = scan(buffer, ' ')
          label = buffer(1:posit)
          buffer = buffer(posit + 1:)

          select case (label)
          case ('L')
            read(buffer, '(I10.1)', iostat=ios) L
            write (*,*) 'System size: ',L
          case ('thermalisation_sweeps')
            read(buffer, '(I10.1)', iostat=ios) therm_sweeps
            write (*,*) 'Thermalisation sweeps: ',therm_sweeps
          case ('measurement_sweeps')
            read(buffer, '(I10.1)', iostat=ios) measurement_sweeps
            write (*,*) 'Measurement sweeps: ',measurement_sweeps
          case ('sample_interval')
            read(buffer, '(I10.1)', iostat=ios) sample_interval
            write (*,*) 'Sample interval: ',sample_interval
          case ('temperature')
            read(buffer, '(F10.1)', iostat=ios) temp
            write (*,*) 'Temperature: ',temp
          case ('lattice_spacing')
            read(buffer, '(F10.1)', iostat=ios) lambda
            write (*,*) 'Lattice spacing: ',lambda
          case ('charge_value')
            read(buffer, '(F10.1)', iostat=ios) q
            write (*,*) 'Charge value: ',q
          case ('delta_max')
            read(buffer, '(F10.1)', iostat=ios) rot_delt
            write (*,*) 'Δ_max for rot. update: ',rot_delt
          case ('charges')
            read(buffer, '(I10.1)', iostat=ios) add_charges
            write (*,*) 'Charges: ',add_charges
          case ('rng_seed')
            read(buffer, '(I10.1)', iostat=ios) seed
            write (*,*) 'RNG seed: ',seed
          case ('hop_ratio')
            read(buffer, '(F10.1)', iostat=ios) hop_ratio
            write (*,*) 'Ratio of hop updates: ',hop_ratio
          case ('rot_ratio')
            read(buffer, '(F10.1)', iostat=ios) rot_ratio
            write (*,*) 'Ratio of rotational updates: ',rot_ratio
          case ('harmonic_ratio')
            read(buffer, '(F10.1)', iostat=ios) g_ratio
            write (*,*) 'Ratio of harmonic updates: ',g_ratio
          case ('bin_size')
            read(buffer, '(F10.1)', iostat=ios) bin_size
            write (*,*) 'Bin size for real space corr.: ',bin_size
          case ('lattice_file')
            read(buffer, '(a)', iostat=ios) lattfile_long
            write (*,*) 'Lattice file name: ',lattfile_long
          case ('energy_file')
            read(buffer, '(a)', iostat=ios) en_long
            write (*,*) 'Energy file name: ',en_long
          case ('squared_energy_file')
            read(buffer, '(a)', iostat=ios) sq_en_long
            write (*,*) 'Squared energy file name: ',sq_en_long
          case ('electric_field_file')
            read(buffer, '(a)', iostat=ios) e_field_long
            write (*,*) 'Electric field file name: ',e_field_long
          case ('charge_structure_factor_file')
            read(buffer, '(a)', iostat=ios) ch_st_l
            write (*,*) 'Charge-charge structure factor file name: ',ch_st_l
          case ('field_structure_factor_file')
            read(buffer, '(a)', iostat=ios) fi_st_l
            write (*,*) 'Field-field structure factor file name: ',fi_st_l
          case ('direct_space_structure_factor_file')
            read(buffer, '(a)', iostat=ios) dir_st_l
            write (*,*) 'Direct space structure factor file name: ',dir_st_l
          case ('direct_space_g(r)_file')
            read(buffer, '(a)', iostat=ios) dir_d_s_l
            write (*,*) 'Direct space g(r) file name: ',dir_d_s_l
          case ('s^(alpha_beta)_file')
            read(buffer, '(a)', iostat=ios) s_ab_l
            write (*,*) 'S^(α β) file name: ',s_ab_l
          case ('s_(perp)_file')
            read(buffer, '(a)', iostat=ios) s_p_l
            write (*,*) 'S_(⊥) file name: ',s_p_l
          case ('field_charge_file')
            read(buffer, '(a)', iostat=ios) fe_ch_l
            write (*,*) 'Charge distribution file name: ',fe_ch_l
          case default
            write (*,*) 'Skipping invalid label at line ',line
          end select
        end if
      end do

      ! Check parameters are physically reasonable
      if (modulo(L,2)==1) then
        write(*,*) "L must be even. Edit input file and try again."
        STOP
      end if
      if (L.le.0) then
        write(*,*) "L must be > 0. Edit input file and try again."
        STOP
      end if
      if (therm_sweeps < 0) then
        write(*,*) "Number of thermalisation sweeps must be > 0.&
          & Edit input file and try again."
        STOP
      end if
      if (measurement_sweeps < 0) then
        write(*,*) "Number of measurement sweeps must be > 0.&
          & Edit input file and try again."
        STOP
      end if
      if (sample_interval < 1) then
        write(*,*) "Sample interval must be ≥ 1.&
          & Edit input file and try again."
        STOP
      end if
      if (temp.le.0.0) then
        write(*,*) "Temperature must be > 0.&
          & Edit input file and try again."
        STOP
      end if
      if (lambda.le.0.0) then
        write(*,*) "Lattice spacing must be > 0.&
          & Edit input file and try again."
        STOP
      end if
      if (q.le.0.0) then
        write(*,*) "Charge magnitude must be > 0.&
          & Edit input file and try again."
        STOP
      end if
      if (rot_delt.le.0.0) then
        write(*,*) "Delta_max for rotational update must be > 0.&
          & Edit input file and try again."
        STOP
      end if
      if (modulo(add_charges,2)==1) then
        write(*,*) "Can't add an odd number of charges; &
          &system must be neutral. Edit input file and try again."
        STOP
      end if
      if (hop_ratio.le.0.0) then
        write(*,*) "Ratio of charge hop updates must be > 0.&
          & Edit input file and try again."
        STOP
      end if
      if (rot_ratio.le.0.0) then
        write(*,*) "Ratio of rotational updates must be > 0.&
          & Edit input file and try again."
        STOP
      end if
      if (g_ratio.le.0.0) then
        write(*,*) "Ratio of harmonic updates must be > 0.&
          & Edit input file and try again."
        STOP
      end if
      if (bin_size.le.0.0) then
        write(*,*) "Bin size for correlation function must be > 0.&
          & Edit input file and try again."
        STOP
      end if

      ! set up other variables, allocations etc.
      volume = lambda**3
      no_measurements = measurement_sweeps / sample_interval
      ! --- NOTE TO SELF ---
      ! is the dimensional analysis sorted out?
      eps_0 = 1.0 / L
      beta = 1.0 / temp

      lattfile = trim(lattfile_long)
      energy_file = trim(en_long)
      sq_energy_file = trim(sq_en_long)
      e_field_file = trim(e_field_long)
      charge_st_file = trim(ch_st_l)
      field_st_file = trim(fi_st_l)
      charge_st_file = trim(ch_st_l)
      field_st_file = trim(fi_st_l)
      dir_st_file = trim(dir_st_l)
      dir_dist_file = trim(dir_d_s_l)
      s_ab_file = trim(s_ab_l)
      s_perp_file = trim(s_p_l)
      field_charge_file = trim(fe_ch_l)

    else
      write (*,'(a)',advance='no') "Can't find an input file at ",arg
      write (*,*) " . Check path and try again."
      stop
    end if

    close(1)

  end subroutine read_input

  subroutine write_output

    ! guess what this one does
    use common
    implicit none
    integer*8 :: i, j, k, N
    real*8 :: avg_e, avg_e2, prefac, sp_he

    avg_e = 0.0
    avg_e2 = 0.0
    sp_he = 0.0
    N = L**3

    open(unit=2, file=energy_file)
    open(unit=3, file=sq_energy_file)
    open(unit=4, file=e_field_file)

    ! write out the parameters in the energy file
    write (2,*) "L",L
    write (2,*) "T",temp
    write (2,*) "number of measurements",no_measurements
    write (2,*) "rot. ratio",rot_ratio
    write (2,*) "ebar ratio",g_ratio

    do i = 1,no_measurements + 1
      write(2,*) i - 1, energy(i)
      write(3,*) i - 1, sq_energy(i)

      avg_e = avg_e + energy(i)
      avg_e2 = avg_e2 + sq_energy(i)
    end do

    do k = 1,L
      do j = 1,L
        do i = 1,L

        write (4,*) i, j, k, e_field(1,i,j,k),&
                  & e_field(2,i,j,k), e_field(3,i,j,k)

        end do
      end do
    end do

    !write (*,*) "<U> = ", avg_e / no_measurements
    !write (*,*) "sum U^2 / iter = ",avg_e2 / no_measurements
    !write (*,*) "N^2 = ",N**2
    !write (*,*) "N T = ",N * temp

    avg_e = avg_e / (no_measurements)
    avg_e2 = avg_e2 / (no_measurements)

    write(*,*)
    write(*,*) " --- averages and specific heat ---"
    write(*,*) "avg_e unnormalised = ",avg_e
    write(*,*) "avg_e^2 unnormalised = ",avg_e2

    ! prefactor troubles
    prefac = 1.0 * N / (temp**2)
    sp_he = prefac * ((avg_e2) - ((avg_e)**2))

    write(*,*) "prefactor = ",prefac
    write(*,*) "sp. heat (C) = ",sp_he
    write(*,*) "C / (N) = ",sp_he / N
    write(*,*)

    avg_e = avg_e / (N)
    avg_e2 = avg_e2 / (N**2)

    write(*,*) "<U> norm. = ",avg_e
    write(*,*) "<U>^2 norm. = ",avg_e**2
    write(*,*) "<U^2> norm. = ",avg_e2
    write(*,*) "prefactor = ",prefac
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
    integer :: i, j, k, kx, ky, kz, dist_bin
    real*8 :: knorm, dist
    real*8, dimension((bz*L)+1,(bz*L)+1,(bz*L)+1) :: s_perp
    !real*8, dimension((bz*L)+1,(bz*L)+1,(bz*L)+1) :: charge_struc,field_struc
    real*8, dimension(ceiling(sqrt(float(3*(((L/2)**2))))*(1 / bin_size))) :: dist_r
    integer, dimension(ceiling(sqrt(float(3*(((L/2)**2))))*(1 / bin_size))) :: bin_count
    character(74) :: struc_format_string
    character(97) :: field_format_string
    character(21) :: dir_format_string
    character(22) :: dir_dist_format_string

    knorm = 0.0; dist_bin = 0; dist = 0.0; dist_r = 0.0; bin_count = 0.0
    s_perp = 0.0
    !charge_struc = 0.0; field_struc = 0.0

    field_format_string = "(ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9)"
    struc_format_string = "(ES18.9, ES18.9, ES18.9, ES18.9)"
    dir_format_string = "(I2, I2, I2, ES18.9)"
    dir_dist_format_string = "(ES18.9, I8, ES18.9)"

    ! normalise
    write (*,*) "ch_ch unnormalised",ch_ch(L,L + 1,L + 1),ch_ch(L + 1,L + 1,L + 1)
    ! delete these! just for testing
    !rho_k_m = rho_k_m * 2
    !rho_k_p = rho_k_p * 2

    rho_k_p = rho_k_p / (no_measurements)
    rho_k_m = rho_k_m / (no_measurements)
    ch_ch = ch_ch / (no_measurements)
    e_kx = e_kx / (no_measurements)
    fe_fe = fe_fe / (no_measurements)
    s_ab = s_ab / (no_measurements)
    !charge_struc = charge_struc / (no_measurements)
    !field_struc = field_struc / (no_measurements)
    !write (*,*) "ch_ch normalised",ch_ch(L,L + 1,L + 1),ch_ch(L + 1,L + 1,L + 1)
    !write (*,*) "rho_k_pm normalised at 0,0,0",rho_k_p(L + 1,L + 1,L + 1),rho_k_m(L + 1,L + 1,L + 1)
    !write (*,*) "charge_struc normalised from main",charge_struc(L,L + 1,L + 1),charge_struc(L + 1,L + 1,L + 1)

    ! what's going on with this?
    ! only difference is this is real space instead of Fourier
    dir_struc = dir_struc / (no_measurements&
              * (0.5 * add_charges / L**3)**2)
    !charge_struc = 0.0; field_struc = 0.0

    open(unit=5, file=charge_st_file)
    open(unit=7, file=field_st_file)
    open(unit=9, file=s_perp_file)
    open(unit=10, file=dir_st_file)
    open(unit=11, file=dir_dist_file)
    open(unit=12, file=s_ab_file)

    write (*,*) "Writing out to files..."

    do k = 1,bz*(L) + 1
      do j = 1,bz*(L) + 1
        do i = 1,bz*(L) + 1

            !field_struc(i,j,k) = field_struc(i,j,k) +&
            !  abs(fe_fe(i,j,k) - e_kx(i,j,k)&
            !  *e_kx(i,j,k))

            !charge_struc(i,j,k) = charge_struc(i,j,k) +&
            !  abs(ch_ch(i,j,k) - rho_k_p(i,j,k)&
            !  *rho_k_m(i,j,k))

            if (k.eq.(L+1).and.j.eq.(L+1).and.i.eq.(L+1)) then
              write (*,*) "charge_struc calculating at 0,0,0",ch_ch(i,j,k),rho_k_p(i,j,k),rho_k_m(i,j,k),&
                ch_ch(i,j,k) - rho_k_p(i,j,k) * rho_k_m(i,j,k),charge_struc(i,j,k)
            end if

            kx = i - 1 - bz*(L/2)
            ky = j - 1 - bz*(L/2)
            kz = k - 1 - bz*(L/2)

            if (kx.eq.0.and.ky.eq.0.and.kz.eq.0) then
              knorm = 0.0
            else
              knorm = 1.0/(kx*kx + ky*ky + kz*kz)
            end if

            if (i.eq.(L-2).and.j.eq.(L+1).and.k.eq.(L+1)) then
              write (*,*) 2*pi*kx/L,2*pi*ky/L,2*pi*kz/L,s_ab(:,:,i,j,k)
            end if

            ! sum over alpha, beta
            s_perp(i,j,k) = (1 - kx*kx*knorm) * s_ab(1,1,i,j,k) +&
                            ((-1)*kx*ky*knorm) * s_ab(1,2,i,j,k)+&
                            ((-1)*kx*kz*knorm) * s_ab(1,3,i,j,k)+&
                            ((-1)*ky*kx*knorm) * s_ab(2,1,i,j,k)+&
                            (1 - ky*ky*knorm) * s_ab(2,2,i,j,k) +&
                            ((-1)*ky*kz*knorm) * s_ab(2,3,i,j,k)+&
                            ((-1)*kz*kx*knorm) * s_ab(3,1,i,j,k)+&
                            ((-1)*kz*ky*knorm) * s_ab(3,2,i,j,k)+&
                            (1 - kz*kz*knorm) * s_ab(3,3,i,j,k)

          ! direct space one
          if (i.le.L/2+1.and.j.le.L/2+1.and.k.le.L/2+1) then

            dist = sqrt(dble((i - 1)**2 + (j - 1)**2 + (k - 1)**2))
            dist_bin = floor(dist / bin_size) + 1
            dist_r(dist_bin) = dist_r(dist_bin) + dir_struc(i,j,k)
            bin_count(dist_bin) = bin_count(dist_bin) + 1

          end if

        end do
      end do
    end do

    do k = 1,bz*(L) + 1
      do i = 1,bz*(L) + 1
        do j = 1,bz*(L) + 1

          ! output is kx, ky, kz, S(\vec{k})

          write (5, struc_format_string)&
          2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
          2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
          2*pi*(k - 1 - bz*(L/2))/(L*lambda),&
          charge_struc(i,j,k)

          ! write direct space one
          if (i.le.L/2+1.and.j.le.L/2+1.and.k.le.L/2+1) then
            write (10, dir_format_string)&
            i - 1,j - 1,k - 1,abs(dir_struc(i,j,k))
          end if

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

          write (12, field_format_string)&
          2*pi*(i - 1 - bz*(l/2))/(l*lambda),&
          2*pi*(j - 1 - bz*(l/2))/(l*lambda),&
          2*pi*(k - 1 - bz*(l/2))/(l*lambda),&
          s_ab(1,1,i,j,k),&
          s_ab(1,2,i,j,k),&
          s_ab(1,3,i,j,k),&
          s_ab(2,1,i,j,k),&
          s_ab(2,2,i,j,k),&
          s_ab(2,3,i,j,k),&
          s_ab(3,1,i,j,k),&
          s_ab(3,2,i,j,k),&
          s_ab(3,3,i,j,k)

        end do
      end do
    end do
    write (*,*) "charge_struc recalculated",charge_struc(L,L + 1,L + 1),charge_struc(L + 1,L + 1,L + 1)

    dist_r = dist_r / (no_measurements)
    do i = 1,ceiling( sqrt(float((3*((L/2)**2)))) * (1 / bin_size) )
        write (11, dir_dist_format_string)&
        i * bin_size, bin_count(i), abs(dist_r(i))
    end do

    close(5)
    close(7)
    close(9)
    close(10)
    close(11)
    close(12)

  end subroutine correlations

end module io
