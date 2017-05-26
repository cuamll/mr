module io
  use common
  implicit none
  logical, private :: start_file_there

  contains

  subroutine read_input

    call getarg(1, arg_long)
    arg = trim(arg_long)

    inquire(file=arg,exist=start_file_there)
    if (start_file_there) then
      open(unit=1, file=arg)
      read(1,'(I10.1)') L
      read(1,'(I10.1)') therm_sweeps
      read(1,'(I10.1)') measurement_sweeps
      read(1,'(I10.1)') sample_interval
      read(1,'(F10.1)') temp
      read(1,'(F10.1)') lambda
      read(1,'(F10.1)') q
      read(1,'(F10.1)') rot_delt
      read(1,'(I10.1)') add_charges
      read(1,'(I10.1)') seed
      read(1,'(F10.1)') hop_ratio
      read(1,'(F10.1)') rot_ratio
      read(1,'(F10.1)') g_ratio
      read(1,'(F10.1)') bin_size
      read(1,'(a)') lattfile_long
      read(1,'(a)') en_long
      read(1,'(a)') sq_en_long
      read(1,'(a)') e_field_long
      read(1,'(a)') ch_st_l
      read(1,'(a)') fi_st_l
      read(1,'(a)') dir_st_l
      read(1,'(a)') dir_d_s_l
      read(1,'(a)') s_ab_l
      read(1,'(a)') s_p_l

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

      write (*,*)
      write (*,*) "--- INPUT PARAMETERS: ---"
      write (*,*) 'L = ',L
      write (*,*) 'thermalisation sweeps = ',therm_sweeps
      write (*,*) 'measurement sweeps = ',measurement_sweeps
      write (*,*) 'sample interval = ',sample_interval
      write (*,*) "no. of measurements = measurement sweeps"&
      &"/ sample interval = ",no_measurements
      write (*,*) 'T = ',temp
      write (*,*) 'lattice spacing = ',lambda
      write (*,*) 'charge value = ',q
      write (*,*) '\Delta_{max} for rotational update = ',rot_delt
      write (*,*) 'number of charges (0 reads in lattice file) = ',add_charges
      write (*,*) 'RNG seed = ',seed
      write (*,*) 'ratio of charge hop updates = ',hop_ratio
      write (*,*) 'ratio of rotational updates = ',rot_ratio
      write (*,*) 'ratio of harmonic mode updates = ',g_ratio
      write (*,*) 'bin size for radial correlation function = ',bin_size
      write (*,*) 'lattice file to read in (if necessary): ',lattfile
      write (*,*) 'energy file: ',energy_file
      write (*,*) 'squared energy file: ',sq_energy_file
      write (*,*) 'E-field file: ',e_field_file
      write (*,*) 'charge struc file: ',charge_st_file
      write (*,*) 'field struc file: ',field_st_file
      write (*,*) 'direct space charge struc file: ',dir_st_file
      write (*,*) 'direct space corr. by distance: ',dir_dist_file
      write (*,*) 'S^{alpha beta} file: ',s_ab_file
      write (*,*) 'S_perp file: ',s_perp_file

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
