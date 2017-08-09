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
          case ('direct_space_structure_factor_file')
            read(buffer, '(a)', iostat=ios) dir_st_l
            write (*,*) 'Direct space structure factor file name: ',dir_st_l
          case ('direct_space_g(r)_file')
            read(buffer, '(a)', iostat=ios) dir_d_s_l
            write (*,*) 'Direct space g(r) file name: ',dir_d_s_l
          case ('charge_structure_factor_file')
            read(buffer, '(a)', iostat=ios) ch_st_l
            write (*,*) 'Charge-charge structure factor file name: ',ch_st_l
          case ('total_field_structure_factor_file')
            read(buffer, '(a)', iostat=ios) fi_st_l
            write (*,*) 'Total field-field structure factor file name: ',fi_st_l
          case ('total_s^(alpha_beta)_file')
            read(buffer, '(a)', iostat=ios) s_ab_l
            write (*,*) 'Total S^(α β) file name: ',s_ab_l
          case ('total_s_(perp)_file')
            read(buffer, '(a)', iostat=ios) s_p_l
            write (*,*) 'Total S_(⊥) file name: ',s_p_l
          case ('total_s_(par)_file')
            read(buffer, '(a)', iostat=ios) spa_l
            write (*,*) 'Total S_(//) file name: ',spa_l
          case ('irrot_field_structure_factor_file')
            read(buffer, '(a)', iostat=ios) ir_fe_l
            write (*,*) 'Irrotational field-field structure factor file name: ',ir_fe_l
          case ('irrot_s^(alpha_beta)_file')
            read(buffer, '(a)', iostat=ios) ir_sab_l
            write (*,*) 'Irrotational S^(α β) file name: ',ir_sab_l
          case ('irrot_s_(perp)_file')
            read(buffer, '(a)', iostat=ios) ir_sp_l
            write (*,*) 'Irrotational S_(⊥) file name: ',ir_sp_l
          case ('irrot_s_(par)_file')
            read(buffer, '(a)', iostat=ios) ir_spa_l
            write (*,*) 'Irrotational S_(//) file name: ',ir_spa_l
          case ('rot_field_structure_factor_file')
            read(buffer, '(a)', iostat=ios) r_fe_l
            write (*,*) 'Rotational field-field structure factor file name: ',r_fe_l
          case ('rot_s^(alpha_beta)_file')
            read(buffer, '(a)', iostat=ios) r_sab_l
            write (*,*) 'Rotational S^(α β) file name: ',r_sab_l
          case ('rot_s_(perp)_file')
            read(buffer, '(a)', iostat=ios) r_sp_l
            write (*,*) 'Rotational S_(⊥) file name: ',r_sp_l
          case ('rot_s_(par)_file')
            read(buffer, '(a)', iostat=ios) r_spa_l
            write (*,*) 'Rotational S_(//) file name: ',r_spa_l
          case ('field_charge_file')
            read(buffer, '(a)', iostat=ios) fe_ch_l
            write (*,*) 'Raw file name for binary fields/charges output: ',fe_ch_l
          case ('average_field_file')
            read(buffer, '(a)', iostat=ios) av_fe_l
            write (*,*) 'Raw file name for average fields/charges: ',av_fe_l
          case ('vertex_ebar_file')
            read(buffer, '(a)', iostat=ios) v_s_l
            write (*,*) 'Raw file name for vertices/ebar: ',v_s_l
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
      volume = lambda**2
      no_measurements = measurement_sweeps / sample_interval
      ! --- NOTE TO SELF ---
      ! is the dimensional analysis sorted out?
      ! eps_0 = 1.0 / L
      eps_0 = 1.0
      q = 2 * pi * q
      write (*,*) "q = ",q
      beta = 1.0 / temp

      lattfile = trim(adjustl(lattfile_long))
      energy_file = trim(adjustl(en_long))
      sq_energy_file = trim(adjustl(sq_en_long))
      e_field_file = trim(adjustl(e_field_long))
      dir_st_file = trim(adjustl(dir_st_l))
      dir_dist_file = trim(adjustl(dir_d_s_l))
      charge_st_file = trim(adjustl(ch_st_l))
      field_st_file = trim(adjustl(fi_st_l))
      s_ab_file = trim(adjustl(s_ab_l))
      s_perp_file = trim(adjustl(s_p_l))
      spar_file = trim(adjustl(spa_l))
      irrot_field_file = trim(adjustl(ir_fe_l))
      irrot_sab_file = trim(adjustl(ir_sab_l))
      irrot_sperp_file = trim(adjustl(ir_sp_l))
      irrot_spar_file = trim(adjustl(ir_spa_l))
      rot_field_file = trim(adjustl(r_fe_l))
      rot_sab_file = trim(adjustl(r_sab_l))
      rot_sperp_file = trim(adjustl(r_sp_l))
      rot_spar_file = trim(adjustl(r_spa_l))
      field_charge_file = trim(adjustl(fe_ch_l))
      avg_field_file = trim(adjustl(av_fe_l))
      vertex_sum_file = trim(adjustl(v_s_l))

    else
      write (*,'(a)',advance='no') "Can't find an input file at ",arg
      write (*,*) " . Check path and try again."
      stop
    end if

    close(10)

  end subroutine read_input

  subroutine write_output

    ! guess what this one does
    use common
    implicit none
    integer*8 :: i, j, k
    real*8 :: avg_e, avg_e2, prefac, sp_he

    avg_e = 0.0
    avg_e2 = 0.0
    sp_he = 0.0

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

      do j = 1,L
        do i = 1,L

        write (4,*) i, j, e_field(1,i,j),&
                  & e_field(2,i,j)

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
    prefac = 1.0 * L**2 / (temp**2)
    sp_he = prefac * ((avg_e2) - ((avg_e)**2))

    write(*,*) "prefactor = ",prefac
    write(*,*) "sp. heat (C) = ",sp_he
    write(*,*) "C / (N) = ",sp_he / L**2
    write(*,*)

    avg_e = avg_e / (L**2)
    avg_e2 = avg_e2 / ((L**2)**2)

    write(*,*) "<U> norm. = ",avg_e
    write(*,*) "<U>^2 norm. = ",avg_e**2
    write(*,*) "<U^2> norm. = ",avg_e2
    write(*,*) "prefactor = ",prefac
    write(*,*)

    sp_he = prefac * ((avg_e2) - ((avg_e)**2))

    write(*,*) "sp. heat (C) = ",sp_he
    write(*,*) "C / (N) = ",sp_he / L**2
    write(*,*) "C / sqrt(L) = ",sp_he / sqrt(float(L))
    write(*,*) "N(<U^2> - <U>^2) = ",L**2 * (avg_e2 - ((avg_e)**2))
    write(*,*)

    close(2)
    close(3)
    close(4)

    call calc_correlations(field_charge_file)

  end subroutine write_output

  subroutine calc_correlations(filename)
    use common
    implicit none
    character(len = *), intent(in) :: filename
    integer :: i,j,k,n,kx,ky,kz,m,p,s,x,y,z,dist_bin,sp
    real(kind=16) :: norm_k,kx_float,ky_float,kz_float,dist
    real(kind=16) :: sp_he_tot, sp_he_rot, sp_he_irrot, prefac
    real(kind=16), dimension(:,:), allocatable :: s_perp, s_perp_irrot, s_perp_rot
    real(kind=16), dimension(:,:), allocatable :: s_par, s_par_irrot, s_par_rot
    real(kind=16), dimension((bz*L)+1,(bz*L)+1) :: charge_struc, field_struc
    real(kind=16), dimension((bz*L)+1,(bz*L)+1) :: field_struc_irrot
    real(kind=16), dimension((bz*L)+1,(bz*L)+1) :: field_struc_rot
    character(100) :: struc_format_string, field_format_string,vertex_format
    character(100) :: avg_field_format,dir_format_string, dir_dist_format_string

    avg_field_format = "(I3, I3, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9)"
    field_format_string = "(ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9)"
    struc_format_string = "(ES18.9, ES18.9, ES18.9, ES18.9)"
    dir_format_string = "(I3, I3, ES18.9)"
    vertex_format = "(I6, I3, I3, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, ES18.9, I3)"
    dir_dist_format_string = "(ES18.9, I9.6, ES18.9)"

    field_struc = 0.0; field_struc_rot = 0.0;
    charge_struc = 0.0; field_struc_irrot = 0.0;

    ebar_sum = (ebar_sum / no_measurements)
    ebar_sq_sum = (ebar_sq_sum / no_measurements)

    rho_k_p = rho_k_p / no_measurements
    rho_k_m = rho_k_m / no_measurements
    ch_ch = ch_ch / no_measurements
    s_ab = s_ab / no_measurements
    s_ab_irrot = s_ab_irrot / no_measurements
    s_ab_rot = s_ab_rot / no_measurements
    dir_struc = dir_struc / no_measurements
    e_tot_avg = e_tot_avg / no_measurements
    e_rot_avg = e_rot_avg / no_measurements
    e_irrot_avg = e_irrot_avg / no_measurements
    v_avg = v_avg / no_measurements

    ener_tot_sum = ener_tot_sum / no_measurements
    ener_rot_sum = ener_rot_sum / no_measurements
    ener_irrot_sum = ener_irrot_sum / no_measurements
    ener_tot_sq_sum = ener_tot_sq_sum / no_measurements
    ener_rot_sq_sum = ener_rot_sq_sum / no_measurements
    ener_irrot_sq_sum = ener_irrot_sq_sum / no_measurements

    !ener_tot_sum = ener_tot_sum / L**2
    !ener_rot_sum = ener_rot_sum / L**2
    !ener_irrot_sum = ener_irrot_sum / L**2
    !ener_tot_sq_sum = ener_tot_sq_sum / ((L**2)**2)
    !ener_rot_sq_sum = ener_rot_sq_sum / ((L**2)**2)
    !ener_irrot_sq_sum = ener_irrot_sq_sum / ((L**2)**2)

    prefac = 1.0 * L**2 / (temp**2)

    sp_he_tot = prefac * (ener_tot_sq_sum - (ener_tot_sum)**2)
    sp_he_rot = prefac * (ener_rot_sq_sum - (ener_rot_sum)**2)
    sp_he_irrot = prefac * (ener_irrot_sq_sum - (ener_irrot_sum)**2)

    write(*,*)
    write (*,'(a)') "   # Specific heat: total, rot., irrot."
    write (*,'(ES18.9, ES18.9, ES18.9, ES18.9)') temp,&
    sp_he_tot, sp_he_rot, sp_he_irrot
    write(*,*) "E^bar averages:"
    write(*,'(a)') "<E^bar_x>, <E^bar_y>, <|E^bar|>"
    write (*,*) ebar_sum(1), ebar_sum(2),&
    sqrt(ebar_sum(1)**2+ebar_sum(2)**2)
    write(*,'(a)') "<(E^bar_x)^2>, <(E^bar_y)^2>, <|(E^bar)^2|>"
    write (*,*) ebar_sq_sum(1), ebar_sq_sum(2),&
    sqrt(ebar_sq_sum(1)**2+ebar_sq_sum(2)**2)

    write(*,'(a)') "<|(E^bar)^2|> - <|E^bar|>^2"
    write (*,*) L**2 * beta * (ebar_sq_sum(1) + ebar_sq_sum(2)&
    - (ebar_sum(1)**2 + ebar_sum(2)**2))

    write (30,'(a)') "   # Specific heat: total, rot., irrot."
    write (30,'(ES18.9, ES18.9, ES18.9, ES18.9)') temp,&
    sp_he_tot, sp_he_rot, sp_he_irrot

    write (30,'(a)') "   # <ebar^2> - <ebar>^2"
    write (30,'(ES18.9, ES18.9)') temp,&
    L**2 * beta * (ebar_sq_sum(1) + ebar_sq_sum(2)&
    - (ebar_sum(1)**2 + ebar_sum(2)**2))

    write (30,'(a)') "   # hop acceptance"
    write (30,'(ES18.9, ES18.9)') temp,&
    (float(accepth) / &
    ((therm_sweeps + measurement_sweeps) * add_charges * hop_ratio))

    close(30)

    ! we can calculate s_perp up to wherever
    sp = 8
    allocate(s_perp((sp*L)+1,(sp*L)+1))
    allocate(s_perp_irrot((sp*L)+1,(sp*L)+1))
    allocate(s_perp_rot((sp*L)+1,(sp*L)+1))
    allocate(s_par((sp*L)+1,(sp*L)+1))
    allocate(s_par_irrot((sp*L)+1,(sp*L)+1))
    allocate(s_par_rot((sp*L)+1,(sp*L)+1))
    s_perp = 0.0; s_perp_rot = 0.0; s_perp_irrot = 0.0;
    s_par = 0.0; s_par_rot = 0.0; s_par_irrot = 0.0;

    do p = (-L/2)*sp,(L/2)*sp
      do m = (-L/2)*sp,(L/2)*sp

        i = m + 1 + sp*(L/2)
        j = p + 1 + sp*(L/2)

        if (j.le.(bz*L + 1).and.i.le.(bz*L + 1)) then
          ! can also subtract e.g. rho_k_p * conjg(rho_k_m)
          charge_struc(i,j) = abs(ch_ch(i,j) - rho_k_p(i,j) * conjg(rho_k_m(i,j)))
          field_struc(i,j) = abs(s_ab(1,1,i,j))
          field_struc_irrot(i,j) = abs(s_ab_irrot(1,1,i,j))
          field_struc_rot(i,j) = abs(s_ab_rot(1,1,i,j))
        end if

        ! use separate variables, we're gonna mess around with values
        kx_float = m * ((2 * pi)/(L * lambda))
        ky_float = p * ((2 * pi)/(L * lambda))

        ! --- NOTE TO SELF ---
        ! need to sort this shit out. this could be a start????
        !gx = floor((kx_float * L) / (2 * pi))
        !gy = floor((ky_float * L) / (2 * pi))
        !qx = mod(kx_float, 2 * pi)
        !qy = mod(ky_float, 2 * pi)

        if (kx_float.eq.0.and.ky_float.eq.0) then
          norm_k = 0.0
        else
          norm_k = 1.0/(kx_float**2 + ky_float**2)
        end if

        ! move back to somewhere we already know about
        ! this is probably not right yet.
        ! but i know something needs doing here
        ! s_ab should be periodic in 2pi so mod should be fine

        if (abs(p).gt.(L/2)*bz) then
          ky = mod(p,(bz*L)/2) + 1 + (bz*L)/2
        else
          ky = p + 1 + (bz*L)/2
        end if
        if (abs(m).gt.(L/2)*bz) then
          kx = mod(m,(bz*L)/2) + 1 + (bz*L)/2
        else
          kx = m + 1 + (bz*L)/2
        end if

        !if (abs(ky_float).gt.(pi)*bz) then
        !  ky_float = mod(ky_float,bz * pi)
        !end if
        !if (abs(kx_float).gt.(pi)*bz) then
        !  kx_float = mod(kx_float,bz * pi)
        !end if

        s_perp(i,j) = (1 - kx_float*kx_float*norm_k) * s_ab(1,1,kx,ky)+&
                        ((-1)*kx_float*ky_float*norm_k) * s_ab(1,2,kx,ky)+&
                        ((-1)*ky_float*kx_float*norm_k) * s_ab(2,1,kx,ky)+&
                        (1 - ky_float*ky_float*norm_k) * s_ab(2,2,kx,ky)

        s_perp_irrot(i,j) = (1 - kx_float*kx_float*norm_k) * s_ab_irrot(1,1,kx,ky)+&
                        ((-1)*kx_float*ky_float*norm_k) * s_ab_irrot(1,2,kx,ky)+&
                        ((-1)*ky_float*kx_float*norm_k) * s_ab_irrot(2,1,kx,ky)+&
                        (1 - ky_float*ky_float*norm_k) * s_ab_irrot(2,2,kx,ky)

        s_perp_rot(i,j) = (1 - kx_float*kx_float*norm_k) * s_ab_rot(1,1,kx,ky)+&
                        ((-1)*kx_float*ky_float*norm_k) * s_ab_rot(1,2,kx,ky)+&
                        ((-1)*ky_float*kx_float*norm_k) * s_ab_rot(2,1,kx,ky)+&
                        (1 - ky_float*ky_float*norm_k) * s_ab_rot(2,2,kx,ky)

        s_par(i,j) = (kx_float*kx_float*norm_k) * s_ab(1,1,kx,ky)+&
                       (kx_float*ky_float*norm_k) * s_ab(1,2,kx,ky)+&
                       (ky_float*kx_float*norm_k) * s_ab(2,1,kx,ky)+&
                       (ky_float*ky_float*norm_k) * s_ab(2,2,kx,ky)

        s_par_irrot(i,j) = (kx_float*kx_float*norm_k)*s_ab_irrot(1,1,kx,ky)+&
                             (kx_float*ky_float*norm_k)*s_ab_irrot(1,2,kx,ky)+&
                             (ky_float*kx_float*norm_k)*s_ab_irrot(2,1,kx,ky)+&
                             (ky_float*ky_float*norm_k)*s_ab_irrot(2,2,kx,ky)

        s_par_rot(i,j) = (kx_float*kx_float*norm_k)*s_ab_rot(1,1,kx,ky)+&
                           (kx_float*ky_float*norm_k)*s_ab_rot(1,2,kx,ky)+&
                           (ky_float*kx_float*norm_k)*s_ab_rot(2,1,kx,ky)+&
                           (ky_float*ky_float*norm_k)*s_ab_rot(2,2,kx,ky)

      end do
    end do ! end p, m loops

    open(unit=10, file=dir_st_file)
    open(unit=11, file=dir_dist_file)
    open(unit=12, file=charge_st_file)
    open(unit=13, file=field_st_file)
    open(unit=14, file=s_ab_file)
    open(unit=15, file=s_perp_file)
    open(unit=16, file=irrot_field_file)
    open(unit=17, file=irrot_sab_file)
    open(unit=18, file=irrot_sperp_file)
    open(unit=19, file=rot_field_file)
    open(unit=20, file=rot_sab_file)
    open(unit=21, file=rot_sperp_file)
    open(unit=22, file=spar_file)
    open(unit=23, file=irrot_spar_file)
    open(unit=24, file=rot_spar_file)
    open(unit=25, file=avg_field_file)

    dist_r = dist_r / (no_measurements)
    do i = 1,ceiling( sqrt(float((3*((L/2)**2)))) * (1 / bin_size) )
        write (11, dir_dist_format_string)&
        i * bin_size, bin_count(i), abs(dist_r(i))
    end do

    close(11)

      do i = 1,sp*(L) + 1
        do j = 1,sp*(L) + 1

          ! output is kx, ky, kz, S(\vec{k})

          if (i.le.L/2+1.and.j.le.L/2+1) then
            write (10, dir_format_string)&
            i - 1,j - 1,abs(dir_struc(i,j))
          end if

          if (i.le.L.and.j.le.L) then
            write (25, avg_field_format)&
            i,j,e_tot_avg(1,i,j),e_tot_avg(2,i,j),&
            e_rot_avg(1,i,j),e_rot_avg(2,i,j),&
            e_irrot_avg(1,i,j),e_irrot_avg(2,i,j),&
            v_avg(i,j)
          end if

          if (j.le.(bz*L + 1).and.i.le.(bz*L + 1)) then

            write (12, struc_format_string)&
            2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
            charge_struc(i,j)

            write (13, struc_format_string)&
            2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
            field_struc(i,j)

            write (14, field_format_string)&
            2*pi*(i - 1 - bz*(l/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(l/2))/(L*lambda),&
            s_ab(1,1,i,j),&
            s_ab(1,2,i,j),&
            s_ab(2,1,i,j),&
            s_ab(2,2,i,j)

            write (16, struc_format_string)&
            2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
            field_struc_irrot(i,j)

            write (17, field_format_string)&
            2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
            s_ab_irrot(1,1,i,j),&
            s_ab_irrot(1,2,i,j),&
            s_ab_irrot(2,1,i,j),&
            s_ab_irrot(2,2,i,j)

            write (19, struc_format_string)&
            2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
            field_struc_rot(i,j)

            write (20, field_format_string)&
            2*pi*(i - 1 - bz*(l/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(l/2))/(L*lambda),&
            s_ab_rot(1,1,i,j),&
            s_ab_rot(1,2,i,j),&
            s_ab_rot(2,1,i,j),&
            s_ab_rot(2,2,i,j)

          end if

            write (15, field_format_string)&
            2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
            s_perp(i,j)

            write (18, field_format_string)&
            2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
            s_perp_irrot(i,j)

            write (21, field_format_string)&
            2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
            s_perp_rot(i,j)

            write (22, field_format_string)&
            2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
            s_par(i,j)

            write (23, field_format_string)&
            2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
            s_par_irrot(i,j)

            write (24, field_format_string)&
            2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
            s_par_rot(i,j)

        end do
      end do

    close(10)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)
    close(22)
    close(23)
    close(24)

    deallocate(s_perp); deallocate(s_perp_rot); deallocate(s_perp_irrot);
    deallocate(s_par); deallocate(s_par_rot); deallocate(s_par_irrot);

  end subroutine calc_correlations

end module io
