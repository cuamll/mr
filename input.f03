module input
  use common
  implicit none
  logical, private :: start_file_there
  character :: corr_char, canon_char
  character(len=200), private :: label, buffer, verb_arg
  integer :: posit
  integer :: ios = 0
  integer :: line = 0

  contains

  subroutine read_input
    use common
    implicit none
    character(len=200) :: lattfile_long, en_long, sq_en_long,&
    e_field_long, arg_long, ch_st_l, fi_st_l, s_ab_l, s_p_l, dir_st_l,&
    dir_d_s_l, fe_ch_l, ir_fe_l, ir_sab_l, ir_sp_l, r_fe_l, r_sab_l,&
    r_sp_l, spa_l, r_spa_l, ir_spa_l, sp_su_l, av_fe_l, eq_l, ch_g, c_ab_l,&
    r_cab_l, ir_cab_l, wf_l, wsq_l, lpl,&
    s_perp_t_l,s_perp_l_l,s_par_t_l,s_par_l_l, s_ab_t_l, s_ab_l_l

    if (command_argument_count().eq.2) then
      call get_command_argument(1, verb_arg)
      if (trim(verb_arg).eq.'-v') then
        verbose = .true.
      else
        write (*,*) "First argument not recognised.&
           & Check and try again."
        stop
      end if
    else if (command_argument_count().gt.2) then
      write (*,*) "Too many arguments given. Check and try again."
      stop
    end if
    ! Final argument should be the input file
    call get_command_argument(command_argument_count(), arg_long)
    arg = trim(arg_long)

    inquire(file=arg,exist=start_file_there)
    if (start_file_there) then
      open(unit=10, file=arg)

      do while (ios == 0)

        if (line.eq.0) then

          if (verbose) then
            write (*,*) "Verbose flag set to true."
            write (*,*)
            write (*,*) "--- INPUT PARAMETERS: ---"
            write (*,*) "Input file: ",arg
          end if

        end if

        read(10, '(A)', iostat=ios) buffer
        if (ios == 0) then
          line = line + 1

          ! Split at whitespace
          posit = scan(buffer, ' ')
          label = buffer(1:posit)
          buffer = buffer(posit + 1:)

          select case (label)
          case ('!')
            cycle
          case ('L')
            read(buffer, '(I10.1)', iostat=ios) L
            if (verbose) then
              write (*,*) 'System size: ',L
            end if
          case ('omp_threads')
            read(buffer, '(I10.1)', iostat=ios) no_threads
            if (verbose) then
              write (*,*) 'OpenMP threads: ',no_threads
            end if
          case ('no_samples')
            read(buffer, '(I10.1)', iostat=ios) no_samples
            if (verbose) then
              write (*,*) 'Number of samples: ',no_samples
            end if
          case ('do_corr')
            read(buffer, '(a)', iostat=ios) corr_char
            if (verbose) then
              write (*,*) 'Calculate correlations? ',corr_char
            end if
          case ('canon')
            read(buffer, '(a)', iostat=ios) canon_char
            if (verbose) then
              write (*,*) 'Canonical? ',canon_char
            end if
          case ('thermalisation_sweeps')
            read(buffer, '(I10.1)', iostat=ios) therm_sweeps
            if (verbose) then
              write (*,*) 'Thermalisation sweeps: ',therm_sweeps
            end if
          case ('measurement_sweeps')
            read(buffer, '(I10.1)', iostat=ios) measurement_sweeps
            if (verbose) then
              write (*,*) 'Measurement sweeps: ',measurement_sweeps
            end if
          case ('sample_interval')
            read(buffer, '(I10.1)', iostat=ios) sample_interval
            if (verbose) then
              write (*,*) 'Sample interval: ',sample_interval
            end if
          case ('temperature')
            read(buffer, '(F10.1)', iostat=ios) temp
            if (verbose) then
              write (*,*) 'Temperature: ',temp
            end if
          case ('lattice_spacing')
            read(buffer, '(F10.1)', iostat=ios) lambda
            if (verbose) then
              write (*,*) 'Lattice spacing: ',lambda
            end if
          case ('charge_value')
            read(buffer, '(F16.1)', iostat=ios) q
            if (verbose) then
              write (*,*) 'Charge value: ',q
            end if
          case ('e_c')
            read(buffer, '(F16.1)', iostat=ios) e_c
            if (verbose) then
              write (*,*) 'Core energy constant: ',e_c
            end if
          case ('delta_max')
            read(buffer, '(F10.1)', iostat=ios) rot_delt
            if (verbose) then
              write (*,*) 'Δ_max for rot. update: ',rot_delt
            end if
          case ('charges')
            read(buffer, '(I10.1)', iostat=ios) add_charges
            if (verbose) then
              write (*,*) 'Charges: ',add_charges
            end if
          case ('rng_seed')
            read(buffer, '(I10.1)', iostat=ios) seed
            if (verbose) then
              write (*,*) 'RNG seed: ',seed
            end if
          case ('hop_ratio')
            read(buffer, '(F10.1)', iostat=ios) hop_ratio
            if (verbose) then
              write (*,*) 'Ratio of hop updates: ',hop_ratio
            end if
          case ('rot_ratio')
            read(buffer, '(F10.1)', iostat=ios) rot_ratio
            if (verbose) then
              write (*,*) 'Ratio of rotational updates: ',rot_ratio
            end if
          case ('harmonic_ratio')
            read(buffer, '(F10.1)', iostat=ios) g_ratio
            if (verbose) then
              write (*,*) 'Ratio of harmonic updates: ',g_ratio
            end if
          case ('bin_size')
            read(buffer, '(F10.1)', iostat=ios) bin_size
            if (verbose) then
              write (*,*) 'Bin size for real space corr.: ',bin_size
            end if
          case ('windings_file')
            read(buffer, '(a)', iostat=ios) wf_l
            if (verbose) then
              write (*,*) 'Windings file name: ',wf_l
            end if
          case ('windings_sq_file')
            read(buffer, '(a)', iostat=ios) wsq_l
            if (verbose) then
              write (*,*) 'Windings file name: ',wsq_l
            end if
          case ('lattice_file')
            read(buffer, '(a)', iostat=ios) lattfile_long
            if (verbose) then
              write (*,*) 'Lattice file name: ',lattfile_long
            end if
          case ('energy_file')
            read(buffer, '(a)', iostat=ios) en_long
            if (verbose) then
              write (*,*) 'Energy file name: ',en_long
            end if
          case ('squared_energy_file')
            read(buffer, '(a)', iostat=ios) sq_en_long
            if (verbose) then
              write (*,*) 'Squared energy file name: ',sq_en_long
            end if
          case ('electric_field_file')
            read(buffer, '(a)', iostat=ios) e_field_long
            if (verbose) then
              write (*,*) 'Electric field file name: ',e_field_long
            end if
          case ('direct_space_structure_factor_file')
            read(buffer, '(a)', iostat=ios) dir_st_l
            if (verbose) then
              write (*,*) 'Direct space structure factor file name: ',dir_st_l
            end if
          case ('direct_space_g(r)_file')
            read(buffer, '(a)', iostat=ios) dir_d_s_l
            if (verbose) then
              write (*,*) 'Direct space g(r) file name: ',dir_d_s_l
            end if
          case ('charge_structure_factor_file')
            read(buffer, '(a)', iostat=ios) ch_st_l
            if (verbose) then
              write (*,*) 'Charge-charge structure factor file name: ',ch_st_l
            end if
          case ('total_field_structure_factor_file')
            read(buffer, '(a)', iostat=ios) fi_st_l
            if (verbose) then
              write (*,*) 'Total field-field structure factor file name: ',fi_st_l
            end if
          case ('total_s^(alpha_beta)_file')
            read(buffer, '(a)', iostat=ios) s_ab_l
            if (verbose) then
              write (*,*) 'Total S^(α β) file name: ',s_ab_l
            end if
          case ('total_chi^(alpha_beta)_file')
            read(buffer, '(a)', iostat=ios) c_ab_l
            if (verbose) then
              write (*,*) 'Total chi^(α β) file name: ',c_ab_l
            end if
          case ('rot_chi^(alpha_beta)_file')
            read(buffer, '(a)', iostat=ios) r_cab_l
            if (verbose) then
              write (*,*) 'Rotational chi^(α β) file name: ',r_cab_l
            end if
          case ('irrot_chi^(alpha_beta)_file')
            read(buffer, '(a)', iostat=ios) ir_cab_l
            if (verbose) then
              write (*,*) 'Irrotational chi^(α β) file name: ',ir_cab_l
            end if
          case ('total_s_(perp)_file')
            read(buffer, '(a)', iostat=ios) s_p_l
            if (verbose) then
              write (*,*) 'Total S_(⊥) file name: ',s_p_l
            end if
          case ('total_s_(par)_file')
            read(buffer, '(a)', iostat=ios) spa_l
            if (verbose) then
              write (*,*) 'Total S_(//) file name: ',spa_l
            end if
          case ('irrot_field_structure_factor_file')
            read(buffer, '(a)', iostat=ios) ir_fe_l
            if (verbose) then
              write (*,*) 'Irrotational field-field structure factor file name: ',ir_fe_l
            end if
          case ('irrot_s^(alpha_beta)_file')
            read(buffer, '(a)', iostat=ios) ir_sab_l
            if (verbose) then
              write (*,*) 'Irrotational S^(α β) file name: ',ir_sab_l
            end if
          case ('irrot_s_(perp)_file')
            read(buffer, '(a)', iostat=ios) ir_sp_l
            if (verbose) then
              write (*,*) 'Irrotational S_(⊥) file name: ',ir_sp_l
            end if
          case ('irrot_s_(par)_file')
            read(buffer, '(a)', iostat=ios) ir_spa_l
            if (verbose) then
              write (*,*) 'Irrotational S_(//) file name: ',ir_spa_l
            end if
          case ('rot_field_structure_factor_file')
            read(buffer, '(a)', iostat=ios) r_fe_l
            if (verbose) then
              write (*,*) 'Rotational field-field structure factor file name: ',r_fe_l
            end if
          case ('rot_s^(alpha_beta)_file')
            read(buffer, '(a)', iostat=ios) r_sab_l
            if (verbose) then
              write (*,*) 'Rotational S^(α β) file name: ',r_sab_l
            end if
          case ('rot_s_(perp)_file')
            read(buffer, '(a)', iostat=ios) r_sp_l
            if (verbose) then
              write (*,*) 'Rotational S_(⊥) file name: ',r_sp_l
            end if
          case ('rot_s_(par)_file')
            read(buffer, '(a)', iostat=ios) r_spa_l
            if (verbose) then
              write (*,*) 'Rotational S_(//) file name: ',r_spa_l
            end if
          case ('field_charge_file')
            read(buffer, '(a)', iostat=ios) fe_ch_l
            if (verbose) then
              write (*,*) 'Raw file name for binary fields/charges output: ',fe_ch_l
            end if
          case ('average_field_file')
            read(buffer, '(a)', iostat=ios) av_fe_l
            if (verbose) then
              write (*,*) 'Raw file name for average fields/charges: ',av_fe_l
            end if
          case ('sphe_sus_file')
            read(buffer, '(a)', iostat=ios) sp_su_l
            if (verbose) then
              write (*,*) 'Raw file name for specific heat/susceptibility: ',sp_su_l
            end if
          case ('lgf_path')
            read(buffer, '(a)', iostat=ios) lpl
            if (verbose) then
              write (*,*) 'Path to LGF binary file:',lpl
            end if
          case ('equil_file')
            read(buffer, '(a)', iostat=ios) eq_l
            if (verbose) then
              write (*,*) 'Raw file name for equilibriation check: ',eq_l
            end if
          case ('s_ab_t_file')
            read(buffer, '(a)', iostat=ios) s_ab_t_l
            if (verbose) then
              write (*,*) 'Raw file name for s_perp transverse file: ',s_ab_t_l
            end if
          case ('s_ab_l_file')
            read(buffer, '(a)', iostat=ios) s_ab_l_l
            if (verbose) then
              write (*,*) 'Raw file name for s_perp long. file: ',s_ab_l_l
            end if
          case ('s_perp_t_file')
            read(buffer, '(a)', iostat=ios) s_perp_t_l
            if (verbose) then
              write (*,*) 'Raw file name for s_perp transverse file: ',s_perp_t_l
            end if
          case ('s_perp_l_file')
            read(buffer, '(a)', iostat=ios) s_perp_l_l
            if (verbose) then
              write (*,*) 'Raw file name for s_perp long. file: ',s_perp_l_l
            end if
          case ('s_par_t_file')
            read(buffer, '(a)', iostat=ios) s_par_t_l
            if (verbose) then
              write (*,*) 'Raw file name for s_perp transverse file: ',s_par_t_l
            end if
          case ('s_par_l_file')
            read(buffer, '(a)', iostat=ios) s_par_l_l
            if (verbose) then
              write (*,*) 'Raw file name for s_perp long. file: ',s_par_l_l
            end if
          case ('charge_generation')
            read(buffer, '(a)', iostat=ios) ch_g
            if (verbose) then
              write (*,*) 'Method of charge generation:',ch_g
            end if
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
      if (rot_delt.lt.0.0) then
        write(*,*) "Delta_max for rotational update must be >= 0.&
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
      if (corr_char.eq.'T'.or.corr_char.eq.'Y') then
        do_corr = .true.
      else
        do_corr = .false.
      end if
      if (canon_char.eq.'T'.or.canon_char.eq.'Y') then
        canon = .true.
      else
        canon = .false.
      end if

      ! --- NOTE TO SELF ---
      ! is the dimensional analysis sorted out?
      eps_0 = 1.0 / (2*pi)
      ! eps_0 = 1.0
      !q = 2 * pi * q
      !write (*,*) "q = ",q
      beta = 1.0 / temp
      g_thr = 1 / real(L)

      if (rot_delt.eq.0) then

        if (temp.lt.0.1) then
          rot_delt = 2 * temp
        else if (temp.gt.10.0) then
          rot_delt = sqrt(temp)
        else
          rot_delt = 1.1 * temp
        end if

        if (verbose) then
          write (*,*) "Delta_max read in as 0; being set to",rot_delt
        end if
      end if

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
      sphe_sus_file = trim(adjustl(sp_su_l))
      equil_file = trim(adjustl(eq_l))
      charge_gen = trim(adjustl(ch_g))
      chi_ab_file = trim(adjustl(c_ab_l))
      rot_chi_ab_file = trim(adjustl(r_cab_l))
      irrot_chi_ab_file = trim(adjustl(ir_cab_l))
      windings_file = trim(adjustl(wf_l))
      windings_sq_file = trim(adjustl(wsq_l))
      s_ab_transverse_file = trim(adjustl(s_ab_t_l))
      s_ab_long_file = trim(adjustl(s_ab_l_l))
      s_perp_transverse_file = trim(adjustl(s_perp_t_l))
      s_perp_long_file = trim(adjustl(s_perp_l_l))
      s_par_transverse_file = trim(adjustl(s_par_t_l))
      s_par_long_file = trim(adjustl(s_par_l_l))
      lgf_path = trim(adjustl(lpl))

      if (charge_gen.ne."RANDOM".and.charge_gen.ne."DIPOLE".and.charge_gen.ne."CRYSTA") then
        write(*,*) "Charge generation method should be either RANDOM,&
          & DIPOLE, or CRYSTAL. Edit input file and try again."
        STOP
      end if


    else
      write (*,'(a)',advance='no') "Can't find an input file at ",arg
      write (*,*) " . Check path and try again."
      stop
    end if

    close(10)

  end subroutine read_input

end module input
