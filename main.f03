! maggs-rossetto algorithm on cubic lattice
! i'll write proper docs at some point
! callum gray, UCL / ENS Lyon, 2016. do wat u like
program mr

  use common
  use linear_solver
  use input
  use output
  use setup
  use update
  use omp_lib
  use mpi
  implicit none
  include 'rev.inc'
  integer(kind=4) :: i, j, k, tot_q, rank, num_procs, mpierr
  real(kind=rk) :: ener_recv
  real(kind=8) :: start_time, end_time
  real(kind=8), dimension(8) :: timings
  integer, dimension(8) :: values

  call MPI_Init(mpierr)
  call MPI_TYPE_CREATE_F90_REAL(prec, expo, MPI_NEW_REAL, mpierr)
  call MPI_TYPE_CREATE_F90_INTEGER(prec, MPI_NEW_INT, mpierr)

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, num_procs, mpierr)

  write (*,'(a,i2.1,a,i2.1)') "Proc. number ", rank, " out of ",num_procs

  if (rank.eq.0) then
    call cpu_time(start_time)

    write (*,*) "Maggs-Rossetto algorithm on the simple cubic lattice."
    write (*,*) "Callum Gray, University College London, 2017."
    write (*,*) "Git revision ",revision
    write (*,*) "Repo at http://github.com/cuamll/mr"

    call date_and_time(VALUES=values)

    write (*,'(a24,I2.2,a1,I2.2,a1,I4.2,a1,I2.2,a1,I2.2,a1,I2.2)')&
      &" Current date and time: ",values(3),"/",values(2),"/"&
      &,values(1)," ",values(5),":",values(6),":",values(7)
    !write (*,'(a,i2.1)') "Max OpenMP threads: ",omp_get_max_threads()
  end if

  call read_input
  if (no_threads.ne.0) then
    call omp_set_num_threads(no_threads)
  end if
  if (rank.eq.0) then
    write (*,'(a,i2.1)') "Max. OpenMP threads: ",omp_get_max_threads()
  end if
  call initial_setup

  do k = 1, no_samples

    call cpu_time(timings(1))

    ! every proc/sample combination should have a different seed,
    ! otherwise we're averaging identical data. i choose to start @ 0
    call setup_wrapper((rank * no_samples) + k - 1)

    call cpu_time(timings(2))

    if (verbose) then
      write (*,'(a,i2.1,a,i2.1,a,i2.1,a,i3.1)') "Proc. ", rank,&
        " out of ",num_procs," starting sample ",k,&
        ". RNG seed = ",((rank * no_samples) + (k - 1))
    end if

    do i = 1,therm_sweeps
      call mc_sweep
    end do

    call cpu_time(timings(3))

    do i = 1,measurement_sweeps
      call mc_sweep

      if (mod(i, sample_interval).eq.0) then
        call measure(i)
      end if

    end do

    if (verbose) then
      write (*,'(a,i2.1,a,i2.1,a,i2.1)') "Proc. ", rank, " out of ",&
        num_procs," finished measurements for sample ",k
    end if

    call cpu_time(timings(4))

    tot_q = sum(v,mask=v.gt.0) ! if this =/= 0 we need to solve Poisson eq.
    if (tot_q.ne.0) then
      call linsol
      tot_q = tot_q * 2 ! previous tot_q def. is the number of charge pairs
    else ! this is probably not needed?
      tot_q = 0
    end if

    call deallocations

    call cpu_time(timings(5))

    if (verbose) then
      write (*,'(a,i2.1,a,i2.1,a,i4.1,a,es18.9)') "Rank: ",rank,&
            " sample: ",k," seed: ",(rank * no_samples) + k - 1,&
            " energy: ",ener_tot_sum
      write (*,'(a,f8.3,a)') "Time taken: ",(timings(5) - timings(1))," seconds."
      write (*,'(a,f8.3,a)') "Setup time: ",(timings(2) - timings(1))," seconds."
      write (*,'(a,f8.5,a)') "Time per therm sweep: "&
        &,(timings(3) - timings(2)) / therm_sweeps," seconds."
      write (*,'(a,f8.5,a)') "Time per measurement sweep: ",&
        &(timings(4) - timings(3)) / measurement_sweeps," seconds."
    end if

  end do ! end k loop: samples

  write (*,*) "Calling MPI barrier from process ",rank
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  call reductions(rank)

  if (rank.eq.0) then

    if (verbose) then
      write (*,*) "--- MPI_Reduce done: unnormalised results ---"
      write (*,*) "Total: ",ener_tot_sum
      write (*,*) "Rotational: ",ener_rot_sum
      write (*,*) "Irrotational: ",ener_irrot_sum
      write (*,*) "Total^2: ",ener_tot_sq_sum
      write (*,*) "Rotational^2: ",ener_rot_sq_sum
      write (*,*) "Irrotational^2: ",ener_irrot_sq_sum
      write (*,*) "Ebar sum: ",ebar_sum(1), ebar_sum(2)
      write (*,*) "Ebar^2 sum: ",ebar_sq_sum(1), ebar_sq_sum(2)
      write (*,*) "Current susc: ",((L**2 * beta) *&
                  (sum(ebar_sq_sum) - sum(ebar_sum * ebar_sum)))&
                  / (no_measurements * no_samples)
      write (*,*) "Avg. charge density: ",rho_avg
      write (*,*)
    end if

    call normalisations(num_procs)

    ! if (do_corr) then
      call write_output
    ! end if

    call cpu_time(end_time)

    write (*,'(a,f10.3,a)') "Simulation finished. Time taken:",&
                            end_time-start_time,"seconds."
    write (*,'(a,i4.1,a,i2.1,a,f10.3,a,a,f10.3,a)') "Total samples: ",&
      no_samples * num_procs, " across ",num_procs,&
      " processes. Time per sample: ",&
      (end_time - start_time) / (no_samples * num_procs)," seconds.",&
      " Time per measurement: ",&
      (end_time - start_time) /&
      (no_samples * num_procs * no_measurements)," seconds."

    ! These probably need a little more thought
    ! and reducing!
    ! write (*,'(a,i12.1,es18.9)') "Charge hops: total, acceptance rate: ",&
    !   accepth, accepth / ((therm_sweeps + measurement_sweeps) * &
    !   hop_ratio * L**2 * (num_procs * no_samples))
    ! write (*,'(a,i12.1,es18.9)') "Plaquette update: total, acceptance rate: ",&
    !   acceptr, acceptr / ((therm_sweeps + measurement_sweeps) * &
    !   rot_ratio * L**2 * (num_procs * no_samples))
    ! write (*,'(a,i12.1,es18.9)') "Harmonic update: total, acceptance rate: ",&
    !   acceptg, acceptg / ((therm_sweeps + measurement_sweeps) * &
    !   g_ratio * L**2 * 2 * (num_procs * no_samples))
    ! write (*,'(a,a,6i12.1)') "Hops, creations, annihilations ",&
    !   "(attempts, accepts): ",&
    !   attempth, accepth, attemptc, acceptc, attempta, accepta
    write (*,'(a,2i12.1,es18.9)') "Hops: total, attempts, rate: ",&
    accepts(1), attempts(1), dble(accepts(1)) / dble(attempts(1))
    write (*,'(a,2i12.1,es18.9)') "Rot.: total, attempts, rate: ",&
    accepts(2), attempts(2), dble(accepts(2)) / dble(attempts(2))
    write (*,'(a,2i12.1,es18.9)') "Harm: total, attempts, rate: ",&
    accepts(3), attempts(3), dble(accepts(3)) / dble(attempts(3))
    write (*,'(a,2i12.1,es18.9)') "Creations: total, attempts, rate: ",&
    accepts(4), attempts(4), dble(accepts(4)) / dble(attempts(4))
    write (*,'(a,2i12.1,es18.9)') "Annihilations: total, attempts, rate: ",&
    accepts(5), attempts(5), dble(accepts(5)) / dble(attempts(5))

  end if

  call MPI_Finalize(mpierr)

  stop

end program mr

subroutine mc_sweep
  use common
  use update
  implicit none

  call hop(int(L**2 * hop_ratio))

  call rot(int(L**2 * rot_ratio))

  if (glob.eq.1) then
    call harm(int(L**2 * g_ratio))
  end if

end subroutine mc_sweep

subroutine reductions(id)
  use common
  use mpi
  implicit none
  integer, intent(in) :: id
  integer :: mpierr

  ! write (*,*) "PRE: rank", id, "e_tot_sum", ener_tot_sum

  if (id.eq.0) then

    call MPI_Reduce(MPI_IN_PLACE, accepts, 5, MPI_NEW_INT,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(MPI_IN_PLACE, attempts, 5, MPI_NEW_INT,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(MPI_IN_PLACE, ener_tot_sum, 1, MPI_NEW_REAL,&
                       MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(MPI_IN_PLACE, ener_rot_sum, 1, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(MPI_IN_PLACE, ener_irrot_sum, 1, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(MPI_IN_PLACE, ener_tot_sq_sum, 1, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(MPI_IN_PLACE, ener_rot_sq_sum, 1, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(MPI_IN_PLACE, ener_irrot_sq_sum, 1, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(MPI_IN_PLACE, ebar_sum, 2, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(MPI_IN_PLACE, ebar_sq_sum, 2, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(MPI_IN_PLACE, ebar_dip_sum, 2, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(MPI_IN_PLACE, ebar_dip_sq_sum, 2, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(MPI_IN_PLACE, ebar_wind_sum, 2, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(MPI_IN_PLACE, ebar_wind_sq_sum, 2, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(MPI_IN_PLACE, rho_avg, 2, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)

  else

    call MPI_Reduce(accepts, accepts, 5, MPI_NEW_INT,&
                       MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(attempts, attempts, 5, MPI_NEW_INT,&
                       MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(ener_tot_sum, ener_tot_sum, 1, MPI_NEW_REAL,&
                       MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(ener_rot_sum, ener_rot_sum, 1, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(ener_irrot_sum, ener_irrot_sum, 1, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(ener_tot_sq_sum, ener_tot_sq_sum, 1, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(ener_rot_sq_sum, ener_rot_sq_sum, 1, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(ener_irrot_sq_sum, ener_irrot_sq_sum, 1, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(ebar_sum, ebar_sum, 2, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(ebar_sq_sum, ebar_sq_sum, 2, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(ebar_dip_sum, ebar_dip_sum, 2, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(ebar_dip_sq_sum, ebar_dip_sq_sum, 2, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(ebar_wind_sum, ebar_wind_sum, 2, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(ebar_wind_sq_sum, ebar_wind_sq_sum, 2, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
    call MPI_Reduce(rho_avg, rho_avg, 2, MPI_NEW_REAL,&
                    MPI_SUM, 0, MPI_COMM_WORLD, mpierr)

  end if

  ! write (6,*) "POST: rank", id, "e_tot_sum", ener_tot_sum

  ! write (6,*) "PRE: rank ",id, "sum(s_ab)",sum(s_ab)

  if (do_corr) then

    if (id.eq.0) then

      call MPI_Reduce(MPI_IN_PLACE, s_ab, size(s_ab), MPI_NEW_REAL,&
                         MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(MPI_IN_PLACE, s_ab_rot, size(s_ab_rot),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(MPI_IN_PLACE, s_ab_irrot, size(s_ab_irrot),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(MPI_IN_PLACE, ch_ch, size(ch_ch),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(MPI_IN_PLACE, rho_k_p, size(rho_k_p),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(MPI_IN_PLACE, rho_k_m, size(rho_k_m),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(MPI_IN_PLACE, dir_struc, size(dir_struc),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(MPI_IN_PLACE, dist_r, size(dist_r),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(MPI_IN_PLACE, v_avg, size(v_avg), MPI_NEW_REAL,&
                         MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(MPI_IN_PLACE, e_tot_avg, size(e_tot_avg),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(MPI_IN_PLACE, e_rot_avg, size(e_rot_avg),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(MPI_IN_PLACE, e_irrot_avg, size(e_irrot_avg),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(MPI_IN_PLACE, bin_count, size(bin_count),&
                         MPI_NEW_INT, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)

    else

      call MPI_Reduce(s_ab, s_ab, size(s_ab), MPI_NEW_REAL,&
                         MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(s_ab_rot, s_ab_rot, size(s_ab_rot),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(s_ab_irrot, s_ab_irrot, size(s_ab_irrot),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(ch_ch, ch_ch, size(ch_ch),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(rho_k_p, rho_k_p, size(rho_k_p),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(rho_k_m, rho_k_m, size(rho_k_m),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(dir_struc, dir_struc, size(dir_struc),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(dist_r, dist_r, size(dist_r),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(v_avg, v_avg, size(v_avg), MPI_NEW_REAL,&
                         MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(e_tot_avg, e_tot_avg, size(e_tot_avg),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(e_rot_avg, e_rot_avg, size(e_rot_avg),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(e_irrot_avg, e_irrot_avg, size(e_irrot_avg),&
                         MPI_NEW_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
      call MPI_Reduce(bin_count, bin_count, size(bin_count),&
                         MPI_NEW_INT, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)

    end if

  end if

  ! write (6,*) "POST: rank ",id,"sum(s_ab)",sum(s_ab)

end subroutine reductions

subroutine normalisations(num_procs)
  use common
  implicit none
  integer, intent(in) :: num_procs
  real(kind=rk) :: denom, sp_he_tot, sp_he_rot, sp_he_irrot,&
  ebar_sus, ebar_dip_sus, ebar_wind_sus
  sp_he_tot = 0.0; sp_he_rot = 0.0; sp_he_irrot = 0.0;
  ebar_sus = 0.0; ebar_dip_sus = 0.0; ebar_wind_sus = 0.0;

  ! between the end of the k loop and MPI_Finalize, we need to
  ! average over measurements/samples, then MPI_Reduce() and
  ! average over procs as well to get a final measurement

  denom = dble(no_samples * num_procs * no_measurements)

  ener_tot_sum = ener_tot_sum / denom
  ener_rot_sum = ener_rot_sum / denom
  ener_irrot_sum = ener_irrot_sum / denom
  ener_tot_sq_sum = ener_tot_sq_sum / denom
  ener_rot_sq_sum = ener_rot_sq_sum / denom
  ener_irrot_sq_sum = ener_irrot_sq_sum / denom
  ebar_sum = ebar_sum / denom
  ebar_sq_sum = ebar_sq_sum / denom
  ebar_dip_sum = ebar_dip_sum / denom
  ebar_dip_sq_sum = ebar_dip_sq_sum / denom
  ebar_wind_sum = ebar_wind_sum / denom
  ebar_wind_sq_sum = ebar_wind_sq_sum / denom
  rho_k_p = rho_k_p / denom
  rho_k_m = rho_k_m / denom
  ch_ch = ch_ch / denom
  s_ab = s_ab / denom
  s_ab_irrot = s_ab_irrot / denom
  s_ab_rot = s_ab_rot / denom
  dir_struc = dir_struc / denom
  dist_r = dist_r / denom
  bin_count = bin_count / denom
  e_tot_avg = e_tot_avg / denom
  e_rot_avg = e_rot_avg / denom
  e_irrot_avg = e_irrot_avg / denom
  v_avg = v_avg / denom
  rho_avg = rho_avg / denom

  sp_he_tot = L**2 * beta**2 * (ener_tot_sq_sum - (ener_tot_sum)**2)
  sp_he_rot = L**2 * beta**2 * (ener_rot_sq_sum - (ener_rot_sum)**2)
  sp_he_irrot = L**2 * beta**2 * (ener_irrot_sq_sum - (ener_irrot_sum)**2)

  ebar_sus = L**2 * beta * (sum(ebar_sq_sum) - sum(ebar_sum * ebar_sum))
  ebar_dip_sus = L**2 * beta *&
  (sum(ebar_dip_sq_sum) - sum(ebar_dip_sum * ebar_dip_sum))
  ebar_wind_sus = L**2 * beta *&
  (sum(ebar_wind_sq_sum) - sum(ebar_wind_sum * ebar_wind_sum))

  write (*,*) "END. MPI_Reduce done, averages calculated."
  write (*,*) "--- Energies: ---"
  write (*,*) "Total: ",ener_tot_sum, ener_tot_sq_sum
  write (*,*) "Rotational: ",ener_rot_sum, ener_rot_sq_sum
  write (*,*) "Irrotational: ",ener_irrot_sum, ener_irrot_sq_sum
  write (*,*) "--- Specific heats: ---"
  write (*,*) "Total: ",sp_he_tot
  write (*,*) "Rotational: ",sp_he_rot
  write (*,*) "Irrotational: ",sp_he_irrot
  write (*,*)
  write (*,*) "Ebar_sum: ",ebar_sum(1),ebar_sum(2)
  write (*,*) "Ebar_sq_sum: ",ebar_sq_sum(1),ebar_sq_sum(2)
  write (*,*) "Ebar susceptibility: ",ebar_sus
  write (*,*) "Ebar_dip_sum: ",ebar_dip_sum(1),ebar_dip_sum(2)
  write (*,*) "Ebar_dip_sq_sum: ",ebar_dip_sq_sum(1),ebar_dip_sq_sum(2)
  write (*,*) "Ebar_dip susceptibility: ",ebar_dip_sus
  write (*,*) "Ebar_wind_sum: ",ebar_wind_sum(1),ebar_wind_sum(2)
  write (*,*) "Ebar_wind_sq_sum: ",ebar_wind_sq_sum(1),ebar_wind_sq_sum(2)
  write (*,*) "Ebar_wind susceptibility: ",ebar_wind_sus
  write (*,*) "Avg. charge density",rho_avg

  open  (30, file=sphe_sus_file // '_MPI_allreduced')
  write (30,'(a)') "# Temp., sp_he^total, sp_he^rot., sp_he^irrot"
  write (30,'(4ES18.9)') temp, sp_he_tot, sp_he_rot, sp_he_irrot
  write (30,'(a)') "# Temp., Chi_Ebar, Chi_{Ebar_dip}, Chi_{Ebar_wind}"
  write (30,'(2ES18.9)') temp, ebar_sus, ebar_dip_sus, ebar_wind_sus
  write (30,'(a)') "# Avg. charge density"
  write (30,'(ES18.9)') rho_avg
  close (30)

end subroutine normalisations
