! maggs-rossetto algorithm on cubic lattice
! i'll write proper docs at some point
! callum gray, UCL / ENS Lyon, 2016. do wat u like
program mr

  use common
  use linear_solver
  use io
  use setup
  implicit none
  include 'revision.inc'
  integer :: i, j, k, tot_q
  real*8 :: start_time, end_time
  real*8, dimension(6) :: timings
  real(kind=16) :: sp_he_tot, sp_he_rot, sp_he_irrot, ebar_sus
  integer, dimension(8) :: values

  write (*,*) "Maggs-Rossetto algorithm on the simple cubic lattice."
  write (*,*) "Callum Gray, University College London, 2017."
  write (*,*) "Git revision ",revision
  write (*,*) "Repo at http://github.com/cuamll/mr"

  e_tot_sum = 0.0; e_rot_sum = 0.0; e_irrot_sum = 0.0;
  e_tot_sq_sum = 0.0; e_rot_sq_sum = 0.0; e_irrot_sq_sum = 0.0;
  ebar_sum = 0.0; ebar_sq_sum = 0.0;
  sp_he_tot = 0.0; sp_he_rot = 0.0; sp_he_irrot = 0.0; ebar_sus = 0.0;

  call cpu_time(start_time)

  call date_and_time(VALUES=values)

  write (*,'(a24,I2.1,a1,I2.1,a1,I4.2,a1,I2.1,a1,I2.1,a1,I2.1)')&
    &" Current date and time: ",values(3),"/",values(2),"/"&
    &,values(1)," ",values(5),":",values(6),":",values(7)

  call read_input

  do k = 1, no_samples

    write (*,'(a, i3.1)') "Sample no. ", k

    if (no_samples.eq.1) then
      call setup_wrapper(seed)
    else
      call setup_wrapper(k - 1)
    end if

    accepth = 0
    acceptr = 0
    acceptg = 0

    write (*,*)
    write (*,'(A10,I5.1,A27)',advance='no') "Beginning ",therm_sweeps," thermalisation sweeps... ["

    call cpu_time(timings(1))

    do i = 1,therm_sweeps
      call mc_sweep
      if (mod(i,therm_sweeps/10).eq.0) then
        write (*,'(a)',advance='no') "*"
      end if
    end do

    call cpu_time(timings(2))

    write (*,'(a)') "]. Completed."
    write (*,'(a,f8.3,a)') "Time taken: ",timings(2)-timings(1)," seconds."
    write (*,*)
    write (*,'(a,I7.1,a,I4.1,a)',advance='no') "Beginning ",measurement_sweeps,&
      &" measurement sweeps with sampling every ",&
      &sample_interval," MC steps... ["

    call cpu_time(timings(3))

    do i = 1,measurement_sweeps
      call mc_sweep

      if (mod(i, sample_interval).eq.0) then
        call measure(i)
      end if

      if (measurement_sweeps.gt.100) then
        if (mod(i,measurement_sweeps / 100).eq.0) then
          write (*,'(a)',advance='no') "*"
        end if
      end if

    end do

    write (*,'(a)') "]"

    call cpu_time(timings(4))

    write (*,*)
    write (*,'(I7.1,a)') measurement_sweeps," measurement sweeps completed."
    write (*,'(a,f8.3,a)') "Time taken: ",timings(4)-timings(3)," seconds."

    if (sum(v,mask=v.gt.0).ne.0) then
      call linsol
    end if

    !write(*,*)
    !write(*,*) "--- END: CHARGE POSITIONS ---"

    tot_q = 0
      do j = 1,L
        do i = 1,L

          tot_q = tot_q + abs(v(i,j))

          !if (v(i,j).eq.1) then
          !  write (*,'(I3.1,I3.1,I3.1,I3.1)') i,j,v(i,j)
          !end if

          !if (v(i,j).eq.-1) then
          !  write (*,'(I3.1,I3.1,I3.1,I3.1)') i,j,v(i,j)
          !end if

        end do
      end do

    !write (*,*) "E_bar:",ebar(1),ebar(2)
    !write(*,*)

    !write (*,*) "Total charges: ",tot_q
    !write(*,*)

    write (*,*)
    write (*,*) '--- MOVE STATS: ---'
    if (tot_q.ne.0) then
      write (*,*) 'hops: moves, acceptance ratio: ',&
        &accepth,float(accepth) / ((therm_sweeps + measurement_sweeps)&
        &* tot_q * hop_ratio)
    end if
    write (*,*) 'rotational updates: moves, acceptance ratio: ',&
      &acceptr,float(acceptr) / ((therm_sweeps + measurement_sweeps)&
      &* L**2 * rot_ratio)
    write (*,*) 'harmonic updates: moves, acceptance ratio: ',&
      &acceptg,float(acceptg) / ((therm_sweeps + measurement_sweeps)&
      &* L**2 * g_ratio * 2)
    write(*,*)

    call cpu_time(timings(5))

    if (do_corr) then
      write (*,*) "Writing output for ",no_measurements," measurements..."
      call write_output
      write (*,*) "Done."
    end if

    call deallocations

    write (*,*)

    call cpu_time(timings(6))

    call cpu_time(end_time)

    write (*,*) "--- Current energies: ---"
    write (*,*) "Total: ",e_tot_sum
    write (*,*) "Rotational: ",e_rot_sum
    write (*,*) "Irrotational: ",e_irrot_sum
    write (*,*) "Ebar sum: ",ebar_sum(1), ebar_sum(2)
    write (*,*) "Ebar^2 sum: ",ebar_sq_sum(1), ebar_sq_sum(2)
    write (*,*) "Current susc: ",((L**2 * beta) * (sum(ebar_sq_sum) - sum(ebar_sum * ebar_sum))) / (no_measurements * k)
    write (*,*)

  end do

  write (*,'(a,f8.3,a)') "Total time taken: ",end_time-start_time," seconds."
  write (*,'(a,f8.3,a)') "Total time per sample: ",(end_time-start_time)/no_samples," seconds."
  write (*,'(a,f8.5,a)') "Time per thermalisation sweep: "&
    &,(timings(2) - timings(1)) / (therm_sweeps)," seconds."
  write (*,'(a,f8.5,a)') "Time per measurement sweep (including&
    & measurements): ",(timings(4) - timings(3)) / (measurement_sweeps)," seconds."

  if (do_corr) then
    write (*,'(a,f8.5,a)') "Time to read in data and calculate correlations&
    &(per step): ",(timings(6) - timings(5)) / (measurement_sweeps)," seconds."
  end if

  write (*,*)


  e_tot_sum = e_tot_sum / (no_samples * no_measurements)
  e_rot_sum = e_rot_sum / (no_samples * no_measurements)
  e_irrot_sum = e_irrot_sum / (no_samples * no_measurements)
  e_tot_sq_sum = e_tot_sq_sum / (no_samples * no_measurements)
  e_rot_sq_sum = e_rot_sq_sum / (no_samples * no_measurements)
  e_irrot_sq_sum = e_irrot_sq_sum / (no_samples * no_measurements)
  ebar_sum = ebar_sum / (no_samples * no_measurements)
  ebar_sq_sum = ebar_sq_sum / (no_samples * no_measurements)

  sp_he_tot = L**2 * beta**2 * (e_tot_sq_sum - (e_tot_sum)**2)
  sp_he_rot = L**2 * beta**2 * (e_rot_sq_sum - (e_rot_sum)**2)
  sp_he_irrot = L**2 * beta**2 * (e_irrot_sq_sum - (e_irrot_sum)**2)

  ebar_sus = L**2 * beta * (ebar_sq_sum(1) + ebar_sq_sum(2) - (ebar_sum(1)**2 + ebar_sum(2)**2))

  write (*,*) "--- Specific heats: ---"
  write (*,*) "Total: ",sp_he_tot
  write (*,*) "Rotational: ",sp_he_rot
  write (*,*) "Irrotational: ",sp_he_irrot
  write (*,*)
  write (*,*) "Ebar susceptibility: ",ebar_sus

  open(30, file=sphe_sus_file)
  write (30,'(a)') "# Temp., sp_he^total, sp_he^rot., sp_he^irrot"
  write (30,'(ES18.9,ES18.9,ES18.9,ES18.9)') temp, sp_he_tot, sp_he_rot, sp_he_irrot
  write (30,'(a)') "# Temp., <Ebar^2> - <Ebar>^2"
  write (30,'(ES18.9, ES18.9)') temp, ebar_sus

  stop

end program mr

subroutine mc_sweep
  use common
  implicit none
  integer :: i,j,k,m,charge,glob,totq,mu1,mu2,pm1
  integer, dimension(2) :: site
  real*8 :: eo1,eo2,eo3,eo4,en1,en2,en3,en4
  real*8 :: u_tot,u_tot_run,increment
  real*8 :: old_e, new_e, delta_e, g_thr

  charge = 0; totq = 0; mu1 = 0; mu2 = 0; pm1 = 0
  site = (/ 0, 0 /)
  eo1 = 0.0; eo2 = 0.0; eo3 = 0.0; eo4 = 0.0
  en1 = 0.0; en2 = 0.0; en3 = 0.0; en4 = 0.0
  u_tot = 0.0; u_tot_run = 0.0; increment = 0.0
  old_e = 0.0; new_e = 0.0; delta_e = 0.0; g_thr = 1 / float(L)

  ! --- START OF UPDATE BLOCKS ---

  ! --- CHARGE HOP UPDATE ---

  do i = 1, int(L**2 * hop_ratio)

    ! NOTE TO SELF: this whole procedure assumes
    ! single-valued charges only.

    ! pick a random site
    site = (/ int(rand() * L) + 1,&
              &int(rand() * L) + 1 /)

    mu1 = floor(2*rand())+1

    ! check we're not doing anything weird with mu1ltiple charges
    if (abs(v(site(1),site(2))).gt.1) then
      write (*,*) "Charge at ",site(1),site(2),&
        " = ",v(site(1),site(2)),". Exiting."
      write (*,*)
      stop
    end if

    ! pick a non-zero charge - canonical!
    if (v(site(1),site(2)).ne.0) then

      charge = v(site(1),site(2))

      ! this takes care of sign issues when hopping
      increment = q * charge / (eps_0 * lambda)

      if (rand().lt.0.5) then ! move it "negative"

        ! get negative in the mu1 direction
        site(mu1) = neg(site(mu1))
        eo1 = e_field(mu1,site(1),site(2))
        en1 = eo1 + increment
        old_e = 0.5 * eps_0 * eo1**2
        new_e = 0.5 * eps_0 * en1**2
        delta_e = new_e - old_e

        if (v(site(1),site(2)).eq.0) then
          if ((delta_e.lt.0.0).or.(exp((-beta)*delta_e).gt.rand())) then

            v(site(1),site(2)) = charge
            e_field(mu1,site(1),site(2)) = en1
            ebar(mu1) = ebar(mu1) + (increment / L**2)

            ! go back to the original site and set the charge to 0
            site(mu1) = pos(site(mu1))
            v(site(1),site(2)) = 0

            accepth = accepth + 1
            u_tot_run = u_tot_run + delta_e

            end if
          end if

      else ! move it "positive"

        eo1 = e_field(mu1,site(1),site(2))
        en1 = eo1 - increment
        old_e = 0.5 * eps_0 * lambda**2 * eo1**2
        new_e = 0.5 * eps_0 * lambda**2 * en1**2
        delta_e = new_e - old_e

        ! get pos in the mu1 direction
        site(mu1) = pos(site(mu1))

        if (v(site(1),site(2)).eq.0) then
          if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

            v(site(1),site(2)) = charge

            ! go back to the original site and set the charge to 0
            site(mu1) = neg(site(mu1))
            v(site(1),site(2)) = 0
            e_field(mu1,site(1),site(2)) = en1
            ebar(mu1) = ebar(mu1) - (increment / L**2)

            accepth = accepth + 1
            u_tot_run = u_tot_run + delta_e

          end if
        end if

      end if ! end "positive" / "negative" choice
    end if ! end charge.ne.0 block

  end do ! end charge hop sweep

  !glob = 1
  if (ebar(mu1).gt.(g_thr).or.ebar(mu1).lt.((-1)*g_thr)) then
    glob = 1
  else
    glob = 0
  end if

  mu1 = 0; increment = 0.0;
  u_tot = 0.0
  u_tot = 0.5 * eps_0 * lambda**2 * sum(e_field * e_field)

  ! --- ROTATIONAL UPDATE ---

  do i = 1,int(L**2 * rot_ratio)

    eo1 = 0.0; eo2 = 0.0; eo3 = 0.0; eo4 = 0.0
    en1 = 0.0; en2 = 0.0; en3 = 0.0; en4 = 0.0
    site = (/ 0, 0 /)

    ! pick at random from interval [-Delta_max, +Delta_max]
    increment = 2 * rot_delt * (rand() - 0.5)

    ! pick xy (1,2), yz (2,3), or zx (3,1) plaquette
    mu1 = floor(2*rand())+1
    mu2 = mod(mu1,2) + 1

    ! give us a coordinate (i,j,k)
    site = (/ int(rand() * L) + 1,&
             &int(rand() * L) + 1 /)
    ! site(mu1) is the coordinate in the first direction
    ! site(mu2) is the coordinate in the second

    ! my convention is that an xy (1,2) plaquette is defined by
    ! field links (1,i,j,k),(2,i,j,k),(1,neg(i),j,k),(2,i,neg(j),k)
    ! so the first two are easy:
    eo1 = e_field(mu1,site(1),site(2))
    eo2 = e_field(mu2,site(1),site(2))

    ! if e.g. we picked xy, the next line does x -> neg(x)
    site(mu1) = neg(site(mu1))
    eo4 = e_field(mu1,site(1),site(2))
    ! and now we need to put it back
    site(mu1) = pos(site(mu1))

    ! this does y -> neg(y)
    site(mu2) = neg(site(mu2))
    eo3 = e_field(mu2,site(1),site(2))
    ! and put it back
    site(mu2) = pos(site(mu2))
    ! so now we have (x,y,z) again

    en1 = eo1 + increment
    en2 = eo2 - increment
    en3 = eo3 + increment
    en4 = eo4 - increment

    old_e = 0.5 * eps_0 * lambda**2 * (eo1**2 + eo2**2 + eo3**2 + eo4**2)
    new_e = 0.5 * eps_0 * lambda**2 * (en1**2 + en2**2 + en3**2 + en4**2)
    delta_e = new_e - old_e

    if ((delta_e.lt.0.0).or.(exp((-beta)*delta_e).gt.rand())) then

      e_field(mu1,site(1),site(2)) = en1
      e_field(mu2,site(1),site(2)) = en2

      site(mu1) = neg(site(mu1))
      e_field(mu1,site(1),site(2)) = en4
      site(mu1) = pos(site(mu1))

      site(mu2) = neg(site(mu2))
      e_field(mu2,site(1),site(2)) = en3
      site(mu2) = pos(site(mu2))

      acceptr = acceptr + 1
      u_tot_run = u_tot_run + delta_e

    end if ! end of Metropolis check

  end do ! end rotational

  increment = 0.0
  u_tot = 0.0
  u_tot = 0.5 * eps_0 * sum(e_field * e_field)

  ! --- HARMONIC UPDATE ---
  ! e bar update
  if (glob.eq.1) then

    ! NOTE TO SELF: again this int cast might need changing
    do i = 1,int(L**2 * g_ratio)

      increment = (q) / (L * eps_0 * lambda)

      do mu1 = 1,2

        if (rand() - 0.5.le.0) then
          pm1 = -1
        else
          pm1 = +1
        end if

        ! NOTE TO SELF - check this little fucker
        !old_e = (u_tot_run + ebar_x)**2
        !new_e = (u_tot_run + ebar_x - g_thr * float(L))**2
        old_e = ebar(mu1)**2
        ! this is not right
        new_e = (ebar(mu1) + pm1 * increment)**2
        delta_e = new_e - old_e
        !delta_e = (0.5 - float(L) * ebar_x)

        if ((delta_e.lt.0.0).or.((exp(-beta*delta_e).gt.rand())&
          .and.(exp(-beta*delta_e).gt.0.00000000001))) then
          ! this block is basically stolen from Michael
          ! not sure what's happening here tbh
          !write (*,*) ebar(1),ebar(2),ebar(3)
          ebar(mu1) = ebar(mu1) + pm1 * increment
          e_field(mu1,:,:) = e_field(mu1,:,:) + pm1 * (increment)

          !write (*,*) ebar(1),ebar(2),ebar(3)
          !do m = 1,L
          !  do k = 1,L
          !    do j = 1,L
          !      e_field(mu1,j,k,m) = e_field(mu1,j,k,m) + pm1 * (increment / L**3)
          !    end do
          !  end do
          !end do

          acceptg = acceptg + 1

        end if ! end weird Metropolis block

      end do ! end mu1 loop

    end do ! end i loop

  end if ! glob.eq.1

  u_tot = 0.0
  u_tot = 0.5 * eps_0 * sum(e_field * e_field)

  ! --- END OF UPDATE BLOCKS --- !

end subroutine mc_sweep

! | --------------- SUBROUTINE MEASURE(STEP_NUMBER) --------------- |
! |                 "MEASURES" RELEVANT QUANTITIES                  |
! |                                                                 |
! | Basically just take an integer and print out the fields.        |
! | Also runs the linear solver and prints out the irrotational     |
! | field, along with the charge distribution. We read this back in |
! | at the end, and do analysis then.                               |
! |                                                                 |
! | --------------------------------------------------------------- |

subroutine measure(step_number)
  use common
  use linear_solver
  implicit none
  logical :: field_ch_exist
  integer,intent(in) :: step_number
  integer :: n
  real*8 :: u_tot_run
  real(kind=8) :: ener_tot, ener_rot, ener_irrot
  real(kind=8), dimension(2,L,L) :: e_rot

  ener_tot = 0.0; ener_rot = 0.0; ener_irrot = 0.0;

  ! this is basically just used for the energy at this point
  ! could move this to the analysis part pretty trivially
  !do_corr = .false.

  !u_tot_run = 0.0
  !u_tot_run = 0.5 * eps_0 * sum(e_field**2)

  !energy(n + 1) = u_tot_run
  !sq_energy(n + 1) = u_tot_run**2

  ! get irrotational part of field
  mnphi = 0.0; e_rot = 0.0
  call linsol

  if (do_corr) then
    n = step_number / sample_interval
    ! if this is the first step, start a new file;
    ! otherwise append to the file that's there
    inquire(file=field_charge_file, exist=field_ch_exist)
    if (n.eq.1) then
      open(15, file=field_charge_file, status="replace", action="write",access="stream",form="unformatted")
    else if (field_ch_exist) then
      open(15, file=field_charge_file, status="old", position="append", action="write", access="stream", form="unformatted")
    else
      write (*,*) "Fields/charges file does not exist? Current step is ",n
      open(15, file=field_charge_file, status="new", action="write",access="stream",form="unformatted")
    end if

    write(15) e_field
    write(15) mnphi
    write(15) v
    !write(15) ebar

    close(15)

  else

    e_rot = e_field - mnphi

    ener_tot = 0.5 * eps_0 * sum(e_field * e_field) / L**2
    ener_rot = 0.5 * eps_0 * sum(e_rot * e_rot) / L**2
    ener_irrot = 0.5 * eps_0 * sum(mnphi * mnphi) / L**2

    e_tot_sum = e_tot_sum + ener_tot
    e_rot_sum = e_rot_sum + ener_rot
    e_irrot_sum = e_irrot_sum + ener_irrot
    e_tot_sq_sum = e_tot_sq_sum + ener_tot**2
    e_rot_sq_sum = e_rot_sq_sum + ener_rot**2
    e_irrot_sq_sum = e_irrot_sq_sum + ener_irrot**2

    ebar(1) = lambda**2 * sum(e_field(1,:,:))
    ebar(2) = lambda**2 * sum(e_field(2,:,:))
    ebar = ebar / L**2
    ebar_sum(1) = ebar_sum(1) + ebar(1)
    ebar_sum(2) = ebar_sum(2) + ebar(2)
    ebar_sq_sum(1) = ebar_sq_sum(1) + (ebar(1)**2)
    ebar_sq_sum(2) = ebar_sq_sum(2) + (ebar(2)**2)

  end if

end subroutine measure
