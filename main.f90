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
  integer :: i, j, k, n, tot_q
  real*8 :: start_time, end_time
  real*8, dimension(4) :: timings
  character(8) :: date
  character(10) :: time
  character(5) :: zone
  integer, dimension(8) :: values

  call cpu_time(start_time)

  write (*,*) "Maggs-Rossetto algorithm on the simple cubic lattice."
  write (*,*) "Callum Gray, University College London, 2017."
  write (*,*) "Git revision ",revision
  write (*,*) "Repo at http://github.com/callumgray/mr"

  call date_and_time(date,time,zone,values)

  write (*,'(a24,I2.1,a1,I2.1,a1,I4.2,a1,I2.1,a1,I2.1,a1,I2.1)')&
    &" Current date and time: ",values(3),"/",values(2),"/"&
    &,values(1)," ",values(5),":",values(6),":",values(7)

  call setup_wrapper

  accepth = 0
  acceptr = 0
  acceptg = 0

  write (*,'(A10,I5.1,A27)',advance='no') "Beginning ",therm_sweeps," thermalisation sweeps... ["

  call cpu_time(timings(1))

  do i = 1,therm_sweeps
    call mc_sweep
    if (mod(i,therm_sweeps/10).eq.0) then
      write (*,'(a)',advance='no') "*"
    end if
  end do

  call cpu_time(timings(2))

  write (*,'(a)') "]"
  write (*,'(I5.1,a)') therm_sweeps," thermalisation sweeps completed."
  write (*,'(a,f8.3,a)') "Time taken: ",timings(2)-timings(1)," seconds."
  write (*,*)
  write (*,'(a,I7.1,a,I4.1,a)') "Beginning ",measurement_sweeps,&
    &" measurement sweeps with sampling every ",&
    &sample_interval," MC steps..."

  call cpu_time(timings(3))

  do i = 1,measurement_sweeps
    call mc_sweep

    if (mod(i, sample_interval).eq.0) then
      call measure(i)
    end if

    if (mod(i, 1000).eq.0) then
      write (*,'(I7.4,a)') i," steps completed..."
    end if

  end do

  call cpu_time(timings(4))

  write (*,*)
  write (*,'(I7.1,a)') measurement_sweeps," measurement sweeps completed."
  write (*,'(a,f8.3,a)') "Time taken: ",timings(4)-timings(3)," seconds."

  if (sum(v,mask=v.gt.0).ne.0) then
    call linsol
  end if

  write(*,*)
  write(*,*) "--- END: CHARGE POSITIONS ---"

  tot_q = 0
  do i = 1,L
    do j = 1,L
      do k = 1,L

        tot_q = tot_q + abs(v(i,j,k))

        if (v(i,j,k).eq.1) then
          write (*,'(I3.1,I3.1,I3.1,I3.1)') i,j,k,v(i,j,k)
        end if

        if (v(i,j,k).eq.-1) then
          write (*,'(I3.1,I3.1,I3.1,I3.1)') i,j,k,v(i,j,k)
        end if

      end do
    end do
  end do

  write (*,*) "Total charges: ",tot_q
  write(*,*)

  write (*,*)
  write (*,*) '--- MOVE STATS: ---'
  write (*,*) 'hops: moves, acceptance ratio: ',&
    &accepth,float(accepth) / ((therm_sweeps + measurement_sweeps)&
    &* tot_q * hop_ratio)
  write (*,*) 'rotational updates: moves, acceptance ratio: ',&
    &acceptr,float(acceptr) / ((therm_sweeps + measurement_sweeps)&
    &* 3*L**3 * rot_ratio)
  write (*,*) 'harmonic updates: moves, acceptance ratio: ',&
    &acceptg,float(acceptg) / ((therm_sweeps + measurement_sweeps)&
    &* L**3 * g_ratio * 3)
  write(*,*)


  write (*,*) "Writing output for ",no_measurements," measurements..."
  call write_output

  call deallocations

  write (*,*) "Done."
  write (*,*)

  call cpu_time(end_time)
  write (*,'(a,f8.3,a)') "Total time taken: ",end_time-start_time," seconds."
  write (*,'(a,f8.5,a)') "Time per thermalisation sweep: "&
    &,(timings(2) - timings(1)) / therm_sweeps," seconds."
  write (*,'(a,f8.5,a)') "Time per measurement sweep (including&
    & measurements): ",(timings(4) - timings(3)) / measurement_sweeps," seconds."
  write (*,*)

  stop

end program mr

subroutine mc_sweep
  use common
  implicit none
  integer :: i,x,y,z,charge,glob,totq
  real*8 :: eo1,eo2,eo3,eo4,en1,en2,en3,en4
  real*8 :: u_tot,u_tot_run,ebar_inc,e_inc
  real*8 :: hop_inc, old_e, new_e, delta_e, g_thr
  real :: chooser, delta

  charge = 0
  glob = 0
  totq = 0
  eo1 = 0.0
  eo2 = 0.0
  eo3 = 0.0
  eo4 = 0.0
  en1 = 0.0
  en2 = 0.0
  en3 = 0.0
  en4 = 0.0
  u_tot = 0.0
  u_tot_run = 0.0
  ebar_inc = 0.0
  e_inc = 0.0
  hop_inc = 0.0
  old_e = 0.0
  new_e = 0.0
  delta_e = 0.0
  g_thr = 1 / float(L)
  ! --- START OF UPDATE BLOCKS ---

  ! --- CHARGE HOP UPDATE ---

  do i = 1, int(L**3 * hop_ratio)

    ! pick a random site
    x = int(rand() * L) + 1
    y = int(rand() * L) + 1
    z = int(rand() * L) + 1

    ! pick a non-zero charge - canonical!
    if (v(x,y,z).ne.0) then

      charge = v(x,y,z)
      ! this takes care of sign issues when hopping
      hop_inc = charge / (eps_0 * lambda**2)
      chooser = rand()

      if (chooser.lt.(1.0 / 3.0)) then
        ! x component
        eo1 = e_x(x,y,z)
        chooser = rand()

        if (chooser.lt.0.5) then ! we try and move the fucker left
          en1 = eo1 + hop_inc
          old_e = 0.5 * eps_0 * eo1**2
          new_e = 0.5 * eps_0 * en1**2
          delta_e = new_e - old_e
          ! i think we can just set v(x,y,z) = 0
          ! if we're enforcing |v| <= 1
          if ((abs(v(x,y,z)).le.1).and.(v(neg(x),y,z).eq.0)) then
            if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

              accepth = accepth + 1
              e_x(x,y,z) = en1
              v(neg(x),y,z) = charge
              v(x,y,z) = 0
              ebar_x = ebar_x - hop_inc
              u_tot_run = u_tot_run + delta_e

            end if
          end if

        else ! move the fucker to the right

          en1 = eo1 - hop_inc
          old_e = 0.5 * eps_0 * eo1**2
          new_e = 0.5 * eps_0 * en1**2
          delta_e = new_e - old_e

          if ((abs(v(x,y,z)).le.1).and.(v(pos(x),y,z).eq.0)) then
            if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

              accepth = accepth + 1
              e_x(x,y,z) = en1
              v(pos(x),y,z) = charge
              v(x,y,z) = 0
              ebar_x = ebar_x + hop_inc
              u_tot_run = u_tot_run + delta_e

            end if
          end if

        end if ! end of left-right movement choice
        if (ebar_x.gt.(1/float(L)).or.ebar_x.lt.((-1)/float(L))) then
          glob = 1
        else
          glob = 0
        end if


      else if (chooser.gt.(2.0 / 3.0)) then ! y component

        eo1 = e_y(x,y,z)
        chooser = rand()

        if (chooser.lt.0.5) then ! we try and move the fucker left
          en1 = eo1 + hop_inc
          old_e = 0.5 * eps_0 * eo1**2
          new_e = 0.5 * eps_0 * en1**2
          delta_e = new_e - old_e
          ! i think we can just set v(x,y,z) = 0
          ! if we're enforcing |v| <= 1
          if ((abs(v(x,y,z)).le.1).and.(v(x,neg(y),z).eq.0)) then
            if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

              accepth = accepth + 1
              e_y(x,y,z) = en1
              v(x,neg(y),z) = charge
              v(x,y,z) = 0
              ebar_y = ebar_y - hop_inc
              u_tot_run = u_tot_run + delta_e

            end if
          end if

        else ! move the fucker to the right

          en1 = eo1 - hop_inc

          old_e = 0.5 * eps_0 * eo1**2
          new_e = 0.5 * eps_0 * en1**2
          delta_e = new_e - old_e

          if ((abs(v(x,y,z)).le.1).and.(v(x,pos(y),z).eq.0)) then
            if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

              accepth = accepth + 1
              e_y(x,y,z) = en1
              v(x,pos(y),z) = charge
              v(x,y,z) = 0
              ebar_y = ebar_y + hop_inc
              u_tot_run = u_tot_run + delta_e

            end if
          end if

        end if ! end of left-right movement choice
        if (ebar_y.gt.(1/float(L)).or.ebar_y.lt.((-1)/float(L))) then
          glob = 1
        else
          glob = 0
        end if

      else ! z component

        eo1 = e_z(x,y,z)
        chooser = rand()

        if (chooser.lt.0.5) then ! we try and move the fucker left
          en1 = eo1 + hop_inc
          old_e = 0.5 * eps_0 * eo1**2
          new_e = 0.5 * eps_0 * en1**2
          delta_e = new_e - old_e
          ! i think we can just set v(x,y,z) = 0
          ! if we're enforcing |v| <= 1
          if ((abs(v(x,y,z)).le.1).and.(v(x,y,neg(z)).eq.0)) then
            if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

              accepth = accepth + 1
              e_z(x,y,z) = en1
              v(x,y,neg(z)) = charge
              v(x,y,z) = 0
              ebar_z = ebar_z - hop_inc
              u_tot_run = u_tot_run + delta_e

            end if
          end if

        else ! move the fucker to the right

          en1 = eo1 - hop_inc
          old_e = 0.5 * eps_0 * eo1**2
          new_e = 0.5 * eps_0 * en1**2
          delta_e = new_e - old_e

          if ((abs(v(x,y,z)).le.1).and.(v(x,y,pos(z)).eq.0)) then
            if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

              accepth = accepth + 1
              e_z(x,y,z) = en1
              v(x,y,pos(z)) = charge
              v(x,y,z) = 0
              ebar_z = ebar_z + hop_inc
              u_tot_run = u_tot_run + delta_e

            end if
          end if

        end if ! end of left-right movement choice
        if (ebar_z.gt.(1/float(L)).or.ebar_z.lt.((-1)/float(L))) then
          glob = 1
        else
          glob = 0
        end if

      end if ! end of component chooser

    end if ! end of v(x,y,z).ne.0 block

  end do ! end charge hop sweep

  u_tot = 0.0
  do x = 1,L
    do y = 1, L
      do z = 1,L
      u_tot = u_tot + 0.5 * eps_0 * (e_x(x,y,z)**2 + e_y(x,y,z)**2 + e_z(x,y,z)**2)
      end do
    end do
  end do

  ! --- ROTATIONAL UPDATE ---

  do i = 1,int(L**3 * rot_ratio)

    eo1 = 0.0
    eo2 = 0.0
    eo3 = 0.0
    eo4 = 0.0
    en1 = 0.0
    en2 = 0.0
    en3 = 0.0
    en4 = 0.0

    x = int(rand() * L) + 1
    y = int(rand() * L) + 1
    z = int(rand() * L) + 1

    ! NOTE TO SELF : next line needs changing
    ! it's for the field increment, not an if statement
    delta = 2 * rot_delt * (rand() - 0.5)

    chooser=rand()
    if (chooser.lt.(1.0/3.0)) then ! xy-plane plaquette

      eo1 = e_x(x,y,z)
      eo2 = e_y(x,y,z)
      eo3 = e_x(x,neg(y),z)
      eo4 = e_y(neg(x),y,z)

      en1 = eo1 + delta
      en2 = eo2 - delta
      en3 = eo3 + delta
      en4 = eo4 - delta

      old_e = 0.5 * eps_0 * (eo1**2 + eo2**2 + eo3**2 + eo4**2)
      new_e = 0.5 * eps_0 * (en1**2 + en2**2 + en3**2 + en4**2)
      delta_e = new_e - old_e

      if ((delta_e.lt.0.0).or.(exp((-beta)*delta_e).gt.rand())) then

        e_x(x,y,z) = en1
        e_y(x,y,z) = en2
        e_x(x,neg(y),z) = en3
        e_y(neg(x),y,z) = en4
        acceptr = acceptr + 1
        u_tot_run = u_tot_run + delta_e

      end if ! end of Metropolis check

    else if (chooser.ge.(2.0/3.0)) then ! xz-plane plaquette

      eo1 = e_x(x,y,z)
      eo2 = e_z(x,y,z)
      eo3 = e_x(x,y,neg(z))
      eo4 = e_z(neg(x),y,z)

      en1 = eo1 + delta
      en2 = eo2 - delta
      en3 = eo3 + delta
      en4 = eo4 - delta

      old_e = 0.5 * eps_0 * (eo1**2 + eo2**2 + eo3**2 + eo4**2)
      new_e = 0.5 * eps_0 * (en1**2 + en2**2 + en3**2 + en4**2)
      delta_e = new_e - old_e

      if ((delta_e.lt.0.0).or.(exp((-beta)*delta_e).gt.rand())) then

        e_x(x,y,z) = en1
        e_z(x,y,z) = en2
        e_x(x,y,neg(z)) = en3
        e_z(neg(x),y,z) = en4
        acceptr = acceptr + 1
        u_tot_run = u_tot_run + delta_e

      end if ! end of Metropolis check

    else ! yz-plane plaquette

      eo1 = e_y(x,y,z)
      eo2 = e_z(x,y,z)
      eo3 = e_y(x,y,neg(z))
      eo4 = e_z(x,neg(y),z)

      en1 = eo1 + delta
      en2 = eo2 - delta
      en3 = eo3 + delta
      en4 = eo4 - delta

      old_e = 0.5 * eps_0 * (eo1**2 + eo2**2 + eo3**2 + eo4**2)
      new_e = 0.5 * eps_0 * (en1**2 + en2**2 + en3**2 + en4**2)
      delta_e = new_e - old_e

      if ((delta_e.lt.0.0).or.(exp((-beta)*delta_e).gt.rand())) then

        !write (*,*) "y-z plane rot:"
        !write (*,*) x, y, z
        !write (*,*) x, y, neg(z)
        !write (*,*) x, neg(y), z

        e_y(x,y,z) = en1
        e_z(x,y,z) = en2
        e_y(x,y,neg(z)) = en3
        e_z(x,neg(y),z) = en4
        acceptr = acceptr + 1
        u_tot_run = u_tot_run + delta_e

      end if ! end of Metropolis check

    end if ! end plane choice block

  end do ! end rotational
  !write (*,*) "utot after rot. = ",u_tot_run

  u_tot = 0.0
  do x = 1,L
    do y = 1, L
      do z = 1,L
      u_tot = u_tot + 0.5 * eps_0 * (e_x(x,y,z)**2 + e_y(x,y,z)**2 + e_z(x,y,z)**2)
      end do
    end do
  end do

  ! --- HARMONIC UPDATE ---
  ! e bar update
  if (glob.eq.1) then

    ! NOTE TO SELF: again this int cast prob needs changing
    do i = 1,int(L**3 * g_ratio)

      ebar_inc = q/(L**2 * eps_0)
      e_inc = ebar_inc / L**3

      ! x-component
      chooser = rand()

      if (chooser.lt.0.5) then
        ! NOTE TO SELF - check this little fucker
        !old_e = (u_tot_run + ebar_x)**2
        !new_e = (u_tot_run + ebar_x - g_thr * float(L))**2
        old_e = ebar_x**2
        ! this is not right
        new_e = (ebar_x - ebar_inc)**2
        delta_e = new_e - old_e
        !delta_e = (0.5 - float(L) * ebar_x)

        if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
          .and.(exp(-beta*delta_e).gt.0.00000000001))) then
          ! this block is basically stolen from Michael
          ! not sure what's happening here tbh
          ebar_x = ebar_x - ebar_inc
          acceptg = acceptg + 1

          do x = 1,L
            do y = 1,L
              do z = 1,L
                e_x(x,y,z) = e_x(x,y,z) - e_inc
              end do
            end do
          end do
        end if ! end weird Metropolis block

      else
        ! NOTE TO SELF - check this little fucker
        !old_e = (u_tot_run + ebar_x)**2
        !new_e = (u_tot_run + ebar_x - g_thr * float(L))**2
        old_e = ebar_x**2
        new_e = (ebar_x + ebar_inc)**2
        delta_e = new_e - old_e
        !delta_e = (0.5 - float(L) * ebar_x)

        if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
          .and.(exp(-beta*delta_e).gt.0.00000000001))) then
          ! this block is basically stolen from Michael
          ! not sure what's happening here tbh
          ebar_x = ebar_x + ebar_inc
          acceptg = acceptg + 1

          do x = 1,L
            do y = 1,L
              do z = 1,L
                e_x(x,y,z) = e_x(x,y,z) + e_inc
              end do
            end do
          end do
        end if ! end weird Metropolis block
      end if

      ! y-component
      chooser = rand()

      if (chooser.lt.0.5) then
        ! NOTE TO SELF - check this little fucker
        !old_e = (u_tot_run + ebar_y)**2
        !new_e = (u_tot_run + ebar_y - g_thr * float(L))**2
        old_e = ebar_y**2
        new_e = (ebar_y - ebar_inc)**2
        delta_e = new_e - old_e
        !delta_e = (0.5 - float(L) * ebar_y)

        if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
          .and.(exp(-beta*delta_e).gt.0.00000000001))) then
          ! this block is basically stolen from Michael
          ! not sure what's happening here tbh
          ebar_y = ebar_y - ebar_inc
          acceptg = acceptg + 1

          do x = 1,L
            do y = 1,L
              do z = 1,L
                e_y(x,y,z) = e_y(x,y,z) - e_inc
              end do
            end do
          end do
        end if ! end weird Metropolis block

      else
        ! NOTE TO SELF - check this little fucker
        !old_e = (u_tot_run + ebar_y)**2
        !new_e = (u_tot_run + ebar_y - g_thr * float(L))**2
        old_e = ebar_y**2
        new_e = (ebar_y + ebar_inc)**2
        delta_e = new_e - old_e
        !delta_e = (0.5 - float(L) * ebar_y)

        if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
          .and.(exp(-beta*delta_e).gt.0.00000000001))) then
          ! this block is basically stolen from Michael
          ! not sure what's happening here tbh
          ebar_y = ebar_y + ebar_inc
          acceptg = acceptg + 1

          do x = 1,L
            do y = 1,L
              do z = 1,L
                e_y(x,y,z) = e_y(x,y,z) + e_inc
              end do
            end do
          end do
        end if ! end weird Metropolis block
      end if

      ! z-component
      chooser = rand()

      if (chooser.lt.0.5) then
        ! NOTE TO SELF - check this little fucker
        !old_e = (u_tot_run + ebar_z)**2
        !new_e = (u_tot_run + ebar_z - g_thr * float(L))**2
        old_e = ebar_z**2
        new_e = (ebar_z - ebar_inc)**2
        delta_e = new_e - old_e
        !delta_e = (0.5 - float(L) * ebar_z)

        if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
          .and.(exp(-beta*delta_e).gt.0.00000000001))) then
          ! this block is basically stolen from Michael
          ! not sure what's happening here tbh
          ebar_z = ebar_z - ebar_inc
          acceptg = acceptg + 1

          do x = 1,L
            do y = 1,L
              do z = 1,L
                e_z(x,y,z) = e_z(x,y,z) - e_inc
              end do
            end do
          end do
        end if ! end weird Metropolis block

      else
        ! NOTE TO SELF - check this little fucker
        !old_e = (u_tot_run + ebar_z)**2
        !new_e = (u_tot_run + ebar_z - g_thr * float(L))**2
        old_e = ebar_z**2
        new_e = (ebar_z + ebar_inc)**2
        delta_e = new_e - old_e
        !write (*,*) ebar_inc,ebar_z,old_e,new_e,delta_e
        !delta_e = (0.5 - float(L) * ebar_z)

        if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
          .and.(exp(-beta*delta_e).gt.0.00000000001))) then
          ! this block is basically stolen from Michael
          ! not sure what's happening here tbh
          ebar_z = ebar_z + ebar_inc
          acceptg = acceptg + 1

          do x = 1,L
            do y = 1,L
              do z = 1,L
                e_z(x,y,z) = e_z(x,y,z) + e_inc
              end do
            end do
          end do
        end if ! end weird Metropolis block
      end if

      end do ! end global update loop

      u_tot_run = 0.0
      do x = 1,L
        do y = 1,L
          do z = 1,L
            u_tot_run = u_tot_run + 0.5 * eps_0 *&
                        (e_x(x,y,z)**2 + e_y(x,y,z)**2 + e_z(x,y,z)**2)
          end do
        end do
      end do

  end if ! end glob.eq.1 block

  ! --- END OF UPDATE BLOCKS ---

end subroutine mc_sweep

subroutine measure(step_number)
  use common
  implicit none
  integer,intent(in) :: step_number
  integer :: x,y,z,i,j,k,m,p,s,kx,ky,kz,n
  real*8 :: norm_k,u_tot_run
  complex*16 :: imag, kdotx
  complex*16 :: e_ky, e_kz

  imag = (0.0,1.0)
  u_tot_run = 0.0

  ! array indexing: we don't sample at each step
  n = step_number / sample_interval

  do x = 1,L
    do y = 1,L
      do z = 1,L
        u_tot_run = u_tot_run + 0.5 * eps_0 *&
                    (e_x(x,y,z)**2 + e_y(x,y,z)**2 + e_z(x,y,z)**2)
      end do
    end do
  end do

  energy(n + 1) = u_tot_run
  sq_energy(n + 1) = u_tot_run**2

  ! --- FOURIER TRANSFORMS ---
  do kz = (-1*L/2)*bz, (L/2)*bz
    do ky = (-1*L/2)*bz, (L/2)*bz
      do kx = (-1*L/2)*bz, (L/2)*bz

        e_ky = (0.0,0.0)
        e_kz = (0.0,0.0)
        kdotx = (0.0,0.0)
        norm_k = 0.0

        ! for array indices
        i = kx + 1 + bz*(L/2)
        j = ky + 1 + bz*(L/2)
        k = kz + 1 + bz*(L/2)

        if (kx.eq.0.and.ky.eq.0.and.kz.eq.0) then
          norm_k = 0.0
        else
          norm_k = 1.0/(((2*pi/(L*lambda))**2)*dble(kx**2 + ky**2 + kz**2))
        end if

        e_kx_t(i,j,k,n) = (0.0,0.0)

        ! m,p,s are the real space coordinates
        do s = 1,L
          do p = 1,L
            do m = 1,L

              ! different offsets for x,y,z
              kdotx = ((-1)*imag*(2*pi/(L*lambda))*((m-(1.0/2))*kx + &
                      ((p-1)*ky) + ((s-1)*kz)))

              ! we want this for every step so we can
              ! average at the end to get field-field struc
              e_kx_t(i,j,k,n) = e_kx_t(i,j,k,n) + exp(kdotx)*e_x(m,p,s)

              kdotx = ((-1)*imag*(2*pi/(L*lambda))*((m-1)*kx + &
                      ((p-(1.0/2))*ky) + ((s-1)*kz)))

              e_ky = e_ky + exp(kdotx)*e_y(m,p,s)

              kdotx = ((-1)*imag*(2*pi/(L*lambda))*((m-1)*kx + &
                      ((p-1)*ky) + ((s-(1.0/2))*kz)))

              e_kz = e_kz + exp(kdotx)*e_z(m,p,s)

              if (v(m,p,s).ne.0) then ! calculate <++ + +->!

                ! FT of charge distribution
                kdotx = ((-1)*imag*2*pi*(((m-1)*kx/(L*lambda)) + &
                        ((p-1)*ky/(L*lambda)) + &
                        ((s-1)*kz/(L*lambda))))

                if (v(m,p,s).eq.-1) then
                  rho_k_m_t(i,j,k,n) = rho_k_m_t(i,j,k,n) + v(m,p,s) * exp(kdotx)
                end if

                if (v(m,p,s).eq.1) then ! take away <++>
                  rho_k_p_t(i,j,k,n) = rho_k_p_t(i,j,k,n) + v(m,p,s)*exp(kdotx)
                end if

                ! --- real space correlation function ---

                if (kx.gt.0.and.kx.le.L&
                  .and.ky.gt.0.and.ky.le.L&
                  .and.kz.gt.0.and.kz.le.L) then

                  ! this should sort it out. kx ky kz here are
                  ! the "r"s, m p s are the "zeros"
                  ! we want to do rho_+(0) rho_-(r)
                  if (v(m,p,s).eq.1) then
                    if (v(kx,ky,kz).eq.-1) then

                      x = abs(m - kx)
                      y = abs(p - ky)
                      z = abs(s - kz)

                      if (x.gt.L/2) then
                        x = L - x
                      end if
                      if (y.gt.L/2) then
                        y = L - y
                      end if
                      if (z.gt.L/2) then
                        z = L - z
                      end if

                      x = x + 1
                      y = y + 1
                      z = z + 1

                      dir_struc_n(x,y,z,n) = dir_struc_n(x,y,z,n) +&
                              v(m,p,s) * v(kx,ky,kz)
                    end if ! neg charge at kx,ky,kz
                  end if ! pos charge at m,p,s

                end if ! kx,ky,kz < L

              end if ! end v != 0 check

            end do ! end s loop
          end do ! end p loop
        end do ! end m loop

        ! normalise, idiot
        rho_k_m_t(i,j,k,n) = rho_k_m_t(i,j,k,n) / float(L**3)
        rho_k_p_t(i,j,k,n) = rho_k_p_t(i,j,k,n) / float(L**3)
        e_kx_t(i,j,k,n) = e_kx_t(i,j,k,n) / float(L**3)
        e_ky = e_ky / float(L**3)
        e_kz = e_kz / float(L**3)

        ch_ch(i,j,k,n) = (rho_k_p_t(i,j,k,n) * conjg(rho_k_m_t(i,j,k,n)))

        fe_fe(i,j,k,n) = (e_kx_t(i,j,k,n)*conjg(e_kx_t(i,j,k,n)))

        s_ab_n(1,1,i,j,k,n) = e_kx_t(i,j,k,n)*conjg(e_kx_t(i,j,k,n))
        s_ab_n(1,2,i,j,k,n) = e_kx_t(i,j,k,n)*conjg(e_ky)
        s_ab_n(1,3,i,j,k,n) = e_kx_t(i,j,k,n)*conjg(e_kz)
        s_ab_n(2,1,i,j,k,n) = e_ky*conjg(e_kx_t(i,j,k,n))
        s_ab_n(2,2,i,j,k,n) = e_ky*conjg(e_ky)
        s_ab_n(2,3,i,j,k,n) = e_ky*conjg(e_kz)
        s_ab_n(3,1,i,j,k,n) = e_kz*conjg(e_kx_t(i,j,k,n))
        s_ab_n(3,2,i,j,k,n) = e_kz*conjg(e_ky)
        s_ab_n(3,3,i,j,k,n) = e_kz*conjg(e_kz)

        ! end kz,ky,kx loops
      end do
    end do
  end do

end subroutine measure
