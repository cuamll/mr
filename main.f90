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
  real*8, dimension(4) :: timings
  integer, dimension(8) :: values

  call cpu_time(start_time)

  write (*,*) "Maggs-Rossetto algorithm on the simple cubic lattice."
  write (*,*) "Callum Gray, University College London, 2017."
  write (*,*) "Git revision ",revision
  write (*,*) "Repo at http://github.com/callumgray/mr"

  call date_and_time(VALUES=values)

  write (*,'(a24,I2.1,a1,I2.1,a1,I4.2,a1,I2.1,a1,I2.1,a1,I2.1)')&
    &" Current date and time: ",values(3),"/",values(2),"/"&
    &,values(1)," ",values(5),":",values(6),":",values(7)

  call setup_wrapper

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

  write(*,*)
  write(*,*) "--- END: CHARGE POSITIONS ---"

  tot_q = 0
  do k = 1,L
    do j = 1,L
      do i = 1,L

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
  integer :: i,j,k,m,charge,glob,totq,mu1,mu2,pm1
  integer, dimension(3) :: site
  real*8 :: eo1,eo2,eo3,eo4,en1,en2,en3,en4
  real*8 :: u_tot,u_tot_run,increment
  real*8 :: old_e, new_e, delta_e, g_thr

  charge = 0; totq = 0; mu1 = 0; mu2 = 0; pm1 = 0
  site =(/ 0, 0, 0 /)
  eo1 = 0.0; eo2 = 0.0; eo3 = 0.0; eo4 = 0.0
  en1 = 0.0; en2 = 0.0; en3 = 0.0; en4 = 0.0
  u_tot = 0.0; u_tot_run = 0.0; increment = 0.0
  old_e = 0.0; new_e = 0.0; delta_e = 0.0; g_thr = 1 / float(L)

  ! --- START OF UPDATE BLOCKS ---

  ! --- CHARGE HOP UPDATE ---

  do i = 1, int(L**3 * hop_ratio)

    ! NOTE TO SELF: this whole procedure assumes
    ! single-valued charges only.

    ! pick a random site
    site = (/ int(rand() * L) + 1,&
               &int(rand() * L) + 1,&
               &int(rand() * L) + 1 /)

    mu1 = floor(3*rand())+1

    ! check we're not doing anything weird with mu1ltiple charges
    if (abs(v(site(1),site(2),site(3))).gt.1) then
      write (*,*) "Charge at ",site(1),site(2),site(3),&
        " = ",v(site(1),site(2),site(3)),". Exiting."
      write (*,*)
      stop
    end if

    ! pick a non-zero charge - canonical!
    if (v(site(1),site(2),site(3)).ne.0) then

      charge = v(site(1),site(2),site(3))
      eo1 = e_field(mu1,site(1),site(2),site(3))

      ! this takes care of sign issues when hopping
      increment = charge / (eps_0 * lambda**2)

      if (rand().lt.0.5) then ! move it "negative"

        en1 = eo1 + increment
        old_e = 0.5 * eps_0 * eo1**2
        new_e = 0.5 * eps_0 * en1**2
        delta_e = new_e - old_e

        ! get negative in the mu1 direction
        site(mu1) = neg(site(mu1))

        if (v(site(1),site(2),site(3)).eq.0) then
          if ((delta_e.lt.0.0).or.(exp((-beta)*delta_e).gt.rand())) then

            v(site(1),site(2),site(3)) = charge

            ! go back to the original site and set the charge to 0
            site(mu1) = pos(site(mu1))
            v(site(1),site(2),site(3)) = 0
            e_field(mu1,site(1),site(2),site(3)) = en1
            ebar(mu1) = ebar(mu1) + increment

            accepth = accepth + 1
            u_tot_run = u_tot_run + delta_e

            end if
          end if

      else ! move it "positive"

        en1 = eo1 - increment
        old_e = 0.5 * eps_0 * eo1**2
        new_e = 0.5 * eps_0 * en1**2
        delta_e = new_e - old_e

        ! get pos in the mu1 direction
        site(mu1) = pos(site(mu1))

        if (v(site(1),site(2),site(3)).eq.0) then
          if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

            v(site(1),site(2),site(3)) = charge

            ! go back to the original site and set the charge to 0
            site(mu1) = neg(site(mu1))
            v(site(1),site(2),site(3)) = 0
            e_field(mu1,site(1),site(2),site(3)) = en1
            ebar(mu1) = ebar(mu1) - increment

            accepth = accepth + 1
            u_tot_run = u_tot_run + delta_e

          end if
        end if

      end if ! end "positive" / "negative" choice
    end if ! end charge.ne.0 block

  end do ! end charge hop sweep

  if (ebar(mu1).gt.(g_thr).or.ebar(mu1).lt.((-1)*g_thr)) then
    glob = 1
  else
    glob = 0
  end if

  mu1 = 0; increment = 0.0; 
  u_tot = 0.0
  u_tot = 0.5 * eps_0 * sum(e_field * e_field)

  ! --- ROTATIONAL UPDATE ---

  do i = 1,int(3 * L**3 * rot_ratio)

    eo1 = 0.0; eo2 = 0.0; eo3 = 0.0; eo4 = 0.0
    en1 = 0.0; en2 = 0.0; en3 = 0.0; en4 = 0.0
    site = (/ 0, 0, 0/)

    ! pick at random from interval [-Delta_max, +Delta_max]
    increment = 2 * rot_delt * (rand() - 0.5)

    ! pick xy (1,2), yz (2,3), or zx (3,1) plaquette
    mu1 = floor(3*rand())+1
    mu2 = mod(mu1,3) + 1

    ! give us a coordinate (i,j,k)
    site = (/ int(rand() * L) + 1,&
               &int(rand() * L) + 1,&
               &int(rand() * L) + 1 /)
    ! site(mu1) is the coordinate in the first direction
    ! site(mu2) is the coordinate in the second

    ! my convention is that an xy (1,2) plaquette is defined by
    ! field links (1,i,j,k),(2,i,j,k),(1,neg(i),j,k),(2,i,neg(j),k)
    ! so the first two are easy:
    eo1 = e_field(mu1,site(1),site(2),site(3))
    eo2 = e_field(mu2,site(1),site(2),site(3))

    ! if e.g. we picked xy, the next line does x -> neg(x)
    site(mu1) = neg(site(mu1))
    eo4 = e_field(mu1,site(1),site(2),site(3))
    ! and now we need to put it back
    site(mu1) = pos(site(mu1))

    ! this does y -> neg(y)
    site(mu2) = neg(site(mu2))
    eo3 = e_field(mu2,site(1),site(2),site(3))
    ! and put it back
    site(mu2) = pos(site(mu2))
    ! so now we have (x,y,z) again

    en1 = eo1 + increment
    en2 = eo2 - increment
    en3 = eo3 + increment
    en4 = eo4 - increment

    old_e = 0.5 * eps_0 * (eo1**2 + eo2**2 + eo3**2 + eo4**2)
    new_e = 0.5 * eps_0 * (en1**2 + en2**2 + en3**2 + en4**2)
    delta_e = new_e - old_e

    if ((delta_e.lt.0.0).or.(exp((-beta)*delta_e).gt.rand())) then

      e_field(mu1,site(1),site(2),site(3)) = en1
      e_field(mu2,site(1),site(2),site(3)) = en2

      site(mu1) = neg(site(mu1))
      e_field(mu1,site(1),site(2),site(3)) = en4
      site(mu1) = pos(site(mu1))

      site(mu2) = neg(site(mu2))
      e_field(mu2,site(1),site(2),site(3)) = en3
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
    do i = 1,int(L**3 * g_ratio)

      increment = q/(L**2 * eps_0)

      do mu1 = 1,3

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

        if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
          .and.(exp(-beta*delta_e).gt.0.00000000001))) then
          ! this block is basically stolen from Michael
          ! not sure what's happening here tbh
          !write (*,*) ebar(1),ebar(2),ebar(3)
          ebar(mu1) = ebar(mu1) + pm1 * increment
          !write (*,*) ebar(1),ebar(2),ebar(3)
          !do m = 1,L
          !  do k = 1,L
          !    do j = 1,L
          !      e_field(mu1,j,k,m) = e_field(mu1,j,k,m) + pm1 * (increment / L**3)
                e_field(mu1,:,:,:) = e_field(mu1,:,:,:) + pm1 * (increment / L**3)
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

  ! --- END OF UPDATE BLOCKS ---

end subroutine mc_sweep

subroutine measure(step_number)
  use common
  use linear_solver
  implicit none
  logical :: field_ch_exist
  character(55) :: f_ch_format
  integer,intent(in) :: step_number
  integer :: x,y,z,i,j,k,m,p,s,kx,ky,kz,n
  real*8 :: norm_k,u_tot_run
  complex*16 :: imag, kdotx
  complex*16 :: e_ky, e_kz

  imag = (0.0,1.0)
  f_ch_format = '(I5.1, I3.1, I3.1, I3.1, ES18.9, ES18.9, ES18.9, I3.1)'

  ! array indexing: we don't sample at each step
  n = step_number / sample_interval

  u_tot_run = 0.0
  u_tot_run = 0.5 * eps_0 * sum(e_field * e_field)

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


  ! get irrotational part of field
  call linsol

  write(15) e_field
  write(15) mnphi
  write(15) v

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

        !if (i.le.L.and.j.le.L.and.k.le.L) then

          ! could do this unformatted if we're just
          ! gonna read it back into fortran anyway
          !write (15,f_ch_format) n, i, j, k, e_field(1,i,j,k),&
          !         e_field(2,i,j,k), e_field(3,i,j,k), v(i,j,k)

        !end if

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
              e_kx_t(i,j,k,n) = e_kx_t(i,j,k,n)&
                              & + exp(kdotx)*e_field(1,m,p,s)

              kdotx = ((-1)*imag*(2*pi/(L*lambda))*((m-1)*kx + &
                      ((p-(1.0/2))*ky) + ((s-1)*kz)))

              e_ky = e_ky + exp(kdotx)*e_field(2,m,p,s)

              kdotx = ((-1)*imag*(2*pi/(L*lambda))*((m-1)*kx + &
                      ((p-1)*ky) + ((s-(1.0/2))*kz)))

              e_kz = e_kz + exp(kdotx)*e_field(3,m,p,s)

              if (v(m,p,s).ne.0) then ! calculate <++ + +->!

                ! FT of charge distribution
                kdotx = ((-1)*imag*2*pi*(((m-1)*kx/(L*lambda)) + &
                        ((p-1)*ky/(L*lambda)) + &
                        ((s-1)*kz/(L*lambda))))

                if (v(m,p,s).eq.-1) then
                  rho_k_m_t(i,j,k,n) = rho_k_m_t(i,j,k,n)&
                                     & + v(m,p,s) * exp(kdotx)
                end if

                if (v(m,p,s).eq.1) then ! take away <++>
                  rho_k_p_t(i,j,k,n) = rho_k_p_t(i,j,k,n)&
                                     & + v(m,p,s)*exp(kdotx)
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

  close(15)

end subroutine measure
