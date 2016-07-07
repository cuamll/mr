! maggs-rossetto algorithm on cubic lattice
! i'll write proper docs at some point
! callum gray, UCL / ENS Lyon, 2016. do wat u like
program mr

  use common
  use linear_solver
  use io
  use setup
  implicit none
  integer :: i, j, k, n, row, col, tot_q

  call do_setup

  call upcan()

  tot_q = 0

  ! best to check again, just in case
  do i = 1,L
    do j = 1,L
      do k = 1,L

        tot_q = tot_q + abs(v(i,j,k))

      end do
    end do
  end do

  if (tot_q.ne.0) then
    call linsol
  end if

  write(*,*)

  call linalg

  call write_output

  call deallocations

  stop

end program mr

! canonical ensemble - MC / Metropolis updates
subroutine upcan()
  use common
  implicit none
  integer :: x,y,z,n,charge,glob,i,j,k,m,p,s,pm1,kx,ky,kz
  real*8 :: eo1,eo2,eo3,eo4,en1,en2,en3,en4
  real*8 :: u_tot,u_tot_run,avg_e,avg_e2,one
  real*8 :: dot,dot_avg,ebar_inc,e_inc,u_diff
  real*8 :: hop_inc, old_e, new_e, delta_e, totq, g_thr
  real :: chooser, delta
  complex*16 :: imag, kdotx

  glob = 0
  totq = 0
  u_tot_run = 0.0
  g_thr = 1 / float(L)
  accepth = 0
  acceptr=0
  acceptg=0
  one=1.0
  imag = (0.0,1.0)

  write(*,*)
  write(*,*) " --- start: charge positions ---"

  do i = 1,L
    do j = 1,L
      do k = 1,L
        u_tot_run = u_tot_run + e_x(i,j,k)**2 + &
                    e_y(i,j,k)**2 + e_z(i,j,k)**2

        totq = totq + abs(v(i,j,k))
        if (v(i,j,k).ne.0) then
          write (*,*) i, j, k, v(i,j,k)
        end if
      end do
    end do
  end do

  write (*,*) "total charges: ",totq
  write(*,*)

  energy(1) = u_tot_run
  sq_energy(1) = u_tot_run**2

  ! charge hop sweep
  do n = 1,iterations

  !write (*,*) "utot at start of step ",n," = ",u_tot_run
    u_tot_run = u_tot

    ! charge hop sweep
    do i = 1, int(L**3 * hop_ratio)

      ! pick a random site
      x = int(rand() * L) + 1
      y = int(rand() * L) + 1
      z = int(rand() * L) + 1

      ! pick a non-zero charge - we want to keep it canonical
      if (v(x,y,z).ne.0) then

        charge = v(x,y,z)
        ! i think this takes care of any sign issues when hopping
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

    !write (*,*) "utot after charge hops = ",u_tot_run

    ! --- ROTATIONAL UPDATE ---

    ! NOTE TO SELF: might need a + 1 next to that int cast
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

      if (rand().gt.0.5) then
        pm1 = 1
      else
        pm1 = -1
      end if

      chooser=rand()
      if (chooser.lt.(1.0/3.0)) then ! xy-plane plaquette

        eo1 = e_x(x,y,z)
        eo2 = e_y(x,y,z)
        eo3 = e_x(x,neg(y),z)
        eo4 = e_y(neg(x),y,z)
        !write(*,*) x,y,z,delta,pm1
        !write(*,*) eo1,eo2,eo3,eo4

        !en1 = eo1 + (delta * sign(one, eo1))
        !write(*,*) (delta * sign(one, eo1))
        !en2 = eo2 - (delta * sign(one, eo2))
        !write(*,*) (delta * sign(one, eo2))
        !en3 = eo3 + (delta * sign(one, eo3))
        !write(*,*) (delta * sign(one, eo3))
        !en4 = eo4 - (delta * sign(one, eo4))
        !write(*,*) (delta * sign(one, eo4))
        en1 = eo1 + delta
        en2 = eo2 - delta
        en3 = eo3 + delta
        en4 = eo4 - delta
        !write(*,*) en1,en2,en3,en4

        old_e = 0.5 * eps_0 * (eo1**2 + eo2**2 + eo3**2 + eo4**2)
        new_e = 0.5 * eps_0 * (en1**2 + en2**2 + en3**2 + en4**2)
        delta_e = new_e - old_e

        if ((delta_e.lt.0.0).or.(exp((-beta)*delta_e).gt.rand())) then

          !write (*,*) "x-y plane rot:"
          !write (*,*) x, y, z
          !write (*,*) x, neg(y), z
          !write (*,*) neg(x), y, z

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

        !en1 = eo1 + (delta * sign(one, eo1))
        !en2 = eo2 - (delta * sign(one, eo2))
        !en3 = eo3 + (delta * sign(one, eo3))
        !en4 = eo4 - (delta * sign(one, eo4))
        en1 = eo1 + delta
        en2 = eo2 - delta
        en3 = eo3 + delta
        en4 = eo4 - delta

        old_e = 0.5 * eps_0 * (eo1**2 + eo2**2 + eo3**2 + eo4**2)
        new_e = 0.5 * eps_0 * (en1**2 + en2**2 + en3**2 + en4**2)
        delta_e = new_e - old_e

        if ((delta_e.lt.0.0).or.(exp((-beta)*delta_e).gt.rand())) then

          !write (*,*) "x-z plane rot:"
          !write (*,*) x, y, z
          !write (*,*) x, y, neg(z)
          !write (*,*) neg(x), y, z

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

        !en1 = eo1 + (delta * sign(one, eo1))
        !en2 = eo2 - (delta * sign(one, eo2))
        !en3 = eo3 + (delta * sign(one, eo3))
        !en4 = eo4 - (delta * sign(one, eo4))

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

    ! --- HARMONIC UPDATE ---
    ! e bar update
    if (glob.eq.1) then
      !write (*,*) "step ",n," start of harmonic update: ",u_tot_run,"total moves accepted ",acceptg

      ! NOTE TO SELF: again this int cast prob needs changing
      do i = 1,int(L**3 * g_ratio)

        ebar_inc = q/(L * eps_0)
        e_inc = ebar_inc / L**3

        ! x-component
        chooser = rand()

        if (chooser.lt.0.5) then
          ! NOTE TO SELF - check this little fucker
          !old_e = (u_tot_run + ebar_x)**2
          !new_e = (u_tot_run + ebar_x - g_thr * float(L))**2
          old_e = u_tot_run**2
          new_e = (u_tot_run - ebar_inc)**2
          delta_e = new_e - old_e
          !delta_e = (0.5 - float(L) * ebar_x)

          if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_x = ebar_x - ebar_inc 
            acceptg = acceptg + 1
            u_tot_run = 0.0 
            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_x(j,k,m) = e_x(j,k,m) - e_inc
                  u_tot_run = u_tot_run + 0.5 * eps_0 *&
                              (e_x(i,j,k)**2 + e_y(i,j,k)**2 + e_z(i,j,k)**2)
                end do
              end do
            end do
            if (n.eq.16) then
              write (*,*) "x - ",old_e,new_e,delta_e,u_tot_run
            end if
          end if ! end weird Metropolis block

        else
          ! NOTE TO SELF - check this little fucker
          !old_e = (u_tot_run + ebar_x)**2
          !new_e = (u_tot_run + ebar_x - g_thr * float(L))**2
          old_e = u_tot_run**2
          new_e = (u_tot_run + ebar_inc)**2
          delta_e = new_e - old_e
          !delta_e = (0.5 - float(L) * ebar_x)

          if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_x = ebar_x + ebar_inc 
            acceptg = acceptg + 1

            u_tot_run = 0.0
            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_x(j,k,m) = e_x(j,k,m) + e_inc
                  u_tot_run = u_tot_run + 0.5 * eps_0 *&
                              (e_x(i,j,k)**2 + e_y(i,j,k)**2 + e_z(i,j,k)**2)
                end do
              end do
            end do
            if (n.eq.16) then
              write (*,*) "x + ",old_e,new_e,delta_e,u_tot_run
            end if
          end if ! end weird Metropolis block
        end if

        ! y-component
        chooser = rand()

        if (chooser.lt.0.5) then
          ! NOTE TO SELF - check this little fucker
          !old_e = (u_tot_run + ebar_y)**2
          !new_e = (u_tot_run + ebar_y - g_thr * float(L))**2
          old_e = u_tot_run**2
          new_e = (u_tot_run - ebar_inc)**2
          delta_e = new_e - old_e
          !delta_e = (0.5 - float(L) * ebar_y)

          if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_y = ebar_y - ebar_inc 
            acceptg = acceptg + 1
            u_tot_run = 0.0
            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_y(j,k,m) = e_y(j,k,m) - e_inc
                  u_tot_run = u_tot_run + 0.5 * eps_0 *&
                              (e_x(i,j,k)**2 + e_y(i,j,k)**2 + e_z(i,j,k)**2)
                end do
              end do
            end do
            if (n.eq.16) then
              write (*,*) "y - ",old_e,new_e,delta_e,u_tot_run
            end if
          end if ! end weird Metropolis block

        else
          ! NOTE TO SELF - check this little fucker
          !old_e = (u_tot_run + ebar_y)**2
          !new_e = (u_tot_run + ebar_y - g_thr * float(L))**2
          old_e = u_tot_run**2
          new_e = (u_tot_run + ebar_inc)**2
          delta_e = new_e - old_e
          !delta_e = (0.5 - float(L) * ebar_y)

          if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_y = ebar_y + ebar_inc 
            acceptg = acceptg + 1

            u_tot_run = 0.0
            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_y(j,k,m) = e_y(j,k,m) + e_inc
                  u_tot_run = u_tot_run + 0.5 * eps_0 *&
                              (e_x(i,j,k)**2 + e_y(i,j,k)**2 + e_z(i,j,k)**2)
                end do
              end do
            end do
            if (n.eq.16) then
              write (*,*) "y + ",old_e,new_e,delta_e,u_tot_run
            end if
          end if ! end weird Metropolis block
        end if

        ! z-component
        chooser = rand()

        if (chooser.lt.0.5) then
          ! NOTE TO SELF - check this little fucker
          !old_e = (u_tot_run + ebar_z)**2
          !new_e = (u_tot_run + ebar_z - g_thr * float(L))**2
          old_e = u_tot_run**2
          new_e = (u_tot_run - ebar_inc)**2
          delta_e = new_e - old_e
          !delta_e = (0.5 - float(L) * ebar_z)

          if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_z = ebar_z - ebar_inc 
            acceptg = acceptg + 1
            u_tot_run = 0.0
            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_z(j,k,m) = e_z(j,k,m) - e_inc
                  u_tot_run = u_tot_run + 0.5 * eps_0 *&
                              (e_x(i,j,k)**2 + e_y(i,j,k)**2 + e_z(i,j,k)**2)
                end do
              end do
            end do
          end if ! end weird Metropolis block

        else
          ! NOTE TO SELF - check this little fucker
          !old_e = (u_tot_run + ebar_z)**2
          !new_e = (u_tot_run + ebar_z - g_thr * float(L))**2
          old_e = u_tot_run**2
          new_e = (u_tot_run - ebar_inc)**2
          delta_e = new_e - old_e
          !delta_e = (0.5 - float(L) * ebar_z)

          if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_z = ebar_z + ebar_inc 
            acceptg = acceptg + 1

            u_tot_run = 0.0
            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_z(j,k,m) = e_z(j,k,m) + e_inc
                  u_tot_run = u_tot_run + 0.5 * eps_0 *&
                              (e_x(i,j,k)**2 + e_y(i,j,k)**2 + e_z(i,j,k)**2)
                end do
              end do
            end do
            if (n.eq.16) then
              write (*,*) "z + ",old_e,new_e,delta_e,u_tot_run
            end if
          end if ! end weird Metropolis block
        end if

        end do ! end global update loop
        !write (*,*) "step ",n," end of harmonic update: ",u_tot_run,"total moves accepted ",acceptg

    end if ! end glob.eq.1 block
    
    u_tot = 0.0
    do i = 1,L
      do j = 1, L
        do k = 1,L
        u_tot = u_tot + 0.5 * eps_0 * (e_x(i,j,k)**2 + e_y(i,j,k)**2 + e_z(i,j,k)**2)
        end do
      end do
    end do

    u_diff = u_tot_run - u_tot
    if (abs(u_diff).ge.0.00001) then
      write (*,*) "u_diff = ",u_diff
    end if

    !write(*,*) n,glob,u_tot_run,u_tot,ebar_x,ebar_y,ebar_z

  ! replace with u_tot_run maybe?
  energy(n + 1) = u_tot
  sq_energy(n + 1) = u_tot**2

  avg_e = 0.0
  avg_e2 = 0.0
  do m=1,n+1
    avg_e = avg_e + energy(m)
    avg_e2 = avg_e2 + sq_energy(m)
  end do

  avg_e = avg_e / (n + 1)
  avg_e2 = avg_e2 / (n + 1)

    ! Fourier transform
    do kx = -L/2, L/2
      do ky = -L/2, L/2
        do kz = -L/2, L/2

          ! for array indices
          i = kx + 1 + L/2
          j = ky + 1 + L/2
          k = kz + 1 + L/2

          e_kx(i,j,k) = 0.0
          e_ky(i,j,k) = 0.0
          e_kz(i,j,k) = 0.0

          ! vec(q) = vec(0) term:
          if (i.eq.0.and.j.eq.0.and.k.eq.0) then
            write (*,*) "q=0 term; need to work this out"
            write (*,*) "e_kx(",i,j,k,") = ",e_kx(i,j,k)
            CYCLE
          end if

          ! m,p,s are the real space coordinates
          do m = 1,L
            do p = 1,L
              do s = 1,L

                kdotx = ((-1)*imag*2*pi*(((m-1)*kx/(L*lambda)) + &
                        ((p-1)*ky/(L*lambda)) + &
                        ((s-1)*kz/(L*lambda))))

                e_kx(i,j,k) = e_kx(i,j,k) + e**(kdotx)*e_x(m,p,s)
                e_ky(i,j,k) = e_ky(i,j,k) + e**(kdotx)*e_y(m,p,s)
                e_kz(i,j,k) = e_kz(i,j,k) + e**(kdotx)*e_z(m,p,s)

                if (v(m,p,s).ne.0) then ! there's a charge

                  ! FT of charge distribution
                  kdotx = ((-1)*imag*2*pi*(((m-1)*kx/(L*lambda)) + &
                          ((p-1)*ky/(L*lambda)) + &
                          ((s-1)*kz/(L*lambda))))
                  rho_k(i,j,k) = rho_k(i,j,k) + v(m,p,s) * e**(kdotx)

                end if

              end do
            end do
          end do

          ! normalise, idiot
          e_kx(i,j,k) = e_kx(i,j,k) / L**3
          e_ky(i,j,k) = e_ky(i,j,k) / L**3
          e_kz(i,j,k) = e_kz(i,j,k) / L**3
          rho_k(i,j,k) = rho_k(i,j,k) / L**3

          ! now we need to do a (?) thermal (?) average to get correlations
          ! this is the dot of the fourier space field with itself at i,j,k
          dot = (e_kx(i,j,k) * CONJG(e_kx(i,j,k))) +&
                (e_ky(i,j,k) * CONJG(e_ky(i,j,k))) +&
                (e_kz(i,j,k) * CONJG(e_kz(i,j,k)))

          !write (*,*) "E \cdot E (",2*pi*i/(L*lambda),2*pi*j/(L*lambda),2*pi*k/(L*lambda),") = ",dot

          dot_avg = dot_avg + dot

          ! second part is weirdly addressed bc of loop structure
          ! ((-1) * k_mu) + 1 + L/2 --> (-\vec{k}) essentially
          ch_ch(i,j,k,n) = ch_ch(i,j,k,n) + (rho_k(i,j,k) * rho_k((-1)*kx + 1 + L/2, (-1)*ky + 1 + L/2, (-1)*kz + 1 + L/2))

        end do
      end do
    end do

  !write(*,*) "<U^2>, <U>^2, ratio = ",avg_e2,avg_e*avg_e,(avg_e2/(avg_e**2))
  !write (*,*) "step = ",n," ebar = ",ebar_x,ebar_y,ebar_z

  end do ! end iteration loop

  write(*,*)
  write(*,*) " --- end: charge positions ---"

  totq = 0
  do j = 1,L
    do k = 1,L
      do m = 1,L
        totq = totq + abs(v(j,k,m))
        if (v(j,k,m).ne.0) then
          write (*,*) j, k, m, v(j,k,m)
        end if
      end do
    end do
  end do
  write (*,*) "total charges: ",totq
  write(*,*)

  write (*,*)
  write (*,*) '----- move stats -----'
  write (*,*) 'hop moves  = ',accepth,float(accepth) / (iterations * totq)
  write (*,*) 'rot moves  = ',acceptr,float(acceptr) / (iterations * L**3 * rot_ratio)
  write (*,*) 'ebar moves = ',acceptg,float(acceptg) / 3 / (iterations * L**3 * g_ratio)

end subroutine upcan
