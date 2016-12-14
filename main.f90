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

  do i = 1,L
    do j = 1,L
      do k = 1,L

        tot_q = tot_q + abs(v(2*i - 1,2*j - 1,2*k - 1))

      end do
    end do
  end do

  if (tot_q.ne.0) then

    call linsol
    call linalg

  end if

  write(*,*)

  call write_output

  call deallocations

  stop

end program mr

! canonical ensemble - MC / Metropolis updates
subroutine upcan()
  use common
  implicit none
  character(2) :: configs
  integer :: x,y,z,n,charge,glob,i,j,k,m,p,s,pm1,kx,ky,kz
  integer, dimension(4) :: dims
  real*8 :: eo1,eo2,eo3,eo4,en1,en2,en3,en4
  real*8 :: u_tot,u_tot_run,avg_e,avg_e2,one
  real*8 :: dot,dot_avg,ebar_inc,e_inc,u_diff
  real*8 :: hop_inc, old_e, new_e, delta_e, totq, g_thr
  real :: chooser, delta
  complex*16 :: imag, kdotx
  type(C_PTR) :: plan_ch, plan_e
  real(C_DOUBLE), dimension(2*L,2*L,2*L) :: fftw_ch_in
  real(C_DOUBLE), dimension(3,2*L,2*L,2*L) :: fftw_e_in
  complex(C_DOUBLE_COMPLEX), dimension(L+1,2*L,2*L) :: fftw_ch_out
  complex(C_DOUBLE_COMPLEX), dimension(L+1,3,2*L,2*L) :: fftw_e_out

  glob = 0
  totq = 0
  u_tot_run = 0.0
  g_thr = 1 / float(L)
  accepth = 0
  acceptr=0
  acceptg=0
  one=1.0
  imag = (0.0,1.0)
  fftw_ch_in = 0.0
  fftw_e_in = 0.0

  ! test electric field structure factor
  configs = "FE"

  ! FFTW stuff
  dims = (/ 3,2*L,2*L,2*L /)

  ! fftw plans for the transforms
  plan_ch = fftw_plan_dft_r2c_3d(2*L,2*L,2*L,fftw_ch_in,fftw_ch_out,FFTW_ESTIMATE)
  plan_e = fftw_plan_dft_r2c(4,dims,fftw_e_in,fftw_e_out,FFTW_ESTIMATE)

  write(*,*)
  write(*,*) "beta = ",beta
  write(*,*) " --- start: charge positions ---"

  do i = 1,L
    do j = 1,L
      do k = 1,L
        u_tot_run = u_tot_run + e_x(2*i,2*j - 1,2*k - 1)**2 + &
                    e_y(2*i - 1,2*j,2*k - 1)**2 + e_z(2*i - 1,2*j - 1,2*k)**2

        totq = totq + abs(v(2*i - 1,2*j - 1,2*k - 1))
        if (v(2*i - 1,2*j - 1,2*k - 1).ne.0) then
          write (*,*) i, j, k, v(2*i - 1,2*j - 1,2*k - 1)
        end if
      end do
    end do
  end do

  write (*,*) "total charges: ",totq
  write(*,*)

  energy(1) = u_tot_run
  sq_energy(1) = u_tot_run**2


  do n = 1,iterations

    u_tot_run = u_tot

    ! --- START OF UPDATE BLOCKS ---

    ! --- CHARGE HOP UPDATE ---

    do i = 1, int(L**3 * hop_ratio)

      ! pick a random site
      x = int(rand() * L) + 1
      y = int(rand() * L) + 1
      z = int(rand() * L) + 1

      ! pick a non-zero charge - canonical!
      if (v(2*x - 1,2*y - 1,2*z - 1).ne.0) then

        charge = v(2*x - 1,2*y - 1,2*z - 1)
        ! this takes care of sign issues when hopping
        hop_inc = charge / (eps_0 * lambda**2)
        chooser = rand()

        if (chooser.lt.(1.0 / 3.0)) then
          ! x component
          eo1 = e_x(2*x,2*y - 1,2*z - 1)
          chooser = rand()

          if (chooser.lt.0.5) then ! we try and move the fucker left
            en1 = eo1 + hop_inc
            old_e = 0.5 * eps_0 * eo1**2
            new_e = 0.5 * eps_0 * en1**2
            delta_e = new_e - old_e
            ! i think we can just set v(2*x - 1,2*y - 1,2*z - 1) = 0
            ! if we're enforcing |v| <= 1
            if ((abs(v(2*x - 1,2*y - 1,2*z - 1)).le.1).and.(v(neg(2*x - 1),2*y - 1,2*z - 1).eq.0)) then
              if ((delta_e.lt.0.0).or.(exp(((-1.0)*(beta)) * delta_e).gt.rand())) then

                accepth = accepth + 1
                e_x(2*x,2*y - 1,2*z - 1) = en1
                v(neg(2*x - 1),2*y - 1,2*z - 1) = charge
                v(2*x - 1,2*y - 1,2*z - 1) = 0
                ebar_x = ebar_x - hop_inc
                u_tot_run = u_tot_run + delta_e

              end if
            end if

          else ! move the fucker to the right

            en1 = eo1 - hop_inc
            old_e = 0.5 * eps_0 * eo1**2
            new_e = 0.5 * eps_0 * en1**2
            delta_e = new_e - old_e

            if ((abs(v(2*x - 1,2*y - 1,2*z - 1)).le.1).and.(v(pos(2*x - 1),2*y - 1,2*z - 1).eq.0)) then
              if ((delta_e.lt.0.0).or.(exp(((-1.0)*(beta)) * delta_e).gt.rand())) then

                accepth = accepth + 1
                e_x(2*x,2*y - 1,2*z - 1) = en1
                v(pos(2*x - 1),2*y - 1,2*z - 1) = charge
                v(2*x - 1,2*y - 1,2*z - 1) = 0
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

          eo1 = e_y(2*x - 1,2*y,2*z - 1)
          chooser = rand()

          if (chooser.lt.0.5) then ! we try and move the fucker left
            en1 = eo1 + hop_inc
            old_e = 0.5 * eps_0 * eo1**2
            new_e = 0.5 * eps_0 * en1**2
            delta_e = new_e - old_e
            ! i think we can just set v(2*x - 1,2*y - 1,2*z - 1) = 0
            ! if we're enforcing |v| <= 1
            if ((abs(v(2*x - 1,2*y - 1,2*z - 1)).le.1).and.(v(2*x - 1,neg(2*y - 1),2*z - 1).eq.0)) then
              if ((delta_e.lt.0.0).or.(exp(((-1.0)*(beta)) * delta_e).gt.rand())) then

                accepth = accepth + 1
                e_y(2*x - 1,2*y,2*z - 1) = en1
                v(2*x - 1,neg(2*y - 1),2*z - 1) = charge
                v(2*x - 1,2*y - 1,2*z - 1) = 0
                ebar_y = ebar_y - hop_inc
                u_tot_run = u_tot_run + delta_e

              end if
            end if

          else ! move the fucker to the right

            en1 = eo1 - hop_inc

            old_e = 0.5 * eps_0 * eo1**2
            new_e = 0.5 * eps_0 * en1**2
            delta_e = new_e - old_e

            if ((abs(v(2*x - 1,2*y - 1,2*z - 1)).le.1).and.(v(2*x - 1,pos(2*y - 1),2*z - 1).eq.0)) then
              if ((delta_e.lt.0.0).or.(exp(((-1.0)*(beta)) * delta_e).gt.rand())) then

                accepth = accepth + 1
                e_y(2*x - 1,2*y,2*z - 1) = en1
                v(2*x - 1,pos(2*y - 1),2*z - 1) = charge
                v(2*x - 1,2*y - 1,2*z - 1) = 0
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

          eo1 = e_z(2*x - 1,2*y - 1,2*z)
          chooser = rand()

          if (chooser.lt.0.5) then ! we try and move the fucker left
            en1 = eo1 + hop_inc
            old_e = 0.5 * eps_0 * eo1**2
            new_e = 0.5 * eps_0 * en1**2
            delta_e = new_e - old_e
            ! i think we can just set v(2*x - 1,2*y - 1,2*z - 1) = 0
            ! if we're enforcing |v| <= 1
            if ((abs(v(2*x - 1,2*y - 1,2*z - 1)).le.1).and.(v(2*x - 1,2*y - 1,neg(2*z - 1)).eq.0)) then
              if ((delta_e.lt.0.0).or.(exp(((-1.0)*(beta)) * delta_e).gt.rand())) then

                accepth = accepth + 1
                e_z(2*x - 1,2*y - 1,2*z) = en1
                v(2*x - 1,2*y - 1,neg(2*z - 1)) = charge
                v(2*x - 1,2*y - 1,2*z - 1) = 0
                ebar_z = ebar_z - hop_inc
                u_tot_run = u_tot_run + delta_e

              end if
            end if

          else ! move the fucker to the right

            en1 = eo1 - hop_inc
            old_e = 0.5 * eps_0 * eo1**2
            new_e = 0.5 * eps_0 * en1**2
            delta_e = new_e - old_e

            if ((abs(v(2*x - 1,2*y - 1,2*z - 1)).le.1).and.(v(2*x - 1,2*y - 1,pos(2*z - 1)).eq.0)) then
              if ((delta_e.lt.0.0).or.(exp(((-1.0)*(beta)) * delta_e).gt.rand())) then

                accepth = accepth + 1
                e_z(2*x - 1,2*y - 1,2*z) = en1
                v(2*x - 1,2*y - 1,pos(2*z - 1)) = charge
                v(2*x - 1,2*y - 1,2*z - 1) = 0
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

      end if ! end of v(2*x - 1,2*y - 1,2*z - 1).ne.0 block

    end do ! end charge hop sweep

    u_tot = 0.0
    do i = 1,L
      do j = 1, L
        do k = 1,L
        u_tot = u_tot + 0.5 * eps_0 * (e_x(2*i,2*j - 1,2*k - 1)**2 + e_y(2*i - 1,2*j,2*k - 1)**2 + e_z(2*i - 1,2*j - 1,2*k)**2)
        end do
      end do
    end do

    !u_diff = u_tot_run - u_tot
    !if (abs(u_diff).ge.0.00001) then
    !  write (*,*) "HOP: u_tot_run = ",u_tot_run," u_tot = ",u_tot," u_diff = ",u_diff
    !end if

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

        eo1 = e_x(2*x,2*y - 1,2*z - 1)
        eo2 = e_y(2*x - 1,2*y,2*z - 1)
        eo3 = e_x(2*x,neg(2*y - 1),2*z - 1)
        eo4 = e_y(neg(2*x - 1),2*y,2*z - 1)
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

        if ((delta_e.lt.0.0).or.(exp(((-1.0)*(beta))*delta_e).gt.rand())) then

          !write (*,*) "x-y plane rot:"
          !write (*,*) x, y, z
          !write (*,*) x, neg(y), z
          !write (*,*) neg(x), y, z

          e_x(2*x,2*y - 1,2*z - 1) = en1
          e_y(2*x - 1,2*y,2*z - 1) = en2
          e_x(2*x,neg(2*y - 1),2*z - 1) = en3
          e_y(neg(2*x - 1),2*y,2*z - 1) = en4
          acceptr = acceptr + 1
          u_tot_run = u_tot_run + delta_e

        end if ! end of Metropolis check

      else if (chooser.ge.(2.0/3.0)) then ! xz-plane plaquette

        eo1 = e_x(2*x,2*y - 1,2*z - 1)
        eo2 = e_z(2*x - 1,2*y - 1,2*z)
        eo3 = e_x(2*x,2*y - 1,neg(2*z - 1))
        eo4 = e_z(neg(2*x - 1),2*y - 1,2*z)

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

        if ((delta_e.lt.0.0).or.(exp((-1.0)*(beta)*delta_e).gt.rand())) then

          !write (*,*) "x-z plane rot:"
          !write (*,*) x, y, z
          !write (*,*) x, y, neg(z)
          !write (*,*) neg(x), y, z

          e_x(2*x,2*y - 1,2*z - 1) = en1
          e_z(2*x - 1,2*y - 1,2*z) = en2
          e_x(2*x,2*y - 1,neg(2*z - 1)) = en3
          e_z(neg(2*x - 1),2*y - 1,2*z) = en4
          acceptr = acceptr + 1
          u_tot_run = u_tot_run + delta_e

        end if ! end of Metropolis check

      else ! yz-plane plaquette

        eo1 = e_y(2*x - 1,2*y,2*z - 1)
        eo2 = e_z(2*x - 1,2*y - 1,2*z)
        eo3 = e_y(2*x - 1,2*y,neg(2*z - 1))
        eo4 = e_z(2*x - 1,neg(2*y - 1),2*z)

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

        if ((delta_e.lt.0.0).or.(exp((-1.0)*(beta)*delta_e).gt.rand())) then

          !write (*,*) "y-z plane rot:"
          !write (*,*) x, y, z
          !write (*,*) x, y, neg(z)
          !write (*,*) x, neg(y), z

          e_y(2*x - 1,2*y,2*z - 1) = en1
          e_z(2*x - 1,2*y - 1,2*z) = en2
          e_y(2*x - 1,2*y,neg(2*z - 1)) = en3
          e_z(2*x - 1,neg(2*y - 1),2*z) = en4
          acceptr = acceptr + 1
          u_tot_run = u_tot_run + delta_e

        end if ! end of Metropolis check

      end if ! end plane choice block

    end do ! end rotational
    !write (*,*) "utot after rot. = ",u_tot_run

    u_tot = 0.0
    do i = 1,L
      do j = 1, L
        do k = 1,L
        u_tot = u_tot + 0.5 * eps_0 * (e_x(2*i,2*j - 1,2*k - 1)**2 + e_y(2*i - 1,2*j,2*k - 1)**2 + e_z(2*i - 1,2*j - 1,2*k)**2)
        end do
      end do
    end do

    !u_diff = u_tot_run - u_tot
    !if (abs(u_diff).ge.0.00001) then
    !  write (*,*) "ROT: u_tot_run = ",u_tot_run," u_tot = ",u_tot," u_diff = ",u_diff
    !end if

    ! --- HARMONIC UPDATE ---
    ! e bar update
    if (glob.eq.1) then
      !write (*,*) "step ",n," start of harmonic update: ",u_tot_run,"total moves accepted ",acceptg

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

          if ((delta_e.lt.0).or.((exp((-1.0)*(beta)*delta_e).gt.rand())&
            .and.(exp((-1.0)*(beta)*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_x = ebar_x - ebar_inc
            acceptg = acceptg + 1

            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_x(2*j,2*k - 1,2*m - 1) = e_x(2*j,2*k - 1,2*m - 1) - e_inc
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

          if ((delta_e.lt.0).or.((exp((-1.0)*(beta)*delta_e).gt.rand())&
            .and.(exp((-1.0)*(beta)*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_x = ebar_x + ebar_inc
            acceptg = acceptg + 1

            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_x(2*j,2*k - 1,2*m - 1) = e_x(2*j,2*k - 1,2*m - 1) + e_inc
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

          if ((delta_e.lt.0).or.((exp((-1.0)*(beta)*delta_e).gt.rand())&
            .and.(exp((-1.0)*(beta)*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_y = ebar_y - ebar_inc
            acceptg = acceptg + 1

            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_y(2*j - 1,2*k,2*m - 1) = e_y(2*j - 1,2*k,2*m - 1) - e_inc
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

          if ((delta_e.lt.0).or.((exp((-1.0)*(beta)*delta_e).gt.rand())&
            .and.(exp((-1.0)*(beta)*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_y = ebar_y + ebar_inc
            acceptg = acceptg + 1

            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_y(2*j - 1,2*k,2*m - 1) = e_y(2*j - 1,2*k,2*m - 1) + e_inc
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

          if ((delta_e.lt.0).or.((exp((-1.0)*(beta)*delta_e).gt.rand())&
            .and.(exp((-1.0)*(beta)*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_z = ebar_z - ebar_inc
            acceptg = acceptg + 1

            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_z(2*j - 1,2*k - 1,2*m) = e_z(2*j - 1,2*k - 1,2*m) - e_inc
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

          if ((delta_e.lt.0).or.((exp((-1.0)*(beta)*delta_e).gt.rand())&
            .and.(exp((-1.0)*(beta)*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_z = ebar_z + ebar_inc
            acceptg = acceptg + 1

            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_z(2*j - 1,2*k - 1,2*m) = e_z(2*j - 1,2*k - 1,2*m) + e_inc
                end do
              end do
            end do
          end if ! end weird Metropolis block
        end if

        end do ! end global update loop

        u_tot_run = 0.0
        do j = 1,L
          do k = 1,L
            do m = 1,L
              u_tot_run = u_tot_run + 0.5 * eps_0 *&
                          (e_x(2*j,2*k - 1,2*m - 1)**2 + e_y(2*j - 1,2*k,2*m - 1)**2 + e_z(2*j - 1,2*k - 1,2*m)**2)
            end do
          end do
        end do

    end if ! end glob.eq.1 block

    ! --- END OF UPDATE BLOCKS ---

    u_tot = 0.0
    do i = 1,2*L
      do j = 1,2*L
        do k = 1,2*L
          u_tot = u_tot + 0.5 * eps_0 * (e_x(i,j,k)**2 + e_y(i,j,k)**2 + e_z(i,j,k)**2)
          ! copy the charge distribution and e field into fftw inputs
          fftw_ch_in(i,j,k) = v(i,j,k)

          ! test field structure factor
          if (configs.eq."NO") then
            fftw_e_in(i,j,k,1) = e_x(i,j,k)
            fftw_e_in(i,j,k,2) = e_y(i,j,k)
            fftw_e_in(i,j,k,3) = e_z(i,j,k)
          else if (configs.eq."AF") then
            fftw_e_in(i,j,k,1) = (-1.0)**((i+j)/2)
            fftw_e_in(i,j,k,2) = (-1.0)**((i+j)/2)
            fftw_e_in(i,j,k,3) = (1.0)
          else if (configs.eq."ST") then
            fftw_e_in(i,j,k,1) = (-1.0)**j
            fftw_e_in(i,j,k,2) = (-1.0)**k
            fftw_e_in(i,j,k,3) = (-1.0)**i
          else if (configs.eq."FE") then
            fftw_e_in(1,i,j,k) = 1.0
            fftw_e_in(2,i,j,k) = 1.0
            fftw_e_in(3,i,j,k) = 1.0
          end if
        end do
      end do
    end do

    !u_diff = u_tot_run - u_tot
    !if (abs(u_diff).ge.0.00001) then
    !  write (*,*) "HARM: u_tot_run = ",u_tot_run," u_tot = ",u_tot," u_diff = ",u_diff
    !end if

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

    ! do the transforms
    call fftw_execute_dft_r2c(plan_ch,fftw_ch_in,fftw_ch_out)
    call fftw_execute_dft_r2c(plan_e,fftw_e_in,fftw_e_out)

    ! according to the fftw docs this is the right normalisation
    fftw_ch_out = fftw_ch_out / sqrt(dble((2*L)**3))
    fftw_e_out = fftw_e_out / sqrt(dble(3*(2*L)**3))

    ! testing the longer expression - didn't work
    !do kx = -L/2, L/2
    !  do ky = -L/2, L/2
    !    do kz = -L/2, L/2

    !      do i=1,L
    !        do j=1,L
    !          do k=1,L

    !            do m=1,L
    !              do p=1,L
    !                do s=1,L

    !                  fe_fe(i,j,k,n) = (fftw_ex_in(i,j,k) * fftw_ex_in(m,p,s) &
    !                  + fftw_ey_in(i,j,k) * fftw_ey_in(m,p,s) &
    !                  + fftw_ez_in(i,j,k) * fftw_ez_in(m,p,s)) &
    !                  * (exp((-1)*imag*(2*pi/L*lambda) * &
    !                  (kx*(i-m) + ky*(j-p) + kz*(s-k))))

    !                end do ! s
    !              end do ! p
    !            end do ! m

    !          end do ! k
    !        end do ! j
    !      end do ! i

    !    end do ! kz
    !  end do ! ky
    !end do ! kx

    do i=1,L + 1
      do j=1,2*L
        do k=1,2*L

          ! products, lovely
          fe_fe(i,j,k,n) = (fftw_e_in(i,1,j,k)*fftw_e_out(i,1,j,k)*conjg(fftw_e_out(i,1,j,k))&
                           + fftw_e_in(i,2,j,k)*fftw_e_out(i,2,j,k)*conjg(fftw_e_out(i,2,j,k))&
                           + fftw_e_in(i,3,j,k)*fftw_e_out(i,3,j,k)*conjg(fftw_e_out(i,3,j,k)))

          ch_ch(i,j,k,n) = v(i,j,k)*fftw_ch_out(i,j,k)*conjg(fftw_ch_out(i,j,k))

          if (n.eq.1) then
            write (*,'(I2.1,I2.1,I2.1,F9.6,F9.6,F9.6)')&
              i,j,k,fftw_e_out(i,2,j,k),fe_fe(i,j,k,n)
          end if

          !if (n.eq.1) then
          !  write (*,*) fftw_e_out(i,1,j,k),&
          !              fftw_e_out(i,2,j,k),fftw_e_out(i,3,j,k)
          !end if

        end do ! k
      end do ! j
    end do ! i

    !! Fourier transform
    !do kx = -L/2, L/2
    !  do ky = -L/2, L/2
    !    do kz = -L/2, L/2

    !      ! for array indices
    !      i = kx + 1 + L/2
    !      j = ky + 1 + L/2
    !      k = kz + 1 + L/2

    !      e_kx(i,j,k) = 0.0
    !      e_ky(i,j,k) = 0.0
    !      e_kz(i,j,k) = 0.0

    !      ! vec(q) = vec(0) term:
    !      if (i.eq.0.and.j.eq.0.and.k.eq.0) then
    !        write (*,*) "q=0 term; need to work this out"
    !        write (*,*) "e_kx(",i,j,k,") = ",e_kx(i,j,k)
    !        CYCLE
    !      end if

    !      ! m,p,s are the real space coordinates
    !      do m = 1,L
    !        do p = 1,L
    !          do s = 1,L

    !            kdotx = ((-1)*imag*2*pi*(((m-1)*kx/(L*lambda)) + &
    !                    ((p-1)*ky/(L*lambda)) + &
    !                    ((s-1)*kz/(L*lambda))))

    !            e_kx(i,j,k) = e_kx(i,j,k) + e**(kdotx)*e_x(m,p,s)
    !            e_ky(i,j,k) = e_ky(i,j,k) + e**(kdotx)*e_y(m,p,s)
    !            e_kz(i,j,k) = e_kz(i,j,k) + e**(kdotx)*e_z(m,p,s)

    !            if (v(m,p,s).eq.1) then ! calculate <++>!

    !              ! FT of charge distribution
    !              kdotx = ((-1)*imag*2*pi*(((m-1)*kx/(L*lambda)) + &
    !                      ((p-1)*ky/(L*lambda)) + &
    !                      ((s-1)*kz/(L*lambda))))
    !              rho_k(i,j,k) = rho_k(i,j,k) + v(m,p,s) * e**(kdotx)

    !            end if

    !          end do
    !        end do
    !      end do

    !      ! normalise, idiot
    !      e_kx(i,j,k) = e_kx(i,j,k) / sqrt(float(L**3))
    !      e_ky(i,j,k) = e_ky(i,j,k) / sqrt(float(L**3))
    !      e_kz(i,j,k) = e_kz(i,j,k) / sqrt(float(L**3))
    !      rho_k(i,j,k) = rho_k(i,j,k) / sqrt(float(L**3))

    !      ! second part is weirdly addressed bc of loop structure
    !      ! ((-1) * k_mu) + 1 + L/2 --> (-\vec{k}) essentially
    !      ch_ch(i,j,k,n) = ch_ch(i,j,k,n) + (rho_k(i,j,k) *&
    !        rho_k((-1)*kx + 1 + L/2, (-1)*ky + 1 + L/2, (-1)*kz + 1 + L/2))

    !      ! loop indices need changing to do this; i dimension is smaller
    !      ch_ch(i,j,k,n) = fftw_output(i,j,k)*conjg(fftw_output(i,j,k))

    !      ! there is probably a better way of doing this
    !      ! folding e_mu into one three-vector would make it easier
    !      fe_fe(1,1,i,j,k,n) = fe_fe(1,1,i,j,k,n) + (e_kx(i,j,k)&
    !        * e_kx((-1)*kx + 1 + L/2, (-1)*ky + 1 + L/2, (-1)*kz + 1 + L/2))
    !      fe_fe(1,2,i,j,k,n) = fe_fe(1,2,i,j,k,n) + (e_kx(i,j,k)&
    !        * e_ky((-1)*kx + 1 + L/2, (-1)*ky + 1 + L/2, (-1)*kz + 1 + L/2))
    !      fe_fe(1,3,i,j,k,n) = fe_fe(1,3,i,j,k,n) + (e_kx(i,j,k)&
    !        * e_kz((-1)*kx + 1 + L/2, (-1)*ky + 1 + L/2, (-1)*kz + 1 + L/2))
    !      fe_fe(2,1,i,j,k,n) = fe_fe(2,1,i,j,k,n) + (e_ky(i,j,k)&
    !        * e_kx((-1)*kx + 1 + L/2, (-1)*ky + 1 + L/2, (-1)*kz + 1 + L/2))
    !      fe_fe(2,2,i,j,k,n) = fe_fe(2,2,i,j,k,n) + (e_ky(i,j,k)&
    !        * e_ky((-1)*kx + 1 + L/2, (-1)*ky + 1 + L/2, (-1)*kz + 1 + L/2))
    !      fe_fe(2,3,i,j,k,n) = fe_fe(2,3,i,j,k,n) + (e_ky(i,j,k)&
    !        * e_kz((-1)*kx + 1 + L/2, (-1)*ky + 1 + L/2, (-1)*kz + 1 + L/2))
    !      fe_fe(3,1,i,j,k,n) = fe_fe(3,1,i,j,k,n) + (e_kz(i,j,k)&
    !        * e_kx((-1)*kx + 1 + L/2, (-1)*ky + 1 + L/2, (-1)*kz + 1 + L/2))
    !      fe_fe(3,2,i,j,k,n) = fe_fe(3,2,i,j,k,n) + (e_kz(i,j,k)&
    !        * e_ky((-1)*kx + 1 + L/2, (-1)*ky + 1 + L/2, (-1)*kz + 1 + L/2))
    !      fe_fe(3,3,i,j,k,n) = fe_fe(3,3,i,j,k,n) + (e_kz(i,j,k)&
    !        * e_kz((-1)*kx + 1 + L/2, (-1)*ky + 1 + L/2, (-1)*kz + 1 + L/2))

    !    end do
    !  end do
    !end do

  end do ! end iteration loop

  write(*,*)
  write(*,*) " --- end: charge positions ---"

  totq = 0
  do j = 1,L
    do k = 1,L
      do m = 1,L
        totq = totq + abs(v(2*j - 1,2*k - 1,2*m - 1))
        if (v(2*j - 1,2*k - 1,2*m - 1).ne.0) then
          write (*,*) 2*j - 1, 2*k - 1, 2*m - 1, v(2*j - 1,2*k - 1,2*m - 1)
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
  write (*,*) 'ebar moves = ',acceptg,float(acceptg) / (iterations * L**3 * g_ratio * 3)

end subroutine upcan
