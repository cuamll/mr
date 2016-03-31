! maggs - rossetto algorithm on cubic lattice
! i'll write proper docs at some point
! callum gray, UCL / ENS Lyon, 2016. do what you like with it
program mr

  use common
  use linear_solver
  use io
  implicit none
  integer :: i, j, k, row, col
  integer, dimension(:,:), allocatable :: v_temp

  ebar_x  =  0.0
  ebar_y  =  0.0
  ebar_z  =  0.0

  call read_input

  allocate(v(L,L,L))
  allocate(pos(L))
  allocate(neg(L))
  allocate(mnphi_x(L,L,L))
  allocate(mnphi_y(L,L,L))
  allocate(mnphi_z(L,L,L))
  allocate(e_rot_x(L,L,L))
  allocate(e_rot_y(L,L,L))
  allocate(e_rot_z(L,L,L))
  allocate(e_x(L,L,L))
  allocate(e_y(L,L,L))
  allocate(e_z(L,L,L))
  allocate(lgf(L,L,L,L,L,L))
  allocate(v_temp(L**2,L))

  write(*,*) 'L = ',L

  call PBCs

  open(unit = 2, file = lattfile)
  read(2,*)((v_temp(row,col),col = 1,L),row = 1,L**2)

  ! 2,3,1 makes x,y,z correspond with what you expect from the file
  ! doesn't actually make any difference so long as you're consistent
  v  =  reshape(v_temp, (/ L,L,L /), ORDER  =  (/ 2,3,1 /))

  deallocate(v_temp) ! we don't need it anymore
  close(2)

  call randinit(seed)

  call linsol

  ! set e_x to irrotational - temporary solution
  do i = 1,L
    do j = 1,L
      do k = 1,L
        e_x(i,j,k) = mnphi_x(i,j,k)
        e_y(i,j,k) = mnphi_y(i,j,k)
        e_z(i,j,k) = mnphi_z(i,j,k)

        e_rot_x(i,j,k) = 0.0
        e_rot_y(i,j,k) = 0.0
        e_rot_z(i,j,k) = 0.0
      end do
    end do
  end do

  call update()

  stop

end program mr

subroutine update()
  use common
  implicit none
  integer :: i,j,k,m,n,x,y,z,v1,v2,glob
  real :: chooser
  real*8 :: old_e,new_e,delta_e,eo1,eo2,eo3,eo4,en1,en2,en3,en4,utotal
  real*8 :: hop_inc,delta,g_thr

  hop_inc = q / (eps_0 * lambda**2)
  g_thr = pi/float(L) ! if |ebar_i| > this, charges are winding
  glob = 0

  do n = 1,iterations

    ! charge hop updates
    accepth = 0

    do i = 1,L**3
      ! one sweep is L**3 attempts
      ! pick a random (x,y,z)
      x = int(rand() * L) + 1
      y = int(rand() * L) + 1
      z = int(rand() * L) + 1

      chooser = rand()
      if (chooser.lt.(1.0 / 3.0)) then ! x component

        eo1 = e_x(x,y,z)
        chooser = rand()

        if (chooser.lt.(0.5)) then ! decrease bond

          en1 = eo1 - hop_inc
          old_e = 0.5 * eo1**2
          new_e = 0.5 * en1**2
          delta_e = new_e - old_e
          v1 = v(x,y,z) + 1
          v2 = v(neg(x),y,z) - 1

          if ((abs(v1).le.1).and.(abs(v2).le.1)) then ! enforce |v| <= 1
            if ((delta_e.lt.0).or.(exp((-beta)*delta_e).gt.rand())) then

              e_x(x,y,z) = en1
              v(x,y,z) = v1
              v(neg(x),y,z) = v2
              ebar_x = ebar_x - (hop_inc / L**3) ! track evolution
              accepth = accepth + 1

              if (ebar_x.gt.g_thr.or.(ebar_x.le.((-1)*g_thr))&
                .and.(glob.eq.0)) then
                glob = 1
              end if ! end global winding check

            end if ! end Metropolis check
          end if ! end |v| < =  1 check

        else ! increase bond

          en1 = eo1 + hop_inc
          old_e = 0.5 * eo1**2
          new_e = 0.5 * en1**2
          delta_e = new_e - old_e
          v1 = v(x,y,z) - 1
          v2 = v(neg(x),y,z) + 1

          if ((abs(v1).le.1).and.(abs(v2).le.1)) then
            if ((delta_e.lt.0).or.(exp((-beta)*delta_e).gt.rand())) then

              e_x(x,y,z) = en1
              v(x,y,z) = v1
              v(neg(x),y,z) = v2
              ebar_x = ebar_x + (hop_inc / L**3)
              accepth = accepth + 1

              if (ebar_x.gt.g_thr.or.(ebar_x.le.((-1)*g_thr))&
                .and.(glob.eq.0)) then
                glob = 1
              end if ! end global winding check

            end if ! end Metropolis check
          end if ! end |v| < =  1 check

        end if ! end of x block

      else if (chooser.ge.(2.0 / 3.0)) then ! y component

        eo1 = e_y(x,y,z)
        chooser = rand()

        if (chooser.lt.(0.5)) then ! decrease bond

          en1 = eo1 - hop_inc
          old_e = 0.5 * eo1**2
          new_e = 0.5 * en1**2
          delta_e = new_e - old_e
          v1 = v(x,y,z) + 1
          v2 = v(x,neg(y),z) - 1

          if ((abs(v1).le.1).and.(abs(v2).le.1)) then ! enforce |v| <= 1
            if ((delta_e.lt.0).or.(exp((-beta)*delta_e).gt.rand())) then

              e_y(x,y,z) = en1
              v(x,y,z) = v1
              v(x,neg(y),z) = v2
              ebar_y = ebar_y - (hop_inc / L**3) ! track evolution
              accepth = accepth + 1

              if (ebar_x.gt.g_thr.or.(ebar_x.le.((-1)*g_thr))&
                .and.(glob.eq.0)) then
                glob = 1
              end if ! end global winding check

            end if ! end Metropolis check
          end if ! end |v| < =  1 check

        else ! increase bond

          en1 = eo1 + hop_inc
          old_e = 0.5 * eo1**2
          new_e = 0.5 * en1**2
          delta_e = new_e - old_e
          v1 = v(x,y,z) - 1
          v2 = v(x,neg(y),z) + 1

          if ((abs(v1).le.1).and.(abs(v2).le.1)) then
            if ((delta_e.lt.0).or.(exp((-beta)*delta_e).gt.rand())) then

              e_y(x,y,z) = en1
              v(x,y,z) = v1
              v(x,neg(y),z) = v2
              ebar_y = ebar_y + (hop_inc / L**3)
              accepth = accepth + 1

              if (ebar_x.gt.g_thr.or.(ebar_x.le.((-1)*g_thr))&
                .and.(glob.eq.0)) then
                glob = 1
              end if ! end global winding check

            end if ! end Metropolis check
          end if ! end |v| < =  1 check

        end if ! end of y block

      else ! z component

        eo1 = e_z(x,y,z)
        chooser = rand()

        if (chooser.lt.(0.5)) then ! decrease bond

          en1 = eo1 - hop_inc
          old_e = 0.5 * eo1**2
          new_e = 0.5 * en1**2
          delta_e = new_e - old_e
          v1 = v(x,y,z) + 1
          v2 = v(x,y,neg(z)) - 1

          if ((abs(v1).le.1).and.(abs(v2).le.1)) then ! enforce |v| <= 1
            if ((delta_e.lt.0).or.(exp((-beta)*delta_e).gt.rand())) then

              e_z(x,y,z) = en1
              v(x,y,z) = v1
              v(x,y,neg(z)) = v2
              ebar_z = ebar_z - (hop_inc / L**3) ! track evolution
              accepth = accepth + 1

              if (ebar_x.gt.g_thr.or.(ebar_x.le.((-1)*g_thr))&
                .and.(glob.eq.0)) then
                glob = 1
              end if ! end global winding check

            end if ! end Metropolis check
          end if ! end |v| < =  1 check

        else ! increase bond

          en1 = eo1 + hop_inc
          old_e = 0.5 * eo1**2
          new_e = 0.5 * en1**2
          delta_e = new_e - old_e
          v1 = v(x,y,z) - 1
          v2 = v(x,y,neg(z)) + 1

          if ((abs(v1).le.1).and.(abs(v2).le.1)) then
            if ((delta_e.lt.0).or.(exp((-beta)*delta_e).gt.rand())) then

              e_z(x,y,z) = en1
              v(x,y,z) = v1
              v(x,y,neg(z)) = v2
              ebar_z = ebar_z + (hop_inc / L**3)
              accepth = accepth + 1

              if (ebar_x.gt.g_thr.or.(ebar_x.le.((-1)*g_thr))&
                .and.(glob.eq.0)) then
                glob = 1
              end if ! end global winding check

            end if ! end Metropolis check
          end if ! end |v| < =  1 check

        end if ! end of z block
      end if ! end of x/y/z block

    end do ! end of charge hop updates????

    ! plaquette rot. update
    acceptr=0

    ! NOTE TO SELF: might need a + 1 next to that int cast
    do i = 1,int(L**3 * rot_ratio)

      x = int(rand() * L) + 1
      y = int(rand() * L) + 1
      z = int(rand() * L) + 1
      ! NOTE TO SELF : next line needs changing
      ! it's for the field increment, not an if statement
      delta = rand()

      chooser=rand()
      if (chooser.lt.(1.0/3.0)) then ! xy-plane plaquette

        eo1 = e_x(x,y,z)
        eo2 = e_y(x,y,z)
        eo3 = e_x(x,neg(y),z)
        eo4 = e_y(neg(x),y,z)

        en1 = eo1 + delta
        en2 = eo2 - delta
        en3 = eo3 - delta
        en4 = eo4 + delta

        old_e = 0.5*(eo1**2 + eo2**2 + eo3**2 + eo4**2)
        new_e = 0.5*(en1**2 + en2**2 + en3**2 + en4**2)
        delta_e = new_e - old_e

        if ((delta_e.lt.0).or.(exp((-beta)*delta_e).gt.rand())) then

          e_x(x,y,z) = en1
          e_y(x,y,z) = en2
          e_x(x,neg(y),z) = en3
          e_y(neg(x),y,z) = en4
          acceptr = acceptr + 1

        end if ! end of Metropolis check

      else if (chooser.ge.(2.0/3.0)) then ! xz-plane plaquette

        eo1 = e_x(x,y,z)
        eo2 = e_z(x,y,z)
        eo3 = e_x(x,y,neg(z))
        eo4 = e_z(neg(x),y,z)

        en1 = eo1 + delta
        en2 = eo2 - delta
        en3 = eo3 - delta
        en4 = eo4 + delta

        old_e = 0.5*(eo1**2 + eo2**2 + eo3**2 + eo4**2)
        new_e = 0.5*(en1**2 + en2**2 + en3**2 + en4**2)
        delta_e = new_e - old_e

        if ((delta_e.lt.0).or.(exp((-beta)*delta_e).gt.rand())) then

          e_x(x,y,z) = en1
          e_z(x,y,z) = en2
          e_x(x,y,neg(z)) = en3
          e_z(neg(x),y,z) = en4
          acceptr = acceptr + 1

        end if ! end of Metropolis check

      else ! yz-plane plaquette

        eo1 = e_y(x,y,z)
        eo2 = e_z(x,y,z)
        eo3 = e_y(x,y,neg(z))
        eo4 = e_z(x,neg(y),z)

        en1 = eo1 + delta
        en2 = eo2 - delta
        en3 = eo3 - delta
        en4 = eo4 + delta

        old_e = 0.5*(eo1**2 + eo2**2 + eo3**2 + eo4**2)
        new_e = 0.5*(en1**2 + en2**2 + en3**2 + en4**2)
        delta_e = new_e - old_e

        if ((delta_e.lt.0).or.(exp((-beta)*delta_e).gt.rand())) then

          e_y(x,y,z) = en1
          e_z(x,y,z) = en2
          e_y(x,y,neg(z)) = en3
          e_z(x,neg(y),z) = en4
          acceptr = acceptr + 1

        end if ! end of Metropolis check

      end if ! end plane choice block

    end do ! end rotational

    ! e bar update
    acceptg = 0
    if (glob.eq.1) then

      ! NOTE TO SELF: again this int cast prob needs changing
      do i = 1,int(L**3 * g_ratio)
        ! x-component
        chooser = rand()

        if (chooser.lt.0.5) then
          ! NOTE TO SELF - check this little fucker
          delta_e = 2 * pi * (pi - float(L) * ebar_x)

          if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_x = ebar_x - 2 * g_thr
            acceptg = acceptg + 1
            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_x(j,k,m) = e_x(j,k,m) - 2 * g_thr
                end do
              end do
            end do
          end if ! end weird Metropolis block

        else
          ! NOTE TO SELF - check this little fucker
          delta_e = 2 * pi * (pi + float(L) * ebar_x)

          if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_x = ebar_x + 2 * g_thr
            acceptg = acceptg + 1

            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_x(j,k,m) = e_x(j,k,m) + 2 * g_thr
                end do
              end do
            end do
          end if ! end weird Metropolis block
        end if

        ! y-component
        chooser = rand()

        if (chooser.lt.0.5) then
          ! NOTE TO SELF - check this little fucker
          delta_e = 2 * pi * (pi - float(L) * ebar_y)

          if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_y = ebar_y - 2 * g_thr
            acceptg = acceptg + 1
            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_y(j,k,m) = e_y(j,k,m) - 2 * g_thr
                end do
              end do
            end do
          end if ! end weird Metropolis block

        else
          ! NOTE TO SELF - check this little fucker
          delta_e = 2 * pi * (pi + float(L) * ebar_y)

          if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_y = ebar_y + 2 * g_thr
            acceptg = acceptg + 1

            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_y(j,k,m) = e_y(j,k,m) + 2 * g_thr
                end do
              end do
            end do
          end if ! end weird Metropolis block
        end if

        ! z-component
        chooser = rand()

        if (chooser.lt.0.5) then
          ! NOTE TO SELF - check this little fucker
          delta_e = 2 * pi * (pi - float(L) * ebar_z)

          if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_z = ebar_z - 2 * g_thr
            acceptg = acceptg + 1
            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_z(j,k,m) = e_z(j,k,m) - 2 * g_thr
                end do
              end do
            end do
          end if ! end weird Metropolis block

        else
          ! NOTE TO SELF - check this little fucker
          delta_e = 2 * pi * (pi + float(L) * ebar_z)

          if ((delta_e.lt.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then
            ! this block is basically stolen from Michael
            ! not sure what's happening here tbh
            ebar_z = ebar_z + 2 * g_thr
            acceptg = acceptg + 1

            do j = 1,L
              do k = 1,L
                do m = 1,L
                  e_z(j,k,m) = e_z(j,k,m) + 2 * g_thr
                end do
              end do
            end do
          end if ! end weird Metropolis block
        end if

        end do ! end global update loop

    end if ! end glob.eq.1 block

  end do ! end n_iteration block

  do j = 1,L
    do k = 1,L
      do m = 1,L
        utotal = utotal + e_x(j,k,m)**2 + e_y(j,k,m)**2 + e_z(j,k,m)**2
      end do
    end do
  end do
  write(*,*) 'U_tot = ',utotal

end subroutine update
