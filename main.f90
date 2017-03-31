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

        tot_q = tot_q + abs(v(i,j,k))

      end do
    end do
  end do

  if (tot_q.ne.0) then
    call linalg
  end if

  !do i = 1,L
  !  do j = 1,L
  !    do k = 1,L

  !    ! does it make sense to do this?
  !      u_rot = u_rot + 0.5 * (e_x(i,j,k) - mnphi_x(i,j,k))**2 +&
  !      (e_y(i,j,k) - mnphi_y(i,j,k))**2 +&
  !      (e_z(i,j,k) - mnphi_z(i,j,k))**2

  !    end do
  !  end do
  !end do

  !write (*,*) 'u_rot. = ',u_rot

  write(*,*)

  call linalg

  call write_output

  call deallocations

  stop

end program mr

! canonical ensemble - MC / Metropolis updates
subroutine upcan()
  use common
  use linear_solver
  implicit none
  integer :: x,y,z,n,charge,glob,i,j,k,m,p,s,pm1,kx,ky,kz
  real*8 :: eo1,eo2,eo3,eo4,en1,en2,en3,en4
  real*8 :: u_tot_run,avg_e,avg_e2,one,norm_k
  real*8 :: dot,dot_avg ! probably won't need the avg
  real*8 :: hop_inc, old_e, new_e, delta_e, totq, g_thr
  real :: chooser, delta
  complex*16 :: imag, kdotx

  glob = 0
  totq = 0
  u_tot_run = 0.0
  old_e = 0.0
  new_e = 0.0
  g_thr = 1 / float(L)
  accepth = 0
  acceptr=0
  acceptg=0
  one=1.0
  imag = (0.0,1.0)

  write(*,*)
  write(*,*) " --- start: charge positions ---"

  call linalg
  u_tot_run = lapack_energy
  do i = 1,L
    do j = 1,L
      do k = 1,L

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

    ! charge hop sweep
    do i = 1, int(L**3 * hop_ratio)

      ! pick a random site
      x = int(rand() * L) + 1
      y = int(rand() * L) + 1
      z = int(rand() * L) + 1

      ! pick a non-zero charge - we want to keep it canonical
      if (v(x,y,z).ne.0) then

        charge = v(x,y,z)

        !choose a direction to move the charge:
        chooser = int(rand() * 6)

        if (chooser.eq.0) then

          ! positive x
          if ((abs(v(x,y,z)).le.1).and.(v(pos(x),y,z).eq.0)) then

            old_e = u_tot_run
            v(x,y,z) = 0
            v(pos(x),y,z) = charge

            call linalg
            ! after calling linalg, the new fields are in e_mu_lapack
            new_e = lapack_energy
            delta_e = new_e - old_e
            !write (*,*) old_e,new_e,delta_e

            if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

              ! move accepted, set the fields/energy to lapack ones
              accepth = accepth + 1
              u_tot_run = lapack_energy
              e_x = e_x_lapack
              e_y = e_y_lapack
              e_z = e_z_lapack

            else

              ! don't change the fields or energy, move the charge back
              v(x,y,z) = charge
              v(pos(x),y,z) = 0
              u_tot_run = old_e

            end if
          end if

        else if (chooser.eq.1) then

          ! negative x
          if ((abs(v(x,y,z)).le.1).and.(v(neg(x),y,z).eq.0)) then

            old_e = u_tot_run
            v(x,y,z) = 0
            v(neg(x),y,z) = charge

            call linalg
            ! after calling linalg, the new fields are in e_mu_lapack
            new_e = lapack_energy
            delta_e = new_e - old_e

            if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

              ! move accepted, set the fields/energy to lapack ones
              accepth = accepth + 1
              u_tot_run = lapack_energy
              e_x = e_x_lapack
              e_y = e_y_lapack
              e_z = e_z_lapack

            else

              ! don't change the fields or energy, move the charge back
              v(x,y,z) = charge
              v(neg(x),y,z) = 0
              u_tot_run = old_e

            end if
          end if

        else if (chooser.eq.2) then

          ! positive y
          if ((abs(v(x,y,z)).le.1).and.(v(x,pos(y),z).eq.0)) then

            old_e = u_tot_run
            v(x,y,z) = 0
            v(x,pos(y),z) = charge

            call linalg
            ! after calling linalg, the new fields are in e_mu_lapack
            new_e = lapack_energy
            delta_e = new_e - old_e

            if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

              ! move accepted, set the fields/energy to lapack ones
              accepth = accepth + 1
              u_tot_run = lapack_energy
              e_x = e_x_lapack
              e_y = e_y_lapack
              e_z = e_z_lapack

            else

              ! don't change the fields or energy, move the charge back
              v(x,y,z) = charge
              v(x,pos(y),z) = 0
              u_tot_run = old_e

            end if
          end if

        else if (chooser.eq.3) then

          ! negative y
          if ((abs(v(x,y,z)).le.1).and.(v(x,neg(y),z).eq.0)) then

            old_e = u_tot_run
            v(x,y,z) = 0
            v(x,neg(y),z) = charge

            call linalg
            ! after calling linalg, the new fields are in e_mu_lapack
            new_e = lapack_energy
            delta_e = new_e - old_e

            if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

              ! move accepted, set the fields/energy to lapack ones
              accepth = accepth + 1
              u_tot_run = lapack_energy
              e_x = e_x_lapack
              e_y = e_y_lapack
              e_z = e_z_lapack

            else

              ! don't change the fields or energy, move the charge back
              v(x,y,z) = charge
              v(x,neg(y),z) = 0
              u_tot_run = old_e

            end if
          end if

        else if (chooser.eq.4) then

          ! positive z
          if ((abs(v(x,y,z)).le.1).and.(v(x,y,pos(z)).eq.0)) then

            old_e = u_tot_run
            v(x,y,z) = 0
            v(x,y,pos(z)) = charge

            call linalg
            ! after calling linalg, the new fields are in e_mu_lapack
            new_e = lapack_energy
            delta_e = new_e - old_e

            if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

              ! move accepted, set the fields/energy to lapack ones
              accepth = accepth + 1
              u_tot_run = lapack_energy
              e_x = e_x_lapack
              e_y = e_y_lapack
              e_z = e_z_lapack

            else

              ! don't change the fields or energy, move the charge back
              v(x,y,z) = charge
              v(x,y,pos(z)) = 0
              u_tot_run = old_e

            end if
          end if

        else if (chooser.eq.5) then

          ! negative z
          if ((abs(v(x,y,z)).le.1).and.(v(x,y,neg(z)).eq.0)) then

            old_e = u_tot_run
            v(x,y,z) = 0
            v(x,y,neg(z)) = charge

            call linalg
            ! after calling linalg, the new fields are in e_mu_lapack
            new_e = lapack_energy
            delta_e = new_e - old_e

            if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

              ! move accepted, set the fields/energy to lapack ones
              accepth = accepth + 1
              u_tot_run = lapack_energy
              e_x = e_x_lapack
              e_y = e_y_lapack
              e_z = e_z_lapack

            else

              ! don't change the fields or energy, move the charge back
              v(x,y,z) = charge
              v(x,y,neg(z)) = 0
              u_tot_run = old_e

            end if
          end if

        end if


      end if ! end of v(x,y,z).ne.0 block

    end do ! end charge hop sweep

  ! replace with u_tot_run maybe?
  energy(n + 1) = u_tot_run
  sq_energy(n + 1) = u_tot_run**2

  avg_e = 0.0
  avg_e2 = 0.0
  do m=1,n+1
    avg_e = avg_e + energy(m)
    avg_e2 = avg_e2 + sq_energy(m)
  end do

  avg_e = avg_e / (n + 1)
  avg_e2 = avg_e2 / (n + 1)

    ! --- FOURIER TRANSFORMS ---
    do kx = (-1*L/2)*bz, (L/2)*bz
      do ky = (-1*L/2)*bz, (L/2)*bz
        do kz = (-1*L/2)*bz, (L/2)*bz

          ! for array indices
          i = kx + 1 + bz*(L/2)
          j = ky + 1 + bz*(L/2)
          k = kz + 1 + bz*(L/2)

          if (kx.eq.0.and.ky.eq.0.and.kz.eq.0) then
            norm_k = 0.0
          else
            norm_k = 1.0/(((2*pi/(L*lambda))**2)*dble(kx**2 + ky**2 + kz**2))
          end if

          !if (n.eq.1) then
          !  write (*,*) kx,ky,kz,(kx*kx + ky*ky + kz*kz)*(norm_k)
          !end if

          e_kx(i,j,k) = 0.0
          e_ky(i,j,k) = 0.0
          e_kz(i,j,k) = 0.0

          ! m,p,s are the real space coordinates
          do m = 1,L
            do p = 1,L
              do s = 1,L

                ! different offsets for x,y,z
                kdotx = ((-1)*imag*(2*pi/(L*lambda))*((m-(1.0/2))*kx + &
                        ((p-1)*ky) + ((s-1)*kz)))

                e_kx(i,j,k) = e_kx(i,j,k) + exp(kdotx)*e_x(m,p,s)

                kdotx = ((-1)*imag*(2*pi/(L*lambda))*((m-1)*kx + &
                        ((p-(1.0/2))*ky) + ((s-1)*kz)))

                e_ky(i,j,k) = e_ky(i,j,k) + exp(kdotx)*e_y(m,p,s)

                kdotx = ((-1)*imag*(2*pi/(L*lambda))*((m-1)*kx + &
                        ((p-1)*ky) + ((s-(1.0/2))*kz)))

                e_kz(i,j,k) = e_kz(i,j,k) + exp(kdotx)*e_z(m,p,s)

                if (v(m,p,s).ne.0) then ! calculate average of <++>,<--><+->!

                  ! FT of charge distribution
                  kdotx = ((-1)*imag*2*pi*(((m-1)*kx/(L*lambda)) + &
                          ((p-1)*ky/(L*lambda)) + &
                          ((s-1)*kz/(L*lambda))))

                  if (v(m,p,s).eq.-1) then
                    rho_k_m(i,j,k) = rho_k_m(i,j,k) + v(m,p,s) * exp(kdotx)
                  end if

                  if (v(m,p,s).eq.1) then
                    rho_k_p(i,j,k) = rho_k_p(i,j,k) + v(m,p,s)*exp(kdotx)
                  end if

                end if

              end do
            end do
          end do

          ! normalise, idiot
          e_kx(i,j,k) = e_kx(i,j,k) /(float(L**3))
          e_ky(i,j,k) = e_ky(i,j,k) /(float(L**3))
          e_kz(i,j,k) = e_kz(i,j,k) /(float(L**3))
          rho_k_m(i,j,k) = rho_k_m(i,j,k) /(float(L**3))
          rho_k_p(i,j,k) = rho_k_p(i,j,k) /(float(L**3))

          e_kx_t(i,j,k,n) = e_kx(i,j,k)
          rho_k_m_t(i,j,k,n) = rho_k_m(i,j,k)
          rho_k_p_t(i,j,k,n) = rho_k_p(i,j,k)

          ch_ch(i,j,k,n) = (rho_k_p(i,j,k)*conjg(rho_k_m(i,j,k)))

          if (kx.eq.((L*lambda/2)).and.ky.eq.(-1*L*lambda/2).and.kz.eq.0) then
            if (n.eq.1) then
              write(*,*) "kx,ky,rho_k_m,rho_k_p,ch_ch(this step),e_kx"
            else
            write (*,*)
            write (*,'(F6.3,F6.3,F12.7,F12.7,F12.7,F12.7,F12.7,F12.7,F12.7)')&
              kx*2*pi/(L*lambda),ky*2*pi/(L*lambda),&
              rho_k_m(i,j,k),rho_k_p(i,j,k),ch_ch(i,j,k,n),e_kx(i,j,k)
            end if
          end if

          !ch_ch(i,j,k,n) = ((rho_k_p(i,j,k) + rho_k_m(i,j,k))&
          !                 *conjg(rho_k_p(i,j,k) + rho_k_m(i,j,k)))
          !ch_ch_pp(i,j,k,n) = (rho_k_pp(i,j,k)*conjg(rho_k_pp(i,j,k)))

          fe_fe(i,j,k,n) = (e_kx(i,j,k)*conjg(e_kx(i,j,k)))

          s_ab_n(1,1,i,j,k,n) = e_kx(i,j,k)*conjg(e_kx(i,j,k))
          s_ab_n(1,2,i,j,k,n) = e_kx(i,j,k)*conjg(e_ky(i,j,k))
          s_ab_n(1,3,i,j,k,n) = e_kx(i,j,k)*conjg(e_kz(i,j,k))
          s_ab_n(2,1,i,j,k,n) = e_ky(i,j,k)*conjg(e_kx(i,j,k))
          s_ab_n(2,2,i,j,k,n) = e_ky(i,j,k)*conjg(e_ky(i,j,k))
          s_ab_n(2,3,i,j,k,n) = e_ky(i,j,k)*conjg(e_kz(i,j,k))
          s_ab_n(3,1,i,j,k,n) = e_kz(i,j,k)*conjg(e_kx(i,j,k))
          s_ab_n(3,2,i,j,k,n) = e_kz(i,j,k)*conjg(e_ky(i,j,k))
          s_ab_n(3,3,i,j,k,n) = e_kz(i,j,k)*conjg(e_kz(i,j,k))

        end do
      end do
    end do

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

end subroutine upcan
