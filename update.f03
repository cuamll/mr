module update
  use common
  ! not sure if linear solver is needed
  use linear_solver
  implicit none

  contains

    subroutine hop(attempts)
      implicit none
      integer, intent(in) :: attempts
      integer :: i,j,mu,charge
      integer, dimension(2) :: site
      real(kind=8) :: eo, en, old_e, new_e, delta_e, increment

      ! --- CHARGE HOP UPDATE ---

      do i = 1, attempts

        ! NOTE TO SELF: this whole procedure assumes
        ! single-valued charges only.

        ! pick a random site
        site = (/ int(rand() * L) + 1,&
                  &int(rand() * L) + 1 /)

        mu = floor(2*rand())+1

        ! check we're not doing anything weird with multiple charges
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

            ! get negative in the mu direction
            site(mu) = neg(site(mu))
            eo = e_field(mu,site(1),site(2))
            en = eo + increment
            old_e = 0.5 * eps_0 * eo**2
            new_e = 0.5 * eps_0 * en**2
            delta_e = new_e - old_e

            if (v(site(1),site(2)).eq.0) then
              if ((delta_e.lt.0.0).or.(exp((-beta)*delta_e).gt.rand())) then

                v(site(1),site(2)) = charge
                e_field(mu,site(1),site(2)) = en
                ebar(mu) = ebar(mu) + (increment / L**2)

                ! go back to the original site and set the charge to 0
                site(mu) = pos(site(mu))
                v(site(1),site(2)) = 0

                accepth = accepth + 1
                u_tot = u_tot + delta_e

                end if
              end if

          else ! move it "positive"

            eo = e_field(mu,site(1),site(2))
            en = eo - increment
            old_e = 0.5 * eps_0 * lambda**2 * eo**2
            new_e = 0.5 * eps_0 * lambda**2 * en**2
            delta_e = new_e - old_e

            ! get pos in the mu direction
            site(mu) = pos(site(mu))

            if (v(site(1),site(2)).eq.0) then
              if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

                v(site(1),site(2)) = charge

                ! go back to the original site and set the charge to 0
                site(mu) = neg(site(mu))
                v(site(1),site(2)) = 0
                e_field(mu,site(1),site(2)) = en
                ebar(mu) = ebar(mu) - (increment / L**2)

                accepth = accepth + 1
                u_tot = u_tot + delta_e

              end if
            end if

          end if ! end "positive" / "negative" choice
        end if ! end charge.ne.0 block

      end do ! end charge hop sweep

      if (ebar(mu).gt.(g_thr).or.ebar(mu).lt.((-1)*g_thr)) then
        glob = 1
      else
        glob = 0
      end if

      mu = 0; increment = 0.0;
      u_tot = 0.0
      u_tot = 0.5 * eps_0 * lambda**2 * sum(e_field * e_field)

    end subroutine hop

    subroutine rot(attempts)
      use common
      implicit none
      integer, intent(in) :: attempts
      integer :: i, mu1, mu2
      integer, dimension(2) :: site
      real(kind=8) :: eo1, eo2, eo3, eo4, en1, en2, en3, en4,&
      increment, old_e, new_e, delta_e

      ! --- ROTATIONAL UPDATE ---

      do i = 1, attempts

        eo1 = 0.0; eo2 = 0.0; eo3 = 0.0; eo4 = 0.0
        en1 = 0.0; en2 = 0.0; en3 = 0.0; en4 = 0.0
        site = (/ 0, 0 /)

        ! pick at random from interval [-Delta_max, +Delta_max]
        increment = 2 * rot_delt * (rand() - 0.5)

        ! this is completely unnecessary in 2d, it's left over
        ! from the 3d version. could be removed.
        mu1 = floor(2*rand())+1
        mu2 = mod(mu1,2) + 1

        site = (/ int(rand() * L) + 1,&
                 &int(rand() * L) + 1 /)
        ! site(mu1) is the coordinate in the first direction
        ! site(mu2) is the coordinate in the second

        ! convention in 2d: a plaquette is defined as
        ! (1,i,j),(2,i,j),(1,neg(i),j),(2,i,neg(j))
        eo1 = e_field(mu1,site(1),site(2))
        eo2 = e_field(mu2,site(1),site(2))

        ! if e.g. we picked xy, the next line does x -> neg(x)
        site(mu1) = neg(site(mu1))
        eo4 = e_field(mu1,site(1),site(2))
        ! and now we need to put it back
        site(mu1) = pos(site(mu1))

        ! same for y
        site(mu2) = neg(site(mu2))
        eo3 = e_field(mu2,site(1),site(2))
        site(mu2) = pos(site(mu2))

        en1 = eo1 + increment
        en2 = eo2 - increment
        en3 = eo3 + increment
        en4 = eo4 - increment

        old_e = 0.5 * eps_0 * lambda**2 *&
                (eo1**2 + eo2**2 + eo3**2 + eo4**2)
        new_e = 0.5 * eps_0 * lambda**2 *&
                (en1**2 + en2**2 + en3**2 + en4**2)
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
          u_tot = u_tot + delta_e

        end if ! end of Metropolis check

      end do ! end rotational

    end subroutine rot

    subroutine harm(attempts)
      use common
      implicit none
      integer, intent(in) :: attempts
      integer :: i, mu, pm1
      real(kind=8) :: old_e, new_e, delta_e, increment

      ! --- HARMONIC UPDATE ---

      do i = 1, attempts

        increment = (q) / (L * eps_0 * lambda)

        do mu = 1,2

          if (rand() - 0.5.le.0) then
            pm1 = -1
          else
            pm1 = +1
          end if

          old_e = ebar(mu)**2
          new_e = (ebar(mu) + pm1 * increment)**2
          delta_e = new_e - old_e

          if ((delta_e.lt.0.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then

            ebar(mu) = ebar(mu) + pm1 * increment
            e_field(mu,:,:) = e_field(mu,:,:) + pm1 * (increment)

            acceptg = acceptg + 1

          end if ! end weird Metropolis block

        end do ! end mu1 loop

      end do ! end i loop

      u_tot = 0.0
      u_tot = 0.5 * eps_0 * sum(e_field * e_field)

    end subroutine harm

    subroutine measure(step_number)
      use common
      use linear_solver
      use omp_lib
      implicit none
      integer,intent(in) :: step_number
      integer :: omp_index,i,j,n,kx,ky,m,p,s,x,y,dist_bin
      real(kind=8) :: norm_k, dist, ener_tot, ener_rot, ener_irrot
      real(kind=8) :: ener_tot_sq, ener_rot_sq, ener_irrot_sq
      real(kind=8), dimension(2,L,L) :: e_rot
      complex(kind=rk) :: rho_k_p_temp, rho_k_m_temp
      complex(kind=rk) :: e_kx_temp, e_ky
      complex(kind=rk) :: mnphi_kx_temp, mnphi_ky
      complex(kind=rk) :: e_rot_kx_temp, e_rot_ky
      complex(kind=rk) :: imag, kdotx

      ! | --------------- SUBROUTINE MEASURE(STEP_NUMBER) --------------- |
      ! |                 "MEASURES" RELEVANT QUANTITIES                  |
      ! |                                                                 |
      ! | Basically just take an integer and print out the fields.        |
      ! | Also runs the linear solver and prints out the irrotational     |
      ! | field, along with the charge distribution. We read this back in |
      ! | at the end, and do analysis then.                               |
      ! |                                                                 |
      ! | --------------------------------------------------------------- |

      norm_k = 0.0; dist = 0.0
      rho_k_p_temp = (0.0,0.0); rho_k_m_temp = (0.0,0.0)
      e_kx_temp = (0.0,0.0); mnphi_kx_temp = (0.0,0.0)
      e_rot_kx_temp = (0.0,0.0)
      e_ky = (0.0,0.0); mnphi_ky = (0.0,0.0); e_rot_ky = (0.0,0.0)
      kdotx = (0.0,0.0); imag = (0.0, 1.0)

      ! get irrotational part of field
      ! this way we can get decomposed parts along with total
      call linsol
      e_rot = e_field - mnphi

      !write (*,*) e_tot_avg(1,1,1), e_field(1,1,1)
      e_tot_avg =       e_tot_avg + e_field
      e_rot_avg =       e_rot_avg + e_rot
      e_irrot_avg =     e_irrot_avg + mnphi
      v_avg =           v_avg + float(v)
      ener_tot =        0.5 * eps_0 * sum(e_field * e_field)
      ener_rot =        0.5 * eps_0 * sum(e_rot * e_rot)
      ener_irrot =      0.5 * eps_0 * sum(mnphi * mnphi)
      ener_tot =        ener_tot / L**2
      ener_rot =        ener_rot / L**2
      ener_irrot =      ener_irrot / L**2
      ener_tot_sq =     ener_tot**2
      ener_rot_sq =     ener_rot**2
      ener_irrot_sq =   ener_irrot**2
      ener_tot_sum =    ener_tot_sum + ener_tot
      ener_rot_sum =    ener_rot_sum + ener_rot
      ener_irrot_sum =  ener_irrot_sum + ener_irrot
      ener_tot_sq_sum =    ener_tot_sq_sum + ener_tot_sq
      ener_rot_sq_sum =    ener_rot_sq_sum + ener_rot_sq
      ener_irrot_sq_sum =  ener_irrot_sq_sum + ener_irrot_sq

      ebar(1) = sum(e_field(1,:,:))
      ebar(2) = sum(e_field(2,:,:))
      ebar = ebar / L**2
      ebar_sum(1) = ebar_sum(1) + ebar(1)
      ebar_sum(2) = ebar_sum(2) + ebar(2)
      ebar_sq_sum(1) = ebar_sq_sum(1) + (ebar(1) * ebar(1))
      ebar_sq_sum(2) = ebar_sq_sum(2) + (ebar(2) * ebar(2))

      if (do_corr) then
        n = step_number / sample_interval

        !$omp parallel do num_threads(2)&
        !$omp& private(i,j,m,p,s,kx,ky,rho_k_p_temp,rho_k_m_temp,e_kx_temp,&
        !$omp& mnphi_kx_temp,e_rot_kx_temp,e_ky,mnphi_ky,e_rot_ky,norm_k,kdotx)&
        !$omp& shared(dir_struc,s_ab,s_ab_rot,s_ab_irrot,dist_r,bin_count)
        do omp_index = 1, ((L*bz)+1)**2

          i = ((omp_index - 1) / ((L*bz)+1)) + 1
          j = mod(omp_index - 1, ((L*bz)+1)) + 1
          kx = i - 1 - bz*(L/2)
          ky = j - 1 - bz*(L/2)

          rho_k_p_temp = (0.0,0.0)
          rho_k_m_temp = (0.0,0.0)
          e_kx_temp = (0.0,0.0)
          mnphi_kx_temp = (0.0,0.0)
          e_rot_kx_temp = (0.0,0.0)
          e_ky = (0.0,0.0)
          mnphi_ky = (0.0,0.0)
          e_rot_ky = (0.0,0.0)
          kdotx = (0.0,0.0)
          norm_k = 0.0

          if (kx.eq.0.and.ky.eq.0) then
            norm_k = 0.0
          else
            norm_k = 1.0/(((2*pi/(L*lambda))**2)*dble(kx**2 + ky**2))
          end if

          do s = 1,L**2
          !do m = 1,L
          !  do p = 1,L
            m = ((s - 1) / L) + 1
            p = mod(s - 1, L) + 1

              ! different offsets for x,y,z
              kdotx = ((-1)*imag*(2*pi/(L*lambda))*((m-(1.0/2))*kx + &
                      ((p-1)*ky)))

              ! we want this for every step so we can
              ! average at the end to get field-field struc
              e_kx_temp = e_kx_temp + exp(kdotx)*e_field(1,m,p)
              mnphi_kx_temp = mnphi_kx_temp + exp(kdotx)*mnphi(1,m,p)
              e_rot_kx_temp = e_rot_kx_temp + exp(kdotx)*e_rot(1,m,p)

              kdotx = ((-1)*imag*(2*pi/(L*lambda))*((m-1)*kx + &
                      ((p-(1.0/2))*ky)))

              e_ky = e_ky + exp(kdotx)*e_field(2,m,p)
              mnphi_ky = mnphi_ky + exp(kdotx)*mnphi(2,m,p)
              e_rot_ky = e_rot_ky + exp(kdotx)*e_rot(2,m,p)

              if (v(m,p).ne.0) then ! calculate <++ + +->!

                ! FT of charge distribution
                kdotx = ((-1)*imag*2*pi*(((m-1)*kx/(L*lambda)) + &
                        ((p-1)*ky/(L*lambda))))

                if (v(m,p).eq.-1) then
                  rho_k_m_temp = rho_k_m_temp + q * v(m,p) * exp(kdotx)
                end if

                if (v(m,p).eq.1) then ! take away <++>
                  rho_k_p_temp = rho_k_p_temp + q * v(m,p) * exp(kdotx)
                end if

                ! --- real space correlation function ---

                if (kx.gt.0.and.kx.le.L.and.&
                    ky.gt.0.and.ky.le.L) then

                  ! this should sort it out. kx ky kz here are
                  ! the "r"s, m p s are the "zeros"
                  ! we want to do rho_+(0) rho_-(r)
                  if (v(m,p).eq.1) then
                    if (v(kx,ky).eq.-1) then

                      x = abs(m - kx)
                      y = abs(p - ky)

                      if (x.gt.L/2) then
                        x = L - x
                      end if
                      if (y.gt.L/2) then
                        y = L - y
                      end if

                      x = x + 1
                      y = y + 1

                      dir_struc(x,y) = dir_struc(x,y) +&
                              v(m,p) * v(kx,ky)
                    end if ! neg charge at kx,ky,kz
                  end if ! pos charge at m,p,s
                end if ! kx, ky, kz < L

              end if ! end v != 0 check

          !  end do
          !end do ! m,p blocks
          end do ! s

          rho_k_p_temp = rho_k_p_temp / float(L**2)
          rho_k_m_temp = rho_k_m_temp / float(L**2)
          e_kx_temp = e_kx_temp / float(L**2)
          e_ky = e_ky / float(L**2)
          mnphi_kx_temp = mnphi_kx_temp / float(L**2)
          mnphi_ky = mnphi_ky / float(L**2)
          e_rot_kx_temp = e_rot_kx_temp / float(L**2)
          e_rot_ky = e_rot_ky / float(L**2)

          rho_k_p(i,j) = rho_k_p(i,j) + rho_k_p_temp
          rho_k_m(i,j) = rho_k_m(i,j) + rho_k_m_temp
          ch_ch(i,j) = ch_ch(i,j) +&
                        (rho_k_p_temp * conjg(rho_k_m_temp))

          s_ab(1,1,i,j) = s_ab(1,1,i,j) +&
          e_kx_temp*conjg(e_kx_temp)
          s_ab(1,2,i,j) = s_ab(1,2,i,j) +&
          e_kx_temp*conjg(e_ky)
          s_ab(2,1,i,j) = s_ab(2,1,i,j) +&
          e_ky*conjg(e_kx_temp)
          s_ab(2,2,i,j) = s_ab(2,2,i,j) +&
          e_ky*conjg(e_ky)

          !write(*,*) i,j,s_ab_rot(1,1,i,j),e_rot_kx_temp
          s_ab_rot(1,1,i,j) = s_ab_rot(1,1,i,j) +&
          e_rot_kx_temp*conjg(e_rot_kx_temp)
          s_ab_rot(1,2,i,j) = s_ab_rot(1,2,i,j) +&
          e_rot_kx_temp*conjg(e_rot_ky)
          s_ab_rot(2,1,i,j) = s_ab_rot(2,1,i,j) +&
          e_rot_ky*conjg(e_rot_kx_temp)
          s_ab_rot(2,2,i,j) = s_ab_rot(2,2,i,j) +&
          e_rot_ky*conjg(e_rot_ky)

          s_ab_irrot(1,1,i,j) = s_ab_irrot(1,1,i,j) +&
          mnphi_kx_temp*conjg(mnphi_kx_temp)
          s_ab_irrot(1,2,i,j) = s_ab_irrot(1,2,i,j) +&
          mnphi_kx_temp*conjg(mnphi_ky)
          s_ab_irrot(2,1,i,j) = s_ab_irrot(2,1,i,j) +&
          mnphi_ky*conjg(mnphi_kx_temp)
          s_ab_irrot(2,2,i,j) = s_ab_irrot(2,2,i,j) +&
          mnphi_ky*conjg(mnphi_ky)

          ! direct space one
          ! write (*,*) i,j,(i.le.L/2+1.and.j.le.L/2+1)
          if (i.le.L/2+1.and.j.le.L/2+1) then

            dist = sqrt(dble((i - 1)**2 + (j - 1)**2))
            dist_bin = floor(dist / bin_size) + 1
            dist_r(dist_bin) = dist_r(dist_bin) + dir_struc(i,j)
            bin_count(dist_bin) = bin_count(dist_bin) + 1

          end if

        end do ! end openmp_index loop
        !$omp end parallel do

      end if ! if do_corr

    end subroutine measure

end module update
