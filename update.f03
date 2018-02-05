
module update
  use common
  ! not sure if linear solver is needed
  use linear_solver
  implicit none

  contains

    subroutine hop(n)
      implicit none
      integer, intent(in) :: n
      integer :: i,j,mu,v1o,v2o,v1n,v2n,pm,charge
      integer, dimension(3) :: site
      real(kind=8) :: eo, en, old_e, new_e, delta_e, increment

      ! --- CHARGE HOP UPDATE ---

      do i = 1, n

        ! NOTE TO SELF: this whole procedure assumes
        ! single-valued charges only.

        ! pick a random site and random component, x or y
        site = (/ int(rand() * L) + 1,&
                  int(rand() * L) + 1,&
                 &int(rand() * L) + 1 /)
        mu = floor(3*rand())+1
        mu_tot = mu_tot + mu

        ! check we're not doing anything weird with multiple charges
        if (abs(v(site(1),site(2),site(3))).gt.1) then
          write (*,*) "Charge at ",site(1),site(2),site(3),&
            " = ",v(site(1),site(2),site(3)),". Exiting."
          write (*,*)
          stop
        end if

        ! pick a non-zero charge - canonical!
        if (v(site(1),site(2),site(3)).ne.0) then

          charge = v(site(1),site(2),site(3))

          ! this takes care of sign issues when hopping
          increment = q * charge / (eps_0 * lambda)

          if (rand().lt.0.5) then ! move it "negative"

            attempts(1) = attempts(1) + 1

            eo = e_field(mu,site(1),site(2),site(3))
            en = eo + increment
            old_e = 0.5 * eps_0 * eo**2
            new_e = 0.5 * eps_0 * en**2
            delta_e = new_e - old_e

            site(mu) = neg(site(mu))

            if (v(site(1),site(2),site(3)).eq.0) then
              if ((delta_e.lt.0.0).or.(exp((-beta)*delta_e).gt.rand())) then

                v(site(1),site(2),site(3)) = charge
                ebar(mu) = ebar(mu) + (increment / L**3)

                ! go back to the original site and set the charge to 0
                site(mu) = pos(site(mu))
                e_field(mu,site(1),site(2),site(3)) = en
                v(site(1),site(2),site(3)) = 0

                accepts(1) = accepts(1) + 1
                u_tot = u_tot + delta_e

                end if
              end if

          else ! move it "positive"

            attempts(1) = attempts(1) + 1

            site(mu) = pos(site(mu))
            eo = e_field(mu,site(1),site(2),site(3))
            en = eo - increment
            old_e = 0.5 * eps_0 * lambda**2 * eo**2
            new_e = 0.5 * eps_0 * lambda**2 * en**2
            delta_e = new_e - old_e

            ! get pos in the mu direction

            if (v(site(1),site(2),site(3)).eq.0) then
              if ((delta_e.lt.0.0).or.(exp((-beta) * delta_e).gt.rand())) then

                v(site(1),site(2),site(3)) = charge
                e_field(mu,site(1),site(2),site(3)) = en
                ebar(mu) = ebar(mu) - (increment / L**3)

                ! go back to the original site and set the charge to 0
                site(mu) = neg(site(mu))
                v(site(1),site(2),site(3)) = 0

                accepts(1) = accepts(1) + 1
                u_tot = u_tot + delta_e

              end if
            end if

          end if ! end "positive" / "negative" choice
        end if ! end charge.ne.0 block

      end do ! end charge hop sweep

      if (abs(ebar(1)).gt.g_thr.or.&
         &abs(ebar(2)).gt.g_thr.or.&
         &abs(ebar(3)).gt.g_thr) then
        glob = 1
      else
        glob = 0
      end if

      mu = 0; increment = 0.0;
      u_tot = 0.0
      u_tot = 0.5 * eps_0 * lambda**3 * sum(e_field * e_field)

    end subroutine hop

    subroutine rot(n)
      use common
      implicit none
      integer, intent(in) :: n
      integer :: i, mu1, mu2
      integer, dimension(3) :: site
      real(kind=8) :: eo1, eo2, eo3, eo4, en1, en2, en3, en4,&
      increment, old_e, new_e, delta_e

      ! --- ROTATIONAL UPDATE ---

      do i = 1, n

        eo1 = 0.0; eo2 = 0.0; eo3 = 0.0; eo4 = 0.0
        en1 = 0.0; en2 = 0.0; en3 = 0.0; en4 = 0.0
        site = (/ 0, 0, 0 /)

        ! pick at random from interval [-Delta_max, +Delta_max]
        ! increment = 2 * (rot_delt / eps_0) * (rand() - 0.5)
        increment = 2 * (rot_delt) * (rand() - 0.5)


        ! this is completely unnecessary in 2d, it's left over
        ! from the 3d version. could be removed.
        mu1 = floor(3 * rand()) + 1
        mu2 = mod(mu1, 3) + 1

        site = (/ int(rand() * L) + 1,&
                  int(rand() * L) + 1,&
                 &int(rand() * L) + 1 /)

        ! convention in 2d: a plaquette is defined as
        ! (1,i,j),(2,i,j),(1,i,neg(j)),(2,neg(i),j)
        eo1 = e_field(mu1,site(1),site(2),site(3))
        eo2 = e_field(mu2,site(1),site(2),site(3))

        ! if e.g. we picked xy, the next line does x -> neg(x)
        site(mu1) = neg(site(mu1))
        eo4 = e_field(mu2,site(1),site(2),site(3))
        ! and now we need to put it back
        site(mu1) = pos(site(mu1))

        ! same for y
        site(mu2) = neg(site(mu2))
        eo3 = e_field(mu1,site(1),site(2),site(3))
        site(mu2) = pos(site(mu2))

        en1 = eo1 + increment
        en2 = eo2 - increment
        en3 = eo3 - increment
        en4 = eo4 + increment

        old_e = 0.5 * lambda**3 *&
                (eo1**2 + eo2**2 + eo3**2 + eo4**2)
        new_e = 0.5 * lambda**3 *&
                (en1**2 + en2**2 + en3**2 + en4**2)
        delta_e = new_e - old_e

        if ((delta_e.lt.0.0).or.(exp((-1.0*beta)*delta_e).gt.rand())) then

          e_field(mu1,site(1),site(2),site(3)) = en1
          e_field(mu2,site(1),site(2),site(3)) = en2

          site(mu1) = neg(site(mu1))
          e_field(mu2,site(1),site(2),site(3)) = en4
          site(mu1) = pos(site(mu1))

          site(mu2) = neg(site(mu2))
          e_field(mu1,site(1),site(2),site(3)) = en3
          site(mu2) = pos(site(mu2))

          accepts(2) = accepts(2) + 1
          u_tot = u_tot + delta_e

        end if ! end of Metropolis check

      end do ! end rotational

    end subroutine rot

    subroutine harm(n)
      use common
      implicit none
      integer, intent(in) :: n
      integer :: i, mu, pm1
      real(kind=8) :: delta_e, increment

      ! --- HARMONIC UPDATE ---

      do i = 1, n

        increment = (q) / (L**2 * eps_0)

        do mu = 1,3

          if ((rand() - 0.5).le.0.0) then
            pm1 = -1
          else
            pm1 = +1
          end if

          ! delta_e = (lambda**3 * (q / L**2)) * ((q * L / (2 * eps_0)) +&
          !           (pm1 * ebar(mu)))
          delta_e = (lambda**3 * (q * L)) * ((q / (2 * L**2 * eps_0)) +&
                    (pm1 * ebar(mu)))

          if ((delta_e.lt.0.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then

            ebar(mu) = ebar(mu) + pm1 * increment
            e_field(mu,:,:,:) = e_field(mu,:,:,:) + pm1 * (increment)

            accepts(3) = accepts(3) + 1

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
      integer :: foo,omp_index,i,j,k,n,kx,ky,kz,m,p,s,x,y,z,dist_bin
      real(kind=8) :: norm_k, dist, ener_tot, ener_rot, ener_irrot
      real(kind=8) :: ener_tot_sq, ener_rot_sq, ener_irrot_sq, dp, np
      real(kind=8), dimension(3,L,L,L) :: e_rot
      complex(kind=rk) :: rho_k_p_temp, rho_k_m_temp
      complex(kind=rk) :: e_kx, e_ky, e_kz
      complex(kind=rk) :: mnphi_kx, mnphi_ky, mnphi_kz
      complex(kind=rk) :: e_rot_kx, e_rot_ky, e_rot_kz
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

      ebar = 0.0; ebar_dip = 0.0; ebar_wind = 0.0; dp = 0.0; np = 0.0
      norm_k = 0.0; dist = 0.0
      rho_k_p_temp = (0.0,0.0); rho_k_m_temp = (0.0,0.0)
      e_kx = (0.0,0.0); mnphi_kx = (0.0,0.0); e_rot_kx = (0.0,0.0)
      e_ky = (0.0,0.0); mnphi_ky = (0.0,0.0); e_rot_ky = (0.0,0.0)
      e_kz = (0.0,0.0); mnphi_kz = (0.0,0.0); e_rot_kz = (0.0,0.0)
      kdotx = (0.0,0.0); imag = (0.0, 1.0)

      ! get irrotational part of field
      ! this way we can get decomposed parts along with total
      call linsol
      e_rot = e_field - mnphi

      e_tot_avg =           e_tot_avg + e_field
      e_rot_avg =           e_rot_avg + e_rot
      e_irrot_avg =         e_irrot_avg + mnphi
      v_avg =               v_avg + float(v)
      rho_avg =             rho_avg + (dble(sum(abs(v))) / L**3)
      ener_tot =            0.5 * lambda**3 * eps_0 * sum(e_field * e_field)
      ener_rot =            0.5 * lambda**3 * eps_0 * sum(e_rot * e_rot)
      ener_irrot =          0.5 * lambda**3 * eps_0 * sum(mnphi * mnphi)
      ener_tot =            ener_tot / L**3
      ener_rot =            ener_rot / L**3
      ener_irrot =          ener_irrot / L**3
      ener_tot_sq =         ener_tot**2
      ener_rot_sq =         ener_rot**2
      ener_irrot_sq =       ener_irrot**2
      ener_tot_sum =        ener_tot_sum + ener_tot
      ener_rot_sum =        ener_rot_sum + ener_rot
      ener_irrot_sum =      ener_irrot_sum + ener_irrot
      ener_tot_sq_sum =     ener_tot_sq_sum + ener_tot_sq
      ener_rot_sq_sum =     ener_rot_sq_sum + ener_rot_sq
      ener_irrot_sq_sum =   ener_irrot_sq_sum + ener_irrot_sq

      do i = 1,3

        ebar(i) = sum(e_field(i,:,:,:))
        dp = ebar(i)
        np = 0

        if (ebar(i).gt.((q * dble(L)) / 2)) then
          dp = dp - (q * dble(L))
          np = ebar(i) - dp
        else if (ebar(i).le.((-1.0 * q * dble(L)) / 2)) then
          dp = dp + (q * dble(L))
          np = ebar(i) - dp
        end if

        ebar_dip(i) = dp
        ebar_wind(i) = np

        avg_field_total(i) = avg_field_total(i) +&
                             sum(abs(real(e_field(i,:,:,:))))
        avg_field_rot(i) = avg_field_rot(i) +&
                           sum(abs(real(e_rot(i,:,:,:))))
        avg_field_irrot(i) = avg_field_irrot(i) +&
                             sum(abs(real(mnphi(i,:,:,:))))
        avg_field_sq_total(i) = avg_field_sq_total(i) + avg_field_total(i)**2
        avg_field_sq_rot(i)   = avg_field_sq_rot(i)   + avg_field_rot(i)**2
        avg_field_sq_irrot(i) = avg_field_sq_irrot(i) + avg_field_irrot(i)**2

      end do

      avg_field_total = avg_field_total / L**3
      avg_field_rot = avg_field_rot / L**3
      avg_field_irrot = avg_field_irrot / L**3
      ebar = ebar / L**3
      ebar_dip = ebar_dip / L**3
      ebar_wind = ebar_wind / L**3

      ebar_sum(1) = ebar_sum(1) + ebar(1)
      ebar_sum(2) = ebar_sum(2) + ebar(2)
      ebar_sum(3) = ebar_sum(3) + ebar(3)
      ebar_sq_sum(1) = ebar_sq_sum(1) + (ebar(1) * ebar(1))
      ebar_sq_sum(2) = ebar_sq_sum(2) + (ebar(2) * ebar(2))
      ebar_sq_sum(3) = ebar_sq_sum(3) + (ebar(3) * ebar(3))

      ebar_dip_sum(1) = ebar_dip_sum(1) + ebar_dip(1)
      ebar_dip_sum(2) = ebar_dip_sum(2) + ebar_dip(2)
      ebar_dip_sum(3) = ebar_dip_sum(3) + ebar_dip(3)
      ebar_dip_sq_sum(1) = ebar_dip_sq_sum(1) + (ebar_dip(1) * ebar_dip(1))
      ebar_dip_sq_sum(2) = ebar_dip_sq_sum(2) + (ebar_dip(2) * ebar_dip(2))
      ebar_dip_sq_sum(3) = ebar_dip_sq_sum(3) + (ebar_dip(3) * ebar_dip(3))

      ebar_wind_sum(1) = ebar_wind_sum(1) + ebar_wind(1)
      ebar_wind_sum(2) = ebar_wind_sum(2) + ebar_wind(2)
      ebar_wind_sum(3) = ebar_wind_sum(3) + ebar_wind(3)
      ebar_wind_sq_sum(1) = ebar_wind_sq_sum(1) + (ebar_wind(1) * ebar_wind(1))
      ebar_wind_sq_sum(2) = ebar_wind_sq_sum(2) + (ebar_wind(2) * ebar_wind(2))
      ebar_wind_sq_sum(3) = ebar_wind_sq_sum(3) + (ebar_wind(3) * ebar_wind(3))

      if (do_corr) then
        n = step_number / sample_interval

        !$omp parallel do &
        !$omp& private(i,j,k,m,p,s,kx,ky,kz,rho_k_p_temp,rho_k_m_temp,e_kx,&
        !$omp& mnphi_kx,e_rot_kx,e_ky,mnphi_ky,e_rot_ky,&
        !$omp& e_kz,mnphi_kz,e_rot_kz,norm_k,kdotx,x,y,z)&
        !$omp& shared(dir_struc,s_ab,s_ab_rot,s_ab_irrot,dist_r,bin_count)
        do omp_index = 1, L**3

          ! need to check i here
          i = ((omp_index - 1) / (L**2)) + 1
          j = mod(((omp_index - 1) / (L)), L) + 1
          k = mod(omp_index - 1, (L)) + 1

          ! repurposing these as the s_ab indices at the end:
          ! we're gonna fill up the positive half of those arrays
          ! and then use symmetry properties to fill the rest at the end.
          ! k_mu = 0 point is at s_ab array index L + 1
          kx = i + L
          ky = j + L
          kz = k + L

          rho_k_p_temp = (0.0,0.0)
          rho_k_m_temp = (0.0,0.0)
          e_kx = (0.0,0.0)
          mnphi_kx = (0.0,0.0)
          e_rot_kx = (0.0,0.0)
          e_ky = (0.0,0.0)
          mnphi_ky = (0.0,0.0)
          e_rot_ky = (0.0,0.0)
          e_kz = (0.0,0.0)
          mnphi_kz = (0.0,0.0)
          e_rot_kz = (0.0,0.0)
          kdotx = (0.0,0.0)
          norm_k = 0.0

          if (kx.eq.0.and.ky.eq.0) then
            norm_k = 0.0
          else
            norm_k = 1.0/(((2*pi/(L*lambda))**3)*dble(kx**2 + ky**2 + kz**2))
          end if

          do foo = 1,L**3
            m = ((foo - 1) / L**2) + 1
            p = mod((foo - 1) / L, L) + 1
            s = mod(foo - 1, L) + 1

              ! kdotx = (-1)*(imag*(2*pi/(L*lambda))*((neg(m)+(1.0/2))*kx + &
              !         ((p)*ky)))
              kdotx = hw(m,i) + fw(p,j) + fw(s,k)

              e_kx = e_kx + exp(kdotx)*e_field(1,m,p,s)
              mnphi_kx = mnphi_kx + exp(kdotx)*mnphi(1,m,p,s)
              e_rot_kx = e_rot_kx + exp(kdotx)*e_rot(1,m,p,s)

              ! kdotx = (-1)*(imag*(2*pi/(L*lambda))*((m)*kx + &
              !         ((neg(p)+(1.0/2))*ky)))
              kdotx = fw(m,i) + hw(p,j) + fw(s,k)

              e_ky =     e_ky     + exp(kdotx)*e_field(2,m,p,s)
              e_rot_ky = e_rot_ky + exp(kdotx)*e_rot(2,m,p,s)
              mnphi_ky = mnphi_ky + exp(kdotx)*mnphi(2,m,p,s)

              ! kdotx = (-1)*(imag*(2*pi/(L*lambda))*((m)*kx + &
              !         ((neg(p)+(1.0/2))*ky)))
              kdotx = fw(m,i) + fw(p,j) + hw(s,k)

              e_kz =     e_kz     + exp(kdotx)*e_field(3,m,p,s)
              e_rot_kz = e_rot_kz + exp(kdotx)*e_rot(3,m,p,s)
              mnphi_kz = mnphi_kz + exp(kdotx)*mnphi(3,m,p,s)

              if (v(m,p,s).ne.0) then ! calculate <++ + +->!

                ! FT of charge distribution
                ! kdotx = (-1)*(imag*2*pi*(((m)*kx/(L*lambda)) + &
                !         ((p)*ky/(L*lambda))))
                kdotx = fw(m,i) + fw(p,j) + fw(s,k)

                if (v(m,p,s).eq.-1) then
                  rho_k_m_temp = rho_k_m_temp + q * v(m,p,s) * exp(kdotx)
                end if

                if (v(m,p,s).eq.1) then ! take away <++>
                  rho_k_p_temp = rho_k_p_temp + q * v(m,p,s) * exp(kdotx)
                end if

                ! --- real space correlation function ---

                if (v(m,p,s).eq.1) then
                  if (v(i,j,k).eq.-1) then

                    x = abs(m - i)
                    y = abs(p - j)
                    z = abs(s - k)

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

                    dir_struc(x,y,z) = dir_struc(x,y,z) +&
                            v(m,p,s) * v(i,j,k)
                  end if ! neg charge at i,j,k
                end if ! pos charge at m,p,s

              end if ! end v != 0 check

          end do ! s

          rho_k_p_temp = rho_k_p_temp / float(L**3)
          rho_k_m_temp = rho_k_m_temp / float(L**3)
          e_kx = e_kx / float(L**3)
          e_ky = e_ky / float(L**3)
          e_kz = e_kz / float(L**3)
          mnphi_kx = mnphi_kx / float(L**3)
          mnphi_ky = mnphi_ky / float(L**3)
          mnphi_kz = mnphi_kz / float(L**3)
          e_rot_kx = e_rot_kx / float(L**3)
          e_rot_ky = e_rot_ky / float(L**3)
          e_rot_kz = e_rot_kz / float(L**3)

          rho_k_p(kx,ky,kz) = rho_k_p(kx,ky,kz) + rho_k_p_temp
          rho_k_m(kx,ky,kz) = rho_k_m(kx,ky,kz) + rho_k_m_temp
          ch_ch(kx,ky,kz) = ch_ch(kx,ky,kz) +&
                        (rho_k_p_temp * conjg(rho_k_m_temp))

          s_ab(1,1,kx,ky,kz) = s_ab(1,1,kx,ky,kz) +&
          e_kx*conjg(e_kx)
          s_ab(1,2,kx,ky,kz) = s_ab(1,2,kx,ky,kz) +&
          e_kx*conjg(e_ky)
          s_ab(1,3,kx,ky,kz) = s_ab(1,3,kx,ky,kz) +&
          e_kx*conjg(e_kz)
          s_ab(2,1,kx,ky,kz) = s_ab(2,1,kx,ky,kz) +&
          e_ky*conjg(e_kx)
          s_ab(2,2,kx,ky,kz) = s_ab(2,2,kx,ky,kz) +&
          e_ky*conjg(e_ky)
          s_ab(2,3,kx,ky,kz) = s_ab(2,3,kx,ky,kz) +&
          e_ky*conjg(e_kz)
          s_ab(3,1,kx,ky,kz) = s_ab(3,1,kx,ky,kz) +&
          e_kz*conjg(e_kx)
          s_ab(3,2,kx,ky,kz) = s_ab(3,2,kx,ky,kz) +&
          e_kz*conjg(e_ky)
          s_ab(3,3,kx,ky,kz) = s_ab(3,3,kx,ky,kz) +&
          e_kz*conjg(e_kz)

          s_ab_rot(1,1,kx,ky,kz) = s_ab_rot(1,1,kx,ky,kz) +&
          e_rot_kx*conjg(e_rot_kx)
          s_ab_rot(1,2,kx,ky,kz) = s_ab_rot(1,2,kx,ky,kz) +&
          e_rot_kx*conjg(e_rot_ky)
          s_ab_rot(1,3,kx,ky,kz) = s_ab_rot(1,3,kx,ky,kz) +&
          e_rot_kx*conjg(e_rot_kz)
          s_ab_rot(2,1,kx,ky,kz) = s_ab_rot(2,1,kx,ky,kz) +&
          e_rot_ky*conjg(e_rot_kx)
          s_ab_rot(2,2,kx,ky,kz) = s_ab_rot(2,2,kx,ky,kz) +&
          e_rot_ky*conjg(e_rot_ky)
          s_ab_rot(2,3,kx,ky,kz) = s_ab_rot(2,3,kx,ky,kz) +&
          e_rot_ky*conjg(e_rot_kz)
          s_ab_rot(3,1,kx,ky,kz) = s_ab_rot(3,1,kx,ky,kz) +&
          e_rot_kz*conjg(e_rot_kx)
          s_ab_rot(3,2,kx,ky,kz) = s_ab_rot(3,2,kx,ky,kz) +&
          e_rot_kz*conjg(e_rot_ky)
          s_ab_rot(3,3,kx,ky,kz) = s_ab_rot(3,3,kx,ky,kz) +&
          e_rot_kz*conjg(e_rot_kz)

          s_ab_irrot(1,1,kx,ky,kz) = s_ab_irrot(1,1,kx,ky,kz) +&
          mnphi_kx*conjg(mnphi_kx)
          s_ab_irrot(1,2,kx,ky,kz) = s_ab_irrot(1,2,kx,ky,kz) +&
          mnphi_kx*conjg(mnphi_ky)
          s_ab_irrot(1,3,kx,ky,kz) = s_ab_irrot(1,3,kx,ky,kz) +&
          mnphi_kx*conjg(mnphi_kz)
          s_ab_irrot(2,1,kx,ky,kz) = s_ab_irrot(2,1,kx,ky,kz) +&
          mnphi_ky*conjg(mnphi_kx)
          s_ab_irrot(2,2,kx,ky,kz) = s_ab_irrot(2,2,kx,ky,kz) +&
          mnphi_ky*conjg(mnphi_ky)
          s_ab_irrot(2,3,kx,ky,kz) = s_ab_irrot(2,3,kx,ky,kz) +&
          mnphi_ky*conjg(mnphi_kz)
          s_ab_irrot(3,1,kx,ky,kz) = s_ab_irrot(3,1,kx,ky,kz) +&
          mnphi_kz*conjg(mnphi_kx)
          s_ab_irrot(3,2,kx,ky,kz) = s_ab_irrot(3,2,kx,ky,kz) +&
          mnphi_kz*conjg(mnphi_ky)
          s_ab_irrot(3,3,kx,ky,kz) = s_ab_irrot(3,3,kx,ky,kz) +&
          mnphi_kz*conjg(mnphi_kz)

          ! direct space one
          if (i.le.L/2+1.and.j.le.L/2+1.and.k.le.L/2+1) then

            dist = sqrt(dble((i - 1)**2 + (j - 1)**2 + (k - 1)**2))
            dist_bin = floor(dist / bin_size) + 1
            dist_r(dist_bin) = dist_r(dist_bin) + dir_struc(i,j,k)
            bin_count(dist_bin) = bin_count(dist_bin) + 1

          end if

        end do ! end openmp_index loop
        !$omp end parallel do

      end if ! if do_corr

    end subroutine measure

end module update
