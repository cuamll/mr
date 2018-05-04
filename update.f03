
module update
  use common
  ! not sure if linear solver is needed
  use linear_solver
  implicit none

  contains

    subroutine markov_chain_HXY
     use common
     implicit none
     integer :: n, i, j
     real(kind=8) :: thetaOld, thetaNew, deltaTheta, Uold, Unew, deltaU
     real(kind=8) :: top1old, top2old, top3old, top4old,&
            top1new, top2new, top3new, top4new

     ! ORIGINAL GRADIENT DIRECTION
     do n = 1, L**2
        i = int(rand() * L) + 1
        j = int(rand() * L) + 1
        deltaU = 0.0

        thetaOld = theta(i,j)
        ! deltaTheta = 2. * proposalInterval * (rand() - 0.5)
        deltaTheta = 2. * rot_delt * (rand() - 0.5)
        thetaNew = thetaOld + deltaTheta
        if (thetaNew .le. -pi) then
           thetaNew = thetaNew + twopi
        else if (thetaNew .gt. pi) then
           thetaNew = thetaNew - twopi
        end if

        ! CALL OLD EMERGENT FIELD

        top1old = top_x(i,j)
        top2old = top_y(i,j)
        top3old = top_x(i,pos(j))
        top4old = top_y(pos(i),j)

        ! PROPOSED EMERGENT FIELD

        top1new = thetanew - theta(i,neg(j))

        ! here is where i should be figuring out the charge logic
        ! in order to do the helmholtz decomposition
        if (top1new .gt. pi) then
           top1new = top1new - twopi
        else if (top1new .le. -pi) then
           top1new = top1new + twopi
        end if

        top2new = - (thetanew - theta(neg(i),j))
        if (top2new .gt. pi) then
           top2new = top2new-twopi
        else if (top2new .le. -pi) then
           top2new = top2new + twopi
        end if

        top3new = theta(i,pos(j)) - thetanew
        if (top3new .gt. pi) then
           top3new = top3new - twopi
        else if (top3new .le. -pi) then
           top3new = top3new + twopi
        end if

        top4new = - (theta(pos(i),j) - thetanew)
        if (top4new .gt. pi) then
           top4new = top4new - twopi
        else if (top4new .le. -pi) then
           top4new = top4new + twopi
        end if

        ! METROPOLIS FILTER

        Uold = 0.5 * eps_0 * (top1old**2 + top2old**2 + top3old**2 + top4old**2)
        Unew = 0.5 * eps_0 * (top1new**2 + top2new**2 + top3new**2 + top4new**2)
        deltaU = Unew - Uold

        attempts(6) = attempts(6) + 1

        if ((deltaU.lt.0.0).or.(exp((-beta)*deltaU).gt.rand())) then
           theta(i,j) = thetaNew
           top_x(i,j) = top1new
           top_y(i,j) = top2new
           top_x(i,pos(j)) = top3new
           top_y(pos(i),j) = top4new
           ! uses same place as harm_fluct from below; not a problem atm
           accepts(6) = accepts(6) + 1
        end if
     end do

     return
    end subroutine markov_chain_HXY

    subroutine harm_fluct(n)
      use common
      implicit none
      integer, intent(in) :: n
      integer :: i, mu
      real(kind=8) :: increment, delta_e, harm_delt

      ! --- HARMONIC FLUCTUATION UPDATE ---

      ! just testing
      harm_delt = 0.5 * rot_delt

      do i = 1, n

        ! pick at random from interval [-Delta_max, +Delta_max]
        increment = 2 * harm_delt * (rand() - 0.5)

        mu = floor(2 * rand()) + 1

        ! expression for an arbitrary increment
        delta_e = ((eps_0 * lambda**2) / 2.0) *&
                  (L**2 * increment) * (increment + 2 * ebar(mu))

        attempts(6) = attempts(6) + 1

        if ((delta_e.lt.0.0).or.(exp((-beta)*delta_e).gt.rand())) then

          accepts(6) = accepts(6) + 1
          ebar(mu) = ebar(mu) + increment
          e_field(mu,:,:) = e_field(mu,:,:) + increment

        end if ! end of Metropolis check

      end do ! i loop

    end subroutine harm_fluct

    subroutine hop(n)
      implicit none
      integer, intent(in) :: n
      integer :: i,j,mu,v1o,v2o,v1n,v2n,pm
      integer, dimension(2) :: site
      real(kind=8) :: eo, en, old_e, new_e, delta_e, increment

      ! --- CHARGE HOP UPDATE ---

      do i = 1, n

        ! NOTE TO SELF: this whole procedure assumes
        ! single-valued charges only.

        ! pick a random site and random component, x or y
        site = (/ int(rand() * L) + 1,&
                  &int(rand() * L) + 1 /)
        mu = floor(2*rand())+1
        mu_tot = mu_tot + mu

        ! check we're not doing anything weird with multiple charges
        if (abs(v(site(1),site(2))).gt.1) then
          write (*,*) "Charge at ",site(1),site(2),&
            " = ",v(site(1),site(2)),". Exiting."
          write (*,*)
          stop
        end if

        ! grand canonical now

        ! charge = v(site(1),site(2))
        increment = q / (eps_0 * lambda)

        if (rand().lt.0.5) then ! increase field bond
          pm = 1
        else ! decrease field bond
          pm = -1
        end if

        eo = e_field(mu,site(1),site(2))
        en = eo + pm * increment
        old_e = 0.5 * eps_0 * eo**2
        new_e = 0.5 * eps_0 * en**2
        delta_e = new_e - old_e

        ! adding to the field bond means charge moving "backwards"
        ! current charges tell us if this is hop/creation/annihilation
        v1o = v(site(1),site(2))
        v1n = v(site(1),site(2)) - pm
        ! we need negative in the mu direction
        site(mu) = neg(site(mu))
        v2o = v(site(1),site(2))
        v2n = v(site(1),site(2)) + pm

        if ((v1o.eq.1.and.v2o.eq.0.and.pm.eq.1).or.&
            (v1o.eq.-1.and.v2o.eq.0.and.pm.eq.-1).or.&
            (v1o.eq.0.and.v2o.eq.1.and.pm.eq.-1).or.&
            (v1o.eq.0.and.v2o.eq.-1.and.pm.eq.1)) then
          ! HOP
          attempts(1) = attempts(1) + 1
        else if ((v1o.eq.0.and.v2o.eq.0)) then
          ! CREATION
          attempts(4) = attempts(4) + 1
        else if (abs(v1o).eq.1.and.abs(v2o).eq.1) then
          ! ANNIHILATION
          attempts(5) = attempts(5) + 1
        end if

        if (abs(v1n).le.1.and.abs(v2n).le.1) then
          if ((delta_e.lt.0.0).or.(exp((-beta)*delta_e).gt.rand())) then

            ! site still pointing at neg(orig. site)
            v(site(1),site(2)) = v2n
            site(mu) = pos(site(mu))
            v(site(1),site(2)) = v1n
            e_field(mu,site(1),site(2)) = en
            ! mnphi(mu,site(1),site(2)) = mnphi(mu,site(1),site(2)) +&
            !                             pm * increment
            ebar(mu) = ebar(mu) + (pm * increment / L**2)

            u_tot = u_tot + delta_e

            if ((v1o.eq.1.and.v2o.eq.0.and.pm.eq.1).or.&
                (v1o.eq.-1.and.v2o.eq.0.and.pm.eq.-1).or.&
                (v1o.eq.0.and.v2o.eq.1.and.pm.eq.-1).or.&
                (v1o.eq.0.and.v2o.eq.-1.and.pm.eq.1)) then
              accepts(1) = accepts(1) + 1
            else if ((v1o.eq.0.and.v2o.eq.0)) then
              accepts(4) = accepts(4) + 1
            else if (abs(v1o).eq.1.and.abs(v2o).eq.1) then
              accepts(5) = accepts(5) + 1
            else
              write (6,'(a,7i3.1)') "SOME WEIRD SHIT GOT ACCEPTED: ",&
              site(1), site(2), v1o, v2o, v1n, v2n, pm
            end if

            end if
          end if

      end do ! end charge hop sweep

      glob = 1

      mu = 0; increment = 0.0;
      u_tot = 0.0
      u_tot = 0.5 * eps_0 * lambda**2 * sum(e_field * e_field)

    end subroutine hop

    subroutine rot(n)
      use common
      implicit none
      integer, intent(in) :: n
      integer :: i, mu1, mu2
      integer, dimension(2) :: site
      real(kind=8) :: eo1, eo2, eo3, eo4, en1, en2, en3, en4,&
      increment, old_e, new_e, delta_e

      ! --- ROTATIONAL UPDATE ---

      do i = 1, n

        eo1 = 0.0; eo2 = 0.0; eo3 = 0.0; eo4 = 0.0
        en1 = 0.0; en2 = 0.0; en3 = 0.0; en4 = 0.0
        site = (/ 0, 0 /)

        ! pick at random from interval [-Delta_max, +Delta_max]
        increment = 2 * rot_delt * (rand() - 0.5)

        ! this is completely unnecessary in 2d, it's left over
        ! from the 3d version. could be removed.
        mu1 = floor(2 * rand()) + 1
        mu2 = mod(mu1, 2) + 1

        site = (/ int(rand() * L) + 1,&
                 &int(rand() * L) + 1 /)

        ! convention in 2d: a plaquette is defined as
        ! (1,i,j),(2,i,j),(1,i,neg(j)),(2,neg(i),j)
        eo1 = e_field(mu1,site(1),site(2))
        eo2 = e_field(mu2,site(1),site(2))

        ! if e.g. we picked xy, the next line does x -> neg(x)
        site(mu1) = neg(site(mu1))
        eo4 = e_field(mu2,site(1),site(2))
        ! and now we need to put it back
        site(mu1) = pos(site(mu1))

        ! same for y
        site(mu2) = neg(site(mu2))
        eo3 = e_field(mu1,site(1),site(2))
        site(mu2) = pos(site(mu2))

        en1 = eo1 + increment
        en2 = eo2 - increment
        en3 = eo3 - increment
        en4 = eo4 + increment

        old_e = 0.5 * eps_0 * lambda**2 *&
                (eo1**2 + eo2**2 + eo3**2 + eo4**2)
        new_e = 0.5 * eps_0 * lambda**2 *&
                (en1**2 + en2**2 + en3**2 + en4**2)
        delta_e = new_e - old_e

        if ((delta_e.lt.0.0).or.(exp((-beta)*delta_e).gt.rand())) then

          e_field(mu1,site(1),site(2)) = en1
          e_field(mu2,site(1),site(2)) = en2

          site(mu1) = neg(site(mu1))
          e_field(mu2,site(1),site(2)) = en4
          site(mu1) = pos(site(mu1))

          site(mu2) = neg(site(mu2))
          e_field(mu1,site(1),site(2)) = en3
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

        increment = (q) / (L * eps_0 * lambda)

        do mu = 1,2

          if (rand() - 0.5.le.0) then
            pm1 = -1
          else
            pm1 = +1
          end if

          delta_e = (lambda**2 * q) *&
                    ((q / (2 * eps_0)) + (pm1 * L * ebar(mu)))

          if ((delta_e.lt.0.0).or.((exp(-beta*delta_e).gt.rand())&
            .and.(exp(-beta*delta_e).gt.0.00000000001))) then

            ebar(mu) = ebar(mu) + pm1 * increment
            e_field(mu,:,:) = e_field(mu,:,:) + pm1 * (increment)

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
      integer :: omp_index,i,j,k,kmax,n,kx,ky,m,p,s,x,y,dist_bin
      real(kind=8) :: norm_k, dist, ener_tot, ener_rot, ener_irrot
      real(kind=8) :: ener_tot_sq, ener_rot_sq, ener_irrot_sq, dp, np
      real(kind=8) :: theta_x, theta_y, vert_sum, diff, cos_k
      real(kind=8), dimension(2,L,L) :: e_rot
      complex(kind=rk) :: rho_k_p_temp, rho_k_m_temp
      complex(kind=rk) :: e_kx, e_ky
      complex(kind=rk) :: theta_kx, theta_ky
      complex(kind=rk) :: mnphi_kx, mnphi_ky
      complex(kind=rk) :: e_rot_kx, e_rot_ky
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
      norm_k = 0.0; dist = 0.0; vert_sum = 0.0
      rho_k_p_temp = (0.0,0.0); rho_k_m_temp = (0.0,0.0)
      e_kx = (0.0,0.0); mnphi_kx = (0.0,0.0)
      e_rot_kx = (0.0,0.0); e_rot = 0.0
      e_ky = (0.0,0.0); mnphi_ky = (0.0,0.0); e_rot_ky = (0.0,0.0)
      kdotx = (0.0,0.0); imag = (0.0, 1.0)
      n = step_number / sample_interval
      kmax = 100

      ! this is all *incredibly ugly* and needs a lot of refactoring
      ! my bad

      ! the sum of emergent fields around a given lattice site
      ! tells us whether there's a vortex there or not:
      ! we need to know this to Helmholtz decompose the field properly
      ! it looks like this:
      !                         |
      !                       - . -
      !                         |
      ! with the charge located at the point (i + 1/2, j + 1/2)
      ! since in the code I take the spins to be located on lattice points.

      do i = 1,L
        do j = 1,L

          ! michael's code measures top_x and top_y again so I do too;
          ! don't think it should be necessary in principle though
          diff = (theta(i,j) - theta(i,neg(j)))
          if (diff.ge.q/2) then
            do while (diff.ge.q/2)
              diff = diff - q
            end do
          else if (diff.le.-1.0*q/2) then
            do while (diff.le.-1.0*q/2)
              diff = diff + q
            end do
          end if
          top_x(i,j) = diff

          diff = -1.0*(theta(i,j) - theta(neg(i),j))
          if (diff.ge.q/2) then
            do while (diff.ge.q/2)
              diff = diff - q
            end do
          else if (diff.le.-1.0*q/2) then
            do while (diff.le.-1.0*q/2)
              diff = diff + q
            end do
          end if
          top_y(i,j) = diff

        end do
      end do

      do i = 1,L
        do j = 1,L

          vert_sum = top_x(i,j) + top_y(i,j) - top_x(neg(i),j) - top_y(i,neg(j))
          vert_sum = vert_sum / q

          ! this is just a check, shouldn't ever happen
          if (vert_sum.gt.2.0001.or.vert_sum.lt.-2.001) then
            write (*,*) n,i,j,vert_sum, top_x(i,j), top_y(i,j), top_x(neg(i),j), top_y(i,neg(j))
          end if

          ! considering the modular operation it's possible
          ! to have vertex sums arbitrarily close to 4 \pi;
          ! they're still just single charges though,
          ! so we have to add/subtract 1 if we find them
          if ((vert_sum.gt.1.999).and.(vert_sum.lt.2.001)) then
            vert_sum = vert_sum - 1.0
          else if ((vert_sum.lt.-1.999).and.(vert_sum.gt.-2.001)) then
            vert_sum = vert_sum + 1.0
          end if

          ! this should cover all bases since we've 
          ! already added/subtracted 1 where necessary
          if (vert_sum.gt.0.999.and.vert_sum.lt.1.001) then
            v(i,j) = 1
          else if (vert_sum.lt.-0.999.and.vert_sum.gt.-1.001) then
            v(i,j) = -1
          else
            v(i,j) = 0
          end if

          ! alternatively could just have something like 
          ! if (vert_sum.gt.0.999.or.vert_sum.lt.-0.999) maybe? 

        end do
      end do

      call linsol
      e_rot(1,:,:) = top_x + mnphi(1,:,:)
      e_rot(2,:,:) = top_y + mnphi(2,:,:)

      ebar(1) = sum(top_x(:,:))
      ebar(2) = sum(top_y(:,:))

      do i = 1,2

        dp = ebar(i)
        np = 0

        do while (abs(dp).gt.((q * dble(L)) / (2 * eps_0)))

          if (dp.gt.((q * dble(L)) / (2 * eps_0))) then
            windings(i,n) = windings(i,n) + 1
            dp = dp - (q * dble(L) / eps_0)
          else if (dp.le.((-1.0 * q * dble(L)) / (2 * eps_0))) then
            windings(i,n) = windings(i,n) - 1
            dp = dp + (q * dble(L) / eps_0)
          end if

        end do

        windings_sq(i,n) = windings(i,n)**2

        np = ebar(i) - dp

        ebar_dip(i) = dp
        ebar_wind(i) = np

      end do

      ebar = ebar / L**2
      ebar_dip = ebar_dip / L**2
      ebar_wind = ebar_wind / L**2

      ebar_sum(1) = ebar_sum(1) + ebar(1)
      ebar_sum(2) = ebar_sum(2) + ebar(2)
      ebar_sq_sum(1) = ebar_sq_sum(1) + (ebar(1) * ebar(1))
      ebar_sq_sum(2) = ebar_sq_sum(2) + (ebar(2) * ebar(2))

      ebar_dip_sum(1) = ebar_dip_sum(1) + ebar_dip(1)
      ebar_dip_sum(2) = ebar_dip_sum(2) + ebar_dip(2)
      ebar_dip_sq_sum(1) = ebar_dip_sq_sum(1) + (ebar_dip(1) * ebar_dip(1))
      ebar_dip_sq_sum(2) = ebar_dip_sq_sum(2) + (ebar_dip(2) * ebar_dip(2))

      ebar_wind_sum(1) = ebar_wind_sum(1) + ebar_wind(1)
      ebar_wind_sum(2) = ebar_wind_sum(2) + ebar_wind(2)
      ebar_wind_sq_sum(1) = ebar_wind_sq_sum(1) + (ebar_wind(1) * ebar_wind(1))
      ebar_wind_sq_sum(2) = ebar_wind_sq_sum(2) + (ebar_wind(2) * ebar_wind(2))

      ! try subtracting ebar here
      ! do i = 1,L
      !   do j = 1,L
      !     e_rot(1,i,j) = e_rot(1,i,j) - ebar(1)
      !     e_rot(2,i,j) = e_rot(2,i,j) - ebar(2)
      !     mnphi(1,i,j) = mnphi(1,i,j) + ebar(1)
      !     mnphi(2,i,j) = mnphi(2,i,j) + ebar(2)
      !   end do
      ! end do

      do i = 1,2

        avg_field_total(i) = avg_field_total(i) + sum(real(e_field(i,:,:)))
        avg_field_rot(i) = avg_field_rot(i) + sum(real(e_rot(i,:,:)))
        avg_field_irrot(i) = avg_field_irrot(i) + sum(real(mnphi(i,:,:)))

        avg_field_sq_total(i) = avg_field_sq_total(i) + avg_field_total(i)**2
        avg_field_sq_rot(i)   = avg_field_sq_rot(i)   + avg_field_rot(i)**2
        avg_field_sq_irrot(i) = avg_field_sq_irrot(i) + avg_field_irrot(i)**2

      end do

      avg_field_total = avg_field_total / L**2
      avg_field_rot = avg_field_rot / L**2
      avg_field_irrot = avg_field_irrot / L**2
      avg_field_sq_total = avg_field_sq_total/ L**2
      avg_field_sq_rot   = avg_field_sq_rot/ L**2
      avg_field_sq_irrot = avg_field_sq_irrot/ L**2

      e_tot_avg(1,:,:) = e_tot_avg(1,:,:) + top_x
      e_tot_avg(2,:,:) = e_tot_avg(2,:,:) + top_y
      e_rot_avg =           e_rot_avg + e_rot
      e_irrot_avg =         e_irrot_avg + mnphi
      v_avg =               v_avg + float(v)
      rho_avg =             rho_avg + (dble(sum(abs(v))) / L**2)
      ener_tot =            0.5 * lambda**2 * eps_0 * sum(top_x**2 + top_y**2)
      ener_rot =            0.5 * lambda**2 * eps_0 * sum(e_rot * e_rot)
      ener_irrot =          0.5 * lambda**2 * eps_0 * sum(mnphi * mnphi)
      ener_tot =            ener_tot / L**2
      ener_rot =            ener_rot / L**2
      ener_irrot =          ener_irrot / L**2
      ener_tot_sq =         ener_tot**2
      ener_rot_sq =         ener_rot**2
      ener_irrot_sq =       ener_irrot**2
      ener_tot_sum =        ener_tot_sum + ener_tot
      ener_rot_sum =        ener_rot_sum + ener_rot
      ener_irrot_sum =      ener_irrot_sum + ener_irrot
      ener_tot_sq_sum =     ener_tot_sq_sum + ener_tot_sq
      ener_rot_sq_sum =     ener_rot_sq_sum + ener_rot_sq
      ener_irrot_sq_sum =   ener_irrot_sq_sum + ener_irrot_sq

      if (do_corr) then

        !$omp parallel do&
        !$omp& private(i,j,k,m,p,s,kx,ky,rho_k_p_temp,rho_k_m_temp,e_kx,&
        !$omp& mnphi_kx,e_rot_kx,theta_x,theta_y,cos_k,&
        !$omp& theta_kx,theta_ky,e_ky,mnphi_ky,e_rot_ky,norm_k,kdotx)&
        !$omp& shared(dir_struc,s_ab,s_ab_rot,s_ab_irrot,s_ab_theta,&
        !$omp& dist_r,bin_count,eps_hxy)
        do omp_index = 1, L**2

          i = ((omp_index - 1) / (L)) + 1
          j = mod(omp_index - 1, (L)) + 1

          ! repurposing these as the s_ab indices at the end:
          ! we're gonna fill up the positive half of those arrays
          ! and then use symmetry properties to fill the rest at the end.
          ! k_mu = 0 point is at s_ab array index L + 1
          kx = i + L
          ky = j + L

          rho_k_p_temp = (0.0,0.0); rho_k_m_temp = (0.0,0.0)
          mnphi_kx = (0.0,0.0); mnphi_ky = (0.0,0.0)
          e_rot_kx = (0.0,0.0); e_rot_ky = (0.0,0.0)
          norm_k = 0.0

          e_kx = (0.0,0.0); e_ky = (0.0,0.0)
          kdotx = (0.0,0.0)
          theta_x = 0.0; theta_y = 0.0
          theta_kx = (0.0,0.0); theta_ky = (0.0,0.0)

          if (kx.eq.0.and.ky.eq.0) then
            norm_k = 0.0
          else
            norm_k = 1.0/(((2*pi/(L*lambda))**2)*dble(kx**2 + ky**2))
          end if

          cos_k = 0.0
          do s = 1,L**2
            m = ((s - 1) / L) + 1
            p = mod(s - 1, L) + 1

            ! we need to calculate this thing to get susceptibilties
            ! NB!!!!! if do_corr is turned off, this will be zero
            ! and the susceptibilities won't mean all that much
            if (omp_index.le.kmax) then
              cos_k = cos_k + cos(omp_index * top_x(m,p)) +&
                              cos(omp_index * top_y(m,p))
              eps_hxy = eps_hxy + (-1)**(omp_index + 1) * cos_k
            end if


              kdotx = fw(m,i) + hw(p,j)

              ! e_kx_temp = e_kx_temp + exp(kdotx)*e_field(1,m,p)
              e_kx = e_kx + exp(kdotx)*top_x(m,p)
              mnphi_kx = mnphi_kx + exp(kdotx)*mnphi(1,m,p)
              e_rot_kx = e_rot_kx + exp(kdotx)*e_rot(1,m,p)

              kdotx = hw(m,i) + fw(p,j)

              ! e_ky =     e_ky     + exp(kdotx)*e_field(2,m,p)
              e_rot_ky = e_rot_ky + exp(kdotx)*e_rot(2,m,p)
              mnphi_ky = mnphi_ky + exp(kdotx)*mnphi(2,m,p)
              e_ky = e_ky + exp(kdotx)*top_y(m,p)

              kdotx = fw(m,i) + fw(p,j)
              theta_x = cos(theta(m,p))
              theta_y = sin(theta(m,p))
              theta_kx = theta_kx + exp(kdotx) * theta_x
              theta_ky = theta_ky + exp(kdotx) * theta_y

              if (v(m,p).ne.0) then ! calculate <++ + +->!

                kdotx = hw(m,i) + hw(p,j)
               
                if (v(m,p).eq.-1) then
                  rho_k_m_temp = rho_k_m_temp + q * v(m,p) * exp(kdotx)
                end if

                if (v(m,p).eq.1) then ! take away <++>
                  rho_k_p_temp = rho_k_p_temp + q * v(m,p) * exp(kdotx)
                end if

              end if

              !   ! --- real space correlation function ---

              !   if (v(m,p).eq.1) then
              !     if (v(i,j).eq.-1) then

              !       x = abs(m - i)
              !       y = abs(p - j)

              !       if (x.gt.L/2) then
              !         x = L - x
              !       end if
              !       if (y.gt.L/2) then
              !         y = L - y
              !       end if

              !       x = x + 1
              !       y = y + 1

              !       dir_struc(x,y) = dir_struc(x,y) +&
              !               v(m,p) * v(i,j)
              !     end if ! neg charge at i,j
              !   end if ! pos charge at m,p

              ! end if ! end v != 0 check

          end do ! s

          rho_k_p_temp = rho_k_p_temp / float(L**2)
          rho_k_m_temp = rho_k_m_temp / float(L**2)
          e_kx         = e_kx / float(2 * L**2)
          e_ky         = e_ky / float(2 * L**2)
          theta_kx     = theta_kx / float(2 * L**2)
          theta_ky     = theta_ky / float(2 * L**2)
          mnphi_kx     = mnphi_kx / float(2 * L**2)
          mnphi_ky     = mnphi_ky / float(2 * L**2)
          e_rot_kx     = e_rot_kx / float(2 * L**2)
          e_rot_ky     = e_rot_ky / float(2 * L**2)

          rho_k_p(kx,ky) = rho_k_p(kx,ky) + rho_k_p_temp
          rho_k_m(kx,ky) = rho_k_m(kx,ky) + rho_k_m_temp
          ch_ch(kx,ky) = ch_ch(kx,ky) +&
                        (rho_k_p_temp * conjg(rho_k_m_temp))

          s_ab(1,1,kx,ky) = s_ab(1,1,kx,ky) +&
          e_kx*conjg(e_kx)
          s_ab(1,2,kx,ky) = s_ab(1,2,kx,ky) +&
          e_kx*conjg(e_ky)
          s_ab(2,1,kx,ky) = s_ab(2,1,kx,ky) +&
          e_ky*conjg(e_kx)
          s_ab(2,2,kx,ky) = s_ab(2,2,kx,ky) +&
          e_ky*conjg(e_ky)

          s_ab_theta(1,1,kx,ky) = s_ab_theta(1,1,kx,ky) +&
          theta_kx*conjg(theta_kx)
          s_ab_theta(1,2,kx,ky) = s_ab_theta(1,2,kx,ky) +&
          theta_kx*conjg(theta_ky)
          s_ab_theta(2,1,kx,ky) = s_ab_theta(2,1,kx,ky) +&
          theta_ky*conjg(theta_kx)
          s_ab_theta(2,2,kx,ky) = s_ab_theta(2,2,kx,ky) +&
          theta_ky*conjg(theta_ky)

          s_ab_rot(1,1,kx,ky) = s_ab_rot(1,1,kx,ky) +&
          e_rot_kx*conjg(e_rot_kx)
          s_ab_rot(1,2,kx,ky) = s_ab_rot(1,2,kx,ky) +&
          e_rot_kx*conjg(e_rot_ky)
          s_ab_rot(2,1,kx,ky) = s_ab_rot(2,1,kx,ky) +&
          e_rot_ky*conjg(e_rot_kx)
          s_ab_rot(2,2,kx,ky) = s_ab_rot(2,2,kx,ky) +&
          e_rot_ky*conjg(e_rot_ky)

          s_ab_irrot(1,1,kx,ky) = s_ab_irrot(1,1,kx,ky) +&
          mnphi_kx*conjg(mnphi_kx)
          s_ab_irrot(1,2,kx,ky) = s_ab_irrot(1,2,kx,ky) +&
          mnphi_kx*conjg(mnphi_ky)
          s_ab_irrot(2,1,kx,ky) = s_ab_irrot(2,1,kx,ky) +&
          mnphi_ky*conjg(mnphi_kx)
          s_ab_irrot(2,2,kx,ky) = s_ab_irrot(2,2,kx,ky) +&
          mnphi_ky*conjg(mnphi_ky)

          ! direct space one
          ! if (i.le.L/2+1.and.j.le.L/2+1) then

          !   dist = sqrt(dble((i - 1)**2 + (j - 1)**2))
          !   dist_bin = floor(dist / bin_size) + 1
          !   dist_r(dist_bin) = dist_r(dist_bin) + dir_struc(i,j)
          !   bin_count(dist_bin) = bin_count(dist_bin) + 1

          ! end if

        end do ! end openmp_index loop
        !$omp end parallel do

      end if ! if do_corr

    end subroutine measure

end module update
