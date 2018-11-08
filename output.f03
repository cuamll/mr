module output
  use common
  implicit none

  contains

  subroutine snapshot(n)
    use common
    use linear_solver
    implicit none
    integer, intent(in) :: n
    integer :: i, j, k
    character(len=20) :: nsl
    character(len=200) :: tfl, ifl, rfl, cfl
    character(:), allocatable :: tot_file, irrot_file, rot_file, charge_file, ns
    real(kind=8), dimension(3,L,L,L) :: e_rot

    ! these snapshots are designed to be plotted using matplotlib's built-in
    ! quiver plots; the documentation for them is available online
    ! take an integer as input; this will be the current step number,
    ! just in case we want to take multiple snapshots

    ! construct file name base: cast step number to string
    write(nsl, '(i20.1)') n
    ns = trim(adjustl(nsl))
    tfl = e_field_file // '_' // ns // '_total.dat'
    ifl = e_field_file // '_' // ns // '_irrot.dat'
    rfl = e_field_file // '_' // ns // '_rot.dat'
    cfl = e_field_file // '_' // ns // '_charges.dat'
    tot_file = trim(adjustl(tfl))
    irrot_file = trim(adjustl(ifl))
    rot_file = trim(adjustl(rfl))
    charge_file = trim(adjustl(cfl))
    ! make sure we have a current Helmholtz decomposition
    call linsol
    e_rot = e_field - mnphi

    open(unit=10, file=tot_file)
    open(unit=11, file=irrot_file)
    open(unit=12, file=rot_file)
    open(unit=13, file=charge_file)

    do i = 1, L
      do j = 1, L
        do k = 1, L

          ! x component
          write (10,'(6f16.12)')  real(i + 0.5), real(j), real(k),&
                                  e_field(1,i,j,k), 0.0, 0.0
          write (11,'(6f16.12)')  real(i + 0.5), real(j), real(k),&
                                  mnphi(1,i,j,k), 0.0, 0.0
          write (12,'(6f16.12)')  real(i + 0.5), real(j), real(k),&
                                  e_rot(1,i,j,k), 0.0, 0.0
          ! y component
          write (10,'(6f16.12)')  real(i), real(j + 0.5), real(k),&
                                  0.0, e_field(2,i,j,k), 0.0
          write (11,'(6f16.12)')  real(i), real(j + 0.5), real(k),&
                                  0.0, mnphi(2,i,j,k), 0.0
          write (12,'(6f16.12)')  real(i), real(j + 0.5), real(k),&
                                  0.0, e_rot(2,i,j,k), 0.0

          ! z-component
          write (10,'(6f16.12)')  real(i), real(j), real(k + 0.5),&
                                  0.0, 0.0, e_field(3,i,j,k)
          write (11,'(6f16.12)')  real(i), real(j), real(k + 0.5),&
                                  0.0, 0.0, mnphi(3,i,j,k)
          write (12,'(6f16.12)')  real(i), real(j), real(k + 0.5),&
                                  0.0, 0.0, e_rot(3,i,j,k)

          ! charges
          write (13,'(3f16.12, i3.1)') real(i), real(j), real(k), v(i,j,k)

        end do
      end do
    end do

    close(10)
    close(11)
    close(12)
    close(13)

  end subroutine snapshot

  subroutine write_output

    call fix_fftw
    call calc_correlations

  end subroutine write_output


  subroutine fix_fftw
    use common
    implicit none
    integer :: i, j, j_eff, k, k_eff, ri, rj, rkk,&
    kx, ky, kz, pmx, pmy, pmz, offdiag

    ! normalisation for correct equipartition result
    s_ab = s_ab * L**3

    ! fftw only does up to pi in the k_x direction: fix that first
    do i = 1, L/2 + 1
      do j = 1, L + 1
        do k = 1, L + 1

          if (j.eq.L+1) then
            j_eff = 1
          else
            j_eff = j
          end if
          if (k.eq.L+1) then
            k_eff = 1
          else
            k_eff = k
          end if

          ri = ((2*L) + 2) - i
          rj = ((2*L) + 2) - j
          rkk = ((2*L) + 2) - k

          s_ab_large(1,1,i+L,j+L,k+L)               = s_ab(i, j_eff, k_eff, 1)
          s_ab_large(1,1,ri,rj,rkk)                 = s_ab(i, j_eff, k_eff, 1)

          s_ab_large(2,2,i+L,j+L,k+L)               = s_ab(i, j_eff, k_eff, 4)
          s_ab_large(2,2,ri,rj,rkk)                 = s_ab(i, j_eff, k_eff, 4)

          s_ab_large(3,3,i+L,j+L,k+L)               = s_ab(i, j_eff, k_eff, 6)
          s_ab_large(3,3,ri,rj,rkk)                 = s_ab(i, j_eff, k_eff, 6)

          ! off diagonals are more annoying
          s_ab_large(1,2,i+L,j+L,k+L)               = s_ab(i, j_eff, k_eff, 2)
          s_ab_large(1,2,ri,rj,rkk)                 = s_ab(i, j_eff, k_eff, 2)

          s_ab_large(1,3,i+L,j+L,k+L)               = s_ab(i, j_eff, k_eff, 3)
          s_ab_large(1,3,ri,rj,rkk)                 = s_ab(i, j_eff, k_eff, 3)

          s_ab_large(2,3,i+L,j+L,k+L)               = s_ab(i, j_eff, k_eff, 5)
          s_ab_large(2,3,ri,rj,rkk)                 = s_ab(i, j_eff, k_eff, 5)

          ! should really check this?
          if (j.eq.1.and.i.eq.1.and.k.eq.1) then

            s_ab_large(1,2,ri,rj,rkk)  = (+1.0)*&
            s_ab_large(1,2,ri,rj,rkk)
            s_ab_large(1,3,ri,rj,rkk)  = (+1.0)*&
            s_ab_large(1,3,ri,rj,rkk)
            s_ab_large(2,3,ri,rj,rkk)  = (+1.0)*&
            s_ab_large(2,3,ri,rj,rkk)

          end if

          if (j.eq.L+1.or.k.eq.L+1) then

            s_ab_large(1,2,i+L,j+L,k+L)  = (+1.0)*&
            s_ab_large(1,2,i+L,j+L,k+L)
            s_ab_large(1,3,i+L,j+L,k+L)  = (+1.0)*&
            s_ab_large(1,3,i+L,j+L,k+L)
            s_ab_large(2,3,i+L,j+L,k+L)  = (+1.0)*&
            s_ab_large(2,3,i+L,j+L,k+L)

          end if
          if (j.eq.L+1.or.k.eq.L+1) then

            s_ab_large(1,2,ri,rj,rkk)  = (+1.0)*&
            s_ab_large(1,2,ri,rj,rkk)
            s_ab_large(1,3,ri,rj,rkk)  = (+1.0)*&
            s_ab_large(1,3,ri,rj,rkk)
            s_ab_large(2,3,ri,rj,rkk)  = (+1.0)*&
            s_ab_large(2,3,ri,rj,rkk)

          end if

        end do
      end do
    end do
    
    ! now out to the other quadrants
    ! we only have the quadrant [kx, ky] \in [0,2pi]
    do kx = -L, L
      do ky = -L, L
        do kz = -L, L

          pmx = 0; pmy = 0; pmz = 0; offdiag = 1
          i = kx + 1 + L
          j = ky + 1 + L
          k = kz + 1 + L

          if (kx.lt.0) then
            pmx = +1
            offdiag = -1 * offdiag
          end if

          if (ky.lt.0) then
            pmy = +1
            offdiag = -1 * offdiag
          end if

          if (kz.lt.0) then
            pmz = +1
            offdiag = -1 * offdiag
          end if

          ! these lines feel a bit off/hacky to me, but they
          ! exactly reproduce what I had before, which in turn
          ! came from directly simulating out to some huge k.
          ! Might be worth rechecking at some point though.
          if ((kx).eq.-L.or.(ky).eq.-L.or.(kz).eq.-L) then
            offdiag = -1
          end if
          if ((kx.eq.L.and.ky.lt.0.and.kz.lt.0).or.&
            (ky.eq.L.and.kx.lt.0.and.kz.lt.0).or.&
            (kz.eq.L.and.ky.lt.0.and.kx.lt.0)) then
            offdiag = +1
          end if
          
          s_ab_large(1,1,i,j,k) =&
          s_ab_large(1,1, i + pmx * L, j + pmy * L, k + pmz * L)
          s_ab_large(2,2,i,j,k) =&
          s_ab_large(2,2, i + pmx * L, j + pmy * L, k + pmz * L)
          s_ab_large(3,3,i,j,k) =&
          s_ab_large(3,3, i + pmx * L, j + pmy * L, k + pmz * L)

          s_ab_large(1,2,i,j,k) = offdiag *&
          s_ab_large(1,2, i + pmx * L, j + pmy * L, k + pmz * L)
          s_ab_large(1,3,i,j,k) = offdiag *&
          s_ab_large(1,3, i + pmx * L, j + pmy * L, k + pmz * L)
          s_ab_large(2,3,i,j,k) = offdiag *&
          s_ab_large(2,3, i + pmx * L, j + pmy * L, k + pmz * L)

        end do
      end do
    end do

    s_ab_large(2,1,:,:,:) = s_ab_large(1,2,:,:,:)
    s_ab_large(3,1,:,:,:) = s_ab_large(1,3,:,:,:)
    s_ab_large(3,2,:,:,:) = s_ab_large(2,3,:,:,:)

  end subroutine fix_fftw

  subroutine calc_correlations
    use common
    implicit none
    integer :: i, j, k, kx, ky, kz, m, p, s, x, y, z, dist_bin, sp
    real(kind=rk) :: norm_k, kx_float, ky_float, kz_float, dist,&
    sp_he_tot, sp_he_rot, sp_he_irrot, prefac,&
    ebar_sus, ebar_dip_sus, ebar_wind_sus
    real(kind=rk), dimension(:,:,:), allocatable :: s_perp,&
    s_perp_irrot, s_perp_rot
    real(kind=rk), dimension(:,:,:), allocatable :: s_par,&
    s_par_irrot, s_par_rot
    real(kind=rk), dimension((bz*L)+1,(bz*L)+1,(bz*L)+1) :: charge_struc,&
    field_struc, field_struc_irrot, field_struc_rot
    real(kind=rk), dimension(3,3,(bz*L)+1,(bz*L)+1,(bz*L)+1) :: chi_ab,&
    chi_ab_rot, chi_ab_irrot
    character(100) :: struc_format_string, field_format_string,&
    vertex_format, avg_field_format,dir_format_string, dir_dist_format_string

    avg_field_format = "(3I3, 7ES18.9)"
    field_format_string = "(12ES18.9)"
    struc_format_string = "(4ES18.9)"
    dir_format_string = "(3I3, ES18.9)"
    vertex_format = "(I6, 3I3, 6ES18.9, I3)"
    dir_dist_format_string = "(ES18.9, I9.6, ES18.9)"

    field_struc = 0.0; field_struc_rot = 0.0;
    charge_struc = 0.0; field_struc_irrot = 0.0;

    prefac = 1.0 * L**3 / (temp**2)

    sp_he_tot = prefac * (ener_tot_sq_sum - (ener_tot_sum)**2)
    sp_he_rot = prefac * (ener_rot_sq_sum - (ener_rot_sum)**2)
    sp_he_irrot = prefac * (ener_irrot_sq_sum - (ener_irrot_sum)**2)
    ebar_sus = L**3 * beta * (sum(ebar_sq_sum) - (sum(ebar_sum))**2)
    ebar_dip_sus = L**3 * beta *&
                   (sum(ebar_dip_sq_sum) - (sum(ebar_dip_sum))**2)
    ebar_wind_sus = L**3 * beta *&
                  (sum(ebar_wind_sq_sum) + (sum(ebar_wind_sum))**2)

    write(*,*)
    write (*,'(a)') "# Specific heat: total, rot., irrot."
    write (*,'(ES18.9, ES18.9, ES18.9, ES18.9)') temp,&
    sp_he_tot, sp_he_rot, sp_he_irrot
    write(*,*) "E^bar averages:"
    write(*,'(a)') "<E^bar_x>, <E^bar_y>, <|E^bar|>"
    write (*,*) ebar_sum(1), ebar_sum(2),ebar_sum(3),&
    sqrt(ebar_sum(1)**2+ebar_sum(2)**2 + ebar_sum(3)**2)
    write(*,'(a)') "<(E^bar_x)^2>, <(E^bar_y)^2>, <|(E^bar)^2|>"
    write (*,*) ebar_sq_sum(1), ebar_sq_sum(2), ebar_sq_sum(3),&
    sqrt(ebar_sq_sum(1)**2+ebar_sq_sum(2)**2+ebar_sq_sum(3)**2)

    write(*,'(a)') "<|(E^bar)^2|> - <|E^bar|>^2"
    write (*,*) L**3 * beta * (sum(ebar_sq_sum) - (sum(ebar_sum))**2)

    ! open  (30, file=sphe_sus_file)
    ! write (30,'(a)') "   # Specific heat: total, rot., irrot."
    ! write (30,'(4ES18.9)') temp,&
    ! sp_he_tot, sp_he_rot, sp_he_irrot

    ! write (30,'(a)') "   # T, Chi_{Ebar}, Chi_{Ebar_dip}, Chi_{Ebar_wind}"
    ! write (30,'(4ES18.9)') temp, ebar_sus, ebar_dip_sus, ebar_wind_sus
    ! write (30,'(a)') "# Avg. x-component: total, rot., irrot.:"
    ! write (30, '(3f18.10)') avg_field_total(1), avg_field_rot(1),&
    !                         avg_field_irrot(1)
    ! write (30,'(a)') "# Avg. y-component: total, rot., irrot.:"
    ! write (30, '(3f18.10)') avg_field_total(2), avg_field_rot(2),&
    !                         avg_field_irrot(2)
    ! write (30,'(a)') "# Avg. z-component: total, rot., irrot.:"
    ! write (30, '(3f18.10)') avg_field_total(3), avg_field_rot(3),&
    !                         avg_field_irrot(3)
    ! write (30,'(a)') "# Avg. x-component^2: total, rot., irrot.:"
    ! write (30, '(3es18.10)') avg_field_sq_total(1), avg_field_sq_rot(1),&
    !                         avg_field_sq_irrot(1)
    ! write (30,'(a)') "# Avg. y-component^2: total, rot., irrot.:"
    ! write (30, '(3es18.10)') avg_field_sq_total(2), avg_field_sq_rot(2),&
    !                         avg_field_sq_irrot(2)
    ! write (30,'(a)') "# Avg. z-component^2: total, rot., irrot.:"
    ! write (30, '(3es18.10)') avg_field_sq_total(3), avg_field_sq_rot(3),&
    !                         avg_field_sq_irrot(3)

    ! write (30,'(a)') "   # hop acceptance"
    ! write (30,'(2ES18.9)') temp,&
    ! (dble(accepts(1)) / dble(attempts(1)))
    ! write (30,'(a,2i12.1,es18.9)') "# Hops: total, attempts, rate: ",&
    ! accepts(1), attempts(1), dble(accepts(1)) / dble(attempts(1))
    ! write (30,'(a,2i12.1,es18.9)') "# Rot.: total, attempts, rate: ",&
    ! accepts(2), attempts(2), dble(accepts(2)) / dble(attempts(2))
    ! write (30,'(a,2i12.1,es18.9)') "# Harm: total, attempts, rate: ",&
    ! accepts(3), attempts(3), dble(accepts(3)) / dble(attempts(3))
    ! write (30,'(a,2i12.1,es18.9)') "# Creations: total, attempts, rate: ",&
    ! accepts(4), attempts(4), dble(accepts(4)) / dble(attempts(4))
    ! write (30,'(a,2i12.1,es18.9)') "# Annihilations: total, attempts, rate: ",&
    ! accepts(5), attempts(5), dble(accepts(5)) / dble(attempts(5))

    ! close(30)

    if (do_corr) then
      sp = 8

      do s = (-L/2)*sp,(L/2)*sp
        do p = (-L/2)*sp,(L/2)*sp
          do m = (-L/2)*sp,(L/2)*sp

            i = m + 1 + sp*(L/2)
            j = p + 1 + sp*(L/2)
            k = s + 1 + sp*(L/2)

            if (k.le.(bz*L + 1).and.j.le.(bz*L + 1).and.i.le.(bz*L + 1)) then
              ! can also subtract e.g. rho_k_p * conjg(rho_k_m)
              charge_struc(i,j,k) = abs(ch_ch(i,j,k) - &
                rho_k_p(i,j,k) * conjg(rho_k_m(i,j,k)))
              ! field_struc(i,j,k) = abs(s_ab_large(1,1,i,j,k))
              ! field_struc_irrot(i,j,k) = abs(s_ab_irrot(1,1,i,j,k))
              ! field_struc_rot(i,j,k) = abs(s_ab_rot(1,1,i,j,k))
            end if

            ! use separate variables, we're gonna mess around with values
            kx_float = m * ((2 * pi)/(L * lambda))
            ky_float = p * ((2 * pi)/(L * lambda))
            kz_float = s * ((2 * pi)/(L * lambda))

            if (kx_float.eq.0.and.ky_float.eq.0.and.kz_float.eq.0) then
              norm_k = 0.0
            else
              norm_k = 1.0/(kx_float**2 + ky_float**2 + kz_float**2)
            end if

            if (abs(s).gt.(L)*bz) then
              kz = modulo(s,(bz*L)) + 1 + (bz*L)
            else
              kz = s + 1 + (bz*L)
            end if
            if (abs(p).gt.(L)*bz) then
              ky = modulo(p,(bz*L)) + 1 + (bz*L)
            else
              ky = p + 1 + (bz*L)
            end if
            if (abs(m).gt.(L)*bz) then
              kx = modulo(m,(bz*L)) + 1 + (bz*L)
            else
              kx = m + 1 + (bz*L)
            end if

            if (s.lt.(-1*(bz*(L/2)))) then
              kz = s
              do while (kz.lt.(-1*(bz*(L/2))))
                kz = kz + bz*(L)
              end do
              ! array index
              kz = kz + 1 + bz*(L/2)
            else if (s.gt.(bz*(L/2))) then
              kz = s
              do while (kz.gt.bz*(L/2))
                kz = kz - bz*(L)
              end do
              kz = kz + 1 + bz*(L/2)
            else
              kz = s + 1 + (bz*(L/2))
            end if

            if (p.lt.(-1*(bz*(L/2)))) then
              ky = p
              do while (ky.lt.(-1*(bz*(L/2))))
                ky = ky + bz*(L)
              end do
              ! array index
              ky = ky + 1 + bz*(L/2)
            else if (p.gt.(bz*(L/2))) then
              ky = p
              do while (ky.gt.bz*(L/2))
                ky = ky - bz*(L)
              end do
              ky = ky + 1 + bz*(L/2)
            else
              ky = p + 1 + (bz*(L/2))
            end if

            if (m.lt.(-1*(bz*(L/2)))) then
              kx = m
              do while (kx.lt.(-1*(bz*(L/2))))
                kx = kx + bz*(L)
              end do
              ! array index
              kx = kx + 1 + bz*(L/2)
            else if (m.gt.(bz*(L/2))) then
              kx = m
              do while (kx.gt.bz*(L/2))
                kx = kx - bz*(L)
              end do
              kx = kx + 1 + bz*(L/2)
            else
              kx = m + 1 + (bz*(L/2))
            end if

          end do
        end do
      end do ! end p, m loops

      open(unit=12, file=charge_st_file)
      open(unit=14, file=s_ab_file)
      open(unit=17, file=irrot_sab_file)
      open(unit=20, file=rot_sab_file)
      open(unit=15, file=s_perp_file)
      open(unit=18, file=irrot_sperp_file)
      open(unit=21, file=rot_sperp_file)
      open(unit=22, file=spar_file)
      open(unit=23, file=irrot_spar_file)
      open(unit=24, file=rot_spar_file)
      open(unit=10, file=dir_st_file)
      open(unit=11, file=dir_dist_file)
      open(unit=25, file=avg_field_file)
      open(unit=26, file=chi_ab_file)
      open(unit=27, file=rot_chi_ab_file)
      open(unit=28, file=irrot_chi_ab_file)
      open(unit=38, file=windings_file)
      open(unit=39, file=windings_sq_file)

      open  (30, file=sphe_sus_file, position='append')
      write (30, '(a)') "# S_ab integrals (* L**2)!"
      write (30, '(a)') "# S_xx: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(s_ab_large(1,1,:,:,:))),&
      sum(real(s_ab_rot_large(1,1,:,:,:))), sum(real(s_ab_irrot_large(1,1,:,:,:)))
      write (30, '(a)') "# S_xy: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(s_ab_large(1,2,:,:,:))),&
      sum(real(s_ab_rot_large(1,2,:,:,:))), sum(real(s_ab_irrot_large(1,2,:,:,:)))
      write (30, '(a)') "# S_xz: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(s_ab_large(1,3,:,:,:))),&
      sum(real(s_ab_rot_large(1,3,:,:,:))), sum(real(s_ab_irrot_large(1,2,:,:,:)))
      write (30, '(a)') "# S_yx: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(s_ab_large(2,1,:,:,:))),&
      sum(real(s_ab_rot_large(2,1,:,:,:))), sum(real(s_ab_irrot_large(2,1,:,:,:)))
      write (30, '(a)') "# S_yy: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(s_ab_large(2,2,:,:,:))),&
      sum(real(s_ab_rot_large(2,2,:,:,:))), sum(real(s_ab_irrot_large(2,2,:,:,:)))
      write (30, '(a)') "# S_yz: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(s_ab_large(2,3,:,:,:))),&
      sum(real(s_ab_rot_large(2,3,:,:,:))), sum(real(s_ab_irrot_large(2,3,:,:,:)))
      write (30, '(a)') "# S_zx: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(s_ab_large(3,1,:,:,:))),&
      sum(real(s_ab_rot_large(3,1,:,:,:))), sum(real(s_ab_irrot_large(2,1,:,:,:)))
      write (30, '(a)') "# S_zy: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(s_ab_large(3,2,:,:,:))),&
      sum(real(s_ab_rot_large(3,2,:,:,:))), sum(real(s_ab_irrot_large(3,2,:,:,:)))
      write (30, '(a)') "# S_zz: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(s_ab_large(3,3,:,:,:))),&
      sum(real(s_ab_rot_large(3,3,:,:,:))), sum(real(s_ab_irrot_large(3,3,:,:,:)))
      write (30, '(a)') "# S_perp integrals: total, rot, irrot"
      write (30, '(3f20.8)') sum(s_perp), sum(s_perp_rot), sum(s_perp_irrot)


      write (*,'(f20.8)') sum(real(chi_ab_rot(1,1,:,:,:)))
      write (30, '(a)') "# Chi_ab integrals (s_ab / temp)!"
      write (30, '(a)') "# Chi_xx: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(chi_ab(1,1,:,:,:))),&
      sum(real(chi_ab_rot(1,1,:,:,:))), sum(real(chi_ab_irrot(1,1,:,:,:)))
      write (30, '(a)') "# Chi_xy: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(chi_ab(1,2,:,:,:))),&
      sum(real(chi_ab_rot(1,2,:,:,:))), sum(real(chi_ab_irrot(1,2,:,:,:)))
      write (30, '(a)') "# Chi_xz: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(chi_ab(1,3,:,:,:))),&
      sum(real(chi_ab_rot(1,3,:,:,:))), sum(real(chi_ab_irrot(1,3,:,:,:)))
      write (30, '(a)') "# Chi_yx: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(chi_ab(2,1,:,:,:))),&
      sum(real(chi_ab_rot(2,1,:,:,:))), sum(real(chi_ab_irrot(2,1,:,:,:)))
      write (30, '(a)') "# Chi_yy: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(chi_ab(2,2,:,:,:))),&
      sum(real(chi_ab_rot(2,2,:,:,:))), sum(real(chi_ab_irrot(2,2,:,:,:)))
      write (30, '(a)') "# Chi_yz: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(chi_ab(2,3,:,:,:))),&
      sum(real(chi_ab_rot(2,3,:,:,:))), sum(real(chi_ab_irrot(2,3,:,:,:)))
      write (30, '(a)') "# Chi_zx: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(chi_ab(3,1,:,:,:))),&
      sum(real(chi_ab_rot(3,1,:,:,:))), sum(real(chi_ab_irrot(3,1,:,:,:)))
      write (30, '(a)') "# Chi_zy: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(chi_ab(3,2,:,:,:))),&
      sum(real(chi_ab_rot(3,2,:,:,:))), sum(real(chi_ab_irrot(3,2,:,:,:)))
      write (30, '(a)') "# Chi_zz: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(chi_ab(3,3,:,:,:))),&
      sum(real(chi_ab_rot(3,3,:,:,:))), sum(real(chi_ab_irrot(3,3,:,:,:)))
      ! write (30, '(a)') "# Chi_perp integrals: total, rot, irrot"
      ! write (30, '(3f20.8)') sum(chi_perp), sum(chi_perp_rot), sum(chi_perp_irrot)
      close (30)

      !dist_r = dist_r / (no_measurements)
      do i = 1,ceiling( sqrt(float((3*((L/2)**3)))) * (1 / bin_size) )
        write (11, dir_dist_format_string)&
        i * bin_size, bin_count(i), abs(dist_r(i))
      end do

      do i = 1,no_measurements
        write (38,*) windings(1,i), windings(2,i), windings(3,i)
        write (39,*) windings_sq(1,i), windings_sq(2,i), windings_sq(3,i)
      end do

      close(11)
      close(38)
      close(39)

      do i = 1,sp*(L) + 1
        do j = 1,sp*(L) + 1
          do k = 1,sp*(L) + 1

            ! output is kx, ky, kz, S(\vec{k})

            if (i.le.L/2+1.and.j.le.L/2+1.and.k.le.L/2+1) then
              write (10, dir_format_string)&
              i - 1,j - 1,k - 1,abs(dir_struc(i,j,k))
            end if

            if (i.le.L.and.j.le.L.and.k.eq.L) then
              write (25, avg_field_format)&
              i,j,k,e_tot_avg(1,i,j,k),e_tot_avg(2,i,j,k),&
              e_rot_avg(1,i,j,k),e_rot_avg(2,i,j,k),&
              e_irrot_avg(1,i,j,k),e_irrot_avg(2,i,j,k),&
              v_avg(i,j,k)
            end if

            if (k.le.(bz*L + 1).and.j.le.(bz*L + 1).and.i.le.(bz*L + 1)) then

              write (12, struc_format_string)&
              2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
              2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
              2*pi*(k - 1 - bz*(L/2))/(L*lambda),&
              charge_struc(i,j,k)

              write (14, field_format_string)&
              2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
              2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
              2*pi*(k - 1 - bz*(L/2))/(L*lambda),&
              real(s_ab_large(1,1,i,j,k)),&
              real(s_ab_large(1,2,i,j,k)),&
              real(s_ab_large(1,3,i,j,k)),&
              real(s_ab_large(2,1,i,j,k)),&
              real(s_ab_large(2,2,i,j,k)),&
              real(s_ab_large(2,3,i,j,k)),&
              real(s_ab_large(3,1,i,j,k)),&
              real(s_ab_large(3,2,i,j,k)),&
              real(s_ab_large(3,3,i,j,k))

            end if

          end do
        end do
      end do

      close(10)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)
      close(19)
      close(20)
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)

    end if ! do_corr

  end subroutine calc_correlations

end module output
