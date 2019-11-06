module output
  use common
  implicit none

  contains

  subroutine snapshot(n)
    use common
    use linear_solver
    implicit none
    integer, intent(in) :: n
    integer :: i, j
    character(len=20) :: nsl
    character(len=200) :: tfl, ifl, rfl, cfl
    character(:), allocatable :: tot_file, irrot_file, rot_file, charge_file, ns
    real(kind=8), dimension(2,L,L) :: e_rot

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
        
        ! x component
        write (10,'(4f18.10)') real(i + 0.5), real(j), e_field(1,i,j), 0.0
        write (11,'(4f18.10)') real(i + 0.5), real(j), mnphi(1,i,j), 0.0
        write (12,'(4f18.10)') real(i + 0.5), real(j), e_rot(1,i,j), 0.0
        ! y component
        write (10,'(4f18.10)') real(i), real(j + 0.5), 0.0, e_field(2,i,j)
        write (11,'(4f18.10)') real(i), real(j + 0.5), 0.0, mnphi(2,i,j)
        write (12,'(4f18.10)') real(i), real(j + 0.5), 0.0, e_rot(2,i,j)

        ! charges
        write (13,'(2f18.12, i3.1)') real(i), real(j), v(i,j)

      end do
    end do

    close(10)
    close(11)
    close(12)
    close(13)

  end subroutine snapshot

  subroutine write_output

    if (do_corr) then
      call fix_fftw
    end if

    call calc_correlations

  end subroutine write_output

  subroutine fix_fftw
    use common
    implicit none
    integer :: i, j, j_eff, ri, rj, kx, ky, pmx, pmy, offdiag

    ! normalisation for correct equipartition result
    s_ab = s_ab * L**2

    ! fftw only does up to pi in the k_x direction: fix that first
    do i = 1, L/2 + 1
      do j = 1, L + 1

        if (j.eq.L+1) then
          j_eff = 1
        else
          j_eff = j
        end if

        ri = ((2*L) + 2) - i
        rj = ((2*L) + 2) - j

        s_ab_large(1,1,i+L,j+L)              = s_ab(i, j_eff, 1)
        s_ab_large(1,1,ri,rj)                = s_ab(i, j_eff, 1)

        s_ab_large(1,2,i+L,j+L)              = s_ab(i,j_eff, 2)
        s_ab_large(1,2,ri,rj)                = s_ab(i,j_eff, 2)

        ! now extend the FFTW ouput out to be symmetrical
        ! to make plotting in gnuplot easier later
        if (j.eq.1.and.i.eq.1) then

          s_ab_large(1,2,ri,rj)  = (-1.0)*&
          s_ab_large(1,2,ri,rj)
          s_ab_large(2,1,ri,rj)  = (-1.0)*&
          s_ab_large(2,1,ri,rj)

        end if

        if (j.eq.L+1) then

          s_ab_large(1,2,i+L,j+L)  = (-1.0)*&
          s_ab_large(1,2,i+L,j+L)
          s_ab_large(2,1,i+L,j+L)  = (-1.0)*&
          s_ab_large(2,1,i+L,j+L)

          s_ab_large(1,2,ri,rj)  = (-1.0)*&
          s_ab_large(1,2,ri,rj)
          s_ab_large(2,1,ri,rj)  = (-1.0)*&
          s_ab_large(2,1,ri,rj)

        end if

        s_ab_large(2,2,i+L,j+L)              = s_ab(i,j_eff, 3)
        s_ab_large(2,2,ri,rj)                = s_ab(i,j_eff, 3)

      end do
    end do
    
    ! now out to the other quadrants
    ! we only have the quadrant [kx, ky] \in [0,2pi]
    do kx = -L, L
      do ky = -L, L

        pmx = 0; pmy = 0; offdiag = 1
        i = kx + 1 + L
        j = ky + 1 + L

        if (kx.lt.0) then
          pmx = +1
          offdiag = -1 * offdiag
        end if

        if (ky.lt.0) then
          pmy = +1
          offdiag = -1 * offdiag
        end if

        ! these lines seem hacky to me but they exactly
        ! reproduce what I had before, which in turn
        ! came from directly simulating out to some huge k.
        ! Might be worth rechecking at some point though.
        if ((kx).eq.-L.or.(ky).eq.-L) then
          offdiag = -1
        end if
        if ((kx.eq.L.and.ky.lt.0).or.(ky.eq.L.and.kx.lt.0)) then
          offdiag = +1
        end if
        
        s_ab_large(1,1,i,j) = s_ab_large(1, 1, i + pmx * L, j + pmy * L)
        s_ab_large(2,2,i,j) = s_ab_large(2, 2, i + pmx * L, j + pmy * L)
        s_ab_large(1,2,i,j) = offdiag * s_ab_large(1, 2, i + pmx * L, j + pmy * L)

      end do
    end do

    s_ab_large(2,1,:,:) = s_ab_large(1,2,:,:)

  end subroutine fix_fftw


  subroutine calc_correlations
    use common
    implicit none
    integer :: i, j, kx, ky, m, p, x, y, dist_bin, sp
    real(kind=rk) :: norm_k, kx_float, ky_float, dist,&
    sp_he_tot, sp_he_rot, sp_he_irrot, prefac,&
    ebar_sus, ebar_dip_sus, ebar_wind_sus
    real(kind=rk), dimension(2,2,(bz*L)+1,(bz*L)+1) :: chi_ab,&
    chi_ab_rot, chi_ab_irrot
    real(kind=rk), dimension((bz*L)+1,(bz*L)+1) :: charge_struc
    character(100) :: struc_format_string, field_format_string,&
    vertex_format, avg_field_format,dir_format_string, dir_dist_format_string

    avg_field_format = "(I3, I3, 7ES18.9)"
    field_format_string = "(12ES18.9)"
    struc_format_string = "(4ES18.9)"
    dir_format_string = "(I3, I3, ES18.9)"
    vertex_format = "(I6, I3, I3, 6ES18.9, I3)"
    dir_dist_format_string = "(ES18.9, I9.6, ES18.9)"

    charge_struc = 0.0

    prefac = 1.0 * L**2 / (temp**2)

    sp_he_tot = prefac * (ener_tot_sq_sum - (ener_tot_sum)**2)
    sp_he_rot = prefac * (ener_rot_sq_sum - (ener_rot_sum)**2)
    sp_he_irrot = prefac * (ener_irrot_sq_sum - (ener_irrot_sum)**2)
    ebar_sus = L**2 * beta * (ebar_sq_sum(1) + ebar_sq_sum(2)&
    - (ebar_sum(1)**2 + ebar_sum(2)**2))
    ebar_dip_sus = L**2 * beta * (ebar_dip_sq_sum(1) + ebar_dip_sq_sum(2)&
    - (ebar_dip_sum(1)**2 + ebar_dip_sum(2)**2))
    ebar_wind_sus = L**2 * beta * (ebar_wind_sq_sum(1) + ebar_wind_sq_sum(2)&
    - (ebar_wind_sum(1)**2 + ebar_wind_sum(2)**2))

    write(*,*)
    write (*,'(a)') "# Specific heat: total, rot., irrot."
    write (*,'(ES18.9, ES18.9, ES18.9, ES18.9)') temp,&
    sp_he_tot, sp_he_rot, sp_he_irrot
    write(*,*) "E^bar averages:"
    write(*,'(a)') "<E^bar_x>, <E^bar_y>, <|E^bar|>"
    write (*,*) ebar_sum(1), ebar_sum(2),&
    sqrt(ebar_sum(1)**2+ebar_sum(2)**2)
    write(*,'(a)') "<(E^bar_x)^2>, <(E^bar_y)^2>, <|(E^bar)^2|>"
    write (*,*) ebar_sq_sum(1), ebar_sq_sum(2),&
    sqrt(ebar_sq_sum(1)**2+ebar_sq_sum(2)**2)

    write(*,'(a)') "<|(E^bar)^2|> - <|E^bar|>^2"
    write (*,*) L**2 * beta * (ebar_sq_sum(1) + ebar_sq_sum(2)&
    - (ebar_sum(1)**2 + ebar_sum(2)**2))

    if (do_corr) then

      chi_ab = s_ab_large / temp

      do p = (-L/2)*sp,(L/2)*sp
        do m = (-L/2)*sp,(L/2)*sp

          i = m + 1 + sp*(L/2)
          j = p + 1 + sp*(L/2)

          ! leaving this here in case I decide to reinstate
          ! charge structure factor calculation
          ! can change loop extents now that s_perp is gone
          if (j.le.(bz*L + 1).and.i.le.(bz*L + 1)) then
            ! can also subtract e.g. rho_k_p * conjg(rho_k_m)
            charge_struc(i,j) = abs(ch_ch(i,j) - &
              rho_k_p(i,j) * conjg(rho_k_m(i,j)))
          end if

        end do
      end do ! end p, m loops

      open(unit=10, file=dir_st_file)
      open(unit=11, file=dir_dist_file)
      open(unit=12, file=charge_st_file)
      open(unit=14, file=s_ab_file)
      open(unit=25, file=avg_field_file)
      open(unit=26, file=chi_ab_file)
      open(unit=38, file=windings_file)
      open(unit=39, file=windings_sq_file)

      open  (30, file=sphe_sus_file, position='append')
      write (30, '(a)') "# S_ab integrals (* L**2)!"
      write (30, '(a)') "# S_xx"
      write (30, '(3f20.8)') sum(real(s_ab_large(1,1,:,:)))
      write (30, '(a)') "# S_xy"
      write (30, '(3f20.8)') sum(real(s_ab_large(1,2,:,:)))
      write (30, '(a)') "# S_yx"
      write (30, '(3f20.8)') sum(real(s_ab_large(2,1,:,:)))
      write (30, '(a)') "# S_yy"
      write (30, '(3f20.8)') sum(real(s_ab_large(2,2,:,:)))
      write (30, '(a)') "# Chi_ab integrals (* L**2)!"
      write (30, '(a)') "# Chi_xx"
      write (30, '(3f20.8)') sum(real(chi_ab(1,1,:,:)))
      write (30, '(a)') "# Chi_xy"
      write (30, '(3f20.8)') sum(real(chi_ab(1,2,:,:)))
      write (30, '(a)') "# Chi_yx"
      write (30, '(3f20.8)') sum(real(chi_ab(2,1,:,:)))
      write (30, '(a)') "# Chi_yy"
      write (30, '(3f20.8)') sum(real(chi_ab(2,2,:,:)))
      close (30)

      !dist_r = dist_r / (no_measurements)
      do i = 1,ceiling( sqrt(float((3*((L/2)**2)))) * (1 / bin_size) )
        write (11, dir_dist_format_string)&
        i * bin_size, bin_count(i), abs(dist_r(i))
      end do

      do i = 1,no_measurements
        write(38,*) windings(1,i), windings(2,i)
        write(39,*) windings_sq(1,i), windings_sq(2,i)
      end do

      close(11)
      close(38)
      close(39)

      do i = 1,sp*(L) + 1
        do j = 1,sp*(L) + 1

          ! output is kx, ky, kz, S(\vec{k})

          if (i.le.L/2+1.and.j.le.L/2+1) then
            write (10, dir_format_string)&
            i - 1,j - 1,abs(dir_struc(i,j))
          end if

          if (i.le.L.and.j.le.L) then
            write (25, avg_field_format)&
            i,j,e_tot_avg(1,i,j),e_tot_avg(2,i,j),&
            v_avg(i,j)
          end if

          if (j.le.(bz*L + 1).and.i.le.(bz*L + 1)) then

            write (12, struc_format_string)&
            2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
            charge_struc(i,j)

            write (14, field_format_string)&
            2*pi*(i - 1 - bz*(l/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(l/2))/(L*lambda),&
            real(s_ab_large(1,1,i,j)),&
            real(s_ab_large(1,2,i,j)),&
            real(s_ab_large(2,1,i,j)),&
            real(s_ab_large(2,2,i,j))

            write (26, field_format_string)&
            2*pi*(i - 1 - bz*(l/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(l/2))/(L*lambda),&
            real(chi_ab(1,1,i,j)),&
            real(chi_ab(1,2,i,j)),&
            real(chi_ab(2,1,i,j)),&
            real(chi_ab(2,2,i,j))

          end if

        end do
      end do

      close(10)
      close(12)
      close(14)
      close(25)
      close(26)

    end if ! do_corr

  end subroutine calc_correlations

end module output
