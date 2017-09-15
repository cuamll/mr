module output
  use common
  implicit none

  contains

  subroutine write_output

    call calc_correlations

  end subroutine write_output

  subroutine calc_correlations
    use common
    implicit none
    integer :: i, j, kx, ky, m, p, x, y, dist_bin, sp
    real(kind=rk) :: norm_k, kx_float, ky_float, dist,&
    sp_he_tot, sp_he_rot, sp_he_irrot, prefac,&
    ebar_sus, ebar_dip_sus, ebar_wind_sus
    real(kind=rk), dimension(:,:), allocatable :: s_perp,&
    s_perp_irrot, s_perp_rot
    real(kind=rk), dimension(:,:), allocatable :: s_par,&
    s_par_irrot, s_par_rot
    real(kind=rk), dimension((bz*L)+1,(bz*L)+1) :: charge_struc,&
    field_struc, field_struc_irrot, field_struc_rot
    character(100) :: struc_format_string, field_format_string,&
    vertex_format, avg_field_format,dir_format_string, dir_dist_format_string

    avg_field_format = "(I3, I3, 7ES18.9)"
    field_format_string = "(12ES18.9)"
    struc_format_string = "(4ES18.9)"
    dir_format_string = "(I3, I3, ES18.9)"
    vertex_format = "(I6, I3, I3, 6ES18.9, I3)"
    dir_dist_format_string = "(ES18.9, I9.6, ES18.9)"

    field_struc = 0.0; field_struc_rot = 0.0;
    charge_struc = 0.0; field_struc_irrot = 0.0;

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

    open  (30, file=sphe_sus_file)
    write (30,'(a)') "   # Specific heat: total, rot., irrot."
    write (30,'(4ES18.9)') temp,&
    sp_he_tot, sp_he_rot, sp_he_irrot

    write (30,'(a)') "   # T, Chi_{Ebar}, Chi_{Ebar_dip}, Chi_{Ebar_wind}"
    write (30,'(4ES18.9)') temp, ebar_sus, ebar_dip_sus, ebar_wind_sus

    ! write (30,'(a)') "   # hop acceptance"
    ! write (30,'(2ES18.9)') temp,&
    ! (dble(accepts(1)) / dble(attempts(1)))
    write (30,'(a,2i12.1,es18.9)') "# Hops: total, attempts, rate: ",&
    accepts(1), attempts(1), dble(accepts(1)) / dble(attempts(1))
    write (30,'(a,2i12.1,es18.9)') "# Rot.: total, attempts, rate: ",&
    accepts(2), attempts(2), dble(accepts(2)) / dble(attempts(2))
    write (30,'(a,2i12.1,es18.9)') "# Harm: total, attempts, rate: ",&
    accepts(3), attempts(3), dble(accepts(3)) / dble(attempts(3))
    write (30,'(a,2i12.1,es18.9)') "# Creations: total, attempts, rate: ",&
    accepts(4), attempts(4), dble(accepts(4)) / dble(attempts(4))
    write (30,'(a,2i12.1,es18.9)') "# Annihilations: total, attempts, rate: ",&
    accepts(5), attempts(5), dble(accepts(5)) / dble(attempts(5))

    close(30)

    if (do_corr) then
      ! we can calculate s_perp up to wherever
      sp = 8
      allocate(s_perp((sp*L)+1,(sp*L)+1))
      allocate(s_perp_irrot((sp*L)+1,(sp*L)+1))
      allocate(s_perp_rot((sp*L)+1,(sp*L)+1))
      allocate(s_par((sp*L)+1,(sp*L)+1))
      allocate(s_par_irrot((sp*L)+1,(sp*L)+1))
      allocate(s_par_rot((sp*L)+1,(sp*L)+1))
      s_perp = 0.0; s_perp_rot = 0.0; s_perp_irrot = 0.0;
      s_par = 0.0; s_par_rot = 0.0; s_par_irrot = 0.0;

      do p = (-L/2)*sp,(L/2)*sp
        do m = (-L/2)*sp,(L/2)*sp

          i = m + 1 + sp*(L/2)
          j = p + 1 + sp*(L/2)

          if (j.le.(bz*L + 1).and.i.le.(bz*L + 1)) then
            ! can also subtract e.g. rho_k_p * conjg(rho_k_m)
            charge_struc(i,j) = abs(ch_ch(i,j) - &
              rho_k_p(i,j) * conjg(rho_k_m(i,j)))
            field_struc(i,j) = abs(s_ab(1,1,i,j))
            field_struc_irrot(i,j) = abs(s_ab_irrot(1,1,i,j))
            field_struc_rot(i,j) = abs(s_ab_rot(1,1,i,j))
          end if

          ! use separate variables, we're gonna mess around with values
          kx_float = m * ((2 * pi)/(L * lambda))
          ky_float = p * ((2 * pi)/(L * lambda))

          if (kx_float.eq.0.and.ky_float.eq.0) then
            norm_k = 0.0
          else
            norm_k = 1.0/(kx_float**2 + ky_float**2)
          end if

          ! move back to somewhere we already know about
          ! this is probably not right yet.
          ! but i know something needs doing here
          ! s_ab should be periodic in 2pi so mod should be fine

          if (abs(p).gt.(L/2)*bz) then
            ky = mod(p,(bz*L)/2) + 1 + (bz*L)/2
          else
            ky = p + 1 + (bz*L)/2
          end if
          if (abs(m).gt.(L/2)*bz) then
            kx = mod(m,(bz*L)/2) + 1 + (bz*L)/2
          else
            kx = m + 1 + (bz*L)/2
          end if

          s_perp(i,j) = (1 - kx_float*kx_float*norm_k) * s_ab(1,1,kx,ky)+&
                          ((-1)*kx_float*ky_float*norm_k) * s_ab(1,2,kx,ky)+&
                          ((-1)*ky_float*kx_float*norm_k) * s_ab(2,1,kx,ky)+&
                          (1 - ky_float*ky_float*norm_k) * s_ab(2,2,kx,ky)

          s_perp_irrot(i,j) = (1 - kx_float*kx_float*norm_k) * s_ab_irrot(1,1,kx,ky)+&
                          ((-1)*kx_float*ky_float*norm_k) * s_ab_irrot(1,2,kx,ky)+&
                          ((-1)*ky_float*kx_float*norm_k) * s_ab_irrot(2,1,kx,ky)+&
                          (1 - ky_float*ky_float*norm_k) * s_ab_irrot(2,2,kx,ky)

          s_perp_rot(i,j) = (1 - kx_float*kx_float*norm_k) * s_ab_rot(1,1,kx,ky)+&
                          ((-1)*kx_float*ky_float*norm_k) * s_ab_rot(1,2,kx,ky)+&
                          ((-1)*ky_float*kx_float*norm_k) * s_ab_rot(2,1,kx,ky)+&
                          (1 - ky_float*ky_float*norm_k) * s_ab_rot(2,2,kx,ky)

          s_par(i,j) = (kx_float*kx_float*norm_k) * s_ab(1,1,kx,ky)+&
                         (kx_float*ky_float*norm_k) * s_ab(1,2,kx,ky)+&
                         (ky_float*kx_float*norm_k) * s_ab(2,1,kx,ky)+&
                         (ky_float*ky_float*norm_k) * s_ab(2,2,kx,ky)

          s_par_irrot(i,j) = (kx_float*kx_float*norm_k)*s_ab_irrot(1,1,kx,ky)+&
                               (kx_float*ky_float*norm_k)*s_ab_irrot(1,2,kx,ky)+&
                               (ky_float*kx_float*norm_k)*s_ab_irrot(2,1,kx,ky)+&
                               (ky_float*ky_float*norm_k)*s_ab_irrot(2,2,kx,ky)

          s_par_rot(i,j) = (kx_float*kx_float*norm_k)*s_ab_rot(1,1,kx,ky)+&
                             (kx_float*ky_float*norm_k)*s_ab_rot(1,2,kx,ky)+&
                             (ky_float*kx_float*norm_k)*s_ab_rot(2,1,kx,ky)+&
                             (ky_float*ky_float*norm_k)*s_ab_rot(2,2,kx,ky)

        end do
      end do ! end p, m loops

      open(unit=10, file=dir_st_file)
      open(unit=11, file=dir_dist_file)
      open(unit=12, file=charge_st_file)
      open(unit=13, file=field_st_file)
      open(unit=14, file=s_ab_file)
      open(unit=15, file=s_perp_file)
      open(unit=16, file=irrot_field_file)
      open(unit=17, file=irrot_sab_file)
      open(unit=18, file=irrot_sperp_file)
      open(unit=19, file=rot_field_file)
      open(unit=20, file=rot_sab_file)
      open(unit=21, file=rot_sperp_file)
      open(unit=22, file=spar_file)
      open(unit=23, file=irrot_spar_file)
      open(unit=24, file=rot_spar_file)
      open(unit=25, file=avg_field_file)

      !dist_r = dist_r / (no_measurements)
      do i = 1,ceiling( sqrt(float((3*((L/2)**2)))) * (1 / bin_size) )
        write (11, dir_dist_format_string)&
        i * bin_size, bin_count(i), abs(dist_r(i))
      end do

      close(11)

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
            e_rot_avg(1,i,j),e_rot_avg(2,i,j),&
            e_irrot_avg(1,i,j),e_irrot_avg(2,i,j),&
            v_avg(i,j)
          end if

          if (j.le.(bz*L + 1).and.i.le.(bz*L + 1)) then

            write (12, struc_format_string)&
            2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
            charge_struc(i,j)

            write (13, struc_format_string)&
            2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
            field_struc(i,j)

            write (14, field_format_string)&
            2*pi*(i - 1 - bz*(l/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(l/2))/(L*lambda),&
            s_ab(1,1,i,j),&
            s_ab(1,2,i,j),&
            s_ab(2,1,i,j),&
            s_ab(2,2,i,j)

            write (16, struc_format_string)&
            2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
            field_struc_irrot(i,j)

            write (17, field_format_string)&
            2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
            s_ab_irrot(1,1,i,j),&
            s_ab_irrot(1,2,i,j),&
            s_ab_irrot(2,1,i,j),&
            s_ab_irrot(2,2,i,j)

            write (19, struc_format_string)&
            2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
            field_struc_rot(i,j)

            write (20, field_format_string)&
            2*pi*(i - 1 - bz*(l/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(l/2))/(L*lambda),&
            s_ab_rot(1,1,i,j),&
            s_ab_rot(1,2,i,j),&
            s_ab_rot(2,1,i,j),&
            s_ab_rot(2,2,i,j)

          end if

          write (15, field_format_string)&
          2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
          2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
          s_perp(i,j)

          write (18, field_format_string)&
          2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
          2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
          s_perp_irrot(i,j)

          write (21, field_format_string)&
          2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
          2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
          s_perp_rot(i,j)

          write (22, field_format_string)&
          2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
          2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
          s_par(i,j)

          write (23, field_format_string)&
          2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
          2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
          s_par_irrot(i,j)

          write (24, field_format_string)&
          2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
          2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
          s_par_rot(i,j)

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

      deallocate(s_perp); deallocate(s_perp_rot); deallocate(s_perp_irrot);
      deallocate(s_par); deallocate(s_par_rot); deallocate(s_par_irrot);
    end if ! do_corr

  end subroutine calc_correlations

end module output
