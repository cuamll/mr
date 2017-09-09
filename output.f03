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
    sp_he_tot, sp_he_rot, sp_he_irrot, prefac
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

    write (30,'(a)') "   # <ebar^2> - <ebar>^2"
    write (30,'(2ES18.9)') temp,&
    L**2 * beta * (ebar_sq_sum(1) + ebar_sq_sum(2)&
    - (ebar_sum(1)**2 + ebar_sum(2)**2))

    write (30,'(a)') "   # hop acceptance"
    write (30,'(2ES18.9)') temp,&
    (dble(accepth) / &
    ((therm_sweeps + measurement_sweeps) * add_charges * hop_ratio))

    close(30)

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


  end subroutine calc_correlations

end module output
