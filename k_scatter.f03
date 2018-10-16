module k_scatter
  use common
  implicit none

  contains

  subroutine write_output

    call read_arrays
    call fix_arrays
    call calc_correlations

  end subroutine write_output


  subroutine read_arrays
    use common
    implicit none
    integer :: i, j, k
    real(kind=rk) :: kx_float, ky_float, kz_float,&
    kx_norm_float, ky_norm_float, kz_norm_float,&
    sxx,sxy,syx,syy
    character(100) :: input_format_string

    input_format_string = "(8F14.10)"

    open(11, file=s_ab_transverse_file, action='read')
    open(12, file=s_ab_long_file, action='read')

    do i = 1, (L) + 1
      do j = 1, (L) + 1

          ! read (11, input_format_string) kx_float, ky_float,&
          read (11, *) kx_float, ky_float,&
          sxx,sxy,syx,syy,&
          kx_norm_float, ky_norm_float

          s_ab_rot_large(1, 1, i + L/2, j + L/2) = cmplx(sxx, 0.0)
          s_ab_rot_large(1, 2, i + L/2, j + L/2) = cmplx(sxy, 0.0)
          s_ab_rot_large(2, 1, i + L/2, j + L/2) = cmplx(syx, 0.0)
          s_ab_rot_large(2, 2, i + L/2, j + L/2) = cmplx(syy, 0.0)

          sxx = 0.0; sxy = 0.0; syx = 0.0; syy = 0.0

          read (12, *) kx_float, ky_float,&
          sxx,sxy,syx,syy,&
          kx_norm_float, ky_norm_float

          ! real(s_ab_irrot(1, 1, i + L/2, j + L/2)),&
          ! real(s_ab_irrot(1, 2, i + L/2, j + L/2)),&
          ! real(s_ab_irrot(2, 1, i + L/2, j + L/2)),&
          ! real(s_ab_irrot(2, 2, i + L/2, j + L/2))

          s_ab_irrot_large(1, 1, i + L/2, j + L/2) = cmplx(sxx, 0.0)
          s_ab_irrot_large(1, 2, i + L/2, j + L/2) = cmplx(sxy, 0.0)
          s_ab_irrot_large(2, 1, i + L/2, j + L/2) = cmplx(syx, 0.0)
          s_ab_irrot_large(2, 2, i + L/2, j + L/2) = cmplx(syy, 0.0)

      end do
    end do

    close(11)
    close(12)

  end subroutine read_arrays


  subroutine fix_arrays
    use common
    implicit none
    integer :: i, j, kx, ky, pmx, pmy, offdiag
    
    do kx = -(bz*L/2), (bz*L/2)
      do ky = -(bz*L/2), (bz*L/2)

        pmx = 0; pmy = 0; offdiag = 1
        i = kx + 1 + L
        j = ky + 1 + L

        if (kx.ge.(L/2)+1) then
          pmx = -1
          offdiag = -1 * offdiag
        end if

        if (kx.lt.-(L/2)+1) then
          pmx = +1
          offdiag = -1 * offdiag
        end if

        if (ky.ge.(L/2)+1) then
          pmy = -1
          offdiag = -1 * offdiag
        end if

        if (ky.lt.-(L/2)+1) then
          pmy = +1
          offdiag = -1 * offdiag
        end if
        
        s_ab_rot_large(1,1,i,j) = s_ab_rot_large(1, 1, i + pmx * L, j + pmy * L)
        s_ab_rot_large(2,2,i,j) = s_ab_rot_large(2, 2, i + pmx * L, j + pmy * L)
        s_ab_rot_large(1,2,i,j) = offdiag * s_ab_rot_large(1, 2, i + pmx * L, j + pmy * L)
        s_ab_rot_large(2,1,i,j) = offdiag * s_ab_rot_large(2, 1, i + pmx * L, j + pmy * L)

        s_ab_irrot_large(1,1,i,j) = s_ab_irrot_large(1, 1, i + pmx * L, j + pmy * L)
        s_ab_irrot_large(2,2,i,j) = s_ab_irrot_large(2, 2, i + pmx * L, j + pmy * L)
        s_ab_irrot_large(1,2,i,j) = offdiag * s_ab_irrot_large(1, 2, i + pmx * L, j + pmy * L)
        s_ab_irrot_large(2,1,i,j) = offdiag * s_ab_irrot_large(2, 1, i + pmx * L, j + pmy * L)

      end do
    end do

  end subroutine fix_arrays


  subroutine calc_correlations
    use common
    implicit none
    integer :: i, j, k, kx, ky, kz, m, p, s, x, y, z, dist_bin, sp
    real(kind=rk) :: norm_k, kx_float, ky_float, kz_float, dist,&
    sp_he_tot, sp_he_rot, sp_he_irrot, prefac,&
    ebar_sus, ebar_dip_sus, ebar_wind_sus
    real(kind=rk), dimension(:,:), allocatable :: s_perp,&
    s_perp_irrot, s_perp_rot
    real(kind=rk), dimension(:,:), allocatable :: s_par,&
    s_par_irrot, s_par_rot
    character(100) :: struc_format_string, field_format_string,&
    vertex_format, avg_field_format,dir_format_string, dir_dist_format_string

    avg_field_format = "(3I3, 7ES18.9)"
    field_format_string = "(12ES18.9)"
    struc_format_string = "(4ES18.9)"
    dir_format_string = "(3I3, ES18.9)"
    vertex_format = "(I6, 3I3, 6ES18.9, I3)"
    dir_dist_format_string = "(ES18.9, I9.6, ES18.9)"

    ! s_ab_rot = s_ab_rot * 0.25
    ! s_ab_irrot = s_ab_irrot * 0.25

    sp = 8
    allocate(s_perp_irrot((sp*L)+1,(sp*L)+1))
    allocate(s_perp_rot((sp*L)+1,(sp*L)+1))
    allocate(s_par_irrot((sp*L)+1,(sp*L)+1))
    allocate(s_par_rot((sp*L)+1,(sp*L)+1))
    s_perp_rot = 0.0; s_perp_irrot = 0.0;
    s_par_rot = 0.0; s_par_irrot = 0.0;

      do p = (-L/2)*sp,(L/2)*sp
        do m = (-L/2)*sp,(L/2)*sp

          i = m + 1 + sp*(L/2)
          j = p + 1 + sp*(L/2)

          ! use separate variables, we're gonna mess around with values
          kx_float = m * ((2 * pi)/(L * lambda))
          ky_float = p * ((2 * pi)/(L * lambda))

          if (kx_float.eq.0.and.ky_float.eq.0) then
            norm_k = 0.0
          else
            norm_k = 1.0/(kx_float**2 + ky_float**2)
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

          if (p.le.(-1*(bz*(L/2)))) then
            ky = p
            do while (ky.le.(-1*(bz*(L/2))))
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

          if (m.le.(-1*(bz*(L/2)))) then
            kx = m
            do while (kx.le.(-1*(bz*(L/2))))
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

          s_perp_irrot(i,j) = (1 - kx_float*kx_float*norm_k) *  real(s_ab_irrot_large(1,1,kx,ky))+&
                          ((-1)*kx_float*ky_float*norm_k) *     real(s_ab_irrot_large(1,2,kx,ky))+&
                          ((-1)*ky_float*kx_float*norm_k) *     real(s_ab_irrot_large(2,1,kx,ky))+&
                          (1 - ky_float*ky_float*norm_k) *      real(s_ab_irrot_large(2,2,kx,ky))

          s_perp_rot(i,j) = (1 - kx_float*kx_float*norm_k) * real(s_ab_rot_large(1,1,kx,ky))+&
                          ((-1)*kx_float*ky_float*norm_k) *  real(s_ab_rot_large(1,2,kx,ky))+&
                          ((-1)*ky_float*kx_float*norm_k) *  real(s_ab_rot_large(2,1,kx,ky))+&
                          (1 - ky_float*ky_float*norm_k) *   real(s_ab_rot_large(2,2,kx,ky))

          s_par_irrot(i,j) = (kx_float*kx_float*norm_k)*  real(s_ab_irrot_large(1,1,kx,ky))+&
                               (kx_float*ky_float*norm_k)*real(s_ab_irrot_large(1,2,kx,ky))+&
                               (ky_float*kx_float*norm_k)*real(s_ab_irrot_large(2,1,kx,ky))+&
                               (ky_float*ky_float*norm_k)*real(s_ab_irrot_large(2,2,kx,ky))

          s_par_rot(i,j) = (kx_float*kx_float*norm_k)*  real(s_ab_rot_large(1,1,kx,ky))+&
                             (kx_float*ky_float*norm_k)*real(s_ab_rot_large(1,2,kx,ky))+&
                             (ky_float*kx_float*norm_k)*real(s_ab_rot_large(2,1,kx,ky))+&
                             (ky_float*ky_float*norm_k)*real(s_ab_rot_large(2,2,kx,ky))

        end do
      end do ! end p, m loops

    open(unit=12, file=s_perp_long_file)
    open(unit=13, file=s_perp_transverse_file)
    open(unit=15, file=s_par_long_file)
    open(unit=16, file=s_par_transverse_file)

    do i = 1,sp*(L) + 1
      do j = 1,sp*(L) + 1

          write (12, field_format_string)&
          2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
          2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
          s_perp_irrot(i,j)

          write (13, field_format_string)&
          2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
          2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
          s_perp_rot(i,j)

          write (15, field_format_string)&
          2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
          2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
          s_par_irrot(i,j)

          write (16, field_format_string)&
          2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
          2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
          s_par_rot(i,j)

        end do
      end do

    close(12)
    close(13)
    close(15)
    close(16)

    deallocate(s_perp_rot); deallocate(s_perp_irrot);
    deallocate(s_par_rot); deallocate(s_par_irrot);

  end subroutine calc_correlations

end module k_scatter
