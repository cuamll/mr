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
        write (10,'(4f16.12)') real(i + 0.5), real(j), e_field(1,i,j), 0.0
        write (11,'(4f16.12)') real(i + 0.5), real(j), mnphi(1,i,j), 0.0
        write (12,'(4f16.12)') real(i + 0.5), real(j), e_rot(1,i,j), 0.0
        ! y component
        write (10,'(4f16.12)') real(i), real(j + 0.5), 0.0, e_field(2,i,j)
        write (11,'(4f16.12)') real(i), real(j + 0.5), 0.0, mnphi(2,i,j)
        write (12,'(4f16.12)') real(i), real(j + 0.5), 0.0, e_rot(2,i,j)

        ! charges
        write (13,'(2f16.12, i3.1)') real(i), real(j), v(i,j)

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
      call fix_arrays
    end if

    call calc_correlations

  end subroutine write_output

  subroutine fix_fftw
    use common
    implicit none
    integer :: i, j, j_eff, ri, rj, kx, ky, pmx, pmy, offdiag

    ! normalisation for correct equipartition result
    sxx = sxx * L**2
    sxy = sxy * L**2
    syy = syy * L**2

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

        fftw_s_ab_total(1,1,i+L,j+L)               = sxx(i,j_eff)
        fftw_s_ab_total(1,1,ri,rj)                 = sxx(i,j_eff)

        fftw_s_ab_total(1,2,i+L,j+L)              = sxy(i,j_eff)
        fftw_s_ab_total(1,2,ri,rj)                = sxy(i,j_eff)

        ! should really check this?
        if (i.eq.1.and.(j.eq.1.or.j_eff.eq.1)) then
          fftw_s_ab_total(1,2,i+L,j+L)  = (-1.0)*&
          fftw_s_ab_total(1,2,i+L,j+L)
          fftw_s_ab_total(2,1,i+L,j+L)  = (-1.0)*&
          fftw_s_ab_total(2,1,i+L,j+L)

          fftw_s_ab_total(1,2,ri,rj)  = (-1.0)*&
          fftw_s_ab_total(1,2,ri,rj)
          fftw_s_ab_total(2,1,ri,rj)  = (-1.0)*&
          fftw_s_ab_total(2,1,ri,rj)
        end if

        fftw_s_ab_total(2,2,i+L,j+L)              = syy(i,j_eff)
        fftw_s_ab_total(2,2,ri,rj)                = syy(i,j_eff)

      end do
    end do
    
    ! now out to the other quadrants
    ! we only have the quadrant [kx, ky] \in [0,2pi]
    do kx = -L, L
      do ky = -L, L

        pmx = 0; pmy = 0; offdiag = 1
        i = kx + 1 + L
        j = ky + 1 + L

        ! if (kx.ge.(L/2)+1) then
        !   pmx = -1
        !   offdiag = -1 * offdiag
        ! end if

        if (kx.le.0) then
          pmx = +1
          offdiag = -1 * offdiag
        end if

        ! if (ky.ge.(L/2)+1) then
        !   pmy = -1
        !   offdiag = -1 * offdiag
        ! end if

        if (ky.le.0) then
          pmy = +1
          offdiag = -1 * offdiag
        end if
        
        fftw_s_ab_total(1,1,i,j) = fftw_s_ab_total(1, 1, i + pmx * L, j + pmy * L)
        fftw_s_ab_total(2,2,i,j) = fftw_s_ab_total(2, 2, i + pmx * L, j + pmy * L)
        fftw_s_ab_total(1,2,i,j) = offdiag * fftw_s_ab_total(1, 2, i + pmx * L, j + pmy * L)
        fftw_s_ab_total(2,1,i,j) = offdiag * fftw_s_ab_total(2, 1, i + pmx * L, j + pmy * L)

        ! s_ab_irrot(1,1,i,j) = s_ab_irrot(1, 1, i + pmx * L, j + pmy * L)
        ! s_ab_irrot(2,2,i,j) = s_ab_irrot(2, 2, i + pmx * L, j + pmy * L)
        ! s_ab_irrot(1,2,i,j) = offdiag * s_ab_irrot(1, 2, i + pmx * L, j + pmy * L)
        ! s_ab_irrot(2,1,i,j) = offdiag * s_ab_irrot(2, 1, i + pmx * L, j + pmy * L)

      end do
    end do

    fftw_s_ab_total(2,1,:,:) = fftw_s_ab_total(1,2,:,:)

  end subroutine fix_fftw


  subroutine fix_arrays
    use common
    implicit none
    integer :: i, j
    
    ! very long and extremely disgusting way of filling up the arrays
    ! so we can then calculate other correlation functions easily.
    ! probably won't work yet. also can probably be simplified a lot

    do i = 1,L
      do j = 1,L

        if (i.eq.1.and.j.eq.1) then

          ! charge-charge ones
          rho_k_p(1,1) = rho_k_p(L+1,L+1)
          rho_k_m(1,1) = rho_k_m(L+1,L+1)
          rho_k_p(L+1,1) = rho_k_p(L+1,L+1)
          rho_k_m(L+1,1) = rho_k_m(L+1,L+1)
          rho_k_p((bz*L)+1,1) = rho_k_p(L+1,L+1)
          rho_k_m((bz*L)+1,1) = rho_k_m(L+1,L+1)
          rho_k_p(1,L+1) = rho_k_p(L+1,L+1)
          rho_k_m(1,L+1) = rho_k_m(L+1,L+1)
          rho_k_p(1,(bz*L)+1) = rho_k_p(L+1,L+1)
          rho_k_m(1,(bz*L)+1) = rho_k_m(L+1,L+1)
          rho_k_p((bz*L)+1,L+1) = rho_k_p(L+1,L+1)
          rho_k_m((bz*L)+1,L+1) = rho_k_m(L+1,L+1)
          rho_k_p(L+1,(bz*L)+1) = rho_k_p(L+1,L+1)
          rho_k_m(L+1,(bz*L)+1) = rho_k_m(L+1,L+1)
          rho_k_p((bz*L)+1,(bz*L)+1) = rho_k_p(L+1,L+1)
          rho_k_m((bz*L)+1,(bz*L)+1) = rho_k_m(L+1,L+1)

          ! s_ab
          s_ab(1,1,1,1)                 = s_ab(1,1,L+1,L+1)
          s_ab(1,1,L+1,1)               = s_ab(1,1,L+1,L+1)
          s_ab(1,1,(bz*L)+1,1)          = s_ab(1,1,L+1,L+1)
          s_ab(1,1,1,L+1)               = s_ab(1,1,L+1,L+1)
          s_ab(1,1,1,(bz*L)+1)          = s_ab(1,1,L+1,L+1)
          s_ab(1,1,(bz*L)+1,L+1)        = s_ab(1,1,L+1,L+1)
          s_ab(1,1,L+1,(bz*L)+1)        = s_ab(1,1,L+1,L+1)
          s_ab(1,1,(bz*L)+1,(bz*L)+1)   = s_ab(1,1,L+1,L+1)

          s_ab(1,2,1,1)                 = (-1) * s_ab(1,2,L+1,L+1)
          s_ab(1,2,L+1,1)               = (-1) * s_ab(1,2,L+1,L+1)
          s_ab(1,2,(bz*L)+1,1)          = (-1) * s_ab(1,2,L+1,L+1)
          s_ab(1,2,1,L+1)               = (-1) * s_ab(1,2,L+1,L+1)
          s_ab(1,2,1,(bz*L)+1)          = (-1) * s_ab(1,2,L+1,L+1)
          s_ab(1,2,(bz*L)+1,L+1)        = (-1) * s_ab(1,2,L+1,L+1)
          s_ab(1,2,L+1,(bz*L)+1)        = (-1) * s_ab(1,2,L+1,L+1)
          s_ab(1,2,(bz*L)+1,(bz*L)+1)   = (-1) * s_ab(1,2,L+1,L+1)

          s_ab(2,1,1,1)                 = (-1) * s_ab(2,1,L+1,L+1)
          s_ab(2,1,L+1,1)               = (-1) * s_ab(2,1,L+1,L+1)
          s_ab(2,1,(bz*L)+1,1)          = (-1) * s_ab(2,1,L+1,L+1)
          s_ab(2,1,1,L+1)               = (-1) * s_ab(2,1,L+1,L+1)
          s_ab(2,1,1,(bz*L)+1)          = (-1) * s_ab(2,1,L+1,L+1)
          s_ab(2,1,(bz*L)+1,L+1)        = (-1) * s_ab(2,1,L+1,L+1)
          s_ab(2,1,L+1,(bz*L)+1)        = (-1) * s_ab(2,1,L+1,L+1)
          s_ab(2,1,(bz*L)+1,(bz*L)+1)   = (-1) * s_ab(2,1,L+1,L+1)

          s_ab(2,2,1,1)                 = s_ab(2,2,L+1,L+1)
          s_ab(2,2,L+1,1)               = s_ab(2,2,L+1,L+1)
          s_ab(2,2,(bz*L)+1,1)          = s_ab(2,2,L+1,L+1)
          s_ab(2,2,1,L+1)               = s_ab(2,2,L+1,L+1)
          s_ab(2,2,1,(bz*L)+1)          = s_ab(2,2,L+1,L+1)
          s_ab(2,2,(bz*L)+1,L+1)        = s_ab(2,2,L+1,L+1)
          s_ab(2,2,L+1,(bz*L)+1)        = s_ab(2,2,L+1,L+1)
          s_ab(2,2,(bz*L)+1,(bz*L)+1)   = s_ab(2,2,L+1,L+1)

          ! s_ab_rot

          s_ab_rot(1,1,1,1)                 = s_ab_rot(1,1,L+1,L+1)
          s_ab_rot(1,1,L+1,1)               = s_ab_rot(1,1,L+1,L+1)
          s_ab_rot(1,1,(bz*L)+1,1)          = s_ab_rot(1,1,L+1,L+1)
          s_ab_rot(1,1,1,L+1)               = s_ab_rot(1,1,L+1,L+1)
          s_ab_rot(1,1,1,(bz*L)+1)          = s_ab_rot(1,1,L+1,L+1)
          s_ab_rot(1,1,(bz*L)+1,L+1)        = s_ab_rot(1,1,L+1,L+1)
          s_ab_rot(1,1,L+1,(bz*L)+1)        = s_ab_rot(1,1,L+1,L+1)
          s_ab_rot(1,1,(bz*L)+1,(bz*L)+1)   = s_ab_rot(1,1,L+1,L+1)

          s_ab_rot(1,2,1,1)                 = (-1) * s_ab_rot(1,2,L+1,L+1)
          s_ab_rot(1,2,L+1,1)               = (-1) * s_ab_rot(1,2,L+1,L+1)
          s_ab_rot(1,2,(bz*L)+1,1)          = (-1) * s_ab_rot(1,2,L+1,L+1)
          s_ab_rot(1,2,1,L+1)               = (-1) * s_ab_rot(1,2,L+1,L+1)
          s_ab_rot(1,2,1,(bz*L)+1)          = (-1) * s_ab_rot(1,2,L+1,L+1)
          s_ab_rot(1,2,(bz*L)+1,L+1)        = (-1) * s_ab_rot(1,2,L+1,L+1)
          s_ab_rot(1,2,L+1,(bz*L)+1)        = (-1) * s_ab_rot(1,2,L+1,L+1)
          s_ab_rot(1,2,(bz*L)+1,(bz*L)+1)   = (-1) * s_ab_rot(1,2,L+1,L+1)

          s_ab_rot(2,1,1,1)                 = (-1) * s_ab_rot(2,1,L+1,L+1)
          s_ab_rot(2,1,L+1,1)               = (-1) * s_ab_rot(2,1,L+1,L+1)
          s_ab_rot(2,1,(bz*L)+1,1)          = (-1) * s_ab_rot(2,1,L+1,L+1)
          s_ab_rot(2,1,1,L+1)               = (-1) * s_ab_rot(2,1,L+1,L+1)
          s_ab_rot(2,1,1,(bz*L)+1)          = (-1) * s_ab_rot(2,1,L+1,L+1)
          s_ab_rot(2,1,(bz*L)+1,L+1)        = (-1) * s_ab_rot(2,1,L+1,L+1)
          s_ab_rot(2,1,L+1,(bz*L)+1)        = (-1) * s_ab_rot(2,1,L+1,L+1)
          s_ab_rot(2,1,(bz*L)+1,(bz*L)+1)   = (-1) * s_ab_rot(2,1,L+1,L+1)

          s_ab_rot(2,2,1,1)                 = s_ab_rot(2,2,L+1,L+1)
          s_ab_rot(2,2,L+1,1)               = s_ab_rot(2,2,L+1,L+1)
          s_ab_rot(2,2,(bz*L)+1,1)          = s_ab_rot(2,2,L+1,L+1)
          s_ab_rot(2,2,1,L+1)               = s_ab_rot(2,2,L+1,L+1)
          s_ab_rot(2,2,1,(bz*L)+1)          = s_ab_rot(2,2,L+1,L+1)
          s_ab_rot(2,2,(bz*L)+1,L+1)        = s_ab_rot(2,2,L+1,L+1)
          s_ab_rot(2,2,L+1,(bz*L)+1)        = s_ab_rot(2,2,L+1,L+1)
          s_ab_rot(2,2,(bz*L)+1,(bz*L)+1)   = s_ab_rot(2,2,L+1,L+1)

          ! s_ab_irrot

          s_ab_irrot(1,1,1,1)                 = s_ab_irrot(1,1,L+1,L+1)
          s_ab_irrot(1,1,L+1,1)               = s_ab_irrot(1,1,L+1,L+1)
          s_ab_irrot(1,1,(bz*L)+1,1)          = s_ab_irrot(1,1,L+1,L+1)
          s_ab_irrot(1,1,1,L+1)               = s_ab_irrot(1,1,L+1,L+1)
          s_ab_irrot(1,1,1,(bz*L)+1)          = s_ab_irrot(1,1,L+1,L+1)
          s_ab_irrot(1,1,(bz*L)+1,L+1)        = s_ab_irrot(1,1,L+1,L+1)
          s_ab_irrot(1,1,L+1,(bz*L)+1)        = s_ab_irrot(1,1,L+1,L+1)
          s_ab_irrot(1,1,(bz*L)+1,(bz*L)+1)   = s_ab_irrot(1,1,L+1,L+1)

          s_ab_irrot(1,2,1,1)                 = (-1) * s_ab_irrot(1,2,L+1,L+1)
          s_ab_irrot(1,2,L+1,1)               = (-1) * s_ab_irrot(1,2,L+1,L+1)
          s_ab_irrot(1,2,(bz*L)+1,1)          = (-1) * s_ab_irrot(1,2,L+1,L+1)
          s_ab_irrot(1,2,1,L+1)               = (-1) * s_ab_irrot(1,2,L+1,L+1)
          s_ab_irrot(1,2,1,(bz*L)+1)          = (-1) * s_ab_irrot(1,2,L+1,L+1)
          s_ab_irrot(1,2,(bz*L)+1,L+1)        = (-1) * s_ab_irrot(1,2,L+1,L+1)
          s_ab_irrot(1,2,L+1,(bz*L)+1)        = (-1) * s_ab_irrot(1,2,L+1,L+1)
          s_ab_irrot(1,2,(bz*L)+1,(bz*L)+1)   = (-1) * s_ab_irrot(1,2,L+1,L+1)

          s_ab_irrot(2,1,1,1)                 = (-1) * s_ab_irrot(2,1,L+1,L+1)
          s_ab_irrot(2,1,L+1,1)               = (-1) * s_ab_irrot(2,1,L+1,L+1)
          s_ab_irrot(2,1,(bz*L)+1,1)          = (-1) * s_ab_irrot(2,1,L+1,L+1)
          s_ab_irrot(2,1,1,L+1)               = (-1) * s_ab_irrot(2,1,L+1,L+1)
          s_ab_irrot(2,1,1,(bz*L)+1)          = (-1) * s_ab_irrot(2,1,L+1,L+1)
          s_ab_irrot(2,1,(bz*L)+1,L+1)        = (-1) * s_ab_irrot(2,1,L+1,L+1)
          s_ab_irrot(2,1,L+1,(bz*L)+1)        = (-1) * s_ab_irrot(2,1,L+1,L+1)
          s_ab_irrot(2,1,(bz*L)+1,(bz*L)+1)   = (-1) * s_ab_irrot(2,1,L+1,L+1)

          s_ab_irrot(2,2,1,1)                 = s_ab_irrot(2,2,L+1,L+1)
          s_ab_irrot(2,2,L+1,1)               = s_ab_irrot(2,2,L+1,L+1)
          s_ab_irrot(2,2,(bz*L)+1,1)          = s_ab_irrot(2,2,L+1,L+1)
          s_ab_irrot(2,2,1,L+1)               = s_ab_irrot(2,2,L+1,L+1)
          s_ab_irrot(2,2,1,(bz*L)+1)          = s_ab_irrot(2,2,L+1,L+1)
          s_ab_irrot(2,2,(bz*L)+1,L+1)        = s_ab_irrot(2,2,L+1,L+1)
          s_ab_irrot(2,2,L+1,(bz*L)+1)        = s_ab_irrot(2,2,L+1,L+1)
          s_ab_irrot(2,2,(bz*L)+1,(bz*L)+1)   = s_ab_irrot(2,2,L+1,L+1)

        else if (i.eq.1.and.j.gt.1) then

          rho_k_p(i,j) = rho_k_p(i+L,j+L)
          rho_k_m(i,j) = rho_k_m(i+L,j+L)
          ch_ch(i,j) = ch_ch(i+L,j+L)

          rho_k_p(i+L,j) = rho_k_p(i+L,j+L)
          rho_k_m(i+L,j) = rho_k_m(i+L,j+L)
          ch_ch(i+L,j) = ch_ch(i+L,j+L)

          rho_k_p(i+2*L,j) = rho_k_p(i+L,j+L)
          rho_k_m(i+2*L,j) = rho_k_m(i+L,j+L)
          ch_ch(i+2*L,j) = ch_ch(i+L,j+L)

          rho_k_p(i,j+L) = rho_k_p(i+L,j+L)
          rho_k_m(i,j+L) = rho_k_m(i+L,j+L)
          ch_ch(i,j+L) = ch_ch(i+L,j+L)

          rho_k_p(i+2*L,j+L) = rho_k_p(i+L,j+L)
          rho_k_m(i+2*L,j+L) = rho_k_m(i+L,j+L)
          ch_ch(i+2*L,j+L) = ch_ch(i+L,j+L)

          ! s_ab
          s_ab(1,1,i,j) = s_ab(1,1,i+L,j+L)
          s_ab(1,2,i,j) = -1 * s_ab(1,2,i+L,j+L)
          s_ab(2,1,i,j) = -1 * s_ab(2,1,i+L,j+L)
          s_ab(2,2,i,j) = s_ab(2,2,i+L,j+L)

          s_ab(1,1,i+L,j) = s_ab(1,1,i+L,j+L)
          s_ab(1,2,i+L,j) = -1 * s_ab(1,2,i+L,j+L)
          s_ab(2,1,i+L,j) = -1 * s_ab(2,1,i+L,j+L)
          s_ab(2,2,i+L,j) = s_ab(2,2,i+L,j+L)

          s_ab(1,1,i+2*L,j) = s_ab(1,1,i+L,j+L)
          s_ab(1,2,i+2*L,j) = -1 * s_ab(1,2,i+L,j+L)
          s_ab(2,1,i+2*L,j) = -1 * s_ab(2,1,i+L,j+L)
          s_ab(2,2,i+2*L,j) = s_ab(2,2,i+L,j+L)

          s_ab(1,1,i,j+L) = s_ab(1,1,i+L,j+L)
          s_ab(1,2,i,j+L) = -1 * s_ab(1,2,i+L,j+L)
          s_ab(2,1,i,j+L) = -1 * s_ab(2,1,i+L,j+L)
          s_ab(2,2,i,j+L) = s_ab(2,2,i+L,j+L)

          s_ab(1,1,i+2*L,j+L) = s_ab(1,1,i+L,j+L)
          s_ab(1,2,i+2*L,j+L) = -1 * s_ab(1,2,i+L,j+L)
          s_ab(2,1,i+2*L,j+L) = -1 * s_ab(2,1,i+L,j+L)
          s_ab(2,2,i+2*L,j+L) = s_ab(2,2,i+L,j+L)

          ! s_ab_rot
          s_ab_rot(1,1,i,j) = s_ab_rot(1,1,i+L,j+L)
          s_ab_rot(1,2,i,j) = -1 * s_ab_rot(1,2,i+L,j+L)
          s_ab_rot(2,1,i,j) = -1 * s_ab_rot(2,1,i+L,j+L)
          s_ab_rot(2,2,i,j) = s_ab_rot(2,2,i+L,j+L)

          s_ab_rot(1,1,i+L,j) = s_ab_rot(1,1,i+L,j+L)
          s_ab_rot(1,2,i+L,j) = -1 * s_ab_rot(1,2,i+L,j+L)
          s_ab_rot(2,1,i+L,j) = -1 * s_ab_rot(2,1,i+L,j+L)
          s_ab_rot(2,2,i+L,j) = s_ab_rot(2,2,i+L,j+L)

          s_ab_rot(1,1,i+2*L,j) = s_ab_rot(1,1,i+L,j+L)
          s_ab_rot(1,2,i+2*L,j) = -1 * s_ab_rot(1,2,i+L,j+L)
          s_ab_rot(2,1,i+2*L,j) = -1 * s_ab_rot(2,1,i+L,j+L)
          s_ab_rot(2,2,i+2*L,j) = s_ab_rot(2,2,i+L,j+L)

          s_ab_rot(1,1,i,j+L) = s_ab_rot(1,1,i+L,j+L)
          s_ab_rot(1,2,i,j+L) = -1 * s_ab_rot(1,2,i+L,j+L)
          s_ab_rot(2,1,i,j+L) = -1 * s_ab_rot(2,1,i+L,j+L)
          s_ab_rot(2,2,i,j+L) = s_ab_rot(2,2,i+L,j+L)

          s_ab_rot(1,1,i+2*L,j+L) = s_ab_rot(1,1,i+L,j+L)
          s_ab_rot(1,2,i+2*L,j+L) = -1 * s_ab_rot(1,2,i+L,j+L)
          s_ab_rot(2,1,i+2*L,j+L) = -1 * s_ab_rot(2,1,i+L,j+L)
          s_ab_rot(2,2,i+2*L,j+L) = s_ab_rot(2,2,i+L,j+L)

          ! s_ab_irrot
          s_ab_irrot(1,1,i,j) = s_ab_irrot(1,1,i+L,j+L)
          s_ab_irrot(1,2,i,j) = -1 * s_ab_irrot(1,2,i+L,j+L)
          s_ab_irrot(2,1,i,j) = -1 * s_ab_irrot(2,1,i+L,j+L)
          s_ab_irrot(2,2,i,j) = s_ab_irrot(2,2,i+L,j+L)

          s_ab_irrot(1,1,i+L,j) = s_ab_irrot(1,1,i+L,j+L)
          s_ab_irrot(1,2,i+L,j) = -1 * s_ab_irrot(1,2,i+L,j+L)
          s_ab_irrot(2,1,i+L,j) = -1 * s_ab_irrot(2,1,i+L,j+L)
          s_ab_irrot(2,2,i+L,j) = s_ab_irrot(2,2,i+L,j+L)

          s_ab_irrot(1,1,i+2*L,j) = s_ab_irrot(1,1,i+L,j+L)
          s_ab_irrot(1,2,i+2*L,j) = -1 * s_ab_irrot(1,2,i+L,j+L)
          s_ab_irrot(2,1,i+2*L,j) = -1 * s_ab_irrot(2,1,i+L,j+L)
          s_ab_irrot(2,2,i+2*L,j) = s_ab_irrot(2,2,i+L,j+L)

          s_ab_irrot(1,1,i,j+L) = s_ab_irrot(1,1,i+L,j+L)
          s_ab_irrot(1,2,i,j+L) = -1 * s_ab_irrot(1,2,i+L,j+L)
          s_ab_irrot(2,1,i,j+L) = -1 * s_ab_irrot(2,1,i+L,j+L)
          s_ab_irrot(2,2,i,j+L) = s_ab_irrot(2,2,i+L,j+L)

          s_ab_irrot(1,1,i+2*L,j+L) = s_ab_irrot(1,1,i+L,j+L)
          s_ab_irrot(1,2,i+2*L,j+L) = -1 * s_ab_irrot(1,2,i+L,j+L)
          s_ab_irrot(2,1,i+2*L,j+L) = -1 * s_ab_irrot(2,1,i+L,j+L)
          s_ab_irrot(2,2,i+2*L,j+L) = s_ab_irrot(2,2,i+L,j+L)

        else if (j.eq.1.and.i.gt.1) then

          rho_k_p(i,j) = rho_k_p(i+L,j+L)
          rho_k_m(i,j) = rho_k_m(i+L,j+L)
          ch_ch(i,j) = ch_ch(i+L,j+L)

          rho_k_p(i,j+L) = rho_k_p(i+L,j+L)
          rho_k_m(i,j+L) = rho_k_m(i+L,j+L)
          ch_ch(i,j+L) = ch_ch(i+L,j+L)

          rho_k_p(i,j+2*L) = rho_k_p(i+L,j+L)
          rho_k_m(i,j+2*L) = rho_k_m(i+L,j+L)
          ch_ch(i,j+2*L) = ch_ch(i+L,j+L)

          rho_k_p(i+L,j) = rho_k_p(i+L,j+L)
          rho_k_m(i+L,j) = rho_k_m(i+L,j+L)
          ch_ch(i+L,j) = ch_ch(i+L,j+L)

          rho_k_p(i+L,j+2*L) = rho_k_p(i+L,j+L)
          rho_k_m(i+L,j+2*L) = rho_k_m(i+L,j+L)
          ch_ch(i+L,j+2*L) = ch_ch(i+L,j+L)

          ! s_ab
          s_ab(1,1,i,j) = s_ab(1,1,i+L,j+L)
          s_ab(1,2,i,j) = -1 * s_ab(1,2,i+L,j+L)
          s_ab(2,1,i,j) = -1 * s_ab(2,1,i+L,j+L)
          s_ab(2,2,i,j) = s_ab(2,2,i+L,j+L)

          s_ab(1,1,i,j+L) = s_ab(1,1,i+L,j+L)
          s_ab(1,2,i,j+L) = -1 * s_ab(1,2,i+L,j+L)
          s_ab(2,1,i,j+L) = -1 * s_ab(2,1,i+L,j+L)
          s_ab(2,2,i,j+L) = s_ab(2,2,i+L,j+L)

          s_ab(1,1,i,j+2*L) = s_ab(1,1,i+L,j+L)
          s_ab(1,2,i,j+2*L) = -1 * s_ab(1,2,i+L,j+L)
          s_ab(2,1,i,j+2*L) = -1 * s_ab(2,1,i+L,j+L)
          s_ab(2,2,i,j+2*L) = s_ab(2,2,i+L,j+L)

          s_ab(1,1,i+L,j) = s_ab(1,1,i+L,j+L)
          s_ab(1,2,i+L,j) = -1 * s_ab(1,2,i+L,j+L)
          s_ab(2,1,i+L,j) = -1 * s_ab(2,1,i+L,j+L)
          s_ab(2,2,i+L,j) = s_ab(2,2,i+L,j+L)

          s_ab(1,1,i+L,j+2*L) = s_ab(1,1,i+L,j+L)
          s_ab(1,2,i+L,j+2*L) = -1 * s_ab(1,2,i+L,j+L)
          s_ab(2,1,i+L,j+2*L) = -1 * s_ab(2,1,i+L,j+L)
          s_ab(2,2,i+L,j+2*L) = s_ab(2,2,i+L,j+L)

          ! s_ab_rot
          s_ab_rot(1,1,i,j) = s_ab_rot(1,1,i+L,j+L)
          s_ab_rot(1,2,i,j) = -1 * s_ab_rot(1,2,i+L,j+L)
          s_ab_rot(2,1,i,j) = -1 * s_ab_rot(2,1,i+L,j+L)
          s_ab_rot(2,2,i,j) = s_ab_rot(2,2,i+L,j+L)

          s_ab_rot(1,1,i,j+L) = s_ab_rot(1,1,i+L,j+L)
          s_ab_rot(1,2,i,j+L) = -1 * s_ab_rot(1,2,i+L,j+L)
          s_ab_rot(2,1,i,j+L) = -1 * s_ab_rot(2,1,i+L,j+L)
          s_ab_rot(2,2,i,j+L) = s_ab_rot(2,2,i+L,j+L)

          s_ab_rot(1,1,i,j+2*L) = s_ab_rot(1,1,i+L,j+L)
          s_ab_rot(1,2,i,j+2*L) = -1 * s_ab_rot(1,2,i+L,j+L)
          s_ab_rot(2,1,i,j+2*L) = -1 * s_ab_rot(2,1,i+L,j+L)
          s_ab_rot(2,2,i,j+2*L) = s_ab_rot(2,2,i+L,j+L)

          s_ab_rot(1,1,i+L,j) = s_ab_rot(1,1,i+L,j+L)
          s_ab_rot(1,2,i+L,j) = -1 * s_ab_rot(1,2,i+L,j+L)
          s_ab_rot(2,1,i+L,j) = -1 * s_ab_rot(2,1,i+L,j+L)
          s_ab_rot(2,2,i+L,j) = s_ab_rot(2,2,i+L,j+L)

          s_ab_rot(1,1,i+L,j+2*L) = s_ab_rot(1,1,i+L,j+L)
          s_ab_rot(1,2,i+L,j+2*L) = -1 * s_ab_rot(1,2,i+L,j+L)
          s_ab_rot(2,1,i+L,j+2*L) = -1 * s_ab_rot(2,1,i+L,j+L)
          s_ab_rot(2,2,i+L,j+2*L) = s_ab_rot(2,2,i+L,j+L)

          ! s_ab_irrot
          s_ab_irrot(1,1,i,j) = s_ab_irrot(1,1,i+L,j+L)
          s_ab_irrot(1,2,i,j) = -1 * s_ab_irrot(1,2,i+L,j+L)
          s_ab_irrot(2,1,i,j) = -1 * s_ab_irrot(2,1,i+L,j+L)
          s_ab_irrot(2,2,i,j) = s_ab_irrot(2,2,i+L,j+L)

          s_ab_irrot(1,1,i,j+L) = s_ab_irrot(1,1,i+L,j+L)
          s_ab_irrot(1,2,i,j+L) = -1 * s_ab_irrot(1,2,i+L,j+L)
          s_ab_irrot(2,1,i,j+L) = -1 * s_ab_irrot(2,1,i+L,j+L)
          s_ab_irrot(2,2,i,j+L) = s_ab_irrot(2,2,i+L,j+L)

          s_ab_irrot(1,1,i,j+2*L) = s_ab_irrot(1,1,i+L,j+L)
          s_ab_irrot(1,2,i,j+2*L) = -1 * s_ab_irrot(1,2,i+L,j+L)
          s_ab_irrot(2,1,i,j+2*L) = -1 * s_ab_irrot(2,1,i+L,j+L)
          s_ab_irrot(2,2,i,j+2*L) = s_ab_irrot(2,2,i+L,j+L)

          s_ab_irrot(1,1,i+L,j) = s_ab_irrot(1,1,i+L,j+L)
          s_ab_irrot(1,2,i+L,j) = -1 * s_ab_irrot(1,2,i+L,j+L)
          s_ab_irrot(2,1,i+L,j) = -1 * s_ab_irrot(2,1,i+L,j+L)
          s_ab_irrot(2,2,i+L,j) = s_ab_irrot(2,2,i+L,j+L)

          s_ab_irrot(1,1,i+L,j+2*L) = s_ab_irrot(1,1,i+L,j+L)
          s_ab_irrot(1,2,i+L,j+2*L) = -1 * s_ab_irrot(1,2,i+L,j+L)
          s_ab_irrot(2,1,i+L,j+2*L) = -1 * s_ab_irrot(2,1,i+L,j+L)
          s_ab_irrot(2,2,i+L,j+2*L) = s_ab_irrot(2,2,i+L,j+L)

        else

          rho_k_p(i,j) = rho_k_p(i+L,j+L)
          rho_k_m(i,j) = rho_k_m(i+L,j+L)
          ch_ch(i,j) = ch_ch(i+L,j+L)

          s_ab(1,1,i,j) = s_ab(1,1,i+L,j+L)
          s_ab(1,2,i,j) = s_ab(1,2,i+L,j+L)
          s_ab(2,1,i,j) = s_ab(2,1,i+L,j+L)
          s_ab(2,2,i,j) = s_ab(2,2,i+L,j+L)

          s_ab(1,1,i+L,j) = s_ab(1,1,i+L,j+L)
          s_ab(1,2,i+L,j) = (-1) * s_ab(1,2,i+L,j+L)
          s_ab(2,1,i+L,j) = (-1) * s_ab(2,1,i+L,j+L)
          s_ab(2,2,i+L,j) = s_ab(2,2,i+L,j+L)

          s_ab(1,1,i,j+L) = s_ab(1,1,i+L,j+L)
          s_ab(1,2,i,j+L) = (-1) * s_ab(1,2,i+L,j+L)
          s_ab(2,1,i,j+L) = (-1) * s_ab(2,1,i+L,j+L)
          s_ab(2,2,i,j+L) = s_ab(2,2,i+L,j+L)

          ! s_ab_rot
          s_ab_rot(1,1,i,j) = s_ab_rot(1,1,i+L,j+L)
          s_ab_rot(1,2,i,j) = s_ab_rot(1,2,i+L,j+L)
          s_ab_rot(2,1,i,j) = s_ab_rot(2,1,i+L,j+L)
          s_ab_rot(2,2,i,j) = s_ab_rot(2,2,i+L,j+L)

          s_ab_rot(1,1,i+L,j) = s_ab_rot(1,1,i+L,j+L)
          s_ab_rot(1,2,i+L,j) = (-1) * s_ab_rot(1,2,i+L,j+L)
          s_ab_rot(2,1,i+L,j) = (-1) * s_ab_rot(2,1,i+L,j+L)
          s_ab_rot(2,2,i+L,j) = s_ab_rot(2,2,i+L,j+L)

          s_ab_rot(1,1,i,j+L) = s_ab_rot(1,1,i+L,j+L)
          s_ab_rot(1,2,i,j+L) = (-1) * s_ab_rot(1,2,i+L,j+L)
          s_ab_rot(2,1,i,j+L) = (-1) * s_ab_rot(2,1,i+L,j+L)
          s_ab_rot(2,2,i,j+L) = s_ab_rot(2,2,i+L,j+L)

          ! s_ab_irrot
          s_ab_irrot(1,1,i,j) = s_ab_irrot(1,1,i+L,j+L)
          s_ab_irrot(1,2,i,j) = s_ab_irrot(1,2,i+L,j+L)
          s_ab_irrot(2,1,i,j) = s_ab_irrot(2,1,i+L,j+L)
          s_ab_irrot(2,2,i,j) = s_ab_irrot(2,2,i+L,j+L)

          s_ab_irrot(1,1,i+L,j) = s_ab_irrot(1,1,i+L,j+L)
          s_ab_irrot(1,2,i+L,j) = (-1) * s_ab_irrot(1,2,i+L,j+L)
          s_ab_irrot(2,1,i+L,j) = (-1) * s_ab_irrot(2,1,i+L,j+L)
          s_ab_irrot(2,2,i+L,j) = s_ab_irrot(2,2,i+L,j+L)

          s_ab_irrot(1,1,i,j+L) = s_ab_irrot(1,1,i+L,j+L)
          s_ab_irrot(1,2,i,j+L) = (-1) * s_ab_irrot(1,2,i+L,j+L)
          s_ab_irrot(2,1,i,j+L) = (-1) * s_ab_irrot(2,1,i+L,j+L)
          s_ab_irrot(2,2,i,j+L) = s_ab_irrot(2,2,i+L,j+L)

        end if

      end do
    end do

  end subroutine fix_arrays

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
    real(kind=rk), dimension(2,2,(bz*L)+1,(bz*L)+1) :: chi_ab,&
    chi_ab_rot, chi_ab_irrot
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

    ! write (30,'(a)') "   # hop acceptance"
    ! write (30,'(2ES18.9)') temp,&
    ! (dble(accepts(1)) / dble(attempts(1)))
    ! if (add_charges.ne.0) then
    !   write (30,'(a,2i12.1,es18.9)') "# Hops: total, attempts, rate: ",&
    !   accepts(1), attempts(1), dble(accepts(1)) / dble(attempts(1))
    ! end if
    ! write (30,'(a,2i12.1,es18.9)') "# Rot.: total, attempts, rate: ",&
    ! accepts(2), attempts(2), dble(accepts(2)) / dble(attempts(2))
    ! write (30,'(a,2i12.1,es18.9)') "# Harm: total, attempts, rate: ",&
    ! accepts(3), attempts(3), dble(accepts(3)) / dble(attempts(3))
    ! ! write (30,'(a,2i12.1,es18.9)') "# Creations: total, attempts, rate: ",&
    ! ! accepts(4), attempts(4), dble(accepts(4)) / dble(attempts(4))
    ! ! write (30,'(a,2i12.1,es18.9)') "# Annihilations: total, attempts, rate: ",&
    ! ! accepts(5), attempts(5), dble(accepts(5)) / dble(attempts(5))
    ! write (30,'(a,2i12.1,es18.9)') "# Harmonic fluctuations: &
    ! &total, attempts, rate: ",&
    ! accepts(6), attempts(6), dble(accepts(6)) / dble(attempts(6))

    ! close(30)

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

      ! renormalise s_ab tensors here: then it propagates through to
      ! s_perp and s_par
      ! s_ab = s_ab * 2 * L**2
      ! s_ab_rot = s_ab_rot * 2 * L**2
      ! s_ab_irrot = s_ab_irrot * 2 * L**2
      s_ab = s_ab * L**2
      s_ab_rot = s_ab_rot * L**2
      s_ab_irrot = s_ab_irrot * L**2

      chi_ab = s_ab / temp
      chi_ab_rot = s_ab_rot / temp
      chi_ab_irrot = s_ab_irrot / temp

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

          s_perp(i,j) = (1 - kx_float*kx_float*norm_k) *   real(s_ab(1,1,kx,ky))+&
                        ((-1)*kx_float*ky_float*norm_k) *  real(s_ab(1,2,kx,ky))+&
                        ((-1)*ky_float*kx_float*norm_k) *  real(s_ab(2,1,kx,ky))+&
                        (1 - ky_float*ky_float*norm_k) *   real(s_ab(2,2,kx,ky))

          s_perp_irrot(i,j) = (1 - kx_float*kx_float*norm_k) *  real(s_ab_irrot(1,1,kx,ky))+&
                          ((-1)*kx_float*ky_float*norm_k) *     real(s_ab_irrot(1,2,kx,ky))+&
                          ((-1)*ky_float*kx_float*norm_k) *     real(s_ab_irrot(2,1,kx,ky))+&
                          (1 - ky_float*ky_float*norm_k) *      real(s_ab_irrot(2,2,kx,ky))

          s_perp_rot(i,j) = (1 - kx_float*kx_float*norm_k) * real(s_ab_rot(1,1,kx,ky))+&
                          ((-1)*kx_float*ky_float*norm_k) *  real(s_ab_rot(1,2,kx,ky))+&
                          ((-1)*ky_float*kx_float*norm_k) *  real(s_ab_rot(2,1,kx,ky))+&
                          (1 - ky_float*ky_float*norm_k) *   real(s_ab_rot(2,2,kx,ky))

          s_par(i,j) = (kx_float*kx_float*norm_k) *   real(s_ab(1,1,kx,ky))+&
                         (kx_float*ky_float*norm_k) * real(s_ab(1,2,kx,ky))+&
                         (ky_float*kx_float*norm_k) * real(s_ab(2,1,kx,ky))+&
                         (ky_float*ky_float*norm_k) * real(s_ab(2,2,kx,ky))

          s_par_irrot(i,j) = (kx_float*kx_float*norm_k)*  real(s_ab_irrot(1,1,kx,ky))+&
                               (kx_float*ky_float*norm_k)*real(s_ab_irrot(1,2,kx,ky))+&
                               (ky_float*kx_float*norm_k)*real(s_ab_irrot(2,1,kx,ky))+&
                               (ky_float*ky_float*norm_k)*real(s_ab_irrot(2,2,kx,ky))

          s_par_rot(i,j) = (kx_float*kx_float*norm_k)*  real(s_ab_rot(1,1,kx,ky))+&
                             (kx_float*ky_float*norm_k)*real(s_ab_rot(1,2,kx,ky))+&
                             (ky_float*kx_float*norm_k)*real(s_ab_rot(2,1,kx,ky))+&
                             (ky_float*ky_float*norm_k)*real(s_ab_rot(2,2,kx,ky))

        end do
      end do ! end p, m loops

      open(unit=10, file=dir_st_file)
      open(unit=11, file=dir_dist_file)
      open(unit=12, file=charge_st_file)
      open(unit=14, file=s_ab_file)
      open(unit=15, file=s_perp_file)
      open(unit=17, file=irrot_sab_file)
      open(unit=18, file=irrot_sperp_file)
      open(unit=20, file=rot_sab_file)
      open(unit=21, file=rot_sperp_file)
      open(unit=22, file=spar_file)
      open(unit=23, file=irrot_spar_file)
      open(unit=24, file=rot_spar_file)
      open(unit=25, file=avg_field_file)
      open(unit=26, file=chi_ab_file)
      open(unit=27, file=irrot_chi_ab_file)
      open(unit=28, file=rot_chi_ab_file)
      open(unit=38, file=windings_file)
      open(unit=39, file=windings_sq_file)

      open  (30, file=sphe_sus_file, position='append')
      write (30, '(a)') "# S_ab integrals (* L**2)!"
      write (30, '(a)') "# S_xx: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(s_ab(1,1,:,:))),&
      sum(real(s_ab_rot(1,1,:,:))), sum(real(s_ab_irrot(1,1,:,:)))
      write (30, '(a)') "# S_xy: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(s_ab(1,2,:,:))),&
      sum(real(s_ab_rot(1,2,:,:))), sum(real(s_ab_irrot(1,2,:,:)))
      write (30, '(a)') "# S_yx: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(s_ab(2,1,:,:))),&
      sum(real(s_ab_rot(2,1,:,:))), sum(real(s_ab_irrot(2,1,:,:)))
      write (30, '(a)') "# S_yy: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(s_ab(2,2,:,:))),&
      sum(real(s_ab_rot(2,2,:,:))), sum(real(s_ab_irrot(2,2,:,:)))
      write (30, '(a)') "# S_perp integrals: total, rot, irrot"
      write (30, '(3f20.8)') sum(s_perp), sum(s_perp_rot), sum(s_perp_irrot)
      write (30, '(a)') "# Chi_ab integrals (* L**2)!"
      write (30, '(a)') "# Chi_xx: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(chi_ab(1,1,:,:))),&
      sum(real(chi_ab_rot(1,1,:,:))), sum(real(chi_ab_irrot(1,1,:,:)))
      write (30, '(a)') "# Chi_xy: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(chi_ab(1,2,:,:))),&
      sum(real(chi_ab_rot(1,2,:,:))), sum(real(chi_ab_irrot(1,2,:,:)))
      write (30, '(a)') "# Chi_yx: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(chi_ab(2,1,:,:))),&
      sum(real(chi_ab_rot(2,1,:,:))), sum(real(chi_ab_irrot(2,1,:,:)))
      write (30, '(a)') "# Chi_yy: total, rot, irrot"
      write (30, '(3f20.8)') sum(real(chi_ab(2,2,:,:))),&
      sum(real(chi_ab_rot(2,2,:,:))), sum(real(chi_ab_irrot(2,2,:,:)))
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

      open(unit=57, file="fftw_sab.dat")

      ! do i = 1, L/2 + 1
      !   do j = 1, L
      do i = 1, bz*L + 1
        do j = 1, bz*L + 1

          write (57, field_format_string)&
          2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
          2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
          real(fftw_s_ab_total(1,1,i,j)),&
          real(fftw_s_ab_total(1,2,i,j)),&
          real(fftw_s_ab_total(2,1,i,j)),&
          real(fftw_s_ab_total(2,2,i,j))

        end do
      end do

      close(57)

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

            write (14, field_format_string)&
            2*pi*(i - 1 - bz*(l/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(l/2))/(L*lambda),&
            real(s_ab(1,1,i,j)),&
            real(s_ab(1,2,i,j)),&
            real(s_ab(2,1,i,j)),&
            real(s_ab(2,2,i,j))

            write (26, field_format_string)&
            2*pi*(i - 1 - bz*(l/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(l/2))/(L*lambda),&
            real(chi_ab(1,1,i,j)),&
            real(chi_ab(1,2,i,j)),&
            real(chi_ab(2,1,i,j)),&
            real(chi_ab(2,2,i,j))

            write (17, field_format_string)&
            2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
            real(s_ab_irrot(1,1,i,j)),&
            real(s_ab_irrot(1,2,i,j)),&
            real(s_ab_irrot(2,1,i,j)),&
            real(s_ab_irrot(2,2,i,j))

            write (27, field_format_string)&
            2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
            real(chi_ab_irrot(1,1,i,j)),&
            real(chi_ab_irrot(1,2,i,j)),&
            real(chi_ab_irrot(2,1,i,j)),&
            real(chi_ab_irrot(2,2,i,j))

            write (20,field_format_string)&
            2*pi*(i - 1 - bz*(l/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(l/2))/(L*lambda),&
            real(s_ab_rot(1,1,i,j)),&
            real(s_ab_rot(1,2,i,j)),&
            real(s_ab_rot(2,1,i,j)),&
            real(s_ab_rot(2,2,i,j))

            write (28,field_format_string)&
            2*pi*(i - 1 - bz*(l/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(l/2))/(L*lambda),&
            real(chi_ab_rot(1,1,i,j)),&
            real(chi_ab_rot(1,2,i,j)),&
            real(chi_ab_rot(2,1,i,j)),&
            real(chi_ab_rot(2,2,i,j))

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
      close(14)
      close(15)
      close(17)
      close(18)
      close(20)
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)

      deallocate(s_perp); deallocate(s_perp_rot); deallocate(s_perp_irrot);
      deallocate(s_par); deallocate(s_par_rot); deallocate(s_par_irrot);
    end if ! do_corr

  end subroutine calc_correlations

end module output
