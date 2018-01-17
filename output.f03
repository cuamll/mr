module output
  use common
  implicit none

  contains

  subroutine write_output

    call fix_arrays
    call calc_correlations

  end subroutine write_output

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
    write (30,'(a)') "# Avg. x-component: total, rot., irrot.:"
    write (30, '(3f18.10)') avg_field_total(1), avg_field_rot(1),&
                            avg_field_irrot(1)
    write (30,'(a)') "# Avg. y-component: total, rot., irrot.:"
    write (30, '(3f18.10)') avg_field_total(2), avg_field_rot(2),&
                            avg_field_irrot(2)

    ! write (30,'(a)') "   # hop acceptance"
    ! write (30,'(2ES18.9)') temp,&
    ! (dble(accepts(1)) / dble(attempts(1)))
    if (add_charges.ne.0) then
      write (30,'(a,2i12.1,es18.9)') "# Hops: total, attempts, rate: ",&
      accepts(1), attempts(1), dble(accepts(1)) / dble(attempts(1))
    end if
    write (30,'(a,2i12.1,es18.9)') "# Rot.: total, attempts, rate: ",&
    accepts(2), attempts(2), dble(accepts(2)) / dble(attempts(2))
    write (30,'(a,2i12.1,es18.9)') "# Harm: total, attempts, rate: ",&
    accepts(3), attempts(3), dble(accepts(3)) / dble(attempts(3))
    ! write (30,'(a,2i12.1,es18.9)') "# Creations: total, attempts, rate: ",&
    ! accepts(4), attempts(4), dble(accepts(4)) / dble(attempts(4))
    ! write (30,'(a,2i12.1,es18.9)') "# Annihilations: total, attempts, rate: ",&
    ! accepts(5), attempts(5), dble(accepts(5)) / dble(attempts(5))
    write (30,'(a,2i12.1,es18.9)') "# Harmonic fluctuations: &
    &total, attempts, rate: ",&
    accepts(6), attempts(6), dble(accepts(6)) / dble(attempts(6))

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

      ! renormalise s_ab tensors here: then it propagates through to
      ! s_perp and s_par
      s_ab = s_ab * L**2
      s_ab_rot = s_ab_rot * L**2
      s_ab_irrot = s_ab_irrot * L**2

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

      open  (30, file=sphe_sus_file, position='append')
      write (30, '(a)') "# S_ab integrals (* L**2)!"
      write (30, '(a)') "# S_xx: total, rot, irrot"
      write (30, '(3f18.10)') sum(real(s_ab(1,1,:,:))),&
      sum(real(s_ab_rot(1,1,:,:))), sum(real(s_ab_irrot(1,1,:,:)))
      write (30, '(a)') "# S_xy: total, rot, irrot"
      write (30, '(3f18.10)') sum(real(s_ab(1,2,:,:))),&
      sum(real(s_ab_rot(1,2,:,:))), sum(real(s_ab_irrot(1,2,:,:)))
      write (30, '(a)') "# S_yx: total, rot, irrot"
      write (30, '(3f18.10)') sum(real(s_ab(2,1,:,:))),&
      sum(real(s_ab_rot(2,1,:,:))), sum(real(s_ab_irrot(2,1,:,:)))
      write (30, '(a)') "# S_yy: total, rot, irrot"
      write (30, '(3f18.10)') sum(real(s_ab(2,2,:,:))),&
      sum(real(s_ab_rot(2,2,:,:))), sum(real(s_ab_irrot(2,2,:,:)))
      write (30, '(a)') "# S_perp integrals: total, rot, irrot"
      write (30, '(3f18.10)') sum(s_perp), sum(s_perp_rot), sum(s_perp_irrot)
      close (30)

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

            write (14, field_format_string)&
            2*pi*(i - 1 - bz*(l/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(l/2))/(L*lambda),&
            real(s_ab(1,1,i,j)),&
            real(s_ab(1,2,i,j)),&
            real(s_ab(2,1,i,j)),&
            real(s_ab(2,2,i,j))

            write (17, field_format_string)&
            2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
            real(s_ab_irrot(1,1,i,j)),&
            real(s_ab_irrot(1,2,i,j)),&
            real(s_ab_irrot(2,1,i,j)),&
            real(s_ab_irrot(2,2,i,j))

            write (20,field_format_string)&
            2*pi*(i - 1 - bz*(l/2))/(L*lambda),&
            2*pi*(j - 1 - bz*(l/2))/(L*lambda),&
            real(s_ab_rot(1,1,i,j)),&
            real(s_ab_rot(1,2,i,j)),&
            real(s_ab_rot(2,1,i,j)),&
            real(s_ab_rot(2,2,i,j))

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

      deallocate(s_perp); deallocate(s_perp_rot); deallocate(s_perp_irrot);
      deallocate(s_par); deallocate(s_par_rot); deallocate(s_par_irrot);
    end if ! do_corr

  end subroutine calc_correlations

end module output
