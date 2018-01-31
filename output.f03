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
    integer :: i, j, k, m, p
    
    ! very long and extremely disgusting way of filling up the arrays
    ! so we can then calculate other correlation functions easily.
    ! probably won't work yet. also can probably be simplified a lot

    do k = 1,L
      do j = 1,L
        do i = 1,L

          if (i.eq.1.and.j.eq.1.and.k.eq.1) then

            ! charge-charge ones
            rho_k_p(1,1,1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(1,1,1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(1,1,L+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(1,1,L+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(1,1,(bz*L)+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(1,1,(bz*L)+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(1,L+1,1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(1,L+1,1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(1,(bz*L)+1,1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(1,(bz*L)+1,1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(1,L+1,L+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(1,L+1,L+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(1,(bz*L)+1,L+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(1,(bz*L)+1,L+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(1,L+1,(bz*L)+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(1,L+1,(bz*L)+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(1,(bz*L)+1,(bz*L)+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(1,(bz*L)+1,(bz*L)+1) = rho_k_m(L+1,L+1,L+1)

            rho_k_p(L+1,1,1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(L+1,1,1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(L+1,1,L+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(L+1,1,L+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(L+1,1,(bz*L)+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(L+1,1,(bz*L)+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(L+1,L+1,1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(L+1,L+1,1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(L+1,(bz*L)+1,1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(L+1,(bz*L)+1,1) = rho_k_m(L+1,L+1,L+1)
            ! rho_k_p(L+1,L+1,L+1) = rho_k_p(L+1,L+1,L+1)
            ! rho_k_m(L+1,L+1,L+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(L+1,(bz*L)+1,L+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(L+1,(bz*L)+1,L+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(L+1,L+1,(bz*L)+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(L+1,L+1,(bz*L)+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(L+1,(bz*L)+1,(bz*L)+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(L+1,(bz*L)+1,(bz*L)+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p(L+1,1,1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m(L+1,1,1) = rho_k_m(L+1,L+1,L+1)

            rho_k_p((bz*L)+1,1,1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m((bz*L)+1,1,1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p((bz*L)+1,1,L+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m((bz*L)+1,1,L+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p((bz*L)+1,1,(bz*L)+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m((bz*L)+1,1,(bz*L)+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p((bz*L)+1,L+1,1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m((bz*L)+1,L+1,1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p((bz*L)+1,(bz*L)+1,1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m((bz*L)+1,(bz*L)+1,1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p((bz*L)+1,L+1,L+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m((bz*L)+1,L+1,L+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p((bz*L)+1,(bz*L)+1,L+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m((bz*L)+1,(bz*L)+1,L+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p((bz*L)+1,L+1,(bz*L)+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m((bz*L)+1,L+1,(bz*L)+1) = rho_k_m(L+1,L+1,L+1)
            rho_k_p((bz*L)+1,(bz*L)+1,(bz*L)+1) = rho_k_p(L+1,L+1,L+1)
            rho_k_m((bz*L)+1,(bz*L)+1,(bz*L)+1) = rho_k_m(L+1,L+1,L+1)

            ! s_ab
            do m=1,3
              do p=1,3
                s_ab(m,p,1,1,1)               = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,1,1,L+1)             = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,1,1,(bz*L)+1)        = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,1,L+1,1)             = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,1,L+1,L+1)           = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,1,L+1,(bz*L)+1)      = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,1,(bz*L)+1,1)        = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,1,(bz*L)+1,L+1)      = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,1,(bz*L)+1,(bz*L)+1) = s_ab(m,p,L+1,L+1,L+1)

                s_ab(m,p,L+1,1,1)               = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,L+1,1,L+1)             = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,L+1,1,(bz*L)+1)        = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,L+1,L+1,1)             = s_ab(m,p,L+1,L+1,L+1)
                !s_ab(m,p,L+1,L+1,L+1)          = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,L+1,L+1,(bz*L)+1)      = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,L+1,(bz*L)+1,1)        = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,L+1,(bz*L)+1,L+1)      = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,L+1,(bz*L)+1,(bz*L)+1) = s_ab(m,p,L+1,L+1,L+1)

                s_ab(m,p,(bz*L)+1,1,1)               = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,(bz*L)+1,1,L+1)             = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,(bz*L)+1,1,(bz*L)+1)        = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,(bz*L)+1,L+1,1)             = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,(bz*L)+1,L+1,L+1)           = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,(bz*L)+1,L+1,(bz*L)+1)      = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,(bz*L)+1,(bz*L)+1,1)        = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,(bz*L)+1,(bz*L)+1,L+1)      = s_ab(m,p,L+1,L+1,L+1)
                s_ab(m,p,(bz*L)+1,(bz*L)+1,(bz*L)+1) = s_ab(m,p,L+1,L+1,L+1)

                ! s_ab_rot
                s_ab_rot(m,p,1,1,1)               = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,1,1,L+1)             = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,1,1,(bz*L)+1)        = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,1,L+1,1)             = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,1,L+1,L+1)           = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,1,L+1,(bz*L)+1)      = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,1,(bz*L)+1,1)        = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,1,(bz*L)+1,L+1)      = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,1,(bz*L)+1,(bz*L)+1) = s_ab_rot(m,p,L+1,L+1,L+1)

                s_ab_rot(m,p,L+1,1,1)               = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,L+1,1,L+1)             = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,L+1,1,(bz*L)+1)        = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,L+1,L+1,1)             = s_ab_rot(m,p,L+1,L+1,L+1)
                !s_ab_rot(m,p,L+1,L+1,L+1)          = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,L+1,L+1,(bz*L)+1)      = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,L+1,(bz*L)+1,1)        = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,L+1,(bz*L)+1,L+1)      = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,L+1,(bz*L)+1,(bz*L)+1) = s_ab_rot(m,p,L+1,L+1,L+1)

                s_ab_rot(m,p,(bz*L)+1,1,1)               = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,(bz*L)+1,1,L+1)             = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,(bz*L)+1,1,(bz*L)+1)        = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,(bz*L)+1,L+1,1)             = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,(bz*L)+1,L+1,L+1)           = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,(bz*L)+1,L+1,(bz*L)+1)      = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,(bz*L)+1,(bz*L)+1,1)        = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,(bz*L)+1,(bz*L)+1,L+1)      = s_ab_rot(m,p,L+1,L+1,L+1)
                s_ab_rot(m,p,(bz*L)+1,(bz*L)+1,(bz*L)+1) = s_ab_rot(m,p,L+1,L+1,L+1)

                ! s_ab_irrot
                s_ab_irrot(m,p,1,1,1)               = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,1,1,L+1)             = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,1,1,(bz*L)+1)        = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,1,L+1,1)             = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,1,L+1,L+1)           = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,1,L+1,(bz*L)+1)      = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,1,(bz*L)+1,1)        = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,1,(bz*L)+1,L+1)      = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,1,(bz*L)+1,(bz*L)+1) = s_ab_irrot(m,p,L+1,L+1,L+1)

                s_ab_irrot(m,p,L+1,1,1)               = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,L+1,1,L+1)             = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,L+1,1,(bz*L)+1)        = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,L+1,L+1,1)             = s_ab_irrot(m,p,L+1,L+1,L+1)
                !s_ab_irrot(m,p,L+1,L+1,L+1)          = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,L+1,L+1,(bz*L)+1)      = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,L+1,(bz*L)+1,1)        = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,L+1,(bz*L)+1,L+1)      = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,L+1,(bz*L)+1,(bz*L)+1) = s_ab_irrot(m,p,L+1,L+1,L+1)

                s_ab_irrot(m,p,(bz*L)+1,1,1)               = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,(bz*L)+1,1,L+1)             = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,(bz*L)+1,1,(bz*L)+1)        = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,(bz*L)+1,L+1,1)             = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,(bz*L)+1,L+1,L+1)           = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,(bz*L)+1,L+1,(bz*L)+1)      = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,(bz*L)+1,(bz*L)+1,1)        = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,(bz*L)+1,(bz*L)+1,L+1)      = s_ab_irrot(m,p,L+1,L+1,L+1)
                s_ab_irrot(m,p,(bz*L)+1,(bz*L)+1,(bz*L)+1) = s_ab_irrot(m,p,L+1,L+1,L+1)

              end do
            end do

          ! rows in the z-direction
          else if (i.eq.1.and.j.eq.1.and.k.gt.1) then

            rho_k_p(i,j,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j+L,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j+L,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j+L,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j+2*L,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j+2*L,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j+2*L,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j+L,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j+L,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j+L,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j+2*L,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j+2*L,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j+2*L,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+2*L,j,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+2*L,j,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+2*L,j,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+2*L,j+L,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+2*L,j+L,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+2*L,j+L,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+2*L,j+2*L,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+2*L,j+2*L,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+2*L,j+2*L,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j+L,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j+L,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j+L,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j+2*L,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j+2*L,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j+2*L,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j+L,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j+L,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j+L,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j+2*L,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j+2*L,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j+2*L,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+2*L,j,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+2*L,j,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+2*L,j,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+2*L,j+L,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+2*L,j+L,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+2*L,j+L,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+2*L,j+2*L,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+2*L,j+2*L,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+2*L,j+2*L,k+L) = ch_ch(i+L,j+L,k+L)

            do m = 1,3
              do p = 1,3

                if (m.eq.p) then

                  s_ab(m,p,i,j,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+2*L,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+2*L,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j+L,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j+2*L,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+2*L,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+2*L,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j+L,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j+2*L,k+L) = s_ab(m,p,i+L,j+L,k+L)

                  s_ab_rot(m,p,i,j,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+2*L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+2*L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j+L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j+2*L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+2*L,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+2*L,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j+L,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j+2*L,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)


                  s_ab_irrot(m,p,i,j,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+2*L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+2*L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j+L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j+2*L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+2*L,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+2*L,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j+L,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j+2*L,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)

                else

                  s_ab(m,p,i,j,k)         = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k)       = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+2*L,k)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k)       = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+2*L,k)   = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j,k)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j+L,k)   = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j+2*L,k) = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j,k+L)         = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k+L)       = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+2*L,k+L)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k+L)       = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k+L)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+2*L,k+L)   = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j,k+L)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j+L,k+L)   = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j+2*L,k+L) = (-1) * s_ab(m,p,i+L,j+L,k+L)

                  s_ab_rot(m,p,i,j,k)         = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k)       = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+2*L,k)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k)       = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+2*L,k)   = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j,k)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j+L,k)   = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j+2*L,k) = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j,k+L)         = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k+L)       = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+2*L,k+L)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k+L)       = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k+L)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+2*L,k+L)   = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j,k+L)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j+L,k+L)   = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j+2*L,k+L) = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)

                  s_ab_irrot(m,p,i,j,k)         = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k)       = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+2*L,k)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k)       = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+2*L,k)   = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j,k)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j+L,k)   = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j+2*L,k) = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j,k+L)         = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k+L)       = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+2*L,k+L)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k+L)       = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k+L)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+2*L,k+L)   = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j,k+L)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j+L,k+L)   = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j+2*L,k+L) = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)

                end if

              end do
            end do

          ! rows in the y direction
          else if (i.eq.1.and.j.gt.1.and.k.eq.1) then

            rho_k_p(i,j,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j,k+2*L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j,k+2*L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j,k+2*L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j,k+2*L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j,k+2*L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j,k+2*L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+2*L,j,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+2*L,j,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+2*L,j,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+2*L,j,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+2*L,j,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+2*L,j,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+2*L,j,k+2*L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+2*L,j,k+2*L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+2*L,j,k+2*L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j+L,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j+L,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j+L,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j+L,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j+L,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j+L,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j+L,k+2*L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j+L,k+2*L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j+L,k+2*L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j+L,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j+L,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j+L,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j+L,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j+L,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j+L,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j+L,k+2*L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j+L,k+2*L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j+L,k+2*L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+2*L,j+L,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+2*L,j+L,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+2*L,j+L,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+2*L,j+L,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+2*L,j+L,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+2*L,j+L,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+2*L,j+L,k+2*L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+2*L,j+L,k+2*L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+2*L,j+L,k+2*L) = ch_ch(i+L,j+L,k+L)

            do m = 1,3
              do p = 1,3

                if (m.eq.p) then

                  s_ab(m,p,i,j,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j,k+2*L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k+2*L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j,k+2*L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k+2*L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k+2*L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j+L,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j+L,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j+L,k+2*L) = s_ab(m,p,i+L,j+L,k+L)

                  s_ab_rot(m,p,i,j,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j,k+2*L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k+2*L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j,k+2*L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k+2*L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k+2*L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j+L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j+L,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j+L,k+2*L) = s_ab_rot(m,p,i+L,j+L,k+L)

                  s_ab_irrot(m,p,i,j,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j,k+2*L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k+2*L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j,k+2*L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k+2*L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k+2*L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j+L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j+L,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j+L,k+2*L) = s_ab_irrot(m,p,i+L,j+L,k+L)

                else


                  s_ab(m,p,i,j,k)           = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j,k+L)         = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j,k+2*L)       = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k)         = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k+L)       = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k+2*L)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j,k)       = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j,k+L)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j,k+2*L)   = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k)         = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k+L)       = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k+2*L)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k)       = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k+L)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k+2*L)   = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j+L,k)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j+L,k+L)   = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+2*L,j+L,k+2*L) = (-1) * s_ab(m,p,i+L,j+L,k+L)

                  s_ab_rot(m,p,i,j,k)           = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j,k+L)         = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j,k+2*L)       = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k)         = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k+L)       = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k+2*L)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j,k)       = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j,k+L)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j,k+2*L)   = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k)         = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k+L)       = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k+2*L)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k)       = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k+L)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k+2*L)   = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j+L,k)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j+L,k+L)   = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+2*L,j+L,k+2*L) = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)

                  s_ab_irrot(m,p,i,j,k)           = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j,k+L)         = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j,k+2*L)       = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k)         = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k+L)       = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k+2*L)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j,k)       = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j,k+L)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j,k+2*L)   = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k)         = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k+L)       = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k+2*L)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k)       = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k+L)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k+2*L)   = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j+L,k)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j+L,k+L)   = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+2*L,j+L,k+2*L) = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)

                end if

              end do
            end do

          ! rows in x direction
          else if (i.gt.1.and.j.eq.1.and.k.eq.1) then

            rho_k_p(i,j,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j,k+2*L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j,k+2*L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j,k+2*L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j+L,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j+L,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j+L,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j+L,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j+L,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j+L,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j+L,k+2*L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j+L,k+2*L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j+L,k+2*L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j+2*L,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j+2*L,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j+2*L,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j+2*L,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j+2*L,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j+2*L,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i,j+2*L,k+2*L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j+2*L,k+2*L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j+2*L,k+2*L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j,k+2*L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j,k+2*L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j,k+2*L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j+L,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j+L,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j+L,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j+L,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j+L,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j+L,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j+L,k+2*L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j+L,k+2*L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j+L,k+2*L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j+2*L,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j+2*L,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j+2*L,k) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j+2*L,k+L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j+2*L,k+L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j+2*L,k+L) = ch_ch(i+L,j+L,k+L)

            rho_k_p(i+L,j+2*L,k+2*L) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i+L,j+2*L,k+2*L) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i+L,j+2*L,k+2*L) = ch_ch(i+L,j+L,k+L)

            do m = 1,3
              do p = 1,3

                if (m.eq.p) then

                  s_ab(m,p,i,j,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j,k+2*L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k+2*L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+2*L,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+2*L,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+2*L,k+2*L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k+2*L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k+2*L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+2*L,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+2*L,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+2*L,k+2*L) = s_ab(m,p,i+L,j+L,k+L)

                  s_ab_rot(m,p,i,j,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j,k+2*L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k+2*L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+2*L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+2*L,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+2*L,k+2*L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k+2*L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k+2*L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+2*L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+2*L,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+2*L,k+2*L) = s_ab_rot(m,p,i+L,j+L,k+L)

                  s_ab_irrot(m,p,i,j,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j,k+2*L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k+2*L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+2*L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+2*L,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+2*L,k+2*L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k+2*L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k+2*L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+2*L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+2*L,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+2*L,k+2*L) = s_ab_irrot(m,p,i+L,j+L,k+L)

                else


                  s_ab(m,p,i,j,k)           = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j,k+L)         = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j,k+2*L)       = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k)         = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k+L)       = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k+2*L)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+2*L,k)       = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+2*L,k+L)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+2*L,k+2*L)   = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k)         = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k+L)       = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k+2*L)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k)       = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k+L)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k+2*L)   = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+2*L,k)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+2*L,k+L)   = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+2*L,k+2*L) = (-1) * s_ab(m,p,i+L,j+L,k+L)

                  s_ab_rot(m,p,i,j,k)           = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j,k+L)         = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j,k+2*L)       = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k)         = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k+L)       = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k+2*L)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+2*L,k)       = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+2*L,k+L)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+2*L,k+2*L)   = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k)         = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k+L)       = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k+2*L)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k)       = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k+L)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k+2*L)   = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+2*L,k)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+2*L,k+L)   = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+2*L,k+2*L) = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)

                  s_ab_irrot(m,p,i,j,k)           = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j,k+L)         = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j,k+2*L)       = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k)         = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k+L)       = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k+2*L)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+2*L,k)       = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+2*L,k+L)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+2*L,k+2*L)   = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k)         = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k+L)       = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k+2*L)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k)       = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k+L)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k+2*L)   = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+2*L,k)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+2*L,k+L)   = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+2*L,k+2*L) = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)

                end if

              end do
            end do

          else

            rho_k_p(i,j,k) = rho_k_p(i+L,j+L,k+L)
            rho_k_m(i,j,k) = rho_k_m(i+L,j+L,k+L)
            ch_ch(i,j,k) = ch_ch(i+L,j+L,k+L)

            do m=1,3
              do p=1,3
                ! these ones have the same amplitude
                s_ab(m,p,i,j,k+L) = s_ab(m,p,i+L,j+L,k+L)
                s_ab(m,p,i,j+L,k) = s_ab(m,p,i+L,j+L,k+L)
                s_ab(m,p,i+L,j,k) = s_ab(m,p,i+L,j+L,k+L)

                s_ab_rot(m,p,i,j,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                s_ab_rot(m,p,i,j+L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                s_ab_rot(m,p,i+L,j,k) = s_ab_rot(m,p,i+L,j+L,k+L)

                s_ab_irrot(m,p,i,j,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                s_ab_irrot(m,p,i,j+L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                s_ab_irrot(m,p,i+L,j,k) = s_ab_irrot(m,p,i+L,j+L,k+L)

                if (m.eq.p) then
                  ! opposite amplitudes on cross terms
                  s_ab(m,p,i,j,k)     = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k+L) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k) = s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k+L) = s_ab(m,p,i+L,j+L,k+L)

                  s_ab_rot(m,p,i,j,k)     = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k) = s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k+L) = s_ab_rot(m,p,i+L,j+L,k+L)

                  s_ab_irrot(m,p,i,j,k)     = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k) = s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k+L) = s_ab_irrot(m,p,i+L,j+L,k+L)

                else

                  s_ab(m,p,i,j,k)     = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i,j+L,k+L) = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j+L,k) = (-1) * s_ab(m,p,i+L,j+L,k+L)
                  s_ab(m,p,i+L,j,k+L) = (-1) * s_ab(m,p,i+L,j+L,k+L)

                  s_ab_rot(m,p,i,j,k)     = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i,j+L,k+L) = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j+L,k) = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)
                  s_ab_rot(m,p,i+L,j,k+L) = (-1) * s_ab_rot(m,p,i+L,j+L,k+L)

                  s_ab_irrot(m,p,i,j,k)     = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i,j+L,k+L) = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j+L,k) = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)
                  s_ab_irrot(m,p,i+L,j,k+L) = (-1) * s_ab_irrot(m,p,i+L,j+L,k+L)

                end if

              end do
            end do

          end if

        end do
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
    real(kind=rk), dimension(:,:,:), allocatable :: s_perp,&
    s_perp_irrot, s_perp_rot
    real(kind=rk), dimension(:,:,:), allocatable :: s_par,&
    s_par_irrot, s_par_rot
    real(kind=rk), dimension((bz*L)+1,(bz*L)+1,(bz*L)+1) :: charge_struc,&
    field_struc, field_struc_irrot, field_struc_rot
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
    write (30,'(a)') "# Avg. z-component: total, rot., irrot.:"
    write (30, '(3f18.10)') avg_field_total(3), avg_field_rot(3),&
                            avg_field_irrot(3)
    write (30,'(a)') "# Avg. x-component^2: total, rot., irrot.:"
    write (30, '(3es18.10)') avg_field_sq_total(1), avg_field_sq_rot(1),&
                            avg_field_sq_irrot(1)
    write (30,'(a)') "# Avg. y-component^2: total, rot., irrot.:"
    write (30, '(3es18.10)') avg_field_sq_total(2), avg_field_sq_rot(2),&
                            avg_field_sq_irrot(2)
    write (30,'(a)') "# Avg. z-component^2: total, rot., irrot.:"
    write (30, '(3es18.10)') avg_field_sq_total(3), avg_field_sq_rot(3),&
                            avg_field_sq_irrot(3)

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
      allocate(s_perp((sp*L)+1,(sp*L)+1,(sp*L)+1))
      allocate(s_perp_irrot((sp*L)+1,(sp*L)+1,(sp*L)+1))
      allocate(s_perp_rot((sp*L)+1,(sp*L)+1,(sp*L)+1))
      allocate(s_par((sp*L)+1,(sp*L)+1,(sp*L)+1))
      allocate(s_par_irrot((sp*L)+1,(sp*L)+1,(sp*L)+1))
      allocate(s_par_rot((sp*L)+1,(sp*L)+1,(sp*L)+1))
      s_perp = 0.0; s_perp_rot = 0.0; s_perp_irrot = 0.0;
      s_par = 0.0; s_par_rot = 0.0; s_par_irrot = 0.0;

      ! renormalise s_ab tensors here: then it propagates through to
      ! s_perp and s_par
      s_ab = s_ab * L**3
      s_ab_rot = s_ab_rot * L**3
      s_ab_irrot = s_ab_irrot * L**3

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
              field_struc(i,j,k) = abs(s_ab(1,1,i,j,k))
              field_struc_irrot(i,j,k) = abs(s_ab_irrot(1,1,i,j,k))
              field_struc_rot(i,j,k) = abs(s_ab_rot(1,1,i,j,k))
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

            ! move back to somewhere we already know about
            ! this is probably not right yet.
            ! but i know something needs doing here
            ! s_ab should be periodic in 2pi so mod should be fine

            ! if (abs(p).gt.(L/2)*bz) then
            !   ky = modulo(p,(bz*L)/2) + 1 + (bz*L)/2
            ! else
            !   ky = p + 1 + (bz*L)/2
            ! end if
            ! if (abs(m).gt.(L/2)*bz) then
            !   kx = modulo(m,(bz*L)/2) + 1 + (bz*L)/2
            ! else
            !   kx = m + 1 + (bz*L)/2
            ! end if

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

            !if (p.lt.(-1*(bz*(L/2)))) then
            !  ky = p
            !  do while (ky.lt.(-1*(bz*(L/2))))
            !    ky = ky + bz*(L/2)
            !  end do
            !  ! array index
            !  ky = ky + 1 + bz*(L/2)
            !else if (p.gt.(bz*(L/2))) then
            !  ky = p
            !  do while (ky.gt.bz*(L/2))
            !    ky = ky - bz*(L/2)
            !  end do
            !  ky = ky + 1 + bz*(L/2)
            !else
            !  ky = p + 1 + (bz*(L/2))
            !end if

            !if (m.lt.(-1*(bz*(L/2)))) then
            !  kx = m
            !  do while (kx.lt.(-1*(bz*(L/2))))
            !    kx = kx + bz*(L/2)
            !  end do
            !  ! array index
            !  kx = kx + 1 + bz*(L/2)
            !else if (m.gt.(bz*(L/2))) then
            !  kx = m
            !  do while (kx.gt.bz*(L/2))
            !    kx = kx - bz*(L/2)
            !  end do
            !  kx = kx + 1 + bz*(L/2)
            !else
            !  kx = m + 1 + (bz*(L/2))
            !end if

            s_perp(i,j,k) = (1 - kx_float*kx_float*norm_k) *&
                            &real(s_ab(1,1,kx,ky,kz))+&
            ((-1)*kx_float*ky_float*norm_k) *  real(s_ab(1,2,kx,ky,kz))+&
            ((-1)*kx_float*kz_float*norm_k) *  real(s_ab(1,3,kx,ky,kz))+&
            ((-1)*ky_float*kx_float*norm_k) *  real(s_ab(2,1,kx,ky,kz))+&
            (1 - ky_float*ky_float*norm_k) *   real(s_ab(2,2,kx,ky,kz))+&
            ((-1)*ky_float*kz_float*norm_k) *  real(s_ab(2,3,kx,ky,kz))+&
            ((-1)*kz_float*kx_float*norm_k) *  real(s_ab(1,3,kx,ky,kz))+&
            ((-1)*kz_float*ky_float*norm_k) *  real(s_ab(3,2,kx,ky,kz))+&
            (1 - kz_float*kz_float*norm_k) *   real(s_ab(3,3,kx,ky,kz))

            s_perp_irrot(i,j,k) = (1 - kx_float*kx_float*norm_k) *&
                                  &real(s_ab_irrot(1,1,kx,ky,kz))+&
            ((-1)*kx_float*ky_float*norm_k) *  real(s_ab_irrot(1,2,kx,ky,kz))+&
            ((-1)*kx_float*kz_float*norm_k) *  real(s_ab_irrot(1,3,kx,ky,kz))+&
            ((-1)*ky_float*kx_float*norm_k) *  real(s_ab_irrot(2,1,kx,ky,kz))+&
            (1 - ky_float*ky_float*norm_k) *   real(s_ab_irrot(2,2,kx,ky,kz))+&
            ((-1)*ky_float*kz_float*norm_k) *  real(s_ab_irrot(2,3,kx,ky,kz))+&
            ((-1)*kz_float*kx_float*norm_k) *  real(s_ab_irrot(1,3,kx,ky,kz))+&
            ((-1)*kz_float*ky_float*norm_k) *  real(s_ab_irrot(3,2,kx,ky,kz))+&
            (1 - kz_float*kz_float*norm_k) *   real(s_ab_irrot(3,3,kx,ky,kz))

            s_perp_rot(i,j,k) = (1 - kx_float*kx_float*norm_k) *&
                                &real(s_ab_rot(1,1,kx,ky,kz))+&
            ((-1)*kx_float*ky_float*norm_k) *  real(s_ab_rot(1,2,kx,ky,kz))+&
            ((-1)*kx_float*kz_float*norm_k) *  real(s_ab_rot(1,3,kx,ky,kz))+&
            ((-1)*ky_float*kx_float*norm_k) *  real(s_ab_rot(2,1,kx,ky,kz))+&
            (1 - ky_float*ky_float*norm_k) *   real(s_ab_rot(2,2,kx,ky,kz))+&
            ((-1)*ky_float*kz_float*norm_k) *  real(s_ab_rot(2,3,kx,ky,kz))+&
            ((-1)*kz_float*kx_float*norm_k) *  real(s_ab_rot(1,3,kx,ky,kz))+&
            ((-1)*kz_float*ky_float*norm_k) *  real(s_ab_rot(3,2,kx,ky,kz))+&
            (1 - kz_float*kz_float*norm_k) *   real(s_ab_rot(3,3,kx,ky,kz))

            s_par(i,j,k) = (kx_float*kx_float*norm_k) *&
                            &real(s_ab(1,1,kx,ky,kz))+&
            (kx_float*ky_float*norm_k) * real(s_ab(1,2,kx,ky,kz))+&
            (kx_float*kz_float*norm_k) * real(s_ab(1,3,kx,ky,kz))+&
            (ky_float*kx_float*norm_k) * real(s_ab(2,1,kx,ky,kz))+&
            (ky_float*ky_float*norm_k) * real(s_ab(2,2,kx,ky,kz))+&
            (ky_float*kz_float*norm_k) * real(s_ab(2,3,kx,ky,kz))+&
            (kz_float*kx_float*norm_k) * real(s_ab(3,1,kx,ky,kz))+&
            (kz_float*ky_float*norm_k) * real(s_ab(3,2,kx,ky,kz))+&
            (kz_float*kz_float*norm_k) * real(s_ab(3,3,kx,ky,kz))

            s_par_irrot(i,j,k) = (kx_float*kx_float*norm_k) *&
                               &real(s_ab_irrot(1,1,kx,ky,kz))+&
            (kx_float*ky_float*norm_k) * real(s_ab_irrot(1,2,kx,ky,kz))+&
            (kx_float*kz_float*norm_k) * real(s_ab_irrot(1,3,kx,ky,kz))+&
            (ky_float*kx_float*norm_k) * real(s_ab_irrot(2,1,kx,ky,kz))+&
            (ky_float*ky_float*norm_k) * real(s_ab_irrot(2,2,kx,ky,kz))+&
            (ky_float*kz_float*norm_k) * real(s_ab_irrot(2,3,kx,ky,kz))+&
            (kz_float*kx_float*norm_k) * real(s_ab_irrot(3,1,kx,ky,kz))+&
            (kz_float*ky_float*norm_k) * real(s_ab_irrot(3,2,kx,ky,kz))+&
            (kz_float*kz_float*norm_k) * real(s_ab_irrot(3,3,kx,ky,kz))

            s_par_rot(i,j,k) = (kx_float*kx_float*norm_k) *&
                               &real(s_ab_rot(1,1,kx,ky,kz))+&
            (kx_float*ky_float*norm_k) * real(s_ab_rot(1,2,kx,ky,kz))+&
            (kx_float*kz_float*norm_k) * real(s_ab_rot(1,3,kx,ky,kz))+&
            (ky_float*kx_float*norm_k) * real(s_ab_rot(2,1,kx,ky,kz))+&
            (ky_float*ky_float*norm_k) * real(s_ab_rot(2,2,kx,ky,kz))+&
            (ky_float*kz_float*norm_k) * real(s_ab_rot(2,3,kx,ky,kz))+&
            (kz_float*kx_float*norm_k) * real(s_ab_rot(3,1,kx,ky,kz))+&
            (kz_float*ky_float*norm_k) * real(s_ab_rot(3,2,kx,ky,kz))+&
            (kz_float*kz_float*norm_k) * real(s_ab_rot(3,3,kx,ky,kz))

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

      ! possible normalisation thing, not sure yet
      ! s_perp = s_perp * L**2
      ! s_perp_rot = s_perp_rot * L**2
      ! s_perp_irrot = s_perp_irrot * L**2
      ! s_par = s_par * L**2
      ! s_par_rot = s_par_rot * L**2
      ! s_par_irrot = s_par_irrot * L**2


      open  (30, file=sphe_sus_file, position='append')
      write (30, '(a)') "# S_ab integrals (* L**2)!"
      write (30, '(a)') "# S_xx: total, rot, irrot"
      write (30, '(3f18.10)') sum(real(s_ab(1,1,:,:,:))),&
      sum(real(s_ab_rot(1,1,:,:,:))), sum(real(s_ab_irrot(1,1,:,:,:)))
      write (30, '(a)') "# S_xy: total, rot, irrot"
      write (30, '(3f18.10)') sum(real(s_ab(1,2,:,:,:))),&
      sum(real(s_ab_rot(1,2,:,:,:))), sum(real(s_ab_irrot(1,2,:,:,:)))
      write (30, '(a)') "# S_xz: total, rot, irrot"
      write (30, '(3f18.10)') sum(real(s_ab(1,3,:,:,:))),&
      sum(real(s_ab_rot(1,2,:,:,:))), sum(real(s_ab_irrot(1,2,:,:,:)))
      write (30, '(a)') "# S_yx: total, rot, irrot"
      write (30, '(3f18.10)') sum(real(s_ab(2,1,:,:,:))),&
      sum(real(s_ab_rot(2,1,:,:,:))), sum(real(s_ab_irrot(2,1,:,:,:)))
      write (30, '(a)') "# S_yy: total, rot, irrot"
      write (30, '(3f18.10)') sum(real(s_ab(2,2,:,:,:))),&
      sum(real(s_ab_rot(2,2,:,:,:))), sum(real(s_ab_irrot(2,2,:,:,:)))
      write (30, '(a)') "# S_yz: total, rot, irrot"
      write (30, '(3f18.10)') sum(real(s_ab(2,3,:,:,:))),&
      sum(real(s_ab_rot(2,3,:,:,:))), sum(real(s_ab_irrot(2,3,:,:,:)))
      write (30, '(a)') "# S_zx: total, rot, irrot"
      write (30, '(3f18.10)') sum(real(s_ab(3,1,:,:,:))),&
      sum(real(s_ab_rot(3,1,:,:,:))), sum(real(s_ab_irrot(2,1,:,:,:)))
      write (30, '(a)') "# S_zy: total, rot, irrot"
      write (30, '(3f18.10)') sum(real(s_ab(3,2,:,:,:))),&
      sum(real(s_ab_rot(3,2,:,:,:))), sum(real(s_ab_irrot(3,2,:,:,:)))
      write (30, '(a)') "# S_zz: total, rot, irrot"
      write (30, '(3f18.10)') sum(real(s_ab(3,3,:,:,:))),&
      sum(real(s_ab_rot(3,3,:,:,:))), sum(real(s_ab_irrot(3,3,:,:,:)))
      write (30, '(a)') "# S_perp integrals: total, rot, irrot"
      write (30, '(3f18.10)') sum(s_perp), sum(s_perp_rot), sum(s_perp_irrot)
      close (30)

      !dist_r = dist_r / (no_measurements)
      do i = 1,ceiling( sqrt(float((3*((L/2)**3)))) * (1 / bin_size) )
        write (11, dir_dist_format_string)&
        i * bin_size, bin_count(i), abs(dist_r(i))
      end do

      close(11)

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

              ! write (13, struc_format_string)&
              ! 2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
              ! 2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
              ! field_struc(i,j,k)

              write (14, field_format_string)&
              2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
              2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
              2*pi*(k - 1 - bz*(L/2))/(L*lambda),&
              real(s_ab(1,1,i,j,k)),&
              real(s_ab(1,2,i,j,k)),&
              real(s_ab(1,3,i,j,k)),&
              real(s_ab(2,1,i,j,k)),&
              real(s_ab(2,2,i,j,k)),&
              real(s_ab(2,3,i,j,k)),&
              real(s_ab(3,1,i,j,k)),&
              real(s_ab(3,2,i,j,k)),&
              real(s_ab(3,3,i,j,k))

              ! write (16, struc_format_string)&
              ! 2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
              ! 2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
              ! field_struc_irrot(i,j,k)

              write (17, field_format_string)&
              2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
              2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
              2*pi*(k - 1 - bz*(L/2))/(L*lambda),&
              real(s_ab_irrot(1,1,i,j,k)),&
              real(s_ab_irrot(1,2,i,j,k)),&
              real(s_ab_irrot(1,3,i,j,k)),&
              real(s_ab_irrot(2,1,i,j,k)),&
              real(s_ab_irrot(2,2,i,j,k)),&
              real(s_ab_irrot(2,3,i,j,k)),&
              real(s_ab_irrot(3,1,i,j,k)),&
              real(s_ab_irrot(3,2,i,j,k)),&
              real(s_ab_irrot(3,3,i,j,k))

              ! write (19, struc_format_string)&
              ! 2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
              ! 2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
              ! field_struc_rot(i,j,k)

              write (20,field_format_string)&
              2*pi*(i - 1 - bz*(L/2))/(L*lambda),&
              2*pi*(j - 1 - bz*(L/2))/(L*lambda),&
              2*pi*(k - 1 - bz*(L/2))/(L*lambda),&
              real(s_ab_rot(1,1,i,j,k)),&
              real(s_ab_rot(1,2,i,j,k)),&
              real(s_ab_rot(1,3,i,j,k)),&
              real(s_ab_rot(2,1,i,j,k)),&
              real(s_ab_rot(2,2,i,j,k)),&
              real(s_ab_rot(2,3,i,j,k)),&
              real(s_ab_rot(3,1,i,j,k)),&
              real(s_ab_rot(3,2,i,j,k)),&
              real(s_ab_rot(3,3,i,j,k))

            end if

            write (15, field_format_string)&
            2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(k - 1 - sp*(L/2))/(L*lambda),&
            s_perp(i,j,k)

            write (18, field_format_string)&
            2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(k - 1 - sp*(L/2))/(L*lambda),&
            s_perp_irrot(i,j,k)

            write (21, field_format_string)&
            2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(k - 1 - sp*(L/2))/(L*lambda),&
            s_perp_rot(i,j,k)

            write (22, field_format_string)&
            2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(k - 1 - sp*(L/2))/(L*lambda),&
            s_par(i,j,k)

            write (23, field_format_string)&
            2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(k - 1 - sp*(L/2))/(L*lambda),&
            s_par_irrot(i,j,k)

            write (24, field_format_string)&
            2*pi*(i - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(j - 1 - sp*(L/2))/(L*lambda),&
            2*pi*(k - 1 - sp*(L/2))/(L*lambda),&
            s_par_rot(i,j,k)

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

      deallocate(s_perp); deallocate(s_perp_rot); deallocate(s_perp_irrot);
      deallocate(s_par); deallocate(s_par_rot); deallocate(s_par_irrot);
    end if ! do_corr

  end subroutine calc_correlations

end module output
