module setup

  use common
  use io
  use linear_solver
  implicit none
  integer, private :: i, j, k, n, row, col, tot_q
  real, public :: u_rot
  integer, dimension(:,:), allocatable, private :: v_temp

  contains

  subroutine allocations

    allocate(v(L,L,L))
    allocate(pos(L))
    allocate(neg(L))
    allocate(mnphi_x(L,L,L))
    allocate(mnphi_y(L,L,L))
    allocate(mnphi_z(L,L,L))
    allocate(e_rot_x(L,L,L))
    allocate(e_rot_y(L,L,L))
    allocate(e_rot_z(L,L,L))
    allocate(e_x(L,L,L))
    allocate(e_y(L,L,L))
    allocate(e_z(L,L,L))
    allocate(e_x_lapack(L,L,L))
    allocate(e_y_lapack(L,L,L))
    allocate(e_z_lapack(L,L,L))
    allocate(phi_lapack(L,L,L))
    allocate(lgf(L,L,L,L,L,L))
    allocate(v_temp(L**2,L))
    allocate(e_kx(bz*(L+1),bz*(L+1),bz*(L+1)))
    allocate(e_ky(bz*(L+1),bz*(L+1),bz*(L+1)))
    allocate(e_kz(bz*(L+1),bz*(L+1),bz*(L+1)))
    allocate(fe_fe(bz*(L+1),bz*(L+1),bz*(L+1),iterations))
    allocate(s_ab(3,3,bz*(L+1),bz*(L+1),bz*(L+1)))
    allocate(s_perp(bz*(L+1),bz*(L+1),bz*(L+1)))
    allocate(s_ab_n(3,3,bz*(L+1),bz*(L+1),bz*(L+1),iterations))
    allocate(rho_k_m(bz*(L+1),bz*(L+1),bz*(L+1)))
    allocate(rho_k_p(bz*(L+1),bz*(L+1),bz*(L+1)))
    allocate(ch_ch(bz*(L+1),bz*(L+1),bz*(L+1),iterations))
    allocate(ch_ch_pp(bz*(L+1),bz*(L+1),bz*(L+1),iterations))
    allocate(field_struc(bz*(L+1),bz*(L+1),bz*(L+1)))
    allocate(charge_struc(bz*(L+1),bz*(L+1),bz*(L+1)))

    v = 0
    pos = 0
    neg = 0
    mnphi_x = 0.0
    mnphi_y = 0.0
    mnphi_z = 0.0
    e_rot_x = 0.0
    e_rot_y = 0.0
    e_rot_z = 0.0
    e_x = 0.0
    e_y = 0.0
    e_z = 0.0
    e_kx = 0.0
    e_ky = 0.0
    e_kz = 0.0
    rho_k_m = 0.0
    rho_k_p = 0.0
    ch_ch = 0.0
    ch_ch_pp = 0.0
    fe_fe = 0.0
    field_struc = 0.0
    charge_struc = 0.0
    s_perp = 0.0
    e_x_lapack = 0.0
    e_y_lapack = 0.0
    e_z_lapack = 0.0
    phi_lapack = 0.0
    lgf = 0.0
    v_temp = 0

  end subroutine allocations

  subroutine deallocations

    deallocate(v)
    deallocate(pos)
    deallocate(neg)
    deallocate(mnphi_x)
    deallocate(mnphi_y)
    deallocate(mnphi_z)
    deallocate(e_rot_x)
    deallocate(e_rot_y)
    deallocate(e_rot_z)
    deallocate(e_x)
    deallocate(e_y)
    deallocate(e_z)
    deallocate(e_x_lapack)
    deallocate(e_y_lapack)
    deallocate(e_z_lapack)
    deallocate(phi_lapack)
    deallocate(lgf)
    deallocate(e_kx)
    deallocate(e_ky)
    deallocate(e_kz)
    deallocate(rho_k_m)
    deallocate(rho_k_p)
    deallocate(ch_ch)
    deallocate(ch_ch_pp)
    deallocate(fe_fe)
    deallocate(charge_struc)
    deallocate(field_struc)
    deallocate(s_ab_n)
    deallocate(s_ab)
    deallocate(s_perp)

  end subroutine deallocations

  subroutine latt_init

    tot_q = 0

    if (add_charges.ne.0) then

      ! zero out the lattice first
      do i = 1,L
        do j = 1,L
          do k = 1,L

            v(i,j,k) = 0

          end do
        end do
      end do

      ! then we place add_charges in our lattice, randomly
      do while (tot_q.lt.add_charges)
        i = int(rand() * L) + 1
        j = int(rand() * L) + 1
        k = int(rand() * L) + 1
        if (v(i,j,k).ne.0) then
          ! need them in different places
          CYCLE
        end if
        ! alternate pos and neg charges
        if (modulo(tot_q,2)==0) then
          v(i,j,k) = 1
        else
          v(i,j,k) = -1
        end if
        ! increment, idiot, otherwise we do this forever
        tot_q = tot_q + 1
        write (*,*) "charge ",tot_q,": ",i,j,k,v(i,j,k)
      end do

    else ! add_charges = 0; read in lattice file

      open(unit = 2, file = lattfile)
      read(2,*)((v_temp(row,col),col = 1,L),row = 1,L**2)

      ! 2,3,1 makes x,y,z correspond with what you expect from the file
      ! doesn't actually make any difference so long as you're consistent
      v  =  reshape(v_temp, (/ L,L,L /), ORDER  =  (/ 2,3,1 /))

      deallocate(v_temp) ! we don't need it anymore
      close(2)

    end if

  end subroutine latt_init

  subroutine arrays_init

    ! if there are zero charges, we don't need to find the LGF
    do i = 1,L
      do j = 1,L
        do k = 1,L

          tot_q = tot_q + abs(v(i,j,k))

        end do
      end do
    end do

    if (tot_q.ne.0) then
      call linsol
    else
      do i = 1,L
        do j = 1,L
          do k = 1,L

            mnphi_x(i,j,k) = 0.0
            mnphi_y(i,j,k) = 0.0
            mnphi_z(i,j,k) = 0.0
            e_x(i,j,k) = 0.0
            e_y(i,j,k) = 0.0
            e_z(i,j,k) = 0.0

          end do
        end do
      end do
    end if

    ! set e_x to irrotational - temporary solution
    do i = 1,L
      do j = 1,L
        do k = 1,L

          e_x(i,j,k) = mnphi_x(i,j,k)
          e_y(i,j,k) = mnphi_y(i,j,k)
          e_z(i,j,k) = mnphi_z(i,j,k)

          e_rot_x(i,j,k) = 0.0
          e_rot_y(i,j,k) = 0.0
          e_rot_z(i,j,k) = 0.0

          charge_struc(i,j,k) = 0.0
          field_struc(i,j,k) = 0.0
          s_perp(i,j,k) = 0.0

          do n = 1,iterations
            ch_ch(i,j,k,n) = 0.0
            fe_fe(i,j,k,n) = 0.0
            s_ab_n(1,1,i,j,k,n) = 0.0
            s_ab_n(1,2,i,j,k,n) = 0.0
            s_ab_n(1,3,i,j,k,n) = 0.0
            s_ab_n(2,1,i,j,k,n) = 0.0
            s_ab_n(2,2,i,j,k,n) = 0.0
            s_ab_n(2,3,i,j,k,n) = 0.0
            s_ab_n(3,1,i,j,k,n) = 0.0
            s_ab_n(3,2,i,j,k,n) = 0.0
            s_ab_n(3,3,i,j,k,n) = 0.0

            s_ab(1,1,i,j,k) = 0.0
            s_ab(1,2,i,j,k) = 0.0
            s_ab(1,3,i,j,k) = 0.0
            s_ab(2,1,i,j,k) = 0.0
            s_ab(2,2,i,j,k) = 0.0
            s_ab(2,3,i,j,k) = 0.0
            s_ab(3,1,i,j,k) = 0.0
            s_ab(3,2,i,j,k) = 0.0
            s_ab(3,3,i,j,k) = 0.0
          end do

        end do
      end do
    end do

    u_rot = 0.0
    ebar_x  =  0.0
    ebar_y  =  0.0
    ebar_z  =  0.0

  end subroutine arrays_init

  subroutine do_setup

    ! wrapper for convenience in main
    call read_input
    call randinit(seed)
    call allocations
    call latt_init
    call PBCs
    call arrays_init

  end subroutine do_setup

end module setup
