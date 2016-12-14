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

    allocate(v(2*L,2*L,2*L))
    allocate(pos(2*L))
    allocate(neg(2*L))
    allocate(ch_sites(L**3,3))
    allocate(ex_sites(L**3,3))
    allocate(ey_sites(L**3,3))
    allocate(ez_sites(L**3,3))
    allocate(mnphi_x(L,L,L))
    allocate(mnphi_y(L,L,L))
    allocate(mnphi_z(L,L,L))
    allocate(e_rot_x(L,L,L))
    allocate(e_rot_y(L,L,L))
    allocate(e_rot_z(L,L,L))
    allocate(e_x(2*L,2*L,2*L))
    allocate(e_y(2*L,2*L,2*L))
    allocate(e_z(2*L,2*L,2*L))
    allocate(e_x_lapack(L,L,L))
    allocate(e_y_lapack(L,L,L))
    allocate(e_z_lapack(L,L,L))
    allocate(phi_lapack(L,L,L))
    allocate(lgf(L,L,L,L,L,L))
    allocate(v_temp(L**2,L))
    allocate(e_kx(L+1,L+1,L+1))
    allocate(e_ky(L+1,L+1,L+1))
    allocate(e_kz(L+1,L+1,L+1))
    allocate(rho_k(L+1,L+1,L+1))
    allocate(ch_ch(L+1,2*L,2*L,iterations))

    !allocate(fe_fe(3,3,L+1,L+1,L+1,iterations))
    !allocate(struc_field(3,3,L+1,L+1,L+1))
    allocate(fe_fe(L+1,2*L,2*L,iterations))
    allocate(struc_field(L+1,2*L,2*L))

    allocate(struc_charge(L+1,2*L,2*L))
    v = 0
    pos = 0
    neg = 0
    ch_sites = 0
    ex_sites = 0
    ey_sites = 0
    ez_sites = 0
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
    rho_k = 0.0
    ch_ch = 0.0
    fe_fe = 0.0
    struc_field = 0.0
    struc_charge = 0.0
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
    deallocate(ch_sites)
    deallocate(ex_sites)
    deallocate(ey_sites)
    deallocate(ez_sites)
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
    deallocate(rho_k)
    deallocate(ch_ch)
    deallocate(fe_fe)
    deallocate(struc_charge)
    deallocate(struc_field)

  end subroutine deallocations

  subroutine latt_init
    implicit none
    integer, dimension(:,:,:), allocatable :: v_temp_3d

    tot_q = 0
    allocate(v_temp_3d(L,L,L))

    ! zero out the lattice first
    do i = 1,2*L
      do j = 1,2*L
        do k = 1,2*L

          v(i,j,k) = 0
          e_x(i,j,k) = 0.0
          e_y(i,j,k) = 0.0
          e_z(i,j,k) = 0.0

        end do
      end do
    end do


    if (add_charges.ne.0) then

      ! then we place add_charges in our lattice, randomly
      do while (tot_q.lt.add_charges)
        i = int(rand() * L) + 1
        j = int(rand() * L) + 1
        k = int(rand() * L) + 1
        if (v(2*i - 1,2*j - 1,2*k - 1).ne.0) then
          ! need them in different places
          CYCLE
        end if
        ! alternate pos and neg charges
        if (modulo(tot_q,2)==0) then
          v(2*i - 1,2*j - 1,2*k - 1) = 1
        else
          v(2*i - 1,2*j - 1,2*k - 1) = -1
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
      v_temp_3d  =  reshape(v_temp, (/ L,L,L /), ORDER  =  (/ 2,3,1 /))
      do i = 1,L
        do j = 1,L
          do k = 1,L
            v(2*i - 1,2*j - 1,2*k - 1) = v_temp_3d(i,j,k)
          end do
        end do
      end do

      deallocate(v_temp) ! we don't need it anymore
      deallocate(v_temp_3d)
      close(2)

    end if

  end subroutine latt_init

  subroutine arrays_init

    ! if there are zero charges, we don't need to find the LGF
    do i = 1,L
      do j = 1,L
        do k = 1,L

          tot_q = tot_q + abs(v(2*i - 1,2*j - 1,2*k - 1))

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

          struc_charge(i,j,k) = 0.0
          struc_field(i,j,k) = 0.0
          !struc_field(1,1,i,j,k) = 0.0
          !struc_field(1,2,i,j,k) = 0.0
          !struc_field(1,3,i,j,k) = 0.0
          !struc_field(2,1,i,j,k) = 0.0
          !struc_field(2,2,i,j,k) = 0.0
          !struc_field(2,3,i,j,k) = 0.0
          !struc_field(3,1,i,j,k) = 0.0
          !struc_field(3,2,i,j,k) = 0.0
          !struc_field(3,3,i,j,k) = 0.0

          do n = 1,iterations
            ch_ch(i,j,k,n) = 0.0
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
