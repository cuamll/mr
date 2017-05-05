module setup

  use common
  use io
  use linear_solver
  implicit none
  integer, private :: i, j, k, n

  contains

  subroutine allocations

    allocate(v(L,L,L))
    allocate(pos(L))
    allocate(neg(L))
    allocate(e_x(L,L,L))
    allocate(e_y(L,L,L))
    allocate(e_z(L,L,L))
    allocate(mnphi_x(L,L,L))
    allocate(mnphi_y(L,L,L))
    allocate(mnphi_z(L,L,L))
    allocate(lgf(L,L,L,L,L,L))
    allocate(ch_ch(bz*(L+1),bz*(L+1),bz*(L+1),no_measurements))
    allocate(rho_k_m_t(bz*(L+1),bz*(L+1),bz*(L+1),no_measurements))
    allocate(rho_k_p_t(bz*(L+1),bz*(L+1),bz*(L+1),no_measurements))
    allocate(dir_struc_n(L/2 + 1,L/2 + 1,L/2 + 1,no_measurements))
    allocate(e_kx_t(bz*(L+1),bz*(L+1),bz*(L+1),no_measurements))
    allocate(fe_fe(bz*(L+1),bz*(L+1),bz*(L+1),no_measurements))
    allocate(s_ab_n(3,3,bz*(L+1),bz*(L+1),bz*(L+1),no_measurements))
    allocate(energy(no_measurements + 1))
    allocate(sq_energy(no_measurements + 1))

    v = 0
    pos = 0
    neg = 0
    e_x = 0.0
    e_y = 0.0
    e_z = 0.0
    mnphi_x = 0.0
    mnphi_y = 0.0
    mnphi_z = 0.0
    lgf = 0.0
    ch_ch = 0.0
    rho_k_m_t = 0.0
    rho_k_p_t = 0.0
    dir_struc_n = 0.0
    e_kx_t = 0.0
    fe_fe = 0.0
    s_ab_n = 0.0
    energy = 0.0
    sq_energy= 0.0

  end subroutine allocations

  subroutine deallocations

    deallocate(v)
    deallocate(pos)
    deallocate(neg)
    deallocate(e_x)
    deallocate(e_y)
    deallocate(e_z)
    deallocate(mnphi_x)
    deallocate(mnphi_y)
    deallocate(mnphi_z)
    deallocate(lgf)
    deallocate(ch_ch)
    deallocate(rho_k_m_t)
    deallocate(rho_k_p_t)
    deallocate(dir_struc_n)
    deallocate(e_kx_t)
    deallocate(fe_fe)
    deallocate(s_ab_n)

  end subroutine deallocations

  subroutine latt_init
    integer, dimension(:,:), allocatable :: v_temp

    n = 0

    if (add_charges.ne.0) then

      write (*,*)
      write (*,'(a,I4.1,a)') "Adding ",add_charges," charges. Charge positions:"

      do while (n.lt.add_charges)

        ! pick a random position, check if there's a charge there
        ! if so, pick again; if not, alternate pos/neg
        i = int(rand() * L) + 1
        j = int(rand() * L) + 1
        k = int(rand() * L) + 1

        if (v(i,j,k).ne.0) then
          CYCLE
        end if

        if (modulo(n,2)==0) then
          v(i,j,k) = 1
        else
          v(i,j,k) = -1
        end if

        n = n + 1

        write (*,'(I4.1,a2,I3.1,I3.1,I3.1,I3.1)') n,": ",i,j,k,v(i,j,k)

      end do

    else ! add_charges = 0; read in lattice file

      allocate(v_temp(L**2,L))

      open(unit = 2, file = lattfile)
      read(2,*)((v_temp(i,j),j = 1,L),i = 1,L**2)

      ! 2,3,1 makes x,y,z correspond with what you expect from the file
      ! doesn't actually make any difference so long as you're consistent
      v  =  reshape(v_temp, (/ L,L,L /), ORDER  =  (/ 2,3,1 /))

      deallocate(v_temp)
      close(2)

    end if

  end subroutine latt_init

  subroutine arrays_init

    ! if there are zero charges, we don't need to find the LGF
    ! should be neutral, so sum(v) is zero; hence mask
    if (sum(v,mask=v.gt.0).ne.0) then
      write (*,*) "Calling Poisson solver..."
      call linsol
    end if

    ! set e_x to irrotational - temporary solution
    do i = 1,L
      do j = 1,L
        do k = 1,L

          e_x(i,j,k) = mnphi_x(i,j,k)
          e_y(i,j,k) = mnphi_y(i,j,k)
          e_z(i,j,k) = mnphi_z(i,j,k)

        end do
      end do
    end do

  end subroutine arrays_init

  subroutine setup_wrapper

    ! wrapper for convenience in main
    call read_input
    call randinit(seed)
    call allocations
    call latt_init
    call PBCs
    call arrays_init

  end subroutine setup_wrapper

end module setup
