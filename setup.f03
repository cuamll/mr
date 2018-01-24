module setup

  use common
  use input
  use linear_solver
  implicit none
  integer, private :: i, j, k, n

  contains

  subroutine initial_setup
    ! various things which can't be zeroed
    ! at the start of every sample
    complex(kind=rk) :: prefac, imag

    allocate(pos(L))
    allocate(neg(L))
    allocate(v_avg(L,L,L))
    allocate(e_tot_avg(3,L,L,L))
    allocate(e_rot_avg(3,L,L,L))
    allocate(e_irrot_avg(3,L,L,L))
    allocate(s_ab(3,3,(bz*L)+1,(bz*L)+1,(bz*L)+1))
    allocate(s_ab_rot(3,3,(bz*L)+1,(bz*L)+1,(bz*L)+1))
    allocate(s_ab_irrot(3,3,(bz*L)+1,(bz*L)+1,(bz*L)+1))
    allocate(ch_ch((bz*L)+1,(bz*L)+1,(bz*L)+1))
    allocate(rho_k_m((bz*L)+1,(bz*L)+1,(bz*L)+1))
    allocate(rho_k_p((bz*L)+1,(bz*L)+1,(bz*L)+1))
    allocate(dir_struc((L/2) + 1,(L/2) + 1,(L/2) + 1))
    allocate(dist_r(ceiling(sqrt(float(3*(((L/2)**3))))*(1 / bin_size))))
    allocate(bin_count(ceiling(sqrt(float(3*(((L/2)**3))))*(1 / bin_size))))
    allocate(lgf(L,L,L,L,L,L))
    allocate(fw(L,L))
    allocate(hw(L,L))

    call PBCs

    v_avg = 0.0; rho_avg = 0.0; runtot = 0.0
    avg_field_total = 0.0; avg_field_rot = 0.0; avg_field_irrot = 0.0
    e_tot_avg = 0.0; e_rot_avg = 0.0; e_irrot_avg = 0.0
    ener_tot_sum = 0.0; ener_rot_sum = 0.0; ener_irrot_sum = 0.0;
    ener_tot_sq_sum = 0.0; ener_rot_sq_sum = 0.0; ener_irrot_sq_sum = 0.0;
    ebar_sum = 0.0; ebar_sq_sum = 0.0; ebar_dip_sum = 0.0;
    ebar_dip_sq_sum = 0.0; ebar_wind_sum = 0.0; ebar_wind_sq_sum = 0.0;
    s_ab = (0.0,0.0); s_ab_rot = (0.0,0.0); s_ab_irrot = (0.0,0.0);
    ch_ch = (0.0,0.0); rho_k_p = (0.0,0.0); rho_k_m = (0.0,0.0); lgf = 0.0
    dir_struc = 0.0; dist_r = 0.0; bin_count = 0.0;
    attempts = 0; accepts = 0
    ! we know in advance how many rot. and harm. attempts we'll make
    attempts(2) = (therm_sweeps + measurement_sweeps) *&
      no_samples * 1 * L**3 * rot_ratio
    attempts(3) = (therm_sweeps + measurement_sweeps) *&
      no_samples * 3 * L**3 * g_ratio

    imag = (0.0, 1.0)
    prefac = (-1)*imag*((2*pi)/(L*lambda))

    do i = 1,L
      do k = 1,L
        fw(i,k) = prefac * i * (k - 1)
        hw(i,k) = prefac * (neg(i) + (1.0/2)) * (k - 1)
      end do
    end do

  end subroutine initial_setup

  subroutine allocations

    allocate(v(L,L,L))
    allocate(e_field(3,L,L,L))
    allocate(mnphi(3,L,L,L))

    v = 0;
    ebar = 0.0;
    e_field = 0.0; mnphi = 0.0;

  end subroutine allocations

  subroutine deallocations

    deallocate(v)
    deallocate(e_field)
    deallocate(mnphi)

  end subroutine deallocations

  subroutine latt_init
    integer, dimension(:,:), allocatable :: v_temp
    logical :: read_lattfile = .false.

    n = 0

    if (add_charges.ne.0) then

      !write (*,*)
      !write (*,'(a,I4.1,a)') "Adding ",add_charges," charges. Charge positions:"

      do while (n.lt.add_charges)

        ! pick a random position, check if there's a charge there
        ! if so, pick again; if not, alternate pos/neg
        i = int(rand() * L) + 1
        j = int(rand() * L) + 1

        if (v(i,j,k).ne.0) then
          CYCLE
        end if

        if (modulo(n,2)==0) then
          v(i,j,k) = 1
        else
          v(i,j,k) = -1
        end if

        n = n + 1

        !write (*,'(I4.1,a2,I3.1,I3.1,I3.1,I3.1)') n,": ",i,j,v(i,j)

      end do

    else ! add_charges = 0; read in lattice file

      if (read_lattfile) then
        allocate(v_temp(L,L**2))

        open(unit = 2, file = lattfile)
        read(2,*)((v_temp(i,j),j = 1,L),i = 1,L**2)

        ! 2,3,1 makes x,y,z correspond with what you expect from the file
        ! doesn't actually make any difference so long as you're consistent
        v  =  reshape(v_temp, (/ L,L,L /), ORDER  =  (/ 2,3,1 /))

        deallocate(v_temp)
        close(2)
      else ! all zeroes
        v = 0
      end if

    end if

  end subroutine latt_init

  subroutine arrays_init

    ! if there are zero charges, we don't need to find the LGF
    ! should be neutral, so sum(v) is zero; hence mask
    if (sum(v,mask=v.gt.0).ne.0) then
      !write (*,*) "Calling Poisson solver..."
      call linsol
    end if

    ! set e_field to irrotational - temporary solution
    e_field = mnphi

  end subroutine arrays_init

  subroutine setup_wrapper(n)
    integer, intent(in) :: n

    ! wrapper for convenience in main
    ! call read_input
    call randinit(n)
    call allocations
    call latt_init
    call arrays_init

  end subroutine setup_wrapper

end module setup
