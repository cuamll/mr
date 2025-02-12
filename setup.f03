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

    allocate(pos(L))
    allocate(neg(L))
    allocate(v_avg(L,L))
    allocate(windings(2,no_measurements))
    allocate(windings_sq(2,no_measurements))
    allocate(e_tot_avg(2,L,L))
    allocate(lgf(0:L/2,0:L/2))
    if (do_corr) then
      allocate(s_ab(L/2+1,L,3))
      allocate(s_ab_large(2,2,(bz*L)+1,(bz*L)+1))
      allocate(ch_ch((bz*L)+1,(bz*L)+1))
      allocate(rho_k_m((bz*L)+1,(bz*L)+1))
      allocate(rho_k_p((bz*L)+1,(bz*L)+1))
      allocate(dir_struc((L/2) + 1,(L/2) + 1))
      allocate(dist_r(ceiling(sqrt(float(3*(((L/2)**2))))*(1 / bin_size))))
      allocate(bin_count(ceiling(sqrt(float(3*(((L/2)**2))))*(1 / bin_size))))
      allocate(ch_in(L,L))
      allocate(e_in(L,L))
      allocate(chk(L/2+1,L))
      allocate(exk(L/2+1,L))
      allocate(eyk(L/2+1,L))
    end if

    call PBCs

    v_avg = 0.0; rho_avg = 0.0; runtot = 0.0
    avg_field_total = 0.0; avg_field_sq_total = 0.0
    e_tot_avg = 0.0
    ener_tot_sum = 0.0
    ener_tot_sq_sum = 0.0
    ebar_sum = 0.0; ebar_sq_sum = 0.0; ebar_dip_sum = 0.0; lgf = 0.0;
    ebar_dip_sq_sum = 0.0; ebar_wind_sum = 0.0; ebar_wind_sq_sum = 0.0;
    div = 0.0; divsq = 0.0

    if (do_corr) then
      s_ab = (0.0,0.0); s_ab_large = (0.0,0.0);
      ch_ch = (0.0,0.0); rho_k_p = (0.0,0.0); rho_k_m = (0.0,0.0)
      dir_struc = 0.0; dist_r = 0.0; bin_count = 0.0;
    end if

    attempts = 0; accepts = 0; windings = 0.0; windings_sq = 0.0
    ! we know in advance how many rot. and harm. attempts we'll make
    attempts(2) = (therm_sweeps + measurement_sweeps) *&
      no_samples * L**2 * rot_ratio
    attempts(3) = (therm_sweeps + measurement_sweeps) *&
      no_samples * L**2 * g_ratio

  ! FFTW plans -- gonna need them a lot
  plan_x = fftw_plan_dft_r2c_2d(L,L,e_in,exk,FFTW_ESTIMATE)
  plan_y = fftw_plan_dft_r2c_2d(L,L,e_in,eyk,FFTW_ESTIMATE)
  plan_ch = fftw_plan_dft_r2c_2d(L,L,ch_in,chk,FFTW_ESTIMATE)

  end subroutine initial_setup

  subroutine allocations

    allocate(v(L,L))
    allocate(e_field(2,L,L))
    allocate(mnphi(2,L,L))

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

      ! this is currently never true, but could come in handy
      if (read_lattfile) then
        allocate(v_temp(L,L))

        open(unit = 2, file = lattfile)
        read(2,*)((v_temp(i,j),j=1,L),i = 1,L)

        ! 2,3,1 makes x,y,z correspond with what you expect from the file
        ! doesn't actually make any difference so long as you're consistent
        v  =  v_temp

        deallocate(v_temp)
        close(2)

      end if

      n = 0

      if (charge_gen.eq."RANDOM") then

        do while (n.lt.add_charges)

          ! pick a random position, check if there's a charge there
          ! if so, pick again; if not, alternate pos/neg
          i = int(rand() * L) + 1
          j = int(rand() * L) + 1

          if (v(i,j).ne.0) then
            CYCLE
          end if

          if (modulo(n,2)==0) then
            v(i,j) = 1
          else
            v(i,j) = -1
          end if

          n = n + 1

        end do

      else if (charge_gen.eq."DIPOLE") then

        do while (n.lt.add_charges)

          i = int(rand() * L) + 1
          j = int(rand() * L) + 1

          if (v(i,j).ne.0) then
            CYCLE
          end if

          ! choose between four orientations of a dipole
          if (rand().lt.0.5) then
            ! x-direction
            if (v(neg(i),j).ne.0) then
              CYCLE
            end if

            if (rand().lt.0.5) then
              ! + -
              v(neg(i),j) = +1
              v(i,j) = -1
            else
              ! - +
              v(neg(i),j) = -1
              v(i,j) = +1
            end if

          else
            ! y-direction
            if (v(i,neg(j)).ne.0) then
              CYCLE
            end if

            if (rand().lt.0.5) then
              ! + -
              v(i,neg(j)) = +1
              v(i,j) = -1
            else
              ! - +
              v(i,neg(j)) = -1
              v(i,j) = +1
            end if

          end if

          n = n + 2

        end do

      ! should change the variable really so the L doesn't get cut off
      else if (charge_gen.eq."CRYSTA") then

        add_charges = L**2
        do i = 1,L**2

          k =((i - 1)/L) + 1
          j = modulo(i - 1, L) + 1
          ! write(6,*) i,k,j, (2*mod(i,2) - 1)
          ! flush(6)
          v(k,j) = 2*mod(k + j, 2) - 1

        end do

      end if

    end subroutine latt_init

  subroutine arrays_init

    ! if there are zero charges, we don't need to find the LGF
    ! should be neutral, so sum(v) is zero; hence mask
    if (sum(v,mask=v.gt.0).ne.0) then
      !write (*,*) "Calling Poisson solver..."
      call linsol
    end if

    ! implicitly set E_rot = 0 at t = 0
    e_field = mnphi

  end subroutine arrays_init

  subroutine setup_wrapper(n)
    integer, intent(in) :: n

    ! wrapper for convenience in main
    call randinit(n)
    call allocations
    call latt_init
    call arrays_init

  end subroutine setup_wrapper

end module setup
