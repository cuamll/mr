! Module containing everything we need for the linear solver
! to get the irrotational E-field part for a cubic lattice
module linear_solver
  use common
  implicit none

  integer, private :: a,b,c,x,y,z,kx,ky,kz,i,nch
  real(kind=8), private :: sum_x,sum_y,sum_z,p1,p2,q1,q2,r1,r2,fkx,fky,fkz
  real(kind=8), private :: m1p1,m1p2,m1q1,m1q2,m1r1,m1r2,charge
  real(kind=8), public :: u_int, u_self, g_zero, g_z_sum
  character(4) :: lgf_path = "lgf/"
  character(2) :: l_char
  character(6) :: lgf_file
  logical :: lgf_there
  save

  contains

  subroutine linsol

  ! some things need initialising
  u_tot=0.0
  u_self = 0.0
  u_int=0.0
  nch=0

  if (have_lgf.eq.0) then

    ! write (l_char,'(i2)') L
    ! lgf_file = lgf_path//l_char

    ! inquire(file=lgf_file, exist = lgf_there)

    ! if (lgf_there) then
    !   open(30, file=lgf_file, status="old", action="read", access="stream", form="unformatted")
    !   read(30) lgf
    !   close(30)
    ! else
      call lgfcalc(lgf_file)
    ! end if

  end if

    do b=1,L
      do a=1,L
        ! Basically, we can think of the {a,b,c} as the sum over \vec{x}
        ! in G(\vec{x},vec{x'}, and the {x,y,z} as the \vec{x'}.

        ! initialise the sums, idiot
        sum_x=0.0
        sum_y=0.0

        !write (*,*) sum_x,sum_y,ebar(1),ebar(2),mnphi(1,a,b),mnphi(2,a,b)
          do y=1,L
            do x=1,L

              if (a.eq.x.and.b.eq.y) then
                g0 = g0 + lgf(a,b,x,y)
              end if

              if (v(x,y).ne.0) then ! non-zero charge at (x,y,z)
                !write (*,*) lgf(a,b,x,y)

                charge = q * v(x,y)
                nch = nch + 1
                sum_x = sum_x + charge * (-1)&
                      * (lgf(a,b,x,y) - lgf(a,b,neg(x),y))
                sum_y = sum_y + charge * (-1)&
                      * (lgf(a,b,x,y) - lgf(a,b,x,neg(y)))
                if (v(a,b).ne.0) then
                  if (a.eq.x.and.b.eq.y) then
                    u_self = u_self + charge**2 * lgf(x,y,x,y)
                  else
                    u_int = u_int + q * v(a,b) * charge * lgf(a,b,x,y)
                  end if
                end if

              end if ! if v(x,y,z).ne.0

          end do ! y do loop
        end do ! x do loop

        mnphi(1,a,b)=sum_x
        mnphi(2,a,b)=sum_y

        ebar(1)=ebar(1)+mnphi(1,a,b)
        ebar(2)=ebar(2)+mnphi(2,a,b)
        !write (*,*) sum_x,sum_y,ebar(1),ebar(2),mnphi(1,a,b),mnphi(2,a,b)

        u_tot=u_tot+0.5*(mnphi(1,a,b)**2&
              +mnphi(2,a,b)**2)

    end do ! b do loop
  end do ! a do loop

  ebar = ebar / L**2
  g0 = g0 / L**2

  nch=nch/L**2 ! bc we count nch once for each abc

  !write (*,*) "ebar = ",ebar(1),ebar(2),ebar(3)
  !write(*,*)
  !write(*,*) "--- LINEAR SOLVER RESULTS: ---"
  !write (*,*) 'sum of irrotational E_ij^2 =',u_tot
  !write(*,*) "self-energy from lgf(0,0) * n charges = ",u_self
  !write(*,*) 'interaction energy = ',u_int
  !write (*,*) 'harmonic term in units of 1/L**3 (V*Ebar^2) = '&
  !  &,L**3*sum(ebar**2)

  end subroutine linsol

  subroutine lgfcalc(filename)
    character(6), intent(in) :: filename

    g0 = 0.0

    do y=1,L
      do x=1,L
        do b=1,L
          do a=1,L

            lgf(a,b,x,y)=0.0

            ! these need to be real, for cos to work later
            p1=float(a-x)
            q1=float(b-y)

            ! these need to be within [L/2,-L/2]. The Green's
            ! function is even though because we use cosines, so it
            ! doesn't matter whether it's positive or negative.

            if (p1.gt.float(L/2)) then
              p1=p1-float(L)
            else if (p1.lt.(-float(L/2))) then
              p1=p1+float(L)
            end if

            if (q1.gt.float(L/2)) then
              q1=q1-float(L)
            else if (q1.lt.(-float(L/2))) then
              q1=q1+float(L)
            end if

              do ky=-(L-1)/2,L/2
                do kx=-(L-1)/2,L/2
                  fky=2*pi*ky/L
                  fkx=2*pi*kx/L
                  if ((kx.eq.0).and.(ky.eq.0)) then

                  else

                    lgf(a,b,x,y)=lgf(a,b,x,y)+(cos(fkx*p1)&
                      *cos(fky*q1))/(2-cos(fkx)-cos(fky))

                  end if ! end of kx=ky=kz=0 block
                end do ! end kx loop
              end do ! end ky loop

            lgf(a,b,x,y)=lgf(a,b,x,y)/(2*L**2)

            if (a.eq.x.and.b.eq.y) then
              g0 = g0 + lgf(a,b,x,y)
            end if

          end do ! end a loop
        end do ! end b loop
      end do ! end x loop
    end do ! end y loop

    g0 = g0 / L**2

    !write (*,*) "g(0) = ",g0
    !write (*,*) "mu = ",-1 * g0 * ((q**2) / (eps_0))

    ! open(30, file=filename, status="new",&
    !      action="write", access="stream", form="unformatted")
    ! write(30) lgf
    ! close(30)

    have_lgf=1

  end subroutine lgfcalc

end module linear_solver
