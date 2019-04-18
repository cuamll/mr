! Module containing everything we need for the linear solver
! to get the irrotational E-field part for a cubic lattice
module linear_solver
  use common
  implicit none

  integer, private :: a,b,c,x,y,z,kx,ky,kz,rx,ry,rpx,rpy,i,nch
  real(kind=8), private :: sum_x,sum_y,sum_z,p1,p2,q1,q2,r1,r2,fkx,fky,fkz
  real(kind=8), private :: m1p1,m1p2,m1q1,m1q2,m1r1,m1r2,charge
  real(kind=8), public :: u_int, u_self, g_zero, g_z_sum
  character(3) :: l_char
  character(8) :: lgf_fn
  character(200) :: lgf_file
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
    call lgfcalc
  else

    do b=1,L
      do a=1,L
        ! Basically, we can think of the {a,b,c} as the sum over \vec{x}
        ! in G(\vec{x},vec{x'}, and the {x,y,z} as the \vec{x'}.

        ! initialise the sums, idiot
        sum_x=0.0
        sum_y=0.0

        do y=1,L
          do x=1,L

            rx = (a - x)
            if (rx.gt.(L/2)) then
              rx = rx - L
            else if (rx.le.(-L/2)) then
              rx = rx + L
            end if
            rx = abs(rx)

            ! get (x + 1), respecting PBCs
            rpx = (a - pos(x))
            if (rpx.gt.(L/2)) then
              rpx = rpx - L
            else if (rpx.le.(-L/2)) then
              rpx = rpx + L
            end if
            rpx = abs(rpx)

            ry = (b - y)
            if (ry.gt.(L/2)) then
              ry = ry - L
            else if (ry.le.(-L/2)) then
              ry = ry + L
            end if
            ry = abs(ry)

            rpy = (b - pos(y))
            if (rpy.gt.(L/2)) then
              rpy = rpy - L
            else if (rpy.le.(-L/2)) then
              rpy = rpy + L
            end if
            rpy = abs(rpy)

            if (a.eq.x.and.b.eq.y) then
              g0 = g0 + lgf(rx,ry)
            end if

            if (v(x,y).ne.0) then ! non-zero charge at (x,y,z)

              charge = q * v(x,y)
              nch = nch + 1

              ! POSITIVE GRAD i.e. - \tilde{\nabla} \phi
              sum_x = sum_x + charge * (+1)&
                    * (lgf(rpx,ry) - lgf(rx,ry))
              sum_y = sum_y + charge * (+1)&
                    * (lgf(rx,rpy) - lgf(rx,ry))

              ! NEGATIVE GRAD i.e. - \hat{\nabla} \phi
              ! old four-index form - left as an example
              ! sum_x = sum_x + charge * (-1) *&
              !       (lgf(a,b,x,y) - lgf(a,b,neg(x),y))
              ! sum_y = sum_y + charge * (-1) *&
              !       (lgf(a,b,x,y) - lgf(a,b,x,neg(y)))
              if (v(a,b).ne.0) then
                if (a.eq.x.and.b.eq.y) then
                  u_self = u_self + charge**2 * lgf(rx,ry)
                else
                  u_int = u_int + q * v(a,b) * charge * lgf(rx,ry)
                end if
              end if

            end if ! if v(x,y,z).ne.0

          end do
        end do ! x,y loops

        ! already multiplied by q
        mnphi(1,a,b)=(1 / eps_0) * sum_x
        mnphi(2,a,b)=(1 / eps_0) * sum_y

      end do
    end do ! a,b loops

  end if

  mnphi = -1.0 * mnphi
  ebar(1) = sum(mnphi(1,:,:))
  ebar(2) = sum(mnphi(2,:,:))
  u_tot = 0.5 * eps_0 * lambda**2 * sum(mnphi*mnphi)

  ebar = ebar / L**2
  g0 = g0 / L**2

  nch=nch/L**2 ! bc we count nch once for each abc

  end subroutine linsol

  subroutine lgfcalc
    implicit none
    real(kind=rk), dimension(0:L) :: fk, cosine

    g0 = 0.0
    lgf = 0.0

    write (l_char,'(i3)') L
    lgf_file = lgf_path//l_char

    inquire(file=lgf_file, exist = lgf_there)

    if (lgf_there) then

      write (*,*) "Reading in LGF from ",lgf_file
      open(30, file=lgf_file, status="old",
           action="read", access="stream", form="unformatted")
      read(30) lgf
      close(30)

    else

      write(*,*) "Calculating LGF"

      do ry=0,L/2
        do rx=0,L/2

          ! these need to be real, for cos to work later
          p1=float(rx)
          q1=float(ry)

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

          do ky=-(L-1)/2, L/2
            do kx=-(L-1)/2, L/2
              fky=2*pi*ky/L
              fkx=2*pi*kx/L

              if ((kx.eq.0).and.(ky.eq.0)) then

              else

                lgf(rx,ry)=lgf(rx,ry)+(cos(fkx*p1)&
                  *cos(fky*q1))/(2-cos(fkx)-cos(fky))

              end if ! end of kx=ky=kz=0 block
            end do ! end kx loop
          end do ! end ky loop

          lgf(rx,ry)=lgf(rx,ry)/(2*L**2)

          if (rx.eq.0.and.ry.eq.0) then
            g0 = g0 + lgf(rx,ry)
          end if

        end do ! end x loop
      end do ! end y loop

      open(30, file=lgf_file, status="new",&
           action="write", access="stream", form="unformatted")
      write(30) lgf
      write (*,*) "Written LGF to ",lgf_file
      close(30)

    end if

    have_lgf=1

  end subroutine lgfcalc

end module linear_solver
