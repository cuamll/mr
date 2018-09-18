! Module containing everything we need for the linear solver
! to get the irrotational E-field part for a cubic lattice
module linear_solver
  use common
  implicit none

  integer, private :: a,b,c,x,y,z,kx,ky,kz,rx,ry,rz,rpx,rpy,rpz,i,nch
  real(kind=8), private :: sum_x,sum_y,sum_z,p1,q1,r1,fkx,fky,fkz
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

    do c=1,L
      do b=1,L
        do a=1,L
          ! Basically, we can think of the {a,b,c} as the sum over \vec{x}
          ! in G(\vec{x},vec{x'}, and the {x,y,z} as the \vec{x'}.

          ! initialise the sums, idiot
          sum_x=0.0
          sum_y=0.0
          sum_z=0.0

          !write (*,*) sum_x,sum_y,ebar(1),ebar(2),mnphi(1,a,b),mnphi(2,a,b)
          do z=1,L
            do y=1,L
              do x=1,L

              rx = (a - x)
              if (rx.gt.(L/2)) then
                rx = rx - L
              else if (rx.le.(-L/2)) then
                rx = rx + L
              end if
              rx = abs(rx)

              ! had a thought
              ! rpx = pos(rx + 1) - 1
              ! rpx = mod(rx + 1, L/2 + 1)

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

              ! rpy = mod(ry + 1, L/2) - 1
              ! rpy = mod(pos(ry), L/2 +  1)

              rpy = (b - pos(y))
              if (rpy.gt.(L/2)) then
                rpy = rpy - L
              else if (rpy.le.(-L/2)) then
                rpy = rpy + L
              end if
              rpy = abs(rpy)

              rz = (c - z)
              if (rz.gt.(L/2)) then
                rz = rz - L
              else if (rz.le.(-L/2)) then
                rz = rz + L
              end if
              rz = abs(rz)

              rpz = (c - pos(z))
              if (rpz.gt.(L/2)) then
                rpz = rpz - L
              else if (rpz.le.(-L/2)) then
                rpz = rpz + L
              end if
              rpz = abs(rpz)

              if (a.eq.x.and.b.eq.y.and.c.eq.z) then
                g0 = g0 + lgf(rx,ry,rz)
              end if

              if (v(x,y,z).ne.0) then ! non-zero charge at (x,y,z)

                charge = q * v(x,y,z)
                ! nch = nch + 1

                sum_x = sum_x + charge * (+1)&
                      * (lgf(rpx,ry,rz) - lgf(rx,ry,rz))
                sum_y = sum_y + charge * (+1)&
                      * (lgf(rx,rpy,rz) - lgf(rx,ry,rz))
                sum_z = sum_z + charge * (+1)&
                      * (lgf(rx,ry,rpz) - lgf(rx,ry,rz))

                if (v(a,b,c).ne.0) then
                  if (a.eq.x.and.b.eq.y.and.c.eq.z) then
                    ! u_self = u_self + (1 / (2 * eps_0)) * charge**2 * lgf(rx,ry,rz)
                    u_self = u_self + charge**2 * lgf(rx,ry,rz)
                  else
                    ! u_int = u_int + (1 / (2 * eps_0)) * q * v(a,b,c) * charge * lgf(rx,ry,rz)
                    u_int = u_int + q * v(a,b,c) * charge * lgf(rx,ry,rz)
                  end if
                end if

              end if ! if v(x,y,z).ne.0

            end do ! y do loop
          end do ! x do loop
        end do

        ! already multiplied by q
        mnphi(1,a,b,c)= (1 / eps_0) * sum_x
        mnphi(2,a,b,c)= (1 / eps_0) * sum_y
        mnphi(3,a,b,c)= (1 / eps_0) * sum_z

      end do
    end do ! b do loop
  end do ! a do loop

  end if

  mnphi = -1.0 * mnphi
  ebar(1) = sum(mnphi(1,:,:,:))
  ebar(2) = sum(mnphi(2,:,:,:))
  ebar(3) = sum(mnphi(3,:,:,:))
  u_tot = 0.5 * eps_0 * lambda**3 * sum(mnphi*mnphi)
  u_self = u_self * lambda**4 * (1 / (1 * eps_0))
  u_int = u_int * lambda**4 * (1 / (1 * eps_0))

  ebar = ebar / L**3
  g0 = g0 / L**3

  nch= sum(abs(v))

  write(*,*)
  write(*,*) "--- LINEAR SOLVER RESULTS: ---"
  write(*,*) "n_ch",nch
  write(*,*) "ebar = ",ebar(1),ebar(2),ebar(3)
  write(*,*) "nn = ",lgf(1,0,0),lgf(0,0,1),lgf(0,1,0)
  write(*,*) 'sum of irrotational E_ij^2 =',u_tot
  write(*,*) "self-energy from lgf(0,0) per charge = ",u_self / nch
  write(*,*) "self-energy from lgf(0,0) = ",u_self
  write(*,*) 'interaction energy = ',u_int
  write(*,*) 'interaction energy * L**3 / 2 = ',u_int * (L**3 / 2)
  write(*,*) '(unnormalised) u_tot = ',sum(mnphi*mnphi)
  write(*,*) 'u_tot - u_self = ',u_tot - u_self
  write(*,*) 'harmonic term in units of 1/L**3 (V*Ebar^2) = '&
    &,L**3*sum(ebar**2)

  end subroutine linsol

  subroutine lgfcalc
    implicit none

    g0 = 0.0
    lgf = 0.0

    write (l_char,'(i3)') L
    lgf_file = lgf_path//l_char

    inquire(file=lgf_file, exist = lgf_there)

    if (lgf_there) then

      write (*,*) "Reading in LGF from ",lgf_file
      open(30, file=lgf_file, status="old", action="read", access="stream", form="unformatted")
      read(30) lgf
      close(30)

    else

    write(*,*) "Calculating LGF"

    do ry=0,L/2
      do rx=0,L/2
        do rz=0,L/2

          ! these need to be real, for cos to work later
          p1=float(rx)
          q1=float(ry)
          r1=float(rz)

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

          if (r1.gt.float(L/2)) then
            r1=r1-float(L)
          else if (r1.lt.(-float(L/2))) then
            r1=r1+float(L)
          end if

          do kz=-(L-1)/2, L/2
            do ky=-(L-1)/2, L/2
              do kx=-(L-1)/2, L/2
            ! do ky=0,L
            !   do kx=0,L
                fkz=2*pi*kz/L
                fky=2*pi*ky/L
                fkx=2*pi*kx/L

                if ((kx.eq.0).and.(ky.eq.0).and.(kz.eq.0)) then
                ! if ((kx.eq.(L/2)+1).and.(ky.eq.(L/2)+1)) then

                else

                  lgf(rx,ry,rz)= lgf(rx,ry,rz) + (cos(fkx*p1)&
                    *cos(fky*q1)*cos(fkz*r1))/(3-cos(fkx)-cos(fky)-cos(fkz))
                  ! lgf(rx,ry)=lgf(rx,ry)+(cos(fk(kx)*p1)&
                  !   *cos(fk(ky)*q1))/(2-cosine(kx)-cosine(ky))

                end if ! end of kx=ky=kz=0 block
              end do ! end kx loop
            end do ! end ky loop
          end do

          lgf(rx,ry,rz)=lgf(rx,ry,rz)/(2*L**3)
          ! lgf(rx,ry,rz)=lgf(rx,ry,rz)/(L**3)

          if (rx.eq.0.and.ry.eq.0.and.rz.eq.0) then
            g0 = g0 + lgf(rx,ry,rz)
          end if

        end do ! end x loop
      end do ! end y loop
    end do ! end z loop

    open(30, file=lgf_file, status="new",&
         action="write", access="stream", form="unformatted")
    write(30) lgf

    write (*,*) "Written LGF to ",lgf_file
    close(30)

    end if

    have_lgf=1

  end subroutine lgfcalc

end module linear_solver
