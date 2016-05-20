! Module containing everything we need for the linear solver
! to get the irrotational E-field part for a cubic lattice
module linear_solver
  use common
  implicit none

  integer, private :: a,b,c,x,y,z,kx,ky,kz,i,nch
  real*8, private :: sum_x,sum_y,sum_z,p1,p2,q1,q2,r1,r2,fkx,fky,fkz
  real*8, private :: m1p1,m1p2,m1q1,m1q2,m1r1,m1r2,charge
  real*8, public :: u_tot, u_int, u_self, g_zero, g_z_sum
  real*8, dimension(:), allocatable, private :: cosine
  save

  contains

  subroutine linsol

  ! some things need initialising
  u_tot=0.0
  u_self = 0.0
  u_int=0.0
  nch=0

  allocate(cosine(L/2-1))

  if (have_lgf.eq.0) then
    call lgfcalc
  end if

  do a=1,L
    do b=1,L
      do c=1,L
        ! Basically, we can think of the {a,b,c} as the sum over \vec{x}
        ! in G(\vec{x},vec{x'}, and the {x,y,z} as the \vec{x'}.

        ! initialise the sums, idiot
        sum_x=0.0
        sum_y=0.0
        sum_z=0.0

        do x=1,L
          do y=1,L
            do z=1,L

              if (v(x,y,z).ne.0) then ! non-zero charge at (x,y,z)

                charge=q*v(x,y,z)
                nch=nch+1
                sum_x=sum_x+charge*(-1)&
                      *(lgf(a,b,c,x,y,z)-lgf(a,b,c,neg(x),y,z))
                sum_y=sum_y+charge*(-1)&
                      *(lgf(a,b,c,x,y,z)-lgf(a,b,c,x,neg(y),z))
                sum_z=sum_z+charge*(-1)&
                      *(lgf(a,b,c,x,y,z)-lgf(a,b,c,x,y,neg(z)))

                if (v(a,b,c).ne.0) then
                  if (a.eq.x.and.b.eq.y.and.c.eq.z) then
                    u_self=u_self+charge**2*lgf(x,y,z,x,y,z)
                  else
                    u_int=u_int+q*v(a,b,c)*charge*lgf(a,b,c,x,y,z)
                  end if
                end if

              end if ! if v(x,y,z).ne.0

            end do ! z do loop
          end do ! y do loop
        end do ! x do loop

        mnphi_x(a,b,c)=sum_x
        mnphi_y(a,b,c)=sum_y
        mnphi_z(a,b,c)=sum_z

        ebar_x=ebar_x+mnphi_x(a,b,c)
        ebar_y=ebar_y+mnphi_x(a,b,c)
        ebar_z=ebar_z+mnphi_x(a,b,c)

        u_tot=u_tot+0.5*(mnphi_x(a,b,c)**2&
              +mnphi_y(a,b,c)**2+mnphi_z(a,b,c)**2)

      end do ! c do loop
    end do ! b do loop
  end do ! a do loop

  nch=nch/L**3 ! bc we count nch once for each abc

  write(*,*)
  write(*,*) " --- linear solver results ---"
  write (*,*) 'irrot E_ij^2 =',u_tot
  write(*,*) "self_e from lgf(0) * n charges = ",u_self
  write(*,*) 'u_int calculated in loop = ',u_int
  write (*,*) 'V*Ebar^2 = ',L**3*ebar_x**2+ebar_y**2+ebar_z**2

  deallocate(cosine)

  end subroutine linsol

  subroutine lgfcalc

  do a=1,L
    do b=1,L
      do c=1,L
        do x=1,L
          do y=1,L
            do z=1,L

              lgf(a,b,c,x,y,z)=0.0

              ! these need to be real, for cos to work later
              p1=float(a-x)
              q1=float(b-y)
              r1=float(c-z)

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

              do kx=-(L-1)/2,L/2
                fkx=2*pi*kx/L
                do ky=-(L-1)/2,L/2
                  fky=2*pi*ky/L
                  do kz=-(L-1)/2,L/2
                    fkz=2*pi*kz/L
                    if ((kx.eq.0).and.(ky.eq.0).and.(kz.eq.0)) then

                    else

                      lgf(a,b,c,x,y,z)=lgf(a,b,c,x,y,z)+(cos(fkx*p1)&
                        *cos(fky*q1)*cos(fkz*r1))&
                        /(3-cos(fkx)-cos(fky)-cos(fkz))

                    end if ! end of kx=ky=kz=0 block
                  end do ! end kz loop
                end do ! end ky loop
              end do ! end kx loop

              lgf(a,b,c,x,y,z)=lgf(a,b,c,x,y,z)/(2*L**3)

            end do ! end z loop
          end do ! end y loop
        end do ! end x loop
      end do ! end c loop
    end do ! end b loop
  end do ! end a loop
  have_lgf=1

  end subroutine lgfcalc
end module linear_solver
