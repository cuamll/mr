! Module containing everything we need for the linear solver
! to get the irrotational E-field part for a cubic lattice
module linear_solver
  use common
  implicit none

  integer, private :: a,b,c,x,y,z,kx,ky,kz,i,nch
  real*8, private :: sum_x,sum_y,sum_z,p1,p2,q1,q2,r1,r2,fkx,fky,fkz,g_x3,g_y3,g_z3
  real*8, private :: m1p1,m1p2,m1q1,m1q2,m1r1,m1r2,charge,g_x1,g_y1,g_z1,g_x2,g_y2,g_z2
  real*8, public :: u_tot, self_e, g_zero, g_z_sum
  real*8, dimension(:), allocatable, private :: cosine
  save

  contains

  subroutine linsol
  ! need to have many things declared already:
  ! L, v(x,y,z), mnabphi_x/y/z, pi, volume, q
!  integer :: a,b,c,x,y,z,kx,ky,kz
!  real*8 :: sum_x,sum_y,sum_z,p1,p2,q1,q2,r1,r2,pi,q,u_tot
!  real*8 :: m1p1,m1p2,m1q1,m1q2,m1r1,m1r2,cosine(L/2-1),charge,volume

  u_tot=0.0
  g_zero=0.0
  allocate(cosine(L/2-1))

  call PBCs
!  call g_zero_calc

  do a=1,L
    do b=1,L
      do c=1,L
        ! Basically, we can think of the {a,b,c} as the sum over \vec{xprime} in
        ! G(\vec{x},vec{xprime}, and the {x,y,z} as the \vec{x}.
        ! {kx,ky,kz} are then the sum over k-space in G(\vec{x},\vec{xprime})

        ! initialise the sums, idiot
        sum_x=0.0
        sum_y=0.0
        sum_z=0.0
        g_x1=0.0
        g_y1=0.0
        g_z1=0.0
        g_x2=0.0
        g_y2=0.0
        g_z2=0.0
        g_x3=0.0
        g_y3=0.0
        g_z3=0.0

        do x=1,L
          do y=1,L
            do z=1,L
              if (v(x,y,z).ne.0) then ! finds a non-zero charge at (x,y,z)
                ! we might want something here depending on charge magnitude
                nch=nch+1
                charge=q*v(x,y,z)

                ! these need to be real numbers for cos etc to work later
                p1=float(a-x)
                p2=float(neg(a)-x)
                q1=float(b-y)
                q2=float(neg(b)-x)
                r1=float(c-z)
                r2=float(neg(c)-z)

                ! these need to be within [L/2,-L/2]. The Green's function is
                ! even though because we use cosines, so it doesn't matter
                ! whether it's positive or negative.

                if (p1.gt.float(L/2)) then
                  p1=p1-float(L)
                else if (p1.lt.(-float(L/2))) then
                  p1=p1+float(L)
                end if

                if (p2.gt.float(L/2)) then
                  p2=p2-float(L)
                else if (p2.lt.(-float(L/2))) then
                  p2=p2+float(L)
                end if

                if (q1.gt.float(L/2)) then
                  q1=q1-float(L)
                else if (q1.lt.(-float(L/2))) then
                  q1=q1+float(L)
                end if

                if (q2.gt.float(L/2)) then
                  q2=q2-float(L)
                else if (q2.lt.(-float(L/2))) then
                  q2=q2+float(L)
                end if

                if (r1.gt.float(L/2)) then
                  r1=r1-float(L)
                else if (r1.lt.(-float(L/2))) then
                  r1=r1+float(L)
                end if

                if (r2.gt.float(L/2)) then
                  r2=r2-float(L)
                else if (r2.lt.(-float(L/2))) then
                  r2=r2+float(L)
                end if


                do kx=-(L-1)/2,L/2
                  fkx=2*pi*kx/L
                  do ky=-(L-1)/2,L/2
                    fky=2*pi*ky/L
                    do kz=-(L-1)/2,L/2
                      fkz=2*pi*kz/L
                      if ((kx.eq.0).and.(ky.eq.0).and.(kz.eq.0)) then
                        g_zero=g_zero ! dunno how to continue or break or whatever
                        sum_x=sum_x
                        sum_y=sum_y
                        sum_z=sum_z
                      else
                        g_zero=g_zero+1/(3-cos(fkx)-cos(fky)-cos(fkz))
                        g_x1=g_x1+cos(fkx)/(3-cos(fkx)-cos(fky)-cos(fkz))
                        g_y1=g_y1+cos(fky)/(3-cos(fkx)-cos(fky)-cos(fkz))
                        g_z1=g_z1+cos(fkz)/(3-cos(fkx)-cos(fky)-cos(fkz))
                        g_x2=g_x2+cos(fkx*2)/(3-cos(fkx)-cos(fky)-cos(fkz))
                        g_y2=g_y2+cos(fky*2)/(3-cos(fkx)-cos(fky)-cos(fkz))
                        g_z2=g_z2+cos(fkz*2)/(3-cos(fkx)-cos(fky)-cos(fkz))
                        g_x3=g_x3+cos(fkx*3)/(3-cos(fkx)-cos(fky)-cos(fkz))
                        g_y3=g_y3+cos(fky*3)/(3-cos(fkx)-cos(fky)-cos(fkz))
                        g_z3=g_z3+cos(fkz*3)/(3-cos(fkx)-cos(fky)-cos(fkz))

                        sum_x=sum_x+charge*((cos(fkx*p2)-cos(fkx*p1))&
                          *(cos(fky*q1)*cos(fkz*r1))/(3-cos(fkx)-cos(fky)-cos(fkz)))
                        sum_y=sum_y+charge*((cos(fkx*q2)-cos(fkx*q1))&
                          *(cos(fky*p1)*cos(fkz*r1))/(3-cos(fkx)-cos(fky)-cos(fkz)))
                        sum_z=sum_z+charge*((cos(fkx*r2)-cos(fkx*r1))&
                        *(cos(fky*q1)*cos(fkz*p1))/(3-cos(fkx)-cos(fky)-cos(fkz)))

                      end if ! end of kx=ky=kz=0 block
                    end do ! end kz loop
                  end do ! end ky loop
                end do ! end kx loop

                ! commented out bc fucked
!                ! These - m1p1 etc - are for simplifying the loops in the case
!                ! where one or more of the (2\pi k_{\mu})/L are (2n+1)\pi, for
!                ! integer n. In this case, the cos will be (-1)^p1 or (-1)^q2,
!                ! etc, as you can see from writing it out in exponentials and
!                ! then factorising. We don't worry about the equivalent for
!                ! 2n\pi, obviously, because 1^n = 1.
!
!                m1p1=(-1)**p1
!                m1q1=(-1)**q1
!                m1r1=(-1)**r1
!                m1p2=(-1)**p2-m1p1
!                m1q2=(-1)**q2-m1q1
!                m1r2=(-1)**r2-m1r1
!
!                ! This vector makes the cos arguments easier to look at
!                ! - add *lambda in denominator if lambda != 1
!                do i=1,((L/2)-1)
!                  cosine(i)=(2*pi*i)/(L)
!                end do
!
!                ! -nabla phi_x(a,b,c) ~ -(G(a,b,c)-G(a-1,b,c))
!                ! The weird terms outside the triple do loop are
!                ! simplifications for when one or more of the
!                ! k_{\mu}s give us (\pm 1)^p1, etc.
!
!                ! if using lambda != 1, need terms like this
!                ! sum_x=sum_x+charge*(cos(pi*r1/lambda)*cos(pi*q1/lambda)&
!                !       *(cos(pi*p2/lambda)-cos(pi*p1/lambda))&
!                !       /(3-3*cos(pi/lambda)))
!
!                ! all three k_mu = L/2:
!                sum_x=sum_x+charge*(m1r1*m1q1*m1p2/6)
!                sum_y=sum_y+charge*(m1r1*m1p1*m1q2/6)
!                sum_z=sum_z+charge*(m1p1*m1q1*m1r2/6)
!                g_zero=g_zero+(1/6)
!
!                ! terms with two L/2 and one 0:
!
!                ! changing the (3/4) (3/2) to (4/4) (4/2) gives the
!                ! right answer to four decimal places. can't be right though
!                sum_x=sum_x+charge*(m1p2*(m1q1+m1r1)/4)
!                sum_y=sum_y+charge*(m1q2*(m1p1+m1r1)/4)
!                sum_z=sum_z+charge*(m1r2*(m1q1+m1p1)/4)
!                g_zero=g_zero+(3/4)
!
!                ! terms with two k_mu = 0 and one k_mu = L/2:
!                sum_x=sum_x+charge*(m1p2/2)
!                sum_y=sum_y+charge*(m1q2/2)
!                sum_z=sum_z+charge*(m1r2/2)
!                g_zero=g_zero+(3/2)
!
!                ! integer division does floor, essentially
!                do kx=1,((L/2)-1)
!                ! terms with two of the k_mu = 0 or L/2 go here
!                sum_x=sum_x+charge*2*(cos(cosine(kx)*p2)-cos(cosine(kx)*p1))*(1+m1q1*m1r1+m1r1+m1q1)&
!                      +charge*2*(m1p2*(cos(cosine(kx)*q1)*(1+m1r1)+(cos(cosine(kx)*r1)*(1+m1q1))))
!                sum_y=sum_y+charge*2*(cos(cosine(kx)*q2)-cos(cosine(kx)*q1))*(1+m1p1*m1r1+m1r1+m1p1)&
!                      +charge*2*(m1q2*(cos(cosine(kx)*p1)*(1+m1r1)+(cos(cosine(kx)*r1)*(1+m1p1))))
!                sum_z=sum_z+charge*2*(cos(cosine(kx)*r2)-cos(cosine(kx)*p1))*(1+m1q1*m1p1+m1p1+m1q1)&
!                      +charge*2*(m1r2*(cos(cosine(kx)*q1)*(1+m1p1)+(cos(cosine(kx)*p1)*(1+m1q1))))
!                g_zero=g_zero+2*((3/(5-cos(cosine(kx))))&
!                       +(3/(1-cos(cosine(kx))))+(6/(3-cos(cosine(kx)))))
!
!                  do ky=1,((L/2)-1)
!                  ! terms with one of the k_mu = 0 or L/2 go here
!                  sum_x=sum_x+charge&
!                        *(4*(cos(cosine(ky)*r1)*(1+m1q1)&
!                        +cos(cosine(ky)*q1)*(1+m1r1)&
!                        *(cos(cosine(kx)*p2)-cos(cosine(kx)*p1)))&
!                        +((cos(cosine(kx)*r1)*cos(cosine(ky)*q1)&
!                          +cos(cosine(kx)*q1)*cos(cosine(ky)*r1))*(m1p2)))
!                  sum_y=sum_y+charge&
!                        *(4*(cos(cosine(ky)*r1)*(1+m1p1)&
!                        +cos(cosine(ky)*p1)*(1+m1r1)&
!                        *(cos(cosine(kx)*q2)-cos(cosine(kx)*q1)))&
!                        +((cos(cosine(kx)*r1)*cos(cosine(ky)*p1)&
!                          +cos(cosine(kx)*p1)*cos(cosine(ky)*r1))*(m1q2)))
!                  sum_z=sum_z+charge&
!                        *(4*(cos(cosine(ky)*p1)*(1+m1q1)&
!                        +cos(cosine(ky)*q1)*(1+m1p1)&
!                        *(cos(cosine(kx)*r2)-cos(cosine(kx)*r1)))&
!                        +((cos(cosine(kx)*p1)*cos(cosine(ky)*q1)&
!                          +cos(cosine(kx)*q1)*cos(cosine(ky)*p1))*(m1r2)))
!                  g_zero=g_zero+4&
!                        *(3/(4-cos(cosine(kx))-cos(cosine(ky)))&
!                        +3/(2-cos(cosine(kx))-cos(cosine(ky))))
!
!                    do kz=1,((L/2)-1)
!
!                      ! There are various ways to make this faster
!                      ! e.g. some components will give the same contribution
!                      ! and also we can pull out some simplifications
!                      ! ...I'm working on it
!
!                      sum_x=sum_x+charge*(-1)*8*(cos(cosine(kx)*p1)&
!                      *cos(cosine(ky)*q1)*cos(cosine(kz)*r1)&
!                      -cos(cosine(kx)*p2)*cos(cosine(ky)*q1)&
!                      *cos(cosine(kz)*r1))/(3-cos(cosine(kx))&
!                      -cos(cosine(ky))-cos(cosine(kz)))
!
!                      sum_y=sum_y+charge*(-1)*8*(cos(cosine(kx)*p1)&
!                      *cos(cosine(ky)*q1)*cos(cosine(kz)*r1)&
!                      -cos(cosine(kx)*p1)*cos(cosine(ky)*q2)&
!                      *cos(cosine(kz)*r1))/(3-cos(cosine(kx))&
!                      -cos(cosine(ky))-cos(cosine(kz)))
!
!                      sum_z=sum_z+charge*(-1)*8*(cos(cosine(kx)*p1)&
!                      *cos(cosine(ky)*q1)*cos(cosine(kz)*r1)&
!                      -cos(cosine(kx)*p1)*cos(cosine(ky)*q1)&
!                      *cos(cosine(kz)*r2))/(3-cos(cosine(kx))&
!                      -cos(cosine(ky))-cos(cosine(kz)))
!
!                    g_zero=g_zero+8*(1/(3-cos(cosine(kx))&
!                      -cos(cosine(ky))-cos(cosine(kz))))
!
!                    end do ! kz do loop
!                  end do ! ky do loop
!                end do ! kx do loop
!
              end if ! if v(x,y,z).ne.0

            end do ! z do loop
          end do ! y do loop
        end do ! x do loop

        ! need to declare these arrays and volume somewhere
        ! had volume here before instead of L**3
        mnphi_x(a,b,c)=sum_x/(2*L**3)
        mnphi_y(a,b,c)=sum_y/(2*L**3)
        mnphi_z(a,b,c)=sum_z/(2*L**3)

        ebar_x=ebar_x+mnphi_x(a,b,c)
        ebar_y=ebar_y+mnphi_x(a,b,c)
        ebar_z=ebar_z+mnphi_x(a,b,c)

        g_zero=g_zero/L**3

        u_tot=u_tot+mnphi_x(a,b,c)**2+mnphi_y(a,b,c)**2+mnphi_z(a,b,c)**2

      end do ! c do loop
    end do ! b do loop
  end do ! a do loop

  ! g_zero is effectively counted once for each charge so the G(0)n_ch
  ! part is taken care of automatically; we just divide by 2.
  self_e=g_zero

  write (*,*) 'sum E_ij^2 =',u_tot

  write (*,*) 'Ebar^2 = ',ebar_x**2+ebar_y**2+ebar_z**2

  write (*,*) 'g(0)*n_charges = ',g_zero
  write (*,*) 'g(a) x = ',g_y1/L**3
  write (*,*) 'g(a) y = ',g_x1/L**3
  write (*,*) 'g(a) z = ',g_z1/L**3
  write (*,*) 'g(2a) x = ',g_y2/L**3
  write (*,*) 'g(2a) y = ',g_x2/L**3
  write (*,*) 'g(2a) z = ',g_z2/L**3
  write (*,*) 'g(3a) x = ',g_y3/L**3
  write (*,*) 'g(3a) y = ',g_x3/L**3
  write (*,*) 'g(3a) z = ',g_z3/L**3

  write (*,*) 'interaction????? = ',u_tot-2*((self_e)-(L**3*(ebar_x**2+ebar_y**2+ebar_z**2)/2))

  end subroutine linsol

  subroutine g_zero_calc

    g_z_sum=0.0

    do i=1,((L/2)-1)
      cosine(i)=(2*pi*i)/(L)
    end do

    g_z_sum=g_z_sum+(1/6)
    g_z_sum=g_z_sum+(3/4)
    g_z_sum=g_z_sum+(3/2)

    do kx=1,((L/2)-1)

    g_z_sum=g_z_sum+2*((3/(5-cos(cosine(kx))))&
           +(3/(1-cos(cosine(kx))))+(6/(3-cos(cosine(kx)))))

      do ky=1,((L/2)-1)

        g_z_sum=g_z_sum+4&
              *(3/(4-cos(cosine(kx))-cos(cosine(ky)))&
              +3/(2-cos(cosine(kx))-cos(cosine(ky))))

        do kz=1,((L/2)-1)
        
          g_z_sum=g_z_sum+8*(1/(3-cos(cosine(kx))&
            -cos(cosine(ky))-cos(cosine(kz))))

        end do
      end do
    end do

    g_zero=g_z_sum/L**3

    write(*,*) g_zero
    write(*,*) "fuck you"

  end subroutine g_zero_calc
end module linear_solver
