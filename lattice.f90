program config_linsolve
integer :: L, i, j, k, row, col
character(len=11) :: filename
integer, dimension(:,:), allocatable :: v_temp
integer, dimension(:,:,:), allocatable :: v,mnphi_x,mnphi_y,mnphi_z

! note: getarg(pos, value) does argc/argv

! these parameters will be put in an input file eventually
filename='lattice.dat'

! read in lattice size (for now) then allocate arrays accordingly
write(*,*) "Lattice dimension L: "
read(*,*) L
allocate(v_temp(L**2,L))
allocate(v(L,L,L))

! irrotational part of E-field
allocate(mnphi_x(L,L,L))
allocate(mnphi_y(L,L,L))
allocate(mnphi_z(L,L,L))

open(unit=1, file=filename)
read(1,*)((v_temp(row,col),col=1,L),row=1,L**2)

do i=1,L**2
  write(*,*) (v_temp(i,j), j=1,L)
enddo

v = reshape(v_temp, (/L,L,L/))
deallocate(v_temp) ! we don't need it anymore

do j=1,L
  do i=1,L
    write(*,*) (v(i,j,k), k=1,L)
  enddo
enddo

call linear_solver

end program config_linsolve

! Subroutine to find the irrotational part of the lattice field
! for a cubic lattice

subroutine linear_solver
implicit none
! need to have many things declared already:
! L, v(x,y,z), mnabphi_x/y/z, pi, volume, q
integer :: a,b,c,x,y,z,kx,ky,kz
real*8 :: sum_x,sum_y,sum_z,p1,p2,q1,q2,r1,r2,pi,q,u_tot
real*8 :: m1p1,m1p2,m1q1,m1q2,m1r1,m1r2,twopibyL(L/2-1),charge,volume

pi=3.141592653589793
u_tot=0
volume=1.0

call PBCs

do a=1,L
  do b=1,L
    do c=1,L
      ! Basically, we can think of the {a,b,c} as the sum over \vec{xprime} in
      ! G(\vec{x},vec{xprime}, and the {x,y,z} as the \vec{x}.
      ! {kx,ky,kz} are then the sum over k-space in G(\vec{x},\vec{xprime})

      ! initialise the sums, idiot
      sum_x=0
      sum_y=0
      sum_z=0

      do x=1,L
        do y=1,L
          do z=1,L
            if (v(x,y,z).ne.0) then ! finds a non-zero charge at (x,y,z)
              ! we might want something here depending on charge magnitude
              q=0.2
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

              if (q2.gt.float(l/2)) then
                q2=q2-float(l)
              else if (q2.lt.(-float(l/2))) then
                q2=q2+float(l)
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

              ! These - m1p1 etc - are for simplifying the loops in the case
              ! where one or more of the (2\pi k_{\mu})/L are (2n+1)\pi, for
              ! integer n. In this case, the cos will be (-1)^p1 or (-1)^q2,
              ! etc, as you can see from writing it out in exponentials and
              ! then factorising. We don't worry about the equivalent for
              ! 2n\pi, obviously, because 1^n = 1.

              m1p1=(-1)**p1
              m1q1=(-1)**q1
              m1r1=(-1)**r1
              m1p2=(-1)**p2-m1p1
              m1q2=(-1)**q2-m1q1
              m1r2=(-1)**r2-m1r1

              ! This vector makes the cos arguments easier to look at
              do i=1,((L/2)-1)
                twopibyL(i)=(2*pi*i)/L
              end do

              ! Now for the actual sum over k-space. Here we define
              ! -nabla phi_x(a,b,c) ~ (G(a,b,c)-G(a-1,b,c)) - i.e. we define
              ! nabla as a backward-difference grad, I find that easier to
              ! explain it. As explained above the weird terms outside the
              ! triple do loop are simplifications for when one or more of the
              ! k_{\mu}s give us (\pm 1)^p1, etc.

              do kx=1,((L/2)-1)
                do ky=1,((L/2)-1)
                  do kz=1,((L/2)-1)

                    ! There are various ways to make this faster
                    ! e.g. some components will give the same contribution
                    ! and also we can pull out some simplifications
                    ! ...I'm working on it

                    sum_x=sum_x+charge*(cos(twopibyL(kx)*p1)&
                    *cos(twopibyL(ky)*q1)*cos(twopibyL(kz)*r1)&
                    -cos(twopibyL(kx)*p2)*cos(twopibyL(ky)*q1)&
                    *cos(twopibyL(kz)*r1))/(3-cos(twopibyL(kx))&
                    -cos(twopibyL(ky))-cos(twopibyL(kz)))

                    sum_y=sum_y+charge*(cos(twopibyL(kx)*p1)&
                    *cos(twopibyL(ky)*q1)*cos(twopibyL(kz)*r1)&
                    -cos(twopibyL(kx)*p1)*cos(twopibyL(ky)*q2)&
                    *cos(twopibyL(kz)*r1))/(3-cos(twopibyL(kx))&
                    -cos(twopibyL(ky))-cos(twopibyL(kz)))

                    sum_z=sum_z+charge*(cos(twopibyL(kx)*p1)&
                    *cos(twopibyL(ky)*q1)*cos(twopibyL(kz)*r1)&
                    -cos(twopibyL(kx)*p1)*cos(twopibyL(ky)*q1)&
                    *cos(twopibyL(kz)*r2))/(3-cos(twopibyL(kx))&
                    -cos(twopibyL(ky))-cos(twopibyL(kz)))

                  end do ! kz do loop
                end do ! ky do loop
              end do ! kx do loop

            end if ! if v(x,y,z).ne.0

          end do ! z do loop
        end do ! y do loop
      end do ! x do loop

      ! need to declare these arrays and volume somewhere
      mnphi_x(a,b,c)=sum_x/volume
      mnphi_y(a,b,c)=sum_y/volume
      mnphi_z(a,b,c)=sum_z/volume

      u_tot=u_tot+mnphi_x(a,b,c)**2+mnphi_y(a,b,c)**2+mnphi_z(a,b,c)**2

    end do ! c do loop
  end do ! b do loop
end do ! a do loop

write (*,*) u_tot

end subroutine linear_solver

subroutine PBCs
integer :: i
do i=1,L
  pos(i)=mod(i,L)+1
  neg(i)=mod(i+L-2,L)+1
end do
return
end subroutine PBCs
