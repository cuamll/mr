!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!        A linear solver for Poisson's equation of          !
!      EM on a two-dimensional, square lattice (FFT)        !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine linear_solver
include "common.h"
integer a,b,i,j,k,l
real*8 sum_x,sum_y,p1,q1,p2,q2
real*8 Factp1,Factp2,Factq1,Factq2,cosine(side/2-1),pivort

! ************* Note that common.h calls stored variables in the full script for the 2D Coulomb gas
! ************* I think if you just make a 2D array v(i,j) (charge configuration)
! ************* and define an E-field array this code should work
! ************* Note that E field is currently Ehat due to historic misunderstandings:
! ************* Ehat=Ebar-nabla phi, but this code finds -nabla phi

call PBC! sets up PBCs

! sweep over lattice to find G(a,b) everywhere

do a=1,side
   do b=1,side

! SUPERPOSITION PRINCIPLE - sweep over lattice to find contribution to field at (a,b) due to charge at (i,j)

      sum_x=0! sets x-component sum to 0
      sum_y=0! sets y-component sum to 0
      do i=1,side
         do j=1,side
            if (v(i,j).ne.0) then! v(i,j) is the array of normalized-value charges, so this finds a non-zero charge at site (i,j)
               pivort=pi*v(i,j)! takes the value q/2 = 2pi/2 = pi (since the elementary charge q is set to 2pi to make comparison with XY physics)

! CALCULATE DISPLACEMENT VECTOR FROM CHARGE TO FIELD SITE
               
               p1=float(a-i)! makes the x component of the separation vector between Green's function site and charge site a real number

! brings separation vector back within the periodicity of the lattice ([-L/2,L/2] - notice that L/2 and -L/2 give the same Green's function contr., so you don't have to pick one)

               if (p1.gt.float(side)/2) then
                  p1=p1-(float(side))
               else if (p1.lt.-float(side)/2) then
                  p1=p1+float(side)
               end if

! you need the site immediately negative to (a,b) as well to find -nabla phi(a,b) ~ -(G(a,b)-G(a-1,b))

               p2=float(neg(a)-i)
               if (p2.gt.float(side)/2) then
                  p2=p2-float(side)
               else if (p2.lt.-float(side)/2) then
                  p2=p2+float(side)
               end if

! same for y component

               q1=float(b-j)
               if (q1.gt.float(side)/2) then
                  q1=q1-(float(side))
               else if (q1.lt.-float(side)/2) then
                  q1=q1+float(side)
               end if

               q2=float(neg(b)-j)
               if (q2.gt.float(side)/2) then
                  q2=q2-float(side)
               else if (q2.lt.-float(side)/2) then
                  q2=q2+float(side)
               end if

! calculate components of the sum over k-space that are used more than once

               Factp1=(-1)**p1
               Factp2=(-1)**p2-Factp1
               Factq1=(-1)**q1
               Factq2=(-1)**q2-Factq1
               do k=1,side/2-1
                  cosine(k)=cos(twopi*k/side)
               end do

! SUM OVER K-SPACE - note that we will define Ehat_x(a,b) ~ -(G(a,b)-G(a-1,b)) in order to agree with top-field/vort def. This should be -nabla phi, not Ehat
! Going through the LGF sum, you find that it splits up like the following

               sum_x=sum_x+pivort*(0.25*Factq1+0.5)*Factp2
               sum_y=sum_y+pivort*(0.25*Factp1+0.5)*Factq2

               do k=1,side/2-1
                  sum_x=sum_x+pivort*(2*(cos(twopi*q1*k/side)*Factp2+Factq1*(cos(twopi*k*p2/side)&
                       -cos(twopi*k*p1/side)))/(3-cosine(k))&
                       +2*(cos(twopi*k*p2/side)-cos(twopi*k*p1/side))/(1-cosine(k)))

                  sum_y=sum_y+pivort*(2*(cos(twopi*p1*k/side)*Factq2+Factp1*(cos(twopi*k*q2/side)&
                       -cos(twopi*k*q1/side)))/(3-cosine(k))&
                       +2*(cos(twopi*k*q2/side)-cos(twopi*k*q1/side))/(1-cosine(k)))
                  do l=1,side/2-1
                     sum_x=sum_x+pivort*(4*cos(twopi*l*q1/side)*(cos(twopi*k*p2/side)-cos(twopi*k*p1/side))&
                          /(2-cosine(k)-cosine(l)))

                     sum_y=sum_y+pivort*(4*cos(twopi*l*p1/side)*(cos(twopi*k*q2/side)-cos(twopi*k*q1/side))&
                          /(2-cosine(k)-cosine(l)))
                  end do
               end do

            end if
         end do
      end do

      Ehat_x(a,b)=sum_x/volume! this was called Ehat before i realized that this doesn't include the k=0 mode
      Ehat_y(a,b)=sum_y/volume
   end do
end do
return
end subroutine


!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!             Sets up our lattice with periodic             !
!                     boundary conditions                   !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine PBC
include "common.h"
integer i
do i=1,side
   pos(i)=mod(i,side)+1
   neg(i)=mod(i+side-2,side)+1
end do
return
end subroutine
