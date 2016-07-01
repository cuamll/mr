!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!  A program to simulate a 2D Maggs Gas (quasiGCE) (|v|<1)  !
!         n.b. q = 2pi * v for ease of comparison           !
!          we attempt to incude global sampling             !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

! quasiGCE refers to our allowing the HXY Model to generate vortex number at each temp., then fixing this from there

program maggs
include "common.h"
integer a,b,i,j,k,l,m,seed,start
open (unit=1,file='initial.in')
open (unit=8,file='sustot.dat')
open (unit=9,file='susdip.dat')
open (unit=100,file='suswind.dat')
call input(seed,start)
call randinit(seed)
write(6,*) rand(seed)
call PBC

! We set field and charge values to zero everywhere

call initial_field
!glob=0! with glob=0, we disallow global updates. As soon as winding becomes non-zero, BKT has broken symmetry so we employ winding updates to improve ergodicity in this new Gibbs' ensemble.

! We will now enter out temperature loop, entering at our min. temp.

T=Tmin
Tincr=(Tmax-Tmin)/Tsteps
do m=0,Tsteps
   write(6,*) T
   beta=1/T
   call initial_measure               ! sets all measurement data to zero
   call sum_field                     ! sums the field to find harmonic contribution
   call fluctuations(therm_sweeps)    ! this thermalizes the system
   quench_fluct=quench_fluct+delta_fluct ! this adds delta(delta phi) to delta phi, the change in the aux. field variable - it is called quench_fluct because that variable was free from an older script. We add delta_fluct at each new temperature in order to keep acceptance rates between 40% and 60%
   if ((T.gt.1.75).and.(T.lt.1.85)) then ! We reset to delta phi = 2.1 in this temperature range as this keeps acceptance rates between 40% and 60%
      quench_fluct=2.1
   end if
   do j=1,measurements
      call fluctuations(sweeps)       !  We allow for thermal fluctuations a 'measurements' amount of times
      call sum_field
      call measure
   end do
   call outputs
   T=T+Tincr
end do
close (1)
close (2)
close (3)
close (4)
stop
end program

!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!            Reads in input values and writes               !
!                information to outfiles                    !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine input(seed, start)
include "common.h"
integer seed,start
read(1,*) side
read(1,*) therm_sweeps
read(1,*) sweeps
read(1,*) ratiorot
read(1,*) ratiowind
read(1,*) bin_size
read(1,*) num_bins
read(1,*) quench_fluct
read(1,*) delta_fluct
read(1,*) Tmin
read(1,*) Tmax
read(1,*) Tsteps
read(1,*) start
read(1,*) seed
read(1,*) outfile1
read(1,*) outfile2
read(1,*) outfile3
read(1,*) outfile4
read(1,*) spinin! files were called spinsin/out in XY, so just leave it
read(1,*) spinout

open(2,file=outfile1)
open(3,file=outfile2)
open(5,file=outfile3)
open(7,file=outfile4)

write(2,*) 'Filename is ',outfile1
write(2,*) 'Graph data is in ',outfile2
write(2,*) 'Cumulant data is in ',outfile3
write(2,*) 'Dielectric data is in ',outfile4
write(2,*) 'Susceptibility data is in ',outfile5
write(2,*) 'Lattice size = ',side

if ((side>max_sides).or.(side<1)) then
   write(2,*) 'Lattice length not allowed. Please make note of max_side and try again with new inputs.'
   stop
end if
sites=side*side
volume=float(sites)
Tlow=Tmin
measurements=bin_size*num_bins

write(2,*) 'Volume: ',sites
write(2,*) 'Number of thermalization sweeps: ',therm_sweeps
write(2,*) 'Number of sweeps per measurement: ',sweeps
write(2,*) 'Bin size: ',bin_size
write(2,*) 'Number of bins: ',num_bins
write(2,*) 'Number of measurements: ',measurements
write(2,*) 'Tmin : ',Tmin
write(2,*) 'Tmax : ',Tmax
write(2,*) 'Tsteps : ',Tsteps
write(2,*) 'Ordered/cold (0) or random/hot (1) start: ',start
write(2,*) 'Seed for the random number generator ',seed
write(2,*) '                                                                         '
write(2,*) '--------------------------------'
write(2,*) '                                                                         '
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

!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!       Initializes our starting spin configuration         !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine initial_spins(start)
include "common.h"
integer start,i,j
real rand
if (start.eq.0) then
   do i=1,side
      do j=1,side
         theta(i,j)=pibytwo
      end do
   end do
else
   do i=1,side
      do j=1,side
         theta(i,j)=rand()*twopi
      end do
   end do
end if
return
end subroutine


!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!       Initializes our starting spin configuration         !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine initial_field
include "common.h"
integer start,i,j
do i=1,side
   do j=1,side
      top_x(i,j)=0
      top_y(i,j)=0
      v(i,j)=0
   end do
end do
return
end subroutine


!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!        Calls stored field and charge configs              !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine field_configuration
include "common.h"
integer i,j
real*8 x
open (10,file=spinin)
do i=1,side
   do j=1,side
      read(10,*) x
      top_x(i,j)=x
      read(10,*) x
      top_y(i,j)=x
      read(10,*) x
      v(i,j)=x
   end do
end do
close (10)
return
end subroutine


!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!                 Thermal fluctuations                      !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine fluctuations(iterations)
include "common.h"
integer i,j,n,m,iterations,v0
real rand
real*8 oldE,newE,deltaE,top1o,top2o,top3o,top4o,top1n,top2n,top3n,top4n,delta,v1,v2
do n=1,iterations

! charge-hop update
   accepth=0

   do m=1,2*side*side!*2! doubled the number of hop attempts by hand
      i=int(rand()*side)+1! picks an x-coordinate at random
      j=int(rand()*side)+1! picks an y-coordinate at random

      delta=rand()!picks x- or y-component at random
      if (delta.lt.0.5) then
         top1o=top_x(i,j)
         delta=rand()
         if (delta.lt.0.5) then! picks an increase or decrease in field bond at random
            top1n=top1o+twopi
            oldE=0.5*top1o**2
            newE=0.5*top1n**2
            deltaE=newE-oldE
            v1=v(i,j)-1
            v2=v(neg(i),j)+1
            if ((abs(v1).le.1).and.(abs(v2).le.1)) then! disallowing |v|>1
               if ((deltaE.lt.0).or.(exp(-beta*deltaE).gt.rand())) then
                  top_x(i,j)=top1n
                  v(i,j)=v1
                  v(neg(i),j)=v2
                  Ebar_x=Ebar_x+twopi/volume! track Ebar evolution
                  accepth=accepth+1
                  if (((Ebar_x.gt.(pi/(float(side)))).or.(Ebar_x.le.(-pi/(float(side))))).and.(glob.eq.0)) then
                     glob=1
                  end if
               end if
            end if
         else!                   decrease field bond
            top1n=top1o-twopi
            oldE=0.5*top1o**2
            newE=0.5*top1n**2
            deltaE=newE-oldE
            v1=v(i,j)+1
            v2=v(neg(i),j)-1
            if ((abs(v1).le.1).and.(abs(v2).le.1)) then! disallowing |v|>1
               if ((deltaE.lt.0).or.(exp(-beta*deltaE).gt.rand())) then
                  top_x(i,j)=top1n
                  v(i,j)=v1
                  v(neg(i),j)=v2
                  Ebar_x=Ebar_x-twopi/volume! track Ebar evolution
                  accepth=accepth+1
                  if (((Ebar_x.gt.(pi/(float(side)))).or.(Ebar_x.le.(-pi/(float(side))))).and.(glob.eq.0)) then
                     glob=1
                  end if
               end if
            end if
         end if
      else!                      y-component
         top1o=top_y(i,j)
         delta=rand()
         if (delta.lt.0.5) then! picks an increase or decrease in field bond at random
            top1n=top1o+twopi
            oldE=0.5*top1o**2
            newE=0.5*top1n**2
            deltaE=newE-oldE
            v1=v(i,j)-1
            v2=v(i,neg(j))+1
            if ((abs(v1).le.1).and.(abs(v2).le.1)) then! disallowing |v|>1
               if ((deltaE.lt.0).or.(exp(-beta*deltaE).gt.rand())) then
                  top_y(i,j)=top1n
                  v(i,j)=v1
                  v(i,neg(j))=v2
                  Ebar_y=Ebar_y+twopi/volume! track Ebar evolution
                  accepth=accepth+1
                  if (((Ebar_y.gt.(pi/(float(side)))).or.(Ebar_y.le.(-pi/(float(side))))).and.(glob.eq.0)) then
                     glob=1
                  end if
               end if
            end if
         else!                   decrease field bond
            top1n=top1o-twopi
            oldE=0.5*top1o**2
            newE=0.5*top1n**2
            deltaE=newE-oldE
            v1=v(i,j)+1
            v2=v(i,neg(j))-1
            if ((abs(v1).le.1).and.(abs(v2).le.1)) then! disallowing |v|>1
               if ((deltaE.lt.0).or.(exp(-beta*deltaE).gt.rand())) then
                  top_y(i,j)=top1n
                  v(i,j)=v1
                  v(i,neg(j))=v2
                  Ebar_y=Ebar_y-twopi/volume! track Ebar evolution
                  accepth=accepth+1
                  if (((Ebar_y.gt.(pi/(float(side)))).or.(Ebar_y.le.(-pi/(float(side))))).and.(glob.eq.0)) then
                     glob=1
                  end if
               end if
            end if
         end if
      end if

   end do

! rot DoF update

   accept=0
   do m=1,side*side*ratiorot
      i=int(rand()*side)+1! picks an x-coordinate at random
      j=int(rand()*side)+1! picks an y-coordinate at random

      delta=2*quench_fluct*(rand()-0.5)! \in [-quench_fluct, quench_fluct) quench_fluct since this was used in the quench subroutine

      top1o=top_x(i,j)
      top2o=top_y(i,j)
      top3o=top_x(i,neg(j))
      top4o=top_y(neg(i),j)

      top1n=top1o+delta
      top2n=top2o-delta
      top3n=top3o-delta
      top4n=top4o+delta

      oldE=0.5*(top1o**2+top2o**2+top3o**2+top4o**2)
      newE=0.5*(top1n**2+top2n**2+top3n**2+top4n**2)
      deltaE=newE-oldE
      if ((deltaE.lt.0).or.(exp(-beta*deltaE).gt.rand())) then
         top_x(i,j)=top1n
         top_y(i,j)=top2n
         top_x(i,neg(j))=top3n
         top_y(neg(i),j)=top4n
         accept=accept+1
      end if
   end do

! global update - we attempt discetized hops equivalent to sampling winding fluctuations - in this version we add the restriction that
!                 Ebar /= 0

   acceptg=0
!   if (glob.eq.1) then
      do m=1,ratiowind
         ! x-component
         delta=rand()! picks a plus or minus global hop
         if (delta.lt.0.5) then
            deltaE=twopi*(pi + float(side)*Ebar_x)
            if ((deltaE.lt.0).or.((exp(-beta*deltaE).gt.rand()).and.(exp(-beta*deltaE).gt.0.000000000001))) then
               Ebar_x=Ebar_x+twopi/(float(side))
               acceptg=acceptg+1
               do i=1,side
                  do j=1,side
                     top_x(i,j)=top_x(i,j)+twopi/(float(side))
                  end do
               end do
            end if
         else
            deltaE=twopi*(pi - float(side)*Ebar_x)
            if ((deltaE.lt.0).or.((exp(-beta*deltaE).gt.rand()).and.(exp(-beta*deltaE).gt.0.000000000001))) then
               Ebar_x=Ebar_x-twopi/(float(side))
               acceptg=acceptg+1
               do i=1,side
                  do j=1,side
                     top_x(i,j)=top_x(i,j)-twopi/(float(side))
                  end do
               end do
            end if
         end if
      
         ! y-component
         delta=rand()! picks a plus or minus global hop
         if (delta.lt.0.5) then
            deltaE=twopi*(pi + float(side)*Ebar_y)
            if ((deltaE.lt.0).or.((exp(-beta*deltaE).gt.rand()).and.(exp(-beta*deltaE).gt.0.000000000001))) then
               Ebar_y=Ebar_y+twopi/(float(side))
               acceptg=acceptg+1
               do i=1,side
                  do j=1,side
                     top_y(i,j)=top_y(i,j)+twopi/(float(side))
                  end do
               end do
            end if
         else
            deltaE=twopi*(pi - float(side)*Ebar_y)
            if ((deltaE.lt.0).or.((exp(-beta*deltaE).gt.rand()).and.(exp(-beta*deltaE).gt.0.000000000001))) then
               Ebar_y=Ebar_y-twopi/(float(side))
               acceptg=acceptg+1
               do i=1,side
                  do j=1,side
                     top_y(i,j)=top_y(i,j)-twopi/(float(side))
                  end do
               end do
            end if
         end if
      end do
!   end if

end do
return
end subroutine


!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!       Writes the final spin configuration of the          !
!   relevant temperature to spins.dat (via initial.in)      !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine record_fields
include "common.h"
integer i,j
open (11,file=spinout)
do i=1,side
   do j=1,side
      write(11,*) top_x(i,j)
      write(11,*) top_y(i,j)
      write(11,*) v(i,j)
   end do
end do
close(11)
return
end subroutine


!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!     Sets all measurements to zero for new temperature     !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine initial_measure
include "common.h"
ibin=0
sum_accept=0.0
sum_accepth=0.0
sum_acceptg=0.0

sum_energy=0.0
sum_enersq=0.0

sum_vort=0.0
sum_vortsq=0.0

sum_Ebar_x=0.0
sum_Ebar_xsq=0.0
sum_Ebar_y=0.0
sum_Ebar_ysq=0.0

sum_Ebardip_x=0.0
sum_Ebardip_xsq=0.0
sum_Ebardip_y=0.0
sum_Ebardip_ysq=0.0

sum_Ebarwind_x=0.0
sum_Ebarwind_xsq=0.0
sum_Ebarwind_y=0.0
sum_Ebarwind_ysq=0.0

bin_energy  = 0.0
bin_enersq=0.0
bin_vort=0.0
bin_Ebar_x=0.0
bin_Ebar_xsq=0.0
bin_Ebar_y=0.0
bin_Ebar_ysq=0.0
bin_Ebardip_x=0.0
bin_Ebardip_xsq=0.0
bin_Ebardip_y=0.0
bin_Ebardip_ysq=0.0
bin_Ebarwind_x=0.0
bin_Ebarwind_xsq=0.0
bin_Ebarwind_y=0.0
bin_Ebarwind_ysq=0.0

avbin_enersq=0.0
avbin_spheat   = 0.0
avbin_spheatsq=0.0
avbin_vortsq=0.0
avbin_totsus=0.0
avbin_totsussq=0.0
avbin_dipsus=0.0
avbin_dipsussq=0.0
avbin_windsus=0.0
avbin_windsussq=0.0

return
end subroutine


!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
! Measures the results of allowing for thermal fluctuations !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine measure
include "common.h"
integer i,j,k,l,m,n,a,x,y,z
real*8 energy,ener
real*8 bin_mean,bin_spheat,bin_meansq,bin_meanf,bin_susp
real*8 top1,top2,hel_e,hel_s
ener=0
vort=0

do i=1,side
   do j=1,side
      top1=top_y(i,j)
      top2=top_x(i,j)
      ener=ener+0.5*(top1**2+top2**2)
      if (v(i,j).ne.0) then
         vort=vort+1
      end if
   end do
end do

! normalize variables - note that Ebar variables are normalized within sum_field

energy=ener/volume
vort=vort/volume
if (ratiorot.ne.0) then
   sum_accept=sum_accept+float(accept)/(volume*ratiorot)
end if
sum_accepth=sum_accepth+float(accepth)/(2*volume)! since this update is over the 2N field bonds
if (ratiowind.ne.0) then
   sum_acceptg=sum_acceptg+float(acceptg)/(2*ratiowind)! since this update is over the 2 Ebar components
end if

! bin data for error bars

ibin=ibin+1

bin_energy = bin_energy + energy
bin_enersq=bin_enersq+energy**2

bin_vort=bin_vort+vort

bin_Ebar_x=bin_Ebar_x+Ebar_x
bin_Ebar_xsq=bin_Ebar_xsq+Ebar_x**2
bin_Ebar_y=bin_Ebar_y+Ebar_y
bin_Ebar_ysq=bin_Ebar_ysq+Ebar_y**2
bin_Ebarsq=bin_Ebarsq+Ebar_x**2+Ebar_y**2

bin_Ebardip_x=bin_Ebardip_x+Ebardip_x
bin_Ebardip_xsq=bin_Ebardip_xsq+Ebardip_x**2
bin_Ebardip_y=bin_Ebardip_y+Ebardip_y
bin_Ebardip_ysq=bin_Ebardip_ysq+Ebardip_y**2
bin_Ebardipsq=bin_Ebardipsq+Ebardip_x**2+Ebardip_y**2

bin_Ebarwind_x=bin_Ebarwind_x+Ebarwind_x
bin_Ebarwind_xsq=bin_Ebarwind_xsq+Ebarwind_x**2
bin_Ebarwind_y=bin_Ebarwind_y+Ebarwind_y
bin_Ebarwind_ysq=bin_Ebarwind_ysq+Ebarwind_y**2
bin_Ebarwindsq=bin_Ebarwindsq+Ebarwind_x**2+Ebarwind_y**2

! open the bin - we bin data into the bin if ibin=bin_size

if (ibin.eq.bin_size) then
   sum_energy  = sum_energy  + bin_energy
   sum_enersq  = sum_enersq  + bin_enersq

   sum_vort=sum_vort+bin_vort

   sum_Ebar_x=sum_Ebar_x+bin_Ebar_x
   sum_Ebar_xsq=sum_Ebar_xsq+bin_Ebar_xsq
   sum_Ebar_y=sum_Ebar_y+bin_Ebar_y
   sum_Ebar_ysq=sum_Ebar_ysq+bin_Ebar_ysq

   sum_Ebardip_x=sum_Ebardip_x+bin_Ebardip_x
   sum_Ebardip_xsq=sum_Ebardip_xsq+bin_Ebardip_xsq
   sum_Ebardip_y=sum_Ebardip_y+bin_Ebardip_y
   sum_Ebardip_ysq=sum_Ebardip_ysq+bin_Ebardip_ysq

   sum_Ebarwind_x=sum_Ebarwind_x+bin_Ebarwind_x
   sum_Ebarwind_xsq=sum_Ebarwind_xsq+bin_Ebarwind_xsq
   sum_Ebarwind_y=sum_Ebarwind_y+bin_Ebarwind_y
   sum_Ebarwind_ysq=sum_Ebarwind_ysq+bin_Ebarwind_ysq

! calculate bin means and mean square

    bin_mean = bin_energy / bin_size
    bin_meansq = bin_mean**2
    avbin_enersq = avbin_enersq + bin_meansq
    bin_spheat = (bin_enersq/bin_size) - bin_meansq
    avbin_spheat = avbin_spheat + bin_spheat
    avbin_spheatsq = avbin_spheatsq + bin_spheat**2

    bin_mean = bin_vort / bin_size
    bin_meansq = bin_mean**2
    avbin_vortsq = avbin_vortsq + bin_meansq

    bin_mean = (bin_Ebar_xsq + bin_Ebar_ysq) / bin_size
    bin_meansq = bin_mean**2
    avbin_Ebarf=avbin_Ebarf+bin_meansq

    bin_mean = (bin_Ebardip_xsq + bin_Ebardip_ysq) / bin_size
    bin_meansq = bin_mean**2
    avbin_Ebardipf=avbin_Ebardipf+bin_meansq

    bin_mean = (bin_Ebarwind_xsq + bin_Ebarwind_ysq) / bin_size
    bin_meansq = bin_mean**2
    avbin_Ebarwindf=avbin_Ebarwindf+bin_meansq


! reset bin sums to 0

    ibin = 0
    bin_energy  = 0.0
    bin_enersq=0.0
    bin_vort=0.0
    bin_Ebar_x=0
    bin_Ebar_xsq=0
    bin_Ebar_y=0
    bin_Ebar_ysq=0
    bin_Ebardip_x=0
    bin_Ebardip_xsq=0
    bin_Ebardip_y=0
    bin_Ebardip_ysq=0
    bin_Ebarwind_x=0
    bin_Ebarwind_xsq=0
    bin_Ebarwind_y=0
    bin_Ebarwind_ysq=0

end if

return
end subroutine


!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!             Write results to output files                 !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine outputs
include "common.h"
real*8 err_ener,err_spheat,err_vort,err_totsus,err_dipsus,err_windsus
real*8 variance,Ebarsus,Ebardipsus,Ebarwindsus
av_energy=sum_energy/measurements
av_enersq=sum_enersq/measurements
av_vort=sum_vort/measurements

av_accept=sum_accept/measurements
av_accepth=sum_accepth/measurements
av_acceptg=sum_acceptg/measurements

av_Ebar_x=sum_Ebar_x/measurements
av_Ebar_y=sum_Ebar_y/measurements
av_Ebar_xsq=sum_Ebar_xsq/measurements
av_Ebar_ysq=sum_Ebar_ysq/measurements

av_Ebardip_x=sum_Ebardip_x/measurements
av_Ebardip_y=sum_Ebardip_y/measurements
av_Ebardip_xsq=sum_Ebardip_xsq/measurements
av_Ebardip_ysq=sum_Ebardip_ysq/measurements

av_Ebarwind_x=sum_Ebarwind_x/measurements
av_Ebarwind_y=sum_Ebarwind_y/measurements
av_Ebarwind_xsq=sum_Ebarwind_xsq/measurements
av_Ebarwind_ysq=sum_Ebarwind_ysq/measurements

avbin_enersq=avbin_enersq/num_bins

avbin_spheat   = avbin_spheat   / num_bins
avbin_spheatsq = avbin_spheatsq / num_bins

avbin_vortsq=avbin_vortsq/num_bins

avbin_Ebarf=avbin_Ebarf/num_bins
avbin_Ebardipf=avbin_Ebardipf/num_bins
avbin_Ebarwindf=avbin_Ebarwindf/num_bins

!  Work out variances and errors.

! Energy

variance = avbin_enersq - av_energy**2
err_ener = sqrt( variance / num_bins )

! Specific heat 

specific_heat  = av_enersq - av_energy**2
specific_heat  = specific_heat * volume * beta**2
variance = avbin_spheatsq - avbin_spheat**2
err_spheat = sqrt( variance / num_bins ) * volume * beta**2

! charge density (labelled vort)

variance = avbin_vortsq - av_vort**2
err_vort = sqrt( variance / num_bins )

! Ebar susceptibilities

Ebarsus=volume * beta * (av_Ebar_xsq + av_Ebar_ysq - av_Ebar_x**2 - av_Ebar_y**2)
variance = avbin_Ebarf - av_Ebar_x**2 - av_Ebar_y**2
err_totsus=volume * beta * sqrt(variance/num_bins)

Ebardipsus=volume * beta * (av_Ebardip_xsq + av_Ebardip_ysq - av_Ebardip_x**2 - av_Ebardip_y**2)
variance = avbin_Ebardipf - av_Ebardip_x**2 - av_Ebardip_y**2
err_dipsus=volume * beta * sqrt(variance/num_bins)

Ebarwindsus=volume * beta * (av_Ebarwind_xsq + av_Ebarwind_ysq - av_Ebarwind_x**2 - av_Ebarwind_y**2)
variance = avbin_Ebarwindf - av_Ebarwind_x**2 - av_Ebarwind_y**2
err_windsus=volume * beta * sqrt(variance/num_bins)

! print results

write(3,400) T,av_energy,av_enersq,specific_heat,av_vort,err_ener,err_spheat,err_vort
write(5,600) T,av_accept,av_accepth,av_acceptg
write(7,700) T,Ebarsus,Ebardipsus,Ebarwindsus,err_totsus,err_dipsus,err_windsus

100  format(A,F12.5)
200  format(A,F10.6,A,F9.6)
300  format(A,F12.5,A,F11.5)
400  format(1X,12F18.8)
500  format(1X,7F18.8)
600  format(1X,8F18.8)
700  format(1X,12F18.8)
return
end subroutine

!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!                Calculates Emergent Field                  !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine top_field
include "common.h"
integer i,j
real*8 theta1
do i=1,side
   do j=1,side
      theta1=theta(i,pos(j))-theta(i,j)!    via our emergent field definition: E_x(x,y)=+(theta(x,y+a)-theta(x,y+a))
      if (theta1.gt.pi) then!                                                  E_y(x,y)=-(theta(x+a,y)-theta(x,y+a))
         theta1=theta1-twopi
      else if (theta1.le.-pi) then
         theta1=theta1+twopi
      end if
      top_x(i,j)=theta1!/twopi
      theta1=theta(i,j)-theta(pos(i),j)
      if (theta1.gt.pi) then
         theta1=theta1-twopi
      else if (theta1.le.-pi) then
         theta1=theta1+twopi
      end if
      top_y(i,j)=theta1!/twopi
   end do
end do
return
end subroutine


!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!          Calculates the Total Emergent Field              !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine sum_field
include "common.h"
integer i,j
Ehatsum_x=0
Ehatsum_y=0
do i=1,side
   do j=1,side
      Ehatsum_x=Ehatsum_x+top_x(i,j)
      Ehatsum_y=Ehatsum_y+top_y(i,j)
   end do
end do
dx=Ehatsum_x
nx=0
dy=Ehatsum_y
ny=0

if (Ehatsum_x.gt.pi*float(side)) then! this finds the Ebar modulo (effectively) 2pi - we don't consider windings whose absolute value is greater than one
   dx=dx-twopi*side!                                 since the harmonic part of the Hamiltonian wants to remove the ±1 contributions anyway
   nx=Ehatsum_x-dx
else if (Ehatsum_x.le.-pi*float(side)) then! we instil a (-pi, pi] restriciton
   dx=dx+twopi*side
   nx=Ehatsum_x-dx
end if

if (Ehatsum_y.gt.pi*float(side)) then
   dy=dy-twopi*side
   ny=Ehatsum_y-dy
else if (Ehatsum_y.le.-pi*float(side)) then
   dy=dy+twopi*side
   ny=Ehatsum_y-dy
end if
Ebar_x=Ehatsum_x/volume
Ebar_y=Ehatsum_y/volume
Ebardip_x=dx/volume
Ebardip_y=dy/volume
Ebarwind_x=nx/volume
Ebarwind_y=ny/volume
return
end subroutine sum_field

!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!                    Records Vortices                       !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine vortices
include "common.h"
integer i,j
real*8 theta1,theta2,theta3,theta4,v0
do i=1,side
   do j=1,side
      theta1=top_x(pos(i),j)!        MUST BE CALLED AFTER top_field
      theta2=top_y(i,pos(j))!        HAVE CHECKED EQUIVALENCE WITH BASIC VORTEX MEASURING CODE
      theta3=-top_x(i,j)
      theta4=-top_y(i,j)
      v0=theta1+theta2+theta3+theta4
      v0=v0/twopi
      if ((v0.lt.1.01).and.(v0.gt.0.99)) then
         v0=1
      else if ((v0.gt.-1.01).and.(v0.lt.-0.99)) then
         v0=-1
      end if
      v(i,j)=v0
   end do
end do
return
end subroutine

!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!                    Records Vortices                       !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

subroutine charges
include "common.h"
integer i,j
real*8 theta1,theta2,theta3,theta4,v0
do i=1,side
   do j=1,side
      theta1=top_x(pos(i),j)!        MUST BE CALLED AFTER top_field
      theta2=top_y(i,pos(j))!        HAVE CHECKED EQUIVALENCE WITH BASIC VORTEX MEASURING CODE
      theta3=-top_x(i,j)
      theta4=-top_y(i,j)
      v0=theta1+theta2+theta3+theta4
      if ((v0.lt.1.01).and.(v0.gt.0.99)) then
         v0=1
      else if ((v0.gt.-1.01).and.(v0.lt.-0.99)) then
         v0=-1
      end if
      v(i,j)=v0
   end do
end do
return
end subroutine



!----------------------------------------------------------------------C
!                                                                      C
!  Lagged Fibonacci random number generator RANMAR.                    C
!  Must be initialized with randinit() before use.                     C
!                                                                      C
!  See F. James, Comp. Phys. Comm. 60, 329 (1990), or                  C
!  G. Marsaglia et al., Stat. Prob. Lett. 9, 35 (1990).                C
!                                                                      C
!----------------------------------------------------------------------C


!----------------------------------------------------------------------C
!                                                                      C
! This is the initialization routine RMARIN for the random number      C
!     generator RANMAR                                                 C
!                                                                      C
! NOTE: The seed variables can have values between:  0 <= IJ <= 31328  C
!                                                    0 <= KL <= 30081  C
!----------------------------------------------------------------------C

      SUBROUTINE randinit(seed)
      IMPLICIT NONE
      INTEGER seed
      INTEGER ij,kl, i,j,k,l, ii,jj, m
      REAL*8 s,t
      INTEGER Maxseed
      PARAMETER (Maxseed = 900000000)
      REAL u(97), c, cd, cm
      INTEGER i97, j97, ivec
      COMMON /raset1/ u, c, cd, cm, i97, j97, ivec 
 
      seed = mod(seed,Maxseed)
      ij = seed / 30082
      kl = seed - (30082 * ij)
      i = mod(ij/177, 177) + 2
      j = mod(ij    , 177) + 2
      k = mod(kl/169, 178) + 1
      l = mod(kl,     169)
      DO 2 ii = 1, 97
        s = 0.0
        t = 0.5
        DO 3 jj = 1, 24
          m = mod(mod(i*j, 179)*k, 179)
          i = j
          j = k
          k = m
          l = mod(53*l+1, 169)
          IF (mod(l*m, 64) .ge. 32) then
            s = s + t
          ENDIF
          t = 0.5 * t
    3   CONTINUE
        u(ii) = s
    2 CONTINUE
      c = 362436.0 / 16777216.0
      cd = 7654321.0 / 16777216.0
      cm = 16777213.0 /16777216.0
      i97 = 97
      j97 = 33
      RETURN
      END



!----------------------------------------------------------------------C
!                                                                      C
!  Lagged Fibonacci random number generator RANMAR().                  C
!                                                                      C
!----------------------------------------------------------------------C

      FUNCTION rand()
      IMPLICIT NONE
      REAL u(97), c, cd, cm, uni, rand
      INTEGER i97, j97, ivec
      COMMON /raset1/ u, c, cd, cm, i97, j97, ivec 

      uni = u(i97) - u(j97)
      IF( uni .LT. 0.0 ) uni = uni + 1.0
      u(i97) = uni
      i97 = i97 - 1
      IF(i97 .EQ. 0) i97 = 97
      j97 = j97 - 1
      IF(j97 .EQ. 0) j97 = 97
      c = c - cd
      IF( c .LT. 0.0 ) c = c + cm
      uni = uni - c
      IF( uni .LT. 0.0 ) uni = uni + 1.0
      rand = uni
      RETURN
      END

