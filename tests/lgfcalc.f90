program lgfcalc

  integer :: a,b,c,x,y,z,kx,ky,kz,L
  real*8 :: p1,q1,r1,fkx,fky,fkz
  real*8,dimension(L,L,L,L,L,L) :: lgf,lgfread
  character, len(11) :: lgf_filename
  logical :: there

  L = 10
  write(lgf_filename,"(A4,I0,A4)") "lgf_",L,".dat"
  write(*,*) "lgf filename is ",lgf_filename

  ! calculate the fucker

  do a=1,L
    do b=1,L
      do c=1,L
        ! Basically, we can think of the {a,b,c} as the sum over \vec{xprime} in
        ! G(\vec{x},vec{xprime}, and the {x,y,z} as the \vec{x}.
        ! {kx,ky,kz} are then the sum over k-space in G(\vec{x},\vec{xprime})

        do x=1,L
          do y=1,L
            do z=1,L

              lgf(a,b,c,x,y,z)=0.0

              ! these need to be real numbers for cos etc to work later
              p1=float(a-x)
              q1=float(b-y)
              r1=float(c-z)

              ! these need to be within [L/2,-L/2]. The Green's function is
              ! even though because we use cosines, so it doesn't matter
              ! whether it's positive or negative.

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

  ! check if file's there, if so read in and compare
  inquire(file=lgf_filename, exist=there)
!  if ( there ) then
!    read(lgf_filename,*)
!  else
!    WRITE
!  end if

end program lgfcalc
