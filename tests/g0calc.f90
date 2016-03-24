program g0calc
  implicit none
  integer :: i,nx,ny,nz,L(15)
  real*8 :: g_zero,kx,ky,kz,pi,a

  pi=3.141592653589793
  L (1:15) = (/ 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 50, 60, 100 /)

  do i=1,15
    g_zero=0.0
    a=1.0

    write(*,*) L(i)

    do nx=-(L(i)-1)/2,L(i)/2
      kx=2*pi*nx/L(i)
      do ny=-(L(i)-1)/2,L(i)/2
        ky=2*pi*ny/L(i)
        do nz=-(L(i)-1)/2,L(i)/2
          kz=2*pi*nz/L(i)
          if ((nx.eq.0).and.(ny.eq.0).and.(nz.eq.0)) then
            g_zero=g_zero ! dunno how to continue or break or whatever
          else
            g_zero=g_zero+1/(3-cos(kx*a)-cos(ky*a)-cos(kz*a))
          end if ! end of kx=ky=kz=0 block
        end do ! end kz loop
      end do ! end ky loop
    end do ! end kx loop

    g_zero=g_zero/(L(i)**3)
    write(*,*) g_zero
  end do ! end loop over L values

end program g0calc
