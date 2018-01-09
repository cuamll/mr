module fft_mod
  use common
  implicit none
contains
 
  ! In place Cooley-Tukey FFT
  recursive subroutine fft_1d(input,output,mult)
    complex(kind=rk), dimension(:), intent(in)  :: input
    complex(kind=rk), dimension(:), intent(out)  :: output
    complex(kind=rk)                               :: t
    integer                                        :: N, i, mult
    complex(kind=rk), dimension(:), allocatable    :: even, odd
 
    N=size(x)
 
    if(N .le. 1) return
 
    allocate(odd((N+1)/2))
    allocate(even(N/2))
 
    ! divide
    odd =input(1:N:2)
    even=input(2:N:2)
 
    ! conquer
    call fft(odd)
    call fft(even)
 
    ! combine
    do i=1,mult*N/2
       t=exp(cmplx(0.0_dp,-2.0_dp*pi*real(i-1,rk)/real(N,rk),kind=rk))*even(i)
       output(i)     = odd(i) + t
       output(i+N/2) = odd(i) - t
    end do
 
    deallocate(odd)
    deallocate(even)
 
  end subroutine fft
 
  subroutine fft_2d(input,output,n1,n2)
    use common
    implicit none



  end subroutine fft_2d

end module fft_mod
