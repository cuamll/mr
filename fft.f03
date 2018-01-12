module fft_mod
  use common
  implicit none
contains
 
  ! In place Cooley-Tukey FFT
  recursive subroutine fft(x,comp)
    complex(kind=rk), dimension(:), intent(inout)  :: x
    complex(kind=rk)                               :: t
    integer                                        :: N, i, comp
    complex(kind=rk), dimension(:), allocatable    :: even, odd
 
    ! need to figure out how to stop this being in-place???

    ! comp tells us which component we're doing: different exponents needed
    ! comp = 1 goes to hw: so rows in x-dir, cols in y-dir
    ! comp = 2 goes to fw: cols in x-dir, rows in y-dir

    N=size(x)
 
    if(N .le. 1) return
 
    allocate(odd((N+1)/2))
    allocate(even(N/2))
 
    ! divide
    odd =x(1:N:2)
    even=x(2:N:2)
 
    ! conquer
    call fft(odd,comp)
    call fft(even,comp)
 
    ! combine
    do i=1,N/2

      if (comp.eq.1) then
        t=exp(cmplx(0.0_rk,-2.0_rk*pi*&
              real(i,rk)/real(N,rk),kind=rk))*even(i)
      else if (comp.eq.2) then
        t=exp(cmplx(0.0_rk,-2.0_rk*pi*&
              real(neg(i)+1.0/2,rk)/real(N,rk),kind=rk))*even(i)
      else
        write (*,*) "FFT problem!!! comp = ",0
      end if

      x(i)     = odd(i) + t
      x(i+N/2) = odd(i) - t

    end do
 
    deallocate(odd)
    deallocate(even)
 
  end subroutine fft
 
  subroutine fft_2d_charges(x,nx,ny)
    use common
    implicit none
    complex(kind=rk), dimension(:,:), intent(inout) :: x
    complex(kind=rk), dimension(:), allocatable     :: strip
    integer, intent(in)                             :: nx,ny
    integer                                         :: i,j

    ! should (should!) perform an in-place FFT on fields:
    ! comp = 2 for both directions here because that corresponds
    ! to the actual lattice point without any offset

    ! do the rows
    allocate(strip(nx))
    strip = (0.0_rk,0.0_rk)
    ! !$omp parallel do&
    ! !$omp& private(i,j)&
    ! !$omp& shared(x,nx,ny,comp,strip)
    do j=1,ny
      ! write (*,*) "CHARGES, ROWS", j

      do i=1,nx
        strip(i) = x(i,j)
      end do

      call fft(strip,1)

      do i=1,nx
        x(i,j) = strip(i)
      end do

    end do
    ! !$omp end parallel do
    deallocate(strip)

    ! do the columns
    allocate(strip(ny))
    strip = (0.0_rk,0.0_rk)
    ! !$omp parallel do&
    ! !$omp& private(i,j)&
    ! !$omp& shared(x,nx,ny,comp,strip)
    do i=1,nx
      ! write (*,*) "CHARGES, COLUMNS", i

      do j=1,ny
        strip(j) = x(i,j)
      end do

      call fft(strip,1)

      do j=1,ny
        x(i,j) = strip(j)
      end do

    end do
    ! !$omp end parallel do
    deallocate(strip)

  end subroutine fft_2d_charges

  subroutine fft_2d_fields(x,nx,ny,comp)
    use common
    implicit none
    complex(kind=rk), dimension(:,:), intent(inout) :: x
    complex(kind=rk), dimension(:), allocatable     :: strip
    integer, intent(in)                             :: nx,ny,comp
    integer                                         :: i,j

    ! should (should!) perform an in-place FFT on fields:
    ! it's awkward because of the different offsets in x and y

    ! do the rows
    allocate(strip(nx))
    ! !$omp parallel do&
    ! !$omp& private(i,j)&
    ! !$omp& shared(x,nx,ny,comp,strip)
    do j=1,ny

      do i=1,nx
        strip(i) = x(i,j)
      end do

      call fft(strip,mod(comp, 2) + 1)

      do i=1,nx
        x(i,j) = strip(i)
      end do

    end do
    ! !$omp end parallel do
    deallocate(strip)

    ! do the columns
    allocate(strip(ny))
    ! !$omp parallel do&
    ! !$omp& private(i,j)&
    ! !$omp& shared(x,nx,ny,comp,strip)
    do i=1,nx

      do j=1,ny
        strip(j) = x(i,j)
      end do

      call fft(strip,mod(comp + 1, 2) + 1)

      do j=1,ny
        x(i,j) = strip(j)
      end do

    end do
    ! !$omp end parallel do
    deallocate(strip)

  end subroutine fft_2d_fields

end module fft_mod
