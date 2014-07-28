program main
    use, intrinsic :: iso_c_binding
    
    implicit none
    
    include 'fftw3.f03'
    
    integer, parameter :: size = 20
    double precision, parameter :: pi = 4d0*atan(1d0)
    type(C_PTR) :: fft_plan, ifft_plan
    complex(C_DOUBLE_COMPLEX), dimension(0:size-1) :: in, out
    double precision, dimension(0:size-1) :: in_real, out_real
    integer :: i
    
    in(:) = 0
    out(:) = 0
    in_real(:) = 0
    out_real(:) = 0
    
    do i=0,size-1
        !in(i) = sin(2*i*pi/size) 
        in_real(i) = sin(2*i*pi/size) 
        !in(i) = sin(i*pi/size)**2
        !in(i) = sin(2*i*pi/size) + sin(4*i*pi/size)
    end do
    in(:) = in_real(:)
    
    do i=0,size-1
        write(*,*) in(i)
    end do
    write(*,*)
    
    !Create the plans (one for forward transform, one for backward one - switches sign in exponent)
    fft_plan = fftw_plan_dft_1d(size, in,out, FFTW_FORWARD,FFTW_ESTIMATE)
    ifft_plan = fftw_plan_dft_1d(size, in,out, FFTW_BACKWARD,FFTW_ESTIMATE)
    
    call fftw_execute_dft(fft_plan, in, out)
    do i=0,size-1
        write(*,*) out(i)
    end do
    write(*,*)
    call fftw_execute_dft(ifft_plan, out, in)
    in_real(:) = realpart(in)
    do i=0,size-1
        !write(*,*) in(i)/size
        write(*,*) in_real(i)/size
    end do
    
    !Destroy the plans since we're done
    call fftw_destroy_plan(fft_plan)
    call fftw_destroy_plan(ifft_plan)
    
end program