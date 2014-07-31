program main
    use, intrinsic :: iso_c_binding
    
    implicit none
    
    include 'fftw3.f03'
    
    integer, parameter :: size = 20
    integer, parameter :: size_spatial = 4*size + 1 !Number of spatial points to transform to
    integer, parameter :: size_dft = 2*(size_spatial) !Number of points to use in the DFT (in order to turn 2*pi*k arguments in the DFT into pi*k arguements)
    double precision, parameter :: dx = 1d0/(size_spatial)
    double precision, parameter :: pi = 4d0*atan(1d0)
    type(C_PTR) :: fft_plan, ifft_plan, dct_plan, idct_plan
    complex(C_DOUBLE_COMPLEX), dimension(0:size) :: in, out, product_spectral_complex
    complex(C_DOUBLE_COMPLEX), dimension(0:size_dft-1) :: in_spatial, out_spatial
    double precision, dimension(0:size) :: in_real, out_real
    double precision, dimension(0:size) :: a_spectral, b_spectral, product_spectral
    double precision, dimension(0:2*size) :: product_spectral_full
    double precision, dimension(0:size_dft-1) :: a_spatial, b_spatial, product_spatial
    double precision :: a0, a1, a6, b1, b2
    integer :: i,n
    
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
    !dct_plan = fftw_plan_dft_1d(size_dft, in_spatial, out_spatial, FFTW_REDFT10, FFTW_ESTIMATE)
    dct_plan = fftw_plan_dft_1d(size_dft, in_spatial, out_spatial, FFTW_FORWARD, FFTW_ESTIMATE)
    !idct_plan = fftw_plan_dft_1d(size_dft, in_spatial, out_spatial, FFTW_REDFT01, FFTW_ESTIMATE)
    idct_plan = fftw_plan_dft_1d(size_dft, in_spatial, out_spatial, FFTW_BACKWARD, FFTW_ESTIMATE)
    fft_plan = fftw_plan_dft_1d(size, in,out, FFTW_FORWARD, FFTW_ESTIMATE)
    ifft_plan = fftw_plan_dft_1d(size, in,out, FFTW_BACKWARD, FFTW_ESTIMATE)
    
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
    write(*,*)
    
    !Next thing to try is extracting the fourier coefficients from a cosine series:
    a0 = 1d0
    a1 = 1d-2
    a6 = 1d-3
    a_spectral(:) = 0d0
    a_spectral(0) = a0
    a_spectral(1) = a1
    a_spectral(6) = a6
    a_spatial(:) = 0d0
    !Convert the spectral components to spatial coordinates:
    do i=0,size_dft-1
        do n=0,size
            a_spatial(i) = a_spatial(i) + a_spectral(n)*cos(2*n*pi*i*dx)
        end do
    end do
    write(*,*) 'In space, we have (direct sum of a components):'
    write(*,*) a_spatial
    write(*,*)
    
    !Do the same with the inverse FFT (frequency -> spatial):
    in_spatial = 0d0
    in_spatial(0:size) = a_spectral
    call fftw_execute_dft(idct_plan, in_spatial, out_spatial)
    a_spatial = out_spatial
    write(*,*) 'In space, we have (iFFT of a):'
    write(*,*) a_spatial
    write(*,*)
    
    !Now, try converting it back to spectral coordinates with the FFT:
    in_spatial = a_spatial(:)
    call fftw_execute_dft(dct_plan, in_spatial, out_spatial)
    a_spectral = realpart(out_spatial(0:size))/(size_dft)
    write(*,*) 'Transformed back to spectral coordinates, we have (a):'
    write(*,*) realpart(out_spatial)/size_dft
    write(*,*)
    write(*,*) imagpart(out_spatial)/size_dft
    write(*,*)
    
    
    !Next thing to try is to reconstruct a product in real space from a spectral transform
    write(*,*) 'Product conversion:'
    a0 = 1d0
    a1 = 1d-2
    a6 = 1d-3
    b1 = 1d-4
    b2 = 3d-4
    a_spectral(:) = 0d0
    a_spectral(0) = a0
    a_spectral(1) = a1
    a_spectral(6) = a6
    b_spectral(:) = 0d0
    b_spectral(1) = b1
    b_spectral(2) = b2
    a_spatial(:) = 0d0
    b_spatial(:) = 0d0
    !Convert the spectral components to spatial coordinates:
    do i=0,size_dft-1
        do n=0,size
            a_spatial(i) = a_spatial(i) + a_spectral(n)*cos(n*pi*i*dx)
            b_spatial(i) = b_spatial(i) + b_spectral(n)*cos(n*pi*i*dx)
        end do
    end do
    
    !Multiply them together in spatial coordinates:
    write(*,*) 'In space, we have (a*b):'
    write(*,*) a_spatial(:)*b_spatial(:)
    write(*,*)
    
    !Now convert to spectral coordinates:
    in_spatial = a_spatial(:)*b_spatial(:)
    call fftw_execute_dft(dct_plan, in_spatial, out_spatial)
    product_spectral_complex = out_spatial(0:size)/(size_dft)
    product_spectral = realpart(out_spatial(0:size))/(size_dft)
    product_spectral(1:size) = 2*product_spectral(1:size)   !not sure why I have to do this yet...
    write(*,*) 'In spectral coordinates, we have:'
    write(*,*) product_spectral
    write(*,*)
    product_spectral_full = 0d0
    product_spectral_full(0) = 5d-1*a1*b1
    product_spectral_full(1) = (a0*b1 + 5d-1*a1*b2)
    product_spectral_full(2) = (a0*b2 + 5d-1*a1*b1)
    product_spectral_full(3) = 5d-1*a1*b2
    product_spectral_full(4) = 5d-1*a6*b2
    product_spectral_full(5) = 5d-1*a6*b1
    product_spectral_full(7) = 5d-1*a6*b1
    product_spectral_full(8) = 5d-1*a6*b2
    write(*,*) 'By hand, the spectral components of the product should be:'
    write(*,*) product_spectral_full
    write(*,*)
    
    !Check by converting back to spatial coordinates:
    product_spatial(:) = 0d0
    do i=0,size_dft-1
        do n=0,size
            product_spatial(i) = product_spatial(i) + product_spectral(n)*cos(n*pi*i*dx)
        end do
    end do
    write(*,*) 'Back in space, we have (a*b):'
    write(*,*) product_spatial
    write(*,*)
    
    !Check the product obtained by hand
    product_spatial(:) = 0d0
    do i=0,size_dft-1
        do n=0,2*size
            product_spatial(i) = product_spatial(i) + product_spectral_full(n)*cos(n*pi*i*dx)
        end do
    end do
    write(*,*) 'Back in space (with product by hand), we have (a*b):'
    write(*,*) product_spatial
    
    !Destroy the plans since we're done
    call fftw_destroy_plan(dct_plan)
    call fftw_destroy_plan(idct_plan)
    call fftw_destroy_plan(fft_plan)
    call fftw_destroy_plan(ifft_plan)
    
end program