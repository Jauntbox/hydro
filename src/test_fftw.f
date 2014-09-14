program main
    use, intrinsic :: iso_c_binding
    
    implicit none
    
    include 'fftw3.f03'
    
    integer, parameter :: size = 20 !Number of spectral components
    integer, parameter :: size_spatial = 2*size + 1 !Number of spatial points to transform to
    integer, parameter :: size_dft = (size_spatial) !Number of points to use in the DFT (in order to turn 2*pi*k arguments in the DFT into pi*k arguements)
    integer, parameter :: size_halfperiod = 2*size
    integer, parameter :: size_dft_halfperiod = 2*size_halfperiod + 1
    integer, parameter :: size_dft_halfperiod_dealias = 3*size_halfperiod + 1
    double precision, parameter :: dx = 1d0/(size_dft)
    double precision, parameter :: pi = 4d0*atan(1d0)
    
    logical, parameter :: use_aligned = .true.
    type(C_PTR) :: fft_plan, ifft_plan, dct_plan, idct_plan, r2c_plan, c2r_plan, &
       dct_plan_halfperiod, idct_plan_halfperiod, dct_plan_halfperiod_dealias, idct_plan_halfperiod_dealias
    double precision, dimension(0:size_dft-1) :: in_r2c
    complex(C_DOUBLE_COMPLEX), dimension(0:size_dft/2) :: out_r2c
    complex(C_DOUBLE_COMPLEX), dimension(0:size-1) :: in, out, product_spectral_complex, &
       in_sin, out_sin, in_cos, out_cos
    complex(C_DOUBLE_COMPLEX), dimension(0:size_dft-1) :: in_spatial, out_spatial
    complex(C_DOUBLE_COMPLEX), pointer :: in_aligned(:), out_aligned(:)
    complex(C_DOUBLE_COMPLEX), pointer :: in_aligned_halfperiod(:), out_aligned_halfperiod(:)
    complex(C_DOUBLE_COMPLEX), pointer :: in_aligned_halfperiod_dealias(:), out_aligned_halfperiod_dealias(:)
    type(C_PTR) :: p_in, p_out, p_in_halfperiod, p_out_halfperiod, &
         p_in_halfperiod_dealias, p_out_halfperiod_dealias
    
    double precision, dimension(0:size-1) :: in_real, out_real
    double precision, dimension(0:size-1) :: a_spectral, b_spectral, product_spectral
    double precision, dimension(0:2*(size-1)) :: product_spectral_full
    double precision, dimension(0:size_dft-1) :: a_spatial, b_spatial, product_spatial
    double precision, dimension(0:size_dft_halfperiod-1) :: a_spatial_halfperiod, &
       b_spatial_halfperiod, product_spatial_halfperiod, product_spectral_halfperiod
    double precision, dimension(0:size_dft_halfperiod_dealias-1) :: a_spatial_halfperiod_dealias, &
       b_spatial_halfperiod_dealias, product_spatial_halfperiod_dealias, &
       product_spectral_halfperiod_dealias
    double precision :: a0, a1, a6, b1, b2
    integer :: i,n
    
    in(:) = 0
    out(:) = 0
    in_real(:) = 0
    out_real(:) = 0
    in_sin(:) = 0
    out_sin(:) = 0
    in_cos(:) = 0
    out_cos(:) = 0
    
    if(use_aligned) then
        p_in = fftw_alloc_complex(int(size_dft, C_SIZE_T))
        call c_f_pointer(p_in, in_aligned, [size_dft])
        p_out = fftw_alloc_complex(int(size_dft, C_SIZE_T))
        call c_f_pointer(p_out, out_aligned, [size_dft])
        
        p_in_halfperiod = fftw_alloc_complex(int(size_dft_halfperiod, C_SIZE_T))
        call c_f_pointer(p_in_halfperiod, in_aligned_halfperiod, [size_dft_halfperiod])
        p_out_halfperiod = fftw_alloc_complex(int(size_dft_halfperiod, C_SIZE_T))
        call c_f_pointer(p_out_halfperiod, out_aligned_halfperiod, [size_dft_halfperiod])
        
        p_in_halfperiod_dealias = fftw_alloc_complex(int(size_dft_halfperiod_dealias, C_SIZE_T))
        call c_f_pointer(p_in_halfperiod_dealias, in_aligned_halfperiod_dealias, [size_dft_halfperiod_dealias])
        p_out_halfperiod_dealias = fftw_alloc_complex(int(size_dft_halfperiod_dealias, C_SIZE_T))
        call c_f_pointer(p_out_halfperiod_dealias, out_aligned_halfperiod_dealias, [size_dft_halfperiod_dealias])
    endif
    
    do i=0,size-1
        !in(i) = sin(2*i*pi/size) 
        !in_real(i) = sin(2*i*pi/size) 
        !in(i) = sin(i*pi/size)**2
        in_sin(i) = 5d-1*sin(2*i*pi/size) + 3d0*sin(4*i*pi/size)
        in_cos(i) = 1d0 + 2d0*cos(2*i*pi/size) + cos(12*i*pi/size)
    end do
    !in(:) = in_real(:)
    
    !do i=0,size-1
    !    write(*,*) in(i)
    !end do
    !write(*,*)
    
    !Create the plans (one for forward transform, one for backward one - switches sign in exponent)
    r2c_plan = fftw_plan_dft_r2c_1d(size_dft,in_r2c,out_r2c,FFTW_ESTIMATE)
    c2r_plan = fftw_plan_dft_c2r_1d(size_dft,out_r2c,in_r2c,FFTW_ESTIMATE)
    if(use_aligned) then
        !dct_plan = fftw_plan_dft_1d(size_dft, in_spatial, out_spatial, FFTW_REDFT10, FFTW_ESTIMATE)
        dct_plan = fftw_plan_dft_1d(size_dft, in_aligned, out_aligned, FFTW_FORWARD, FFTW_ESTIMATE)
        !idct_plan = fftw_plan_dft_1d(size_dft, in_spatial, out_spatial, FFTW_REDFT01, FFTW_ESTIMATE)
        idct_plan = fftw_plan_dft_1d(size_dft, in_aligned, out_aligned, FFTW_BACKWARD, FFTW_ESTIMATE)
        
        dct_plan_halfperiod = fftw_plan_dft_1d(size_dft_halfperiod, in_aligned_halfperiod, &
            out_aligned_halfperiod, FFTW_FORWARD, FFTW_ESTIMATE)
        idct_plan_halfperiod = fftw_plan_dft_1d(size_dft_halfperiod, in_aligned_halfperiod, &
            out_aligned_halfperiod, FFTW_BACKWARD, FFTW_ESTIMATE)
        
        dct_plan_halfperiod_dealias = fftw_plan_dft_1d(size_dft_halfperiod_dealias, &
            in_aligned_halfperiod_dealias, out_aligned_halfperiod, FFTW_FORWARD, FFTW_ESTIMATE)
        idct_plan_halfperiod_dealias = fftw_plan_dft_1d(size_dft_halfperiod_dealias, &
            in_aligned_halfperiod_dealias, out_aligned_halfperiod, FFTW_BACKWARD, FFTW_ESTIMATE)
    else
        !dct_plan = fftw_plan_dft_1d(size_dft, in_spatial, out_spatial, FFTW_REDFT10, FFTW_ESTIMATE)
        dct_plan = fftw_plan_dft_1d(size_dft, in_spatial, out_spatial, FFTW_FORWARD, FFTW_ESTIMATE)
        !idct_plan = fftw_plan_dft_1d(size_dft, in_spatial, out_spatial, FFTW_REDFT01, FFTW_ESTIMATE)
        idct_plan = fftw_plan_dft_1d(size_dft, in_spatial, out_spatial, FFTW_BACKWARD, FFTW_ESTIMATE)
    endif
    fft_plan = fftw_plan_dft_1d(size, in,out, FFTW_FORWARD, FFTW_ESTIMATE)
    ifft_plan = fftw_plan_dft_1d(size, in,out, FFTW_BACKWARD, FFTW_ESTIMATE)
    
    !Try to replace this with some sort of r2c or c2r transform
    !call fftw_execute_dft(fft_plan, in, out)
    !do i=0,size-1
    !    write(*,*) out(i)
    !end do
    !write(*,*)
    !call fftw_execute_dft(ifft_plan, out, in)
    !in_real(:) = realpart(in)
    !do i=0,size-1
        !write(*,*) in(i)/size
    !    write(*,*) in_real(i)/size
    !end do
    !write(*,*)
    
    !Next thing to try is extracting the fourier coefficients from a cosine series:
    a0 = 1d0
    a1 = 2d0
    a6 = 1d0
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
    if(use_aligned) then
        in_aligned = 0d0
        !Convert the cosine spectral components to complex exponential components:
        !Note: the aligned memory arrays are indexed from 1:N (ugh...)
        in_aligned(1) = a_spectral(0)
        do i=2,size
           in_aligned(i) = a_spectral(i-1)/2d0
           in_aligned(size_dft-(i-2)) = a_spectral(i-1)/2d0
        end do
        !in_aligned(0:size-1) = a_spectral
        call fftw_execute_dft(idct_plan, in_aligned, out_aligned)
        a_spatial = out_aligned
    else
        in_spatial = 0d0
        in_spatial(0:size-1) = a_spectral
        call fftw_execute_dft(idct_plan, in_spatial, out_spatial)
        a_spatial = out_spatial
    endif
    !out_r2c = 0d0
    !out_r2c(0:size) = a_spectral
    !call fftw_execute_dft_c2r(c2r_plan, out_r2c, in_r2c)
    !a_spatial = in_r2c
    write(*,*) 'In space, we have (iFFT of a):'
    !write(*,*) out_aligned
    write(*,*) a_spatial
    write(*,*)
    
    !Now, try converting it back to spectral coordinates with the FFT:
    if(use_aligned) then
        in_aligned = a_spatial(:)
        call fftw_execute_dft(dct_plan, in_aligned, out_aligned)
        !a_spectral = realpart(out_aligned(0:size-1))/(size_dft)
        a_spectral = out_aligned(0:size-1)/size_dft
        write(*,*) 'Transformed back to spectral coordinates, we have (a):'
        write(*,*) realpart(out_aligned)/size_dft
        write(*,*)
        write(*,*) imagpart(out_aligned)/size_dft
        write(*,*)
    else
        in_spatial = a_spatial(:)
        call fftw_execute_dft(dct_plan, in_spatial, out_spatial)
        a_spectral = realpart(out_spatial(0:size-1))/(size_dft)
        write(*,*) 'Transformed back to spectral coordinates, we have (a):'
        write(*,*) realpart(out_spatial)/size_dft
        write(*,*)
        write(*,*) imagpart(out_spatial)/size_dft
        write(*,*)
    endif
    out_r2c = 0d0
    out_r2c(0:size-1) = a_spectral
    call fftw_execute_dft_c2r(c2r_plan, out_r2c, in_r2c)
    a_spatial = in_r2c
    
    
    !Next thing to try is to reconstruct a product in real space from a spectral transform
    write(*,*) 'Product conversion:'
    a0 = 1d0
    a1 = 2d0
    a6 = 1d0
    b1 = 5d-1
    b2 = 3d0
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
        do n=0,size-1
            a_spatial(i) = a_spatial(i) + a_spectral(n)*cos(2*n*pi*i*dx)
            b_spatial(i) = b_spatial(i) + b_spectral(n)*cos(2*n*pi*i*dx)
        end do
    end do
    
    !Multiply them together in spatial coordinates:
    write(*,*) 'In space, we have (a*b):'
    write(*,*) a_spatial(:)*b_spatial(:)
    write(*,*)
    
    !Now convert to spectral coordinates:
    in_spatial = a_spatial(:)*b_spatial(:)
    call fftw_execute_dft(dct_plan, in_spatial, out_spatial)
    product_spectral_complex = out_spatial(0:size-1)/(size_dft)
    product_spectral = realpart(out_spatial(0:size-1))/(size_dft)
    product_spectral(1:size-1) = 2*product_spectral(1:size-1)   !not sure why I have to do this yet...
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
        do n=0,size-1
            product_spatial(i) = product_spatial(i) + product_spectral(n)*cos(2*n*pi*i*dx)
        end do
    end do
    write(*,*) 'Back in space, we have (a*b):'
    write(*,*) product_spatial
    write(*,*)
    
    !Check the product obtained by hand
    product_spatial(:) = 0d0
    do i=0,size_dft-1
        do n=0,size-1
            product_spatial(i) = product_spatial(i) + product_spectral_full(n)*cos(2*n*pi*i*dx)
        end do
    end do
    write(*,*) 'Back in space (with product by hand), we have (a*b):'
    write(*,*) product_spatial
    write(*,*)
    
    !Finally, compare to the spectral transform method:
    write(*,*) 'Now, compare to the product we get from a spectral transform.'
    write(*,*) 'The spectral components of a and b are:'
    write(*,*) 'a_spectral:'
    write(*,*) a_spectral
    write(*,*) 'b_spectral:'
    write(*,*) b_spectral
    write(*,*)
    
    write(*,*) 'Transform each into space, multiply together, and transform back.'
    write(*,*) 'The spectral components of the product are:'
    if(use_aligned) then
       in_aligned = 0d0
       in_aligned(1) = a_spectral(0)
       do i=2,size
          in_aligned(i) = a_spectral(i-1)/2d0
          in_aligned(size_dft-(i-2)) = a_spectral(i-1)/2d0
       end do
       call fftw_execute_dft(idct_plan, in_aligned, out_aligned)
       a_spatial = out_aligned
       in_aligned = 0d0
       in_aligned(1) = b_spectral(0)
       do i=2,size
          in_aligned(i) = b_spectral(i-1)/2d0
          in_aligned(size_dft-(i-2)) = b_spectral(i-1)/2d0
       end do
       call fftw_execute_dft(idct_plan, in_aligned, out_aligned)
       b_spatial = out_aligned
       in_aligned = a_spatial(:)*b_spatial(:)
       call fftw_execute_dft(dct_plan, in_aligned, out_aligned)
       product_spectral = out_aligned/size_dft
       !Convert from complex exponential coeffs to pure cosine coeffs:
       product_spectral(1:size-1) = 2*product_spectral(1:size-1)
    endif
    !write(*,*) a_spatial
    !write(*,*)
    !write(*,*) b_spatial
    !write(*,*)
    write(*,*) product_spectral
    write(*,*)
    
    write(*,*) 'Now, if we want to use n*pi*i*dx as our argument instead of 2*n*pi*i*dx, we'
    write(*,*) 'just double the DFT size:'
    !a_spectral_halfperiod = 0d0
    !a_spectral_halfperiod(0:size-1) = a_spectral(0:size-1)
    !b_spectral_halfperiod = 0d0
    !b_spectral_halfperiod(0:size-1) = b_spectral(0:size-1)
    a_spatial_halfperiod = 0d0
    b_spatial_halfperiod = 0d0
    do i=0,size_dft_halfperiod-1
        do n=0,size-1
            a_spatial_halfperiod(i) = a_spatial_halfperiod(i) + a_spectral(n)*cos(n*pi*i*dx)
            b_spatial_halfperiod(i) = b_spatial_halfperiod(i) + b_spectral(n)*cos(n*pi*i*dx)
        end do
    end do
    write(*,*) a_spatial_halfperiod
    write(*,*)
    write(*,*) b_spatial_halfperiod
    write(*,*)
    if(use_aligned) then
       in_aligned_halfperiod = 0d0
       in_aligned_halfperiod(1) = a_spectral(0)
       do i=2,size
          in_aligned_halfperiod(i) = a_spectral(i-1)/2d0
          in_aligned_halfperiod(size_dft_halfperiod-(i-2)) = a_spectral(i-1)/2d0
       end do
       call fftw_execute_dft(idct_plan_halfperiod, in_aligned_halfperiod, out_aligned_halfperiod)
       a_spatial_halfperiod = out_aligned_halfperiod
       in_aligned_halfperiod = 0d0
       do i=2,size
          in_aligned_halfperiod(i) = b_spectral(i-1)/2d0
          in_aligned_halfperiod(size_dft_halfperiod-(i-2)) = b_spectral(i-1)/2d0
       end do
       call fftw_execute_dft(idct_plan_halfperiod, in_aligned_halfperiod, out_aligned_halfperiod)
       b_spatial_halfperiod = out_aligned_halfperiod
       in_aligned_halfperiod = a_spatial_halfperiod(:)*b_spatial_halfperiod(:)
       write(*,*) a_spatial_halfperiod
       write(*,*)
       write(*,*) b_spatial_halfperiod
       write(*,*)
       write(*,*) in_aligned_halfperiod
       write(*,*)
       call fftw_execute_dft(dct_plan_halfperiod, in_aligned_halfperiod, out_aligned_halfperiod)
       product_spectral_halfperiod = out_aligned_halfperiod/size_dft_halfperiod
       !Convert from complex exponential coeffs to pure cosine coeffs:
       product_spectral_halfperiod(1:size-1) = 2*product_spectral_halfperiod(1:size-1)
    endif
    write(*,*) product_spectral_halfperiod
    write(*,*)
    
    write(*,*) 'Finally, if we want to throw out aliased terms so as not to make up information'
    write(*,*) 'on resolutions smaller than our grid size, then we up the spatial_size to'
    write(*,*) '(3*size + 1) instead of the Nyquist limit (2*size + 1):'
    a_spatial_halfperiod_dealias = 0d0
    b_spatial_halfperiod_dealias = 0d0
    do i=0,size_dft_halfperiod_dealias-1
        do n=0,size-1
            a_spatial_halfperiod_dealias(i) = a_spatial_halfperiod_dealias(i) + a_spectral(n)*cos(n*pi*i*dx)
            b_spatial_halfperiod_dealias(i) = b_spatial_halfperiod_dealias(i) + b_spectral(n)*cos(n*pi*i*dx)
        end do
    end do
    if(use_aligned) then
       in_aligned_halfperiod_dealias = 0d0
       in_aligned_halfperiod_dealias(1) = a_spectral(0)
       do i=2,size
          in_aligned_halfperiod_dealias(i) = a_spectral(i-1)/2d0
          in_aligned_halfperiod_dealias(size_dft_halfperiod_dealias-(i-2)) = a_spectral(i-1)/2d0
       end do
       call fftw_execute_dft(idct_plan_halfperiod_dealias, in_aligned_halfperiod_dealias, out_aligned_halfperiod_dealias)
       a_spatial_halfperiod_dealias = out_aligned_halfperiod_dealias
       in_aligned_halfperiod_dealias = 0d0
       do i=2,size
          in_aligned_halfperiod_dealias(i) = b_spectral(i-1)/2d0
          in_aligned_halfperiod_dealias(size_dft_halfperiod_dealias-(i-2)) = b_spectral(i-1)/2d0
       end do
       call fftw_execute_dft(idct_plan_halfperiod_dealias, in_aligned_halfperiod_dealias, out_aligned_halfperiod_dealias)
       b_spatial_halfperiod_dealias = out_aligned_halfperiod_dealias
       in_aligned_halfperiod_dealias = a_spatial_halfperiod_dealias(:)*b_spatial_halfperiod_dealias(:)
       call fftw_execute_dft(dct_plan_halfperiod_dealias, in_aligned_halfperiod_dealias, out_aligned_halfperiod_dealias)
       product_spectral_halfperiod_dealias = out_aligned_halfperiod_dealias/size_dft_halfperiod_dealias
       !Convert from complex exponential coeffs to pure cosine coeffs:
       product_spectral_halfperiod_dealias(1:size-1) = 2*product_spectral_halfperiod_dealias(1:size-1)
    endif
    write(*,*) product_spectral_halfperiod_dealias
    write(*,*)
    
    !Destroy the plans since we're done
    call fftw_destroy_plan(dct_plan)
    call fftw_destroy_plan(idct_plan)
    call fftw_destroy_plan(fft_plan)
    call fftw_destroy_plan(ifft_plan)
    
    !Deallocate arrays as well:
    call fftw_free(p_in)
    call fftw_free(p_out)
    call fftw_free(p_in_halfperiod)
    call fftw_free(p_out_halfperiod)
    call fftw_free(p_in_halfperiod_dealias)
    call fftw_free(p_out_halfperiod_dealias)
    
end program