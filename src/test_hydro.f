!Code to solve the Boussinesq equations for Rayleigh-Bernard convection in 2D using
!spectral decomposition in the horizontal, and finite-difference in the vertical.
!-
!Written by: Kevin Moore
!6/14/14

module hydro
    
    use, intrinsic :: iso_c_binding

	implicit none
    
    include 'fftw3.f03'
	
	!Overall simulation parameters:
	!Only solve the linear problem (eg. equations )
	logical, parameter :: do_linear = .false.
    logical, parameter :: do_spectral_transform = .false.
	logical, parameter :: use_pgplot = .true.
	logical, parameter :: dbg = .false.		!If true, print debugging output to the screen
	
	!Domain size parameters
	integer, parameter :: nz_vertical = 101	    !Number of spatial vertical zones in simulation
	integer, parameter :: nn_horizontal = 50	!Number of Fourier modes in the horizontal
    !Number of spatial zones in horizontal (for spectral transform):
    integer, parameter :: nz_horizontal = 3*nn_horizontal + 1
    !Number of zones to use in the DFT in order to get odd modes (pi*n*x/L), rather than just the normal (2*pi*n*x/L) modes:
    integer, parameter :: n_dft = 2*nz_horizontal
	double precision, parameter :: dz = 1d0/(nz_vertical-1)	!spacing in the z-direction
    double precision, parameter :: dx = 1d0/(nz_horizontal)	!spacing in the x-direction
	integer, parameter :: a = 1.0	!horizontal/vertical length scale (integer for FFT'ing)
	
	!Timestep parameters:
	double precision :: dt	!Gets set in initialize()
    double precision :: dt_old  !For recalculating coefficients in Adams-Bashforth
	double precision, parameter :: fixed_dt = 3d-6	!If >0, then use this as the fixed dt value
	!Number of old time steps to keep track of when time-stepping:
	integer, parameter :: num_time_differences = 2
	integer, parameter :: num_steps = 1000000		!Number of time steps to take
	character (len=16) :: solver	!'RK2' or 'AB2' for now (not implemented yet)
    integer, parameter :: output_tick = 1000 !Print screen output every output_tick time steps
    integer, parameter :: cfl_tick = 10 !Check the CFL condition every cfl_tick time steps
    double precision :: avg_cfl_x, avg_cfl_z    !Avg CFL RHS values (to guess at increasing timestep)
    
	
	!Physical parameters:
	double precision, parameter :: rayleigh = 1d6	!Rayleigh number (dimensionless)
	double precision, parameter :: prandtl = 0.5	!Prandtl number (dimensionless)
	
	!Numerical constants:
	double precision, parameter :: pi = 4d0*atan(1d0)
	
	!Fluid variables:
	double precision, dimension(1:nz_vertical) :: z, sub, dia, sup, wk1, wk2
	double precision, dimension(1:nz_vertical, 0:nn_horizontal) :: psi, omega, temp, prev_temp
	!Nonlinear terms: [(v*del)variable]
	double precision, dimension(1:nz_vertical, 0:nn_horizontal) :: nonlinear_temp, nonlinear_omega
	!Arrays of first and second derivatives in the z-direction:
	double precision, dimension(1:nz_vertical, 0:nn_horizontal) :: dtemp_dz2, domega_dz2, &
		dtemp_dz, dpsi_dz, domega_dz
	double precision, dimension(1:nz_vertical, 0:nn_horizontal, 1:num_time_differences) :: &
		domega_dt, dtemp_dt
    double precision, dimension(1:nz_vertical, 1:nz_vertical) :: ux_cfl, uz_cfl
    !Purely spatial veriables for spectral transform method:
    double precision, dimension(1:nz_vertical, 0:n_dft-1) :: ux_spatial, uz_spatial, & 
        omega_spatial, temp_spatial, dtempdx_spatial, domegadx_spatial, &
        domegadz_spatial, dtempdz_spatial, galerkin_test
        
    !FFTW variables:
    !FFT, discrete cosine, and discrete sine transform plans, real to complex transforms:
    type(C_PTR) :: fft_plan, ifft_plan, dct_plan, dst_plan, r2c_plan, c2r_plan
    double precision, dimension(0:n_dft-1) :: in_r2c
    complex(C_DOUBLE_COMPLEX), dimension(0:n_dft/2) :: out_r2c
    complex(C_DOUBLE_COMPLEX), dimension(0:n_dft-1) :: in_fftw, out_fftw
    !complex(C_DOUBLE_COMPLEX), pointer :: in_fftw(:), out_fftw(:)
    type(C_PTR) :: p_in, p_out
		
	!PGPLOT variables:
	integer :: pgplot_id		!PGPLOT handle
	character (len=16) :: dev	!PGPLOT device name (eg. '/xwin')
	double precision, dimension(1:nz_vertical, 0:nz_vertical) :: temp_display, omega_display
    integer, parameter :: plot_tick = 20    !Update the plot window every plot_tick timesteps
	
	contains
	
	!Sets up all the variables, initial conditions, pgplot window, etc.
	subroutine initialize
	
		integer :: k,n	!Loop indicies
		integer, dimension(:), pointer :: itime
		integer :: pgopen, isize, c1, c2
		real :: x1, x2, y1, y2	!PGPLOT window bounds
      real :: rmin, rmax, gmin, gmax, bmin, bmax  !PGPLOT colors
      integer :: omp_get_max_threads
	
		write(*,*) 'Initializing hydro module!'
		write(*,*)
      
      !Check max number of threads (should be set by $OMP_NUM_THREADS)
      if(omp_get_max_threads().gt.1) then
         write(*,*) 'omp_get_num_threads():', omp_get_max_threads()
         write(*,*) 'Sorry, code is currently not working with multiple threads'
         write(*,*)
         stop 1
      endif
		
		!Make sure we satisfy the numerical diffusion limit:
		if(prandtl.lt.1d0) then
			dt = (dz**2)/4d0
		else
			dt = (dz**2)/(4d0*prandtl)
		endif
		
		dt = 0.8*dt
		
		if(fixed_dt > 0) dt = fixed_dt 
        dt_old = dt
		
		write(*,*) 'Domain size parameters:'
      write(*,'(a25,I25)') 'nn_horizontal:', nn_horizontal
      write(*,'(a25,I25)') 'nz_vertical:', nz_vertical
      !write(*,'(a25,I25)') 'a (aspect ratio):', a
      write(*,'(a25,ES25.10)') 'dz:', dz
      write(*,*)
      write(*,*) 'Timestep parameters:'
      write(*,'(a25,ES25.10)') 'dt:', dt
      write(*,'(a25,I25)') 'num_time_differences:', num_time_differences
      write(*,'(a25,I25)') 'num_steps:', num_steps
      write(*,*)
      write(*,*) 'Physical parameters:'
      write(*,'(a25,ES25.10)') 'rayleigh:', rayleigh
      write(*,'(a25,ES25.10)') 'prandtl:', prandtl
      write(*,*)
      write(*,*) 'Numerical constants:'
      write(*,'(a25,ES25.10)') 'pi:', pi
      write(*,*)
		
		!Initialize the random number generator
		call random_seed(size=isize)
		allocate(itime(isize))
		itime = time()
		call random_seed(put=itime)
		deallocate(itime)
		
		!Initialize the physical variables:
		omega(:,:) = 0d0
		psi(:,:) = 0d0
		!call random_number(omega)
		!omega = 10d0*omega - 5d0
		!call random_number(psi)
		!psi = 10d0*psi - 5d0
		temp(:,:) = 0d0
		prev_temp(:,:) = 0d0
		temp_display(:,:) = 0d0
		omega_display(:,:) = 0d0
		nonlinear_temp(:,:) = 0d0
		nonlinear_omega(:,:) = 0d0
      ux_spatial(:,:) = 0d0
      uz_spatial(:,:) = 0d0
      temp_spatial(:,:) = 0d0
      omega_spatial(:,:) = 0d0
		dtemp_dt(:,:,:) = 0d0
		domega_dt(:,:,:) = 0d0
		dtemp_dz(:,:) = 0d0 
		dpsi_dz(:,:) = 0d0
		domega_dz(:,:) = 0d0
		dtemp_dz2(:,:) = 0d0
		domega_dz2(:,:) = 0d0
		do k=1,nz_vertical
			z(k) = (k - 1)*dz
			!Pick a temperature that satisfies the boundary conditions:
			temp(k,0) = 1d0 - z(k)
			temp(k,1) = 0.01*sin(pi*z(k))
			temp(k,8) = 0.01*sin(pi*z(k))
			!temp(k,1:nn_horizontal) = sin(pi*z(k))	!Works for all spatial rows
		end do
		do k=1,nz_vertical
			do n=0,nn_horizontal
				temp_display(:,k) = temp_display(:,k) + temp(:,n)*cos(n*pi*dz*k)
			end do
		end do
		
		!Set the subdiagonal and superdiagonal of the Poisson solver matrix:
		sub(:) = -1d0/(dz**2)
		sup(:) = -1d0/(dz**2)
		sub(nz_vertical) = 0d0
		sup(1) = 0d0
		dia(:) = 1d0	!Only for the boundary conditions 2:nz_vertical-1 will be overwritten
        
      avg_cfl_x = 0d0
      avg_cfl_z = 0d0

      !Create the FFTW plans if we're using the spectral transform method:
      if(do_spectral_transform) then
         !Allocate aligned arrays to speed things up (doesn't crash, but isn't agreeing with normal way):
         !p_in = fftw_alloc_complex(int(n_dft, C_SIZE_T))
         !call c_f_pointer(p_in, in_fftw, [n_dft])
         !p_out = fftw_alloc_complex(int(n_dft, C_SIZE_T))
         !call c_f_pointer(p_out, out_fftw, [n_dft])
   
         !Old (C) way:
         !fft_plan = fftw_plan_dft_1d(n_dft, in_fftw, out_fftw, FFTW_FORWARD, FFTW_PATIENT)
         !ifft_plan = fftw_plan_dft_1d(n_dft, in_fftw, out_fftw, FFTW_BACKWARD, FFTW_PATIENT)
         !dct_plan = fftw_plan_dft_1d(n_dft, in_fftw, out_fftw, FFTW_REDFT10, FFTW_ESTIMATE)
         !dst_plan = fftw_plan_dft_1d(n_dft, in_fftw, out_fftw, FFTW_RODFT10, FFTW_ESTIMATE)
   
         !New (Fortran) way:
         call dfftw_plan_dft_r2c_1d(r2c_plan, n_dft, in_r2c, out_r2c, FFTW_PATIENT)
         call dfftw_plan_dft_c2r_1d(c2r_plan, n_dft, out_r2c, in_r2c, FFTW_PATIENT)
         call dfftw_plan_dft_1d(fft_plan, n_dft, in_fftw, out_fftw, FFTW_FORWARD, FFTW_PATIENT)
         call dfftw_plan_dft_1d(ifft_plan, n_dft, in_fftw, out_fftw, FFTW_BACKWARD, FFTW_PATIENT)
         call dfftw_plan_dft_1d(dct_plan, n_dft, in_fftw, out_fftw, FFTW_REDFT10, FFTW_ESTIMATE)
         call dfftw_plan_dft_1d(dst_plan, n_dft, in_fftw, out_fftw, FFTW_RODFT10, FFTW_ESTIMATE)
      endif
		
		!Start up PGPLOT window:
		if(use_pgplot) then
			dev = '/xwin'
			pgplot_id = pgopen(trim(dev))
			call PGQWIN (X1, X2, Y1, Y2)
			write(*,*) X1, X2, Y1, Y2
			X2 = nz_vertical*1.0
			Y2 = nz_vertical*1.0
			call PGWNAD (1.0, X2, 1.0, Y2)
			call PGQWIN (X1, X2, Y1, Y2)
			write(*,*) X1, X2, Y1, Y2
            
            !Set up the colors through interpolation
            call PGQCIR(c1, c2)
            rmin = 0.0
            rmax = 1.0
            gmin = 0.5
            gmax = 0.5
            bmin = 1.0
            bmax = 0.0
            do k = c1, c2
                call PGSCR (k, rmin + (rmax-rmin)*(k-c1)/(c2-c1), &
                    gmin + (gmax-gmin)*(k-c1)/(c2-c1), &
                    bmin + (bmax-bmin)*(k-c1)/(c2-c1))
            end do
		endif
		
	end subroutine initialize
	
	!Simple pgplot test to make sure it's installed and linked correctly
	subroutine pgplot_test
		INTEGER I, IER, PGBEG
		REAL XR(100), YR(100)
		REAL XS(5), YS(5)
		DATA XS/1.,2.,3.,4.,5./
		DATA YS/1.,4.,9.,16.,25./
		IER = PGBEG(0,'?',1,1)
		IF (IER.NE.1) STOP
		CALL PGENV(0.,10.,0.,20.,0,1)
		CALL PGLAB('(x)', '(y)', 'A Simple Graph')
		CALL PGPT(5,XS,YS,9)
		DO 10 I=1,60
		  XR(I) = 0.1*I
		  YR(I) = XR(I)**2
		10 CONTINUE
		CALL PGLINE(60,XR,YR)
		CALL PGEND()
	end subroutine pgplot_test
	
	!Generates the tridiagonal matrix used in the Poisson solve to obtain psi (the
	!velocity streamfunction) from omega (the vorticity) in the y-direction
	!Parameters:
	!n - integer of the current Fourier mode
	subroutine gen_tridiag_matrix(n)
		integer, intent(in) :: n
		integer :: k	!Loop index
		
		do k=2,nz_vertical-1
			dia(k) = (n*pi/a)**2 + 2d0/(dz**2)
		end do
	end subroutine gen_tridiag_matrix
	
	!Gary's solver from Appendix A of his book:
	!This is a tridiagonal matrix solver written in Fortran based on LINPACK and LAPACK 
	!routines. Note, in this routine the first index of the arrays is 1 and the last is Nz; 
	!if one wishes to start the arrays with index 0 and end with Nz − 1, all index 
	!references in this routine would need to be reduced by 1.
	subroutine tridi(nz,rhs,sol,sub,dia,sup,wk1,wk2)
		double precision, dimension (1:nz), intent(inout) :: rhs,sol,sub,dia,sup,wk1,wk2 
		integer, intent(in) :: nz
		integer :: i
	
		wk1(1)=1./dia(1) 
		wk2(1)=sup(1)*wk1(1) 
		do i=2,nz-1
			wk1(i)=1./(dia(i)-sub(i)*wk2(i-1))
			wk2(i)=sup(i)*wk1(i) 
		enddo
		wk1(nz)=1./(dia(nz)-sub(nz)*wk2(nz-1))
		sol(1)=rhs(1)*wk1(1) 
		do i=2,nz
			sol(i)=(rhs(i)-sub(i)*sol(i-1))*wk1(i) 
		enddo
		do i=nz-1,1,-1 
			sol(i)=sol(i)-wk2(i)*sol(i+1)
		enddo
	end subroutine tridi
	
	!This subroutine calculates the values of the time derivaties of temp & omega, used in
	!generalizing the timestepping. 
	subroutine eval_rhs
	end subroutine eval_rhs
	
	!Calculates all the values for the next time step from the current values
	subroutine evolve_step
	
		integer :: k,n,i	!Loop indicies
		
		!Initialize the derivatives we compute in this subroutine:
		dtemp_dt(:,:,2) = 0d0
		domega_dt(:,:,2) = 0d0
	
		!First step is to calculate the necessary z-derivatives of variables using the
		!finite difference method:
		!Note, we don't need to solve the equations on the boundaries since our solutions
		!are constructed to always satisfy them.
		!$OMP PARALLEL DO
		do k=2,nz_vertical-1
			dtemp_dz2(k,:) = (temp(k+1,:) - 2d0*temp(k,:) + temp(k-1,:))/(dz**2)
			domega_dz2(k,:) = (omega(k+1,:) - 2d0*omega(k,:) + omega(k-1,:))/(dz**2)
		end do
		!$OMP END PARALLEL DO
		
		!We also need the first derivatives of psi and temp when calculating the nonlinear
		!terms:
		if(.not.do_linear) then
			!Need to reset the nonlinear terms since we're recomputing them from a sum:
			nonlinear_temp(:,:) = 0d0
			nonlinear_omega(:,:) = 0d0
			!$OMP PARALLEL DO
			do k=2,nz_vertical-1
				dtemp_dz(k,:) = (temp(k+1,:) - temp(k-1,:))/(2*dz)
				domega_dz(k,:) = (omega(k+1,:) - omega(k-1,:))/(2*dz)
				dpsi_dz(k,:) = (psi(k+1,:) - psi(k-1,:))/(2*dz)
				do i=1,nn_horizontal
					nonlinear_temp(k,0) = nonlinear_temp(k,0) - &
						pi/(2*a)*(i*dpsi_dz(k,i)*temp(k,i) + &
						i*psi(k,i)*dtemp_dz(k,i))
				end do
			end do
			!$OMP END PARALLEL DO
            
		!	!$OMP PARALLEL DO
			!do i=1,nn_horizontal
			!	nonlinear_temp(2:nz_vertical-1,0) = nonlinear_temp(2:nz_vertical-1,0) - &
			!		pi/(2*a)*(i*dpsi_dz(2:nz_vertical-1,i)*temp(2:nz_vertical-1,i) + &
			!		i*psi(2:nz_vertical-1,i)*dtemp_dz(2:nz_vertical-1,i))
			!end do
		!	!$OMP END PARALLEL DO
		endif
		
		!Second step is to calculate the time derivatives of omega and temp. These come
		!from eqs. (3.3 & 3.4) in Gary's book:
		
		!Add in the nonlinear terms using a Galerkin method (if necessary):
		!First is the n=0 term for the cosine series of [(v*del)T]_n:
		!Which n,m values do we want? Those where n-m = 0	
		!First, the n=0 terms:
		dtemp_dt(2:nz_vertical-1,0,2) = dtemp_dz2(2:nz_vertical-1,0)
		!domega_dt(2:nz_vertical-1,0,2) = prandtl*domega_dz2(2:nz_vertical-1,0) !No n=0 term for omega
        
        if(.not.do_spectral_transform) then
            !$OMP PARALLEL DO
    		do n=1,nn_horizontal
    			!Always need to compute the linear terms:
    			if(do_linear) then
    				dtemp_dt(2:nz_vertical-1,n,2) = (n*pi/a)*psi(2:nz_vertical-1,n) + &
    					(dtemp_dz2(2:nz_vertical-1,n) - temp(2:nz_vertical-1,n)*(n*pi/a)**2)
    				domega_dt(2:nz_vertical-1,n,2) = rayleigh*prandtl*(n*pi/a)*temp(2:nz_vertical-1,n) + &
    					prandtl*(domega_dz2(2:nz_vertical-1,n) - omega(2:nz_vertical-1,n)*(n*pi/a)**2)
    			else
    				dtemp_dt(2:nz_vertical-1,n,2) = &	!Remove the linear approximation to the advection term
    					(dtemp_dz2(2:nz_vertical-1,n) - temp(2:nz_vertical-1,n)*(n*pi/a)**2)
    				domega_dt(2:nz_vertical-1,n,2) = rayleigh*prandtl*(n*pi/a)*temp(2:nz_vertical-1,n) + &
    					prandtl*(domega_dz2(2:nz_vertical-1,n) - omega(2:nz_vertical-1,n)*(n*pi/a)**2)
    				!Add in the nonlinear terms using a Galerkin method (figure out which modes i
    				!contribute to the mode n being considered in the outer loop)
    				do i=1,nn_horizontal
    					!i+j = nn:
    					if((1.le.(n-i)).and.((n-i).le.nn_horizontal)) then
    						!write(*,*) 'Entering 1 < (n-i) < nn_horizontal branch:', (n-i)
    						nonlinear_omega(2:nz_vertical-1,n) = nonlinear_omega(2:nz_vertical-1,n) - &
    							pi/(2*a)*((n-i)*psi(2:nz_vertical-1,n-i)*domega_dz(2:nz_vertical-1,i) - &
    							i*dpsi_dz(2:nz_vertical-1,n-i)*omega(2:nz_vertical-1,i))
    						nonlinear_temp(2:nz_vertical-1,n) = nonlinear_temp(2:nz_vertical-1,n) - &
    							pi/(2*a)*(-i*dpsi_dz(2:nz_vertical-1,n-i)*temp(2:nz_vertical-1,i) + &
    							(n-i)*psi(2:nz_vertical-1,n-i)*dtemp_dz(2:nz_vertical-1,i))
    					endif
    					!i-j = nn:
    					if((1.le.(i-n)).and.((i-n).le.nn_horizontal)) then
    						!write(*,*) 'Entering 1 < (i-n) < nn_horizontal branch:', (i-n)
    						nonlinear_omega(2:nz_vertical-1,n) = nonlinear_omega(2:nz_vertical-1,n) - &
    							pi/(2*a)*((i-n)*psi(2:nz_vertical-1,i-n)*domega_dz(2:nz_vertical-1,i) + &
    							i*dpsi_dz(2:nz_vertical-1,i-n)*omega(2:nz_vertical-1,i))
    						nonlinear_temp(2:nz_vertical-1,n) = nonlinear_temp(2:nz_vertical-1,n) - &
    							pi/(2*a)*(i*dpsi_dz(2:nz_vertical-1,i-n)*temp(2:nz_vertical-1,i) + &
    							(i-n)*psi(2:nz_vertical-1,i-n)*dtemp_dz(2:nz_vertical-1,i))
                            !write(*,*) nonlinear_temp(3, n)
    					endif
    					!j-i = nn:
    					if((1.le.(i+n)).and.((i+n).le.nn_horizontal)) then
    						!write(*,*) 'Entering 1 < (i+n) < nn_horizontal branch:', (i+n)
    						nonlinear_omega(2:nz_vertical-1,n) = nonlinear_omega(2:nz_vertical-1,n) + & 
    							pi/(2*a)*((i+n)*psi(2:nz_vertical-1,i+n)*domega_dz(2:nz_vertical-1,i) + &
    							i*dpsi_dz(2:nz_vertical-1,i+n)*omega(2:nz_vertical-1,i))
    						nonlinear_temp(2:nz_vertical-1,n) = nonlinear_temp(2:nz_vertical-1,n) - &
    							pi/(2*a)*(i*dpsi_dz(2:nz_vertical-1,i+n)*temp(2:nz_vertical-1,i) + &
    							(i+n)*psi(2:nz_vertical-1,i+n)*dtemp_dz(2:nz_vertical-1,i))
                            !write(*,*) nonlinear_temp(3, n)
    					endif
					
    					!read(*,*)
    				end do
				
    				!lastly, the term outside the double sum with dT0/dz (only depends on n):
    				nonlinear_temp(2:nz_vertical-1,n) = nonlinear_temp(2:nz_vertical-1,n) - &
    					(n*pi/a)*dtemp_dz(2:nz_vertical-1,0)*psi(2:nz_vertical-1,n)
					
    				if(dbg) then
    					write(*,*) 'Final nonlinear terms:'
    					write(*,*) n, nonlinear_temp(nz_vertical/3,n), nonlinear_omega(nz_vertical/3,n)
    				endif
    			endif !do_linear 
            end do
            !$OMP END PARALLEL DO
        endif !do_spectral transform
        
        if(do_spectral_transform) then
            do n=1,nn_horizontal
    			dtemp_dt(2:nz_vertical-1,n,2) = &	!Remove the linear approximation to the advection term
    				(dtemp_dz2(2:nz_vertical-1,n) - temp(2:nz_vertical-1,n)*(n*pi/a)**2)
    			domega_dt(2:nz_vertical-1,n,2) = rayleigh*prandtl*(n*pi/a)*temp(2:nz_vertical-1,n) + &
    				prandtl*(domega_dz2(2:nz_vertical-1,n) - omega(2:nz_vertical-1,n)*(n*pi/a)**2)
            end do
            !write(*,*) "Engaging spectral transform!"
            ux_spatial(:,:) = 0d0
            uz_spatial(:,:) = 0d0
            omega_spatial(:,:) = 0d0
            temp_spatial(:,:) = 0d0
            domegadz_spatial(:,:) = 0d0
            dtempdz_spatial(:,:) = 0d0
            dtempdx_spatial(:,:) = 0d0
            domegadx_spatial(:,:) = 0d0
            !First step is to convert to spatial coordinates:
            !Convert the horizontal Fourier modes into real space
        	!do i=0,n_dft-1
        	!	do n=0,nn_horizontal
        	!		ux_spatial(:,i) = ux_spatial(:,i) - dpsi_dz(:,n)*sin(n*pi*dx*i/a)
        	!		uz_spatial(:,i) = uz_spatial(:,i) + (n*pi/a)*psi(:,n)*cos(n*pi*dx*i/a)
            !       omega_spatial(:,i) = omega_spatial(:,i) + omega(:,n)*sin(n*pi*dx*i/a)
            !       temp_spatial(:,i) = temp_spatial(:,i) + temp(:,n)*cos(n*pi*dx*i/a)
            !       domegadz_spatial(:,i) = domegadz_spatial(:,i) + domega_dz(:,n)*sin(n*pi*dx*i/a)
            !       dtempdz_spatial(:,i) = dtempdz_spatial(:,i) + dtemp_dz(:,n)*cos(n*pi*dx*i/a)
            !       dtempdx_spatial(:,i) = dtempdx_spatial(:,i) - (n*pi/a)*temp(:,n)*sin(n*pi*dx*i/a)
            !       domegadx_spatial(:,i) = domegadx_spatial(:,i) + (n*pi/a)*omega(:,n)*cos(n*pi*dx*i/a)
        	!	end do
            !end do
            !galerkin_test = dtempdx_spatial
            !write(*,*) 'ux_spatial (nz_vertical/3,:) direct sum:'
            !write(*,*) ux_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'uz_spatial (nz_vertical/3,:) direct sum:'
            !write(*,*) uz_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'omega_spatial (nz_vertical/3,:) direct sum:'
            !write(*,*) omega_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'domegadz_spatial (nz_vertical/3,:) direct sum:'
            !write(*,*) domegadz_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'temp_spatial (nz_vertical/3,:) direct sum:'
            !write(*,*) temp_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'dtempdz_spatial (nz_vertical/3,:) direct sum:'
            !write(*,*) dtempdz_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'dtempdz_spatial (nz_vertical/3,:) direct sum:'
            !write(*,*) dtempdz_spatial(nz_vertical/3,:)
            !write(*,*)
            
            !Faster to do the frequency -> spatial conversion with an inverse FFT:
            do k=1, nz_vertical
                in_fftw = 0d0
                in_fftw(0:nn_horizontal) = -dpsi_dz(k,0:nn_horizontal)
                call fftw_execute_dft(ifft_plan, in_fftw, out_fftw)
                ux_spatial(k,:) = imagpart(out_fftw)
                
                in_fftw = 0d0
                do n=0, nn_horizontal
                    in_fftw(n) = (n*pi/a)*psi(k,n)
                end do
                call fftw_execute_dft(ifft_plan, in_fftw, out_fftw)
                uz_spatial(k,:) = realpart(out_fftw)
                
                in_fftw = 0d0
                in_fftw(0:nn_horizontal) = omega(k,0:nn_horizontal)
                call fftw_execute_dft(ifft_plan, in_fftw, out_fftw)
                omega_spatial(k,:) = imagpart(out_fftw)
                
                in_fftw = 0d0
                in_fftw(0:nn_horizontal) = domega_dz(k,0:nn_horizontal)
                call fftw_execute_dft(ifft_plan, in_fftw, out_fftw)
                domegadz_spatial(k,:) = imagpart(out_fftw)
                
                in_fftw = 0d0
                in_fftw(0:nn_horizontal) = temp(k,0:nn_horizontal)
                call fftw_execute_dft(ifft_plan, in_fftw, out_fftw)
                temp_spatial(k,:) = realpart(out_fftw)
                
                in_fftw = 0d0
                in_fftw(0:nn_horizontal) = dtemp_dz(k,0:nn_horizontal)
                call fftw_execute_dft(ifft_plan, in_fftw, out_fftw)
                dtempdz_spatial(k,:) = realpart(out_fftw)
                
                in_fftw = 0d0
                do n=0, nn_horizontal
                    in_fftw(n) = -(n*pi/a)*temp(k,n)
                end do
                call fftw_execute_dft(ifft_plan, in_fftw, out_fftw)
                dtempdx_spatial(k,:) = imagpart(out_fftw)
                
                in_fftw = 0d0
                do n=0, nn_horizontal
                    in_fftw(n) = (n*pi/a)*omega(k,n)
                end do
                call fftw_execute_dft(ifft_plan, in_fftw, out_fftw)
                domegadx_spatial(k,:) = realpart(out_fftw)
            end do
            !write(*,*) 'ux_spatial (nz_vertical/3,:) iFFT:'
            !write(*,*) ux_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'uz_spatial (nz_vertical/3,:) iFFT:'
            !write(*,*) uz_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'omega_spatial (nz_vertical/3,:) iFFT:'
            !write(*,*) omega_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'domegadz_spatial (nz_vertical/3,:) iFFT:'
            !write(*,*) domegadz_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'temp_spatial (nz_vertical/3,:) iFFT:'
            !write(*,*) temp_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'dtempdz_spatial (nz_vertical/3,:) iFFT:'
            !write(*,*) dtempdz_spatial(nz_vertical/3,:)
            !write(*,*)
            
            !write(*,*) 'absolute difference: ux_spatial (nz_vertical/3,:)'
            !write(*,*) galerkin_test(nz_vertical/3,:) - ux_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'absolute difference: uz_spatial (nz_vertical/3,:)'
            !write(*,*) galerkin_test(nz_vertical/3,:) - uz_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'absolute difference: omega_spatial (nz_vertical/3,:)'
            !write(*,*) galerkin_test(nz_vertical/3,:) - omega_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'absolute difference: domegadz_spatial (nz_vertical/3,:)'
            !write(*,*) galerkin_test(nz_vertical/3,:) - domegadz_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'absolute difference: temp_spatial (nz_vertical/3,:)'
            !write(*,*) galerkin_test(nz_vertical/3,:) - temp_spatial(nz_vertical/3,:)
            !write(*,*)
            !write(*,*) 'absolute difference: dtempdz_spatial (nz_vertical/3,:)'
            !write(*,*) galerkin_test(nz_vertical/3,:) - dtempdz_spatial(nz_vertical/3,:)
            !write(*,*)
            
            !Ok - all absolute differences have been checked to be ~1d-16 (machine precision) (07/31/14)
            
            !Next, we multiply the nonlinear terms together and convert them back to 
            !spectral space (each row at a time). One question - should we differentiate in the 
            !horiontal direction before or after we conver to spatial coordinates?
            
            !Differentiating after the inverse FFT:
            !(perhaps I should differentiate before the iFFT to agree with Galerkin?)
            !dtempdx_spatial(:,:) = 0d0
            !domegadx_spatial(:,:) = 0d0
            !do i=1,n_dft-2
            !    dtempdx_spatial(:,i) = (temp_spatial(:,i+1) - temp_spatial(:,i-1))/(2*dx)
            !    domegadx_spatial(:,i) = (omega_spatial(:,i+1) - omega_spatial(:,i-1))/(2*dx)
            !end do
            
            !write(*,*) 'absolute difference: dtempdx_spatial (nz_vertical/3,:)'
            !write(*,*) galerkin_test(nz_vertical/3,:) - dtempdx_spatial(nz_vertical/3,:)
            !write(*,*)
            
            !Output Galerkin method results:
            !write(*,*) 'Galerkin:'
            !write(*,'(a5,2a25)') 'n','nonlinear_temp','nonlinear_omega'
    	    !do i=0,nn_horizontal
            !    write(*,'(i5,3ES25.10)') i, nonlinear_temp(nz_vertical/3,i), nonlinear_omega(nz_vertical/3,i)
    		!end do 
            !write(*,*)
            
            do k=1,nz_vertical
                in_fftw = 0d0
                in_fftw = -ux_spatial(k,:)*dtempdx_spatial(k,:) - uz_spatial(k,:)*dtempdz_spatial(k,:)
                call fftw_execute_dft(fft_plan, in_fftw, out_fftw)
                nonlinear_temp(k,:) = realpart(out_fftw(0:nn_horizontal))/n_dft
                nonlinear_temp(k,1:nn_horizontal) = 2*nonlinear_temp(k,1:nn_horizontal) !not sure why yet
                
                in_fftw = 0d0
                in_fftw = -ux_spatial(k,:)*domegadx_spatial(k,:) - uz_spatial(k,:)*domegadz_spatial(k,:)
                call fftw_execute_dft(fft_plan, in_fftw, out_fftw)
                nonlinear_omega(k,:) = -imagpart(out_fftw(0:nn_horizontal))/n_dft !WHY IS A NEGATIVE SIGN NEEDED HERE?? Perhaps something to do with a complex conjugate somewhere...
                nonlinear_omega(k,1:nn_horizontal) = 2*nonlinear_omega(k,1:nn_horizontal) !not sure why yet
            end do
            
            !Output spectral transform results:
            !write(*,*) 'Spectral transform:'
            !write(*,'(a5,2a25)') 'n','nonlinear_temp','nonlinear_omega'
    		!do i=0,nn_horizontal
            !    write(*,'(i5,3ES25.10)') i, nonlinear_temp(nz_vertical/3,i), nonlinear_omega(nz_vertical/3,i)
    		!end do 
            !write(*,*)
            
            !Galerkin & spectral transform calculations of nonlinear terms verified to be within machine
            !precision of each other! (07/31/14)
            
        endif !do_spectral_transform

		!Add in all the nonlinear terms (only nonzero if do_linear = .false.)
		dtemp_dt(2:nz_vertical-1,:,2) = dtemp_dt(2:nz_vertical-1,:,2) + &
			nonlinear_temp(2:nz_vertical-1,:)
		domega_dt(2:nz_vertical-1,:,2) = domega_dt(2:nz_vertical-1,:,2) + &
			nonlinear_omega(2:nz_vertical-1,:)
		
		!Keep the old value of temp to examine the growth rate
		prev_temp(:,:) = temp(:,:)
		
		!Next, we use the current time derivatives along with those from the previous
		!step to calculate the fluid quantities at the next time step.
		
		!For the first step, we only have initial condition information so can't go back the
		!required two steps in the Adams-Bashforth scheme. Instead, the first step will be computed
		!using a 2nd order Runge-Kutta solver (to be implemented)
		!if(nstep.eq.1) then
		!	temp(:,:) = temp(:,:) + dt*
		!else
			!Here, we're using the Adams-Bashforth scheme (explicit):
			!temp(:,:) = temp(:,:) + dt/2d0*(3*dtemp_dt(:,:,2) -  dtemp_dt(:,:,1))
			!omega(:,:) = omega(:,:) + dt/2d0*(3*domega_dt(:,:,2) -  domega_dt(:,:,1))
			temp(:,:) = temp(:,:) + dt*((1d0 + dt/(2*dt_old))*dtemp_dt(:,:,2) -  &
                (dt/(2*dt_old))*dtemp_dt(:,:,1))
			omega(:,:) = omega(:,:) + dt*((1d0 + dt/(2*dt_old))*domega_dt(:,:,2) -  &
                (dt/(2*dt_old))*domega_dt(:,:,1))
		!endif
		
		!Check for NaN values:
		! !$OMP PARALLEL DO
		!do k=1,nz_vertical
		!	do n=0,nn_horizontal
		!		if (isnan(temp(k,n))) then
		!			write(*,*) '"temp" is a NaN', k, n
		!			stop 20
		!		endif
		!		if (isnan(omega(k,n))) then
		!			write(*,*) '"omega" is a NaN', k, n
		!			stop 20
		!		endif
		!	end do
		!end do
		! !$OMP END PARALLEL DO
		
		!Then find the streamfunction psi by doing a Poisson solve with our tridiagonal
		!matrix routine.
		!$OMP PARALLEL DO
		do n=1,nn_horizontal
			call gen_tridiag_matrix(n)
			call tridi(nz_vertical,omega(:,n),psi(:,n),sub,dia,sup,wk1,wk2)
		end do
		!$OMP END PARALLEL DO
		
		!Now update the old timestep derivatives.
		dtemp_dt(:,:,1) = dtemp_dt(:,:,2)
		domega_dt(:,:,1) = domega_dt(:,:,2)
		
	end subroutine evolve_step
    
    !Checks the CFL condition (||
    !ux_cfl = -dpsi_dz, uz_cfl = dpsi_dx
    subroutine check_cfl(step_num)
        
        integer, intent(in) :: step_num
        integer :: k,n
        double precision :: cfl_x, cfl_z
        !Reduce timestep if the CFL limit gets within cfl_timestep_limit_factor*dt (cfl_timestep_limit_factor > 1)
        double precision :: cfl_timestep_limit_factor
        !If reducing timestep due to CFL condition, then make the new timestep dt*cfl_reduction_factor
        double precision :: cfl_reduction_factor
        
        cfl_reduction_factor = 0.8
        cfl_timestep_limit_factor = 1.2
        
        if(cfl_timestep_limit_factor.le.1.0) then
            write(*,*) 'cfl_timestep_limit_factor must be > 1.0'
            stop 1
        endif
        
        ux_cfl(:,:) = 0d0
        uz_cfl(:,:) = 0d0
        dt_old = dt
        
        !Convert the horizontal Fourier modes into real space
		do k=1,nz_vertical
			do n=0,nn_horizontal
				ux_cfl(:,k) = ux_cfl(:,k) - dpsi_dz(:,n)*sin(n*pi*dz*k)
				uz_cfl(:,k) = uz_cfl(:,k) + (n*pi/a)*psi(:,n)*cos(n*pi*dz*k)
			end do
		end do
        
        cfl_x = (a*1d0/nn_horizontal)/maxval(abs(ux_cfl))
        cfl_z = dz/maxval(abs(uz_cfl))
        !write(*,*) 'max ux_cfl, dx/ux_cfl:', (a/nn_horizontal)/cfl_x, cfl_x
        !write(*,*) 'max uz_cfl, dz/uz_cfl:', dz/cfl_z, cfl_z
        !Compute running averages:
        avg_cfl_x = (avg_cfl_x*(step_num/cfl_tick - 1) + cfl_x)/(step_num/cfl_tick)
        avg_cfl_z = (avg_cfl_z*(step_num/cfl_tick - 1) + cfl_z)/(step_num/cfl_tick)
        
        !Reduce the timestep, if necessary:
        if((a*1d0/nn_horizontal)/maxval(ux_cfl).lt.dt*cfl_timestep_limit_factor) then
            dt = dt*cfl_reduction_factor
            write(*,*) 'Reducing timestep due to CFL condition (x-direction)', dt
        endif
        if(dz/maxval(uz_cfl).lt.dt*cfl_timestep_limit_factor) then
            dt = dt*cfl_reduction_factor
            write(*,*) 'Reducing timestep due to CFL condition (z-direction)', dt
        endif
        
        !Increase the timestep, if necessary:
        !maybe try this later with a running average of the last ~1000 timesteps?
        !if(((step_num > 10000).and.(avg_cfl_x > 2*dt*cfl_timestep_limit_factor)).and.&
        !    (avg_cfl_z > 2*dt*cfl_timestep_limit_factor)) then
        !    dt = dt/cfl_reduction_factor
        !endif
        
    end subroutine check_cfl
	
	!Outputs selected fluid variables to the screen
	subroutine output_vals
		integer :: i
		
		!do i=1,nz_vertical
		!	write(*,'(a5,i3,a5,99ES25.10)') 'temp(',i,',1:2)', temp(i,1:2), dtemp_dz2(i,1:2)
		!end do

		!write(*,'(a25,99ES25.10)') 'temp(1,:):', temp(1,:)
		!write(*,'(a25,99ES25.10)') 'temp(nz/3,:):',temp(nz_vertical/3,:)
		!write(*,'(a25,99ES25.10)') 'temp(2*nz/3,:):',temp(2*nz_vertical/3,:)
		!write(*,'(a25,99ES25.10)') 'temp(nz,:):',temp(nz_vertical,:)
		!write(*,'(a25,99ES25.10)') 'omega(1,:):', omega(1,:)
		!write(*,'(a25,99ES25.10)') 'omega(nz/3,:):', omega(nz_vertical/3,:)
		!write(*,'(a25,99ES25.10)') 'omega(nz,:):', omega(nz_vertical,:)
		!write(*,'(a25,99ES25.10)') 'psi(nz/3,:):', psi(nz_vertical/3,:)
		!write(*,'(a25,99ES25.10)') 'dtemp_dt(nz/3,:,1):', dtemp_dt(nz_vertical/3,:,1)
		!write(*,'(a25,99ES25.10)') 'domega_dt(nz/3,:,1):', domega_dt(nz_vertical/3,:,1)
		
		!write(*,'(a40)') 'non-zero mode temp growth(nz/3,:):'
		!do i=0,nn_horizontal
		!	if((temp(nz_vertical/3,i).ne.0d0).and.(prev_temp(nz_vertical/3,i).ne.0d0)) then
		!		write(*,'(a17,i2,a3,ES25.10)') 'temp growth(nz/3,',i,'):', &
		!			log(abs(temp(nz_vertical/3,i))) - log(abs(prev_temp(nz_vertical/3,i)))
		!	endif
			!write(*,'(ES25.10)',advance='no') 'temp growth(nz/3,:):', &
			!	abs(log(temp(nz_vertical/3,i)) - log(prev_temp(nz_vertical/3,i)))
            !end do
		!write(*,*)
        
        write(*,'(a5,3a25)') 'n','temp','omega','psi'
		do i=0,nn_horizontal
            write(*,'(i5,3ES25.10)') i, temp(nz_vertical/3,i), omega(nz_vertical/3,i), psi(nz_vertical/3,i)
			!write(*,'(a17,i2,a3,ES25.10)') 'temp(nz/3,',i,'):', temp(nz_vertical/3,i)
            !write(*,'(a17,i2,a3,ES25.10)') 'omega(nz/3,',i,'):', omega(nz_vertical/3,i)
            !write(*,'(a17,i2,a3,ES25.10)') 'psi(nz/3,',i,'):', psi(nz_vertical/3,i)
		end do
        
		!do i=0,nn_horizontal
		!	write(*,'(a17,i2,a3,999ES25.10)') 'dtemp_dt(z(i),',i,'):', dtemp_dt(1:nz_vertical,i,2)
        !    write(*,'(a17,i2,a3,999ES25.10)') 'domega_dt(z(i),',i,'):', domega_dt(1:nz_vertical,i,2)

			!write(*,'(a17,i2,a3,99ES25.10)') 'temp(z(i),',i,'):', temp(1:nz_vertical,i)
            !write(*,'(a17,i2,a3,99ES25.10)') 'omega(z(i),',i,'):', omega(1:nz_vertical,i)
        !    write(*,'(a17,i2,a3,999ES25.10)') 'psi(z(i),',i,'):', psi(1:nz_vertical,i)
        !end do
		
		write(*,*)
	end subroutine output_vals
	
	!Calls a pgplot window to update with the new fluid quantities
	subroutine plot_vals
		integer :: k,n
		integer :: i1,i2,j1,j2	!Array bounds for PGPLOT
		real :: a1,a2	!Colormap bounds for the plotting
		integer :: c1,c2,itf
		double precision :: pgplot_theta	!Rotation angle
		real, dimension(6) :: tr	!Transformation matrix (array grid -> plot)
		double precision, dimension(5) :: hl, hr, hg, hb, alev

		i1 = 1
		i2 = nz_vertical
		j1 = 1
		j2 = nz_vertical
		!
		!a1 = -0.1e1
		!a2 = 0.1e1
        a1 = 0.0
        a2 = 1.0
		
		!Possible transformation matrix for rotated coordinate system:
		pgplot_theta = 0d0
		tr(1) = -0.5
		tr(2) = cos(pgplot_theta)
        tr(3) = -sin(pgplot_theta)
        tr(4) = -0.5
        tr(5) = sin(pgplot_theta)
        tr(6) = cos(pgplot_theta)
		
		!Transformation matrix:
		tr = (/0.0, 0.0, 1.0, 0.0, 1.0, 0.0/)	!First index is vertical, second is horizontal
        
        !Query and set color indicies:
        !call PGQCIR(c1, c2)
        !call PGQITF(itf)
        !write(*,*) c1, c2, itf
        !itf = 1
        !call PGSITF(itf)
        !read(*,*)

		!hl = (/0.0, 0.2, 0.4, 0.6, 1.0/)
		!hr = (/0.0, 0.5, 1.0, 1.0, 1.0/)
		!hg = (/0.0, 0.0, 0.5, 1.0, 1.0/) 
		!hb = (/0.0, 0.0, 0.0, 0.3, 1.0/)
		!hl = (/0.0, 0.25, 0.5, 0.75, 1.0/)
		!hr = (/0.0, 0.1, 0.8, 1.0, 1.0/)
		!hg = (/0.0, 0.0, 0.0, 0.0, 0.0/) 
		!hb = (/1.0, 0.9, 0.2, 0.0, 0.0/)
		!call pgctab(hl, hr, hg, hb, 5, 1.0, 0.5)
        
        !Transform the horizontal components into real space
		temp_display(:,:) = 0d0
		omega_display(:,:) = 0d0
		!$OMP PARALLEL DO
		do k=0,nz_vertical-1
			do n=0,nn_horizontal
				temp_display(:,k) = temp_display(:,k) + temp(:,n)*cos(n*pi*dz*k)
				omega_display(:,k) = omega_display(:,k) + omega(:,n)*sin(n*pi*dz*k)
			end do
		end do
		!$OMP END PARALLEL DO
		
		!call pgimag (real(temp_display), nz_vertical, nn_horizontal, i1, i2, j1, j2, a1, a2, tr)
		!call pgpage()
		call pgimag (real(temp_display), nz_vertical, nz_vertical, i1, i2, j1, j2, a1, a2, tr)
	end subroutine plot_vals
	
	!Deallocates and cleans up all data, file IO streams, etc.
	subroutine shutdown
        if(do_spectral_transform) then
            call fftw_free(p_in)
            call fftw_free(p_out)
            call fftw_destroy_plan(fft_plan)
            call fftw_destroy_plan(ifft_plan)
            call fftw_destroy_plan(dst_plan)
            call fftw_destroy_plan(dct_plan)
            call fftw_destroy_plan(r2c_plan)
            call fftw_destroy_plan(c2r_plan)
        endif
		call pgclos()
	end subroutine shutdown

end module hydro


!Main driver program that loops through time steps
program main

	use hydro
	
	implicit none
	
	integer :: i	!Loop index
	
	call initialize()
	!call pgplot_test()
	call output_vals()
	if(use_pgplot) call plot_vals()
	
	do i=1,num_steps
		call evolve_step()
		if(use_pgplot.and.(mod(i,plot_tick).eq.2)) call plot_vals()
        
        !Remove timestep increaser for now, it seems to be breaking things
        !if((mod(i,100000).eq.1).and.(fixed_dt.gt.0)) then
        !    dt_old = dt
        !    dt = fixed_dt    !Try to reset the timestep a bit larger every 50k timesteps (offset from CFL check to not interfere)
        !endif
        
        if(mod(i,output_tick).eq.0) then
            write(*,*) 'Evolving step',i
            call output_vals()
        endif
        
        if(mod(i,cfl_tick).eq.0) then
            !write(*,*) 'Evolving step',i
            call check_cfl(i)
        endif
        
        if(dbg) read(*,*)
	end do
	
	call shutdown()
	
end program