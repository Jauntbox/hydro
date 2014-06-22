!Code to solve the Boussinesq equations for Rayleigh-Bernard convection in 2D using
!spectral decomposition in the horizontal, and finite-difference in the vertical.
!-
!Written by: Kevin Moore
!6/14/14

module hydro

	implicit none
	
	!Overall simulation parameters:
	!Only solve the linear problem (eg. equations )
	logical, parameter :: do_linear = .false.
	logical, parameter :: use_pgplot = .true.
	
	!Domain size parameters
	integer, parameter :: nz_vertical = 500		!Number of vertical zones in simulation
	integer, parameter :: nz_horizontal = 10	!Number of Fourier modes in the horizontal
	double precision, parameter :: dz = 1d0/(nz_vertical-1)	!spacing in the z-direction
	double precision, parameter :: a = sqrt(1.0)	!horizontal/vertical length scale
	
	!Timestep parameters:
	double precision :: dt	!Gets set in initialize()
	!Number of old time steps to keep track of when time-stepping:
	integer, parameter :: num_time_differences = 2
	integer, parameter :: num_steps = 10000		!Number of time steps to take
	
	!Physical parameters:
	double precision, parameter :: rayleigh = 4000.0	!Rayleigh number (dimensionless)
	double precision, parameter :: prandtl = 0.1	!Prandtl number (dimensionless)
	
	!Numerical constants:
	double precision, parameter :: pi = 4d0*atan(1d0)
	
	!Fluid variables:
	double precision, dimension(1:nz_vertical) :: u_x, u_z, z, sub, dia, sup, wk1, wk2
	double precision, dimension(1:nz_vertical, 0:nz_horizontal) :: psi, omega, temp, prev_temp
	!Nonlinear terms: [(v*del)variable]
	double precision, dimension(1:nz_vertical, 0:nz_horizontal) :: nonlinear_temp, nonlinear_omega
	!Arrays of first and second derivatives in the z-direction:
	double precision, dimension(1:nz_vertical, 0:nz_horizontal) :: dtemp_dz2, domega_dz2, &
		dtemp_dz, dpsi_dz, domega_dz
	double precision, dimension(1:nz_vertical, 0:nz_horizontal, 1:num_time_differences) :: &
		domega_dt, dtemp_dt
		
	!PGPLOT variables:
	integer :: pgplot_id		!PGPLOT handle
	character (len=16) :: dev	!PGPLOT device name (eg. '/xwin')
	double precision, dimension(1:nz_vertical, 0:nz_vertical) :: temp_display, omega_display
	
	contains
	
	!Sets up all the variables, initial conditions, pgplot window, etc.
	subroutine initialize
	
		integer :: k,n	!Loop indicies
		integer, dimension(:), pointer :: itime
		integer :: pgopen, isize
		real :: x1, x2, y1, y2	!PGPLOT window bounds
	
		write(*,*) 'Initializing hydro module!'
		write(*,*)
		
		!Make sure we satisfy the CFL limit:
		if(prandtl.lt.1d0) then
			dt = (dz**2)/4d0
		else
			dt = (dz**2)/(4d0*prandtl)
		endif
		
		dt = 0.8*dt
		
		write(*,*) 'Domain size parameters:'
		write(*,'(a25,I25)') 'nz_horizontal:', nz_horizontal
		write(*,'(a25,I25)') 'nz_vertical:', nz_vertical
		write(*,'(a25,ES25.10)') 'a (aspect ratio):', a
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
		call random_number(omega)
		omega = 10d0*omega - 5d0
		!call random_number(psi)
		!psi = 10d0*psi - 5d0
		temp(:,:) = 0d0
		prev_temp(:,:) = 0d0
		temp_display(:,:) = 0d0
		omega_display(:,:) = 0d0
		nonlinear_temp(:,:) = 0d0
		nonlinear_omega(:,:) = 0d0
		dtemp_dt(:,:,:) = 0d0
		domega_dt(:,:,:) = 0d0
		dtemp_dz(:,:) = 0d0 
		dpsi_dz(:,:) = 0d0
		domega_dz(:,:) = 0d0
		dtemp_dz2(:,:) = 0d0
		domega_dz2(:,:) = 0d0
		do k=1,nz_vertical
			z(k) = (k - 1)*dz
			!Satisfies the boundary conditions:
			temp(k,0) = 1d0 - z(k)
			temp(k,3:4) = sin(pi*z(k))
			!temp(k,1:nz_horizontal) = sin(pi*z(k))
		end do
		do k=1,nz_vertical
			do n=1,nz_horizontal
				temp_display(:,k) = temp_display(:,k) + temp(:,n)*cos(n*pi*dz*(k-1))
			end do
		end do
		
		!Set the subdiagonal and superdiagonal of the Poisson solver matrix:
		sub(:) = -1d0/(dz**2)
		sup(:) = -1d0/(dz**2)
		sub(nz_vertical) = 0d0
		sup(1) = 0d0
		dia(:) = 1d0	!Only for the boundary conditions 2:nz_vertical-1 will be overwritten
		
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
		CALL PGEND
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
	
	!Calculates all the values for the next time step from the current values
	subroutine evolve_step
	
		integer :: k,n,i	!Loop indicies
	
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
			end do
			!$OMP END PARALLEL DO
		endif
		
		!Second step is to calculate the time derivatives of omega and temp. These come
		!from eqs. (3.3 & 3.4) in Gary's book:
		
		!Add in the nonlinear terms using a Galerkin method (if necessary):
		!First is the n=0 term for the cosine series of [(v*del)T]_n:
		!Which n,m values do we want? Those where n-m = 0	
		!$OMP PARALLEL DO
		do n=1,nz_horizontal
			!Always need to compute the linear terms:
			dtemp_dt(2:nz_vertical-1,n,2) = (n*pi/a)*psi(2:nz_vertical-1,n) + &
				(dtemp_dz2(2:nz_vertical-1,n) - temp(2:nz_vertical-1,n)*(n*pi/a)**2)
			domega_dt(2:nz_vertical-1,n,2) = rayleigh*prandtl*(n*pi/a)*temp(2:nz_vertical-1,n) + &
				prandtl*(domega_dz2(2:nz_vertical-1,n) - omega(2:nz_vertical-1,n)*(n*pi/a)**2)
			if(.not.do_linear) then
				!Add in the nonlinear terms using a Galerkin method:
				do i=1,nz_horizontal
					!i+j = n:	(can combine this term with the next using abs(n-i))
					if((1.le.(n-i)).and.((n-i).le.nz_horizontal)) then
						nonlinear_omega(2:nz_vertical-1,n) = nonlinear_omega(2:nz_vertical-1,n) + &
							pi/(2*a)*(i*psi(2:nz_vertical-1,i)*domega_dz(2:nz_vertical-1,n-i) - &
							(n-i)*dpsi_dz(2:nz_vertical-1,i)*omega(2:nz_vertical-1,n-i))
						nonlinear_temp(2:nz_vertical-1,n) = nonlinear_temp(2:nz_vertical-1,n) + &
							pi/(2*a)*(-(n-i)*dpsi_dz(2:nz_vertical-1,i)*temp(2:nz_vertical-1,n-i) + &
							i*psi(2:nz_vertical-1,i)*dtemp_dz(2:nz_vertical-1,n-i))
					endif
					!i+j = n:
					if((1.le.(i-n)).and.((i-n).le.nz_horizontal)) then
						nonlinear_omega(2:nz_vertical-1,n) = nonlinear_omega(2:nz_vertical-1,n) + &
							pi/(2*a)*(i*psi(2:nz_vertical-1,i)*domega_dz(2:nz_vertical-1,i-n) - &
							(i-n)*dpsi_dz(2:nz_vertical-1,i)*omega(2:nz_vertical-1,i-n))
						nonlinear_temp(2:nz_vertical-1,n) = nonlinear_temp(2:nz_vertical-1,n) + &
							pi/(2*a)*(-(i-n)*dpsi_dz(2:nz_vertical-1,i)*temp(2:nz_vertical-1,i-n) + &
							i*psi(2:nz_vertical-1,i)*dtemp_dz(2:nz_vertical-1,i-n))
					endif
					!i-j = n:
					if((1.le.(i-n)).and.((i-n).le.nz_horizontal)) then
						nonlinear_omega(2:nz_vertical-1,n) = nonlinear_omega(2:nz_vertical-1,n) - &
							pi/(2*a)*(i*psi(2:nz_vertical-1,i)*domega_dz(2:nz_vertical-1,i-n) + &
							(i-n)*dpsi_dz(2:nz_vertical-1,i)*omega(2:nz_vertical-1,i-n))
						nonlinear_temp(2:nz_vertical-1,n) = nonlinear_temp(2:nz_vertical-1,n) + &
							pi/(2*a)*((i-n)*dpsi_dz(2:nz_vertical-1,i)*temp(2:nz_vertical-1,i-n) + &
							i*psi(2:nz_vertical-1,i)*dtemp_dz(2:nz_vertical-1,i-n))
					endif
					!j-i = n:
					if((1.le.(i+n)).and.((i+n).le.nz_horizontal)) then
						nonlinear_omega(2:nz_vertical-1,n) = nonlinear_omega(2:nz_vertical-1,n) + & 
							pi/(2*a)*(i*psi(2:nz_vertical-1,i)*domega_dz(2:nz_vertical-1,n+i) + &
							(n+i)*dpsi_dz(2:nz_vertical-1,i)*omega(2:nz_vertical-1,n+i))
						nonlinear_temp(2:nz_vertical-1,n) = nonlinear_temp(2:nz_vertical-1,n) + &
							pi/(2*a)*((i+n)*dpsi_dz(2:nz_vertical-1,i)*temp(2:nz_vertical-1,i+n) + &
							i*psi(2:nz_vertical-1,i)*dtemp_dz(2:nz_vertical-1,i+n))
					endif
				end do
				!Now that the nonlinear terms are computed, add them to the time derivatives:
				dtemp_dt(2:nz_vertical-1,n,2) = dtemp_dt(2:nz_vertical-1,n,2) + &
					nonlinear_temp(2:nz_vertical-1,n)
				domega_dt(2:nz_vertical-1,n,2) = domega_dt(2:nz_vertical-1,n,2) + &
					nonlinear_omega(2:nz_vertical-1,n)
			endif
		end do
		!$OMP END PARALLEL DO
		
		!Keep the old value of temp to examine the growth rate
		prev_temp(:,:) = temp(:,:)
		
		!Next, we use the current time derivatives along with those from the previous
		!step to calculate the fluid quantities at the next time step.
		!Here, we're using the Adams-Bashforth scheme (explicit):
		temp(:,:) = temp(:,:) + dt/2d0*(3*dtemp_dt(:,:,2) -  dtemp_dt(:,:,1))
		omega(:,:) = omega(:,:) + dt/2d0*(3*domega_dt(:,:,2) -  domega_dt(:,:,1))
		temp_display(:,:) = 0d0
		omega_display(:,:) = 0d0
		!$OMP PARALLEL DO
		do k=1,nz_vertical
			do n=0,nz_horizontal
				temp_display(:,k) = temp_display(:,k) + temp(:,n)*cos(n*pi*dz*k)
				omega_display(:,k) = omega_display(:,k) + omega(:,n)*sin(n*pi*dz*k)
			end do
		end do
		!$OMP END PARALLEL DO
		
		!Then find the streamfunction psi by doing a Poisson solve with our tridiagonal
		!matrix routine.
		!$OMP PARALLEL DO
		do n=1,nz_horizontal
			call gen_tridiag_matrix(n)
			call tridi(nz_vertical,omega(:,n),psi(:,n),sub,dia,sup,wk1,wk2)
		end do
		!$OMP END PARALLEL DO
		
		!Now update the old timestep derivatives.
		dtemp_dt(:,:,1) = dtemp_dt(:,:,2)
		domega_dt(:,:,1) = domega_dt(:,:,2)
		
	end subroutine evolve_step
	
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
		
		write(*,'(a25)') 'temp growth(nz/3,:):'
		do i=0,nz_horizontal
			if((temp(nz_vertical/3,i).ne.0d0).and.(prev_temp(nz_vertical/3,i).ne.0d0)) then
				write(*,'(a17,i2,a3,ES25.10)') 'temp growth(nz/3,',i,'):', &
					log(abs(temp(nz_vertical/3,i))) - log(abs(prev_temp(nz_vertical/3,i)))
			endif
			!write(*,'(ES25.10)',advance='no') 'temp growth(nz/3,:):', &
			!	abs(log(temp(nz_vertical/3,i)) - log(prev_temp(nz_vertical/3,i)))
		end do
		write(*,*)
		
		write(*,*)
	end subroutine output_vals
	
	!Calls a pgplot window to update with the new fluid quantities
	subroutine plot_vals
		integer :: i
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
		a1 = -5e1
		a2 = 5e1
		
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

		hl = (/0.0, 0.2, 0.4, 0.6, 1.0/)
		hr = (/0.0, 0.5, 1.0, 1.0, 1.0/)
		hg = (/0.0, 0.0, 0.5, 1.0, 1.0/) 
		hb = (/0.0, 0.0, 0.0, 0.3, 1.0/)
		call pgctab(hl, hr, hg, hb, 5, 1.0, 0.5)
		
		!call pgimag (real(temp_display), nz_vertical, nz_horizontal, i1, i2, j1, j2, a1, a2, tr)
		!call pgpage()
		call pgimag (real(omega_display), nz_vertical, nz_vertical, i1, i2, j1, j2, a1, a2, tr)
	end subroutine plot_vals
	
	subroutine shutdown
		call pgclos()
	end subroutine shutdown

end module hydro

program main

	use hydro
	
	implicit none
	
	integer :: i	!Loop index
	
	call initialize()
	!call pgplot_test()
	call output_vals()
	if(use_pgplot) call plot_vals()
	
	do i=1,num_steps
		write(*,*) 'Evolving step',i
		call evolve_step()
		if(use_pgplot) call plot_vals()
		!call output_vals()
	end do
	
	call shutdown()
	
end program