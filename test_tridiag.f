!Code to test Gary's tridiagonal matrix solver and compare to the ones in LAPACK (3.5.0)
!-----------------------
!Written by: Kevin Moore
!6/6/14

module tridiag

	implicit none

	contains
	
	!Generates a tridiagonal matrix with random values, but dominant diagonal entries
	subroutine gen_tridiag_matrix(nz,sub,dia,sup)
		integer, intent(in) :: nz
		double precision, dimension (1:nz), intent(inout) :: sub,dia,sup
		double precision :: x
		integer :: i
		
		do i=1,nz
			call random_number(x)
			sub(i) = x
			call random_number(x)
			sup(i) = x
			call random_number(x)
			dia(i) = 5*x + 1
		end do
	end subroutine gen_tridiag_matrix

	!Gary's solver from Appendix A of his book:
	!This is a tridiagonal matrix solver written in Fortran based on LINPACK and LAPACK 
	!routines. Note, in this routine the first index of the arrays is 1 and the last is Nz; 
	!if one wishes to start the arrays with index 0 and end with Nz − 1, all index 
	!references in this routine would need to be reduced by 1.
	subroutine tridi(nz,rhs,sol,sub,dia,sup,wk1,wk2)
		double precision, dimension (1:nz) :: rhs,sol,sub,dia,sup,wk1,wk2 
		integer :: i,nz
		
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

end module tridiag

program test_tridiag
	
	use tridiag
	
	implicit none
	
	integer, parameter :: nz = 1000000	!Going to 1e7 overflows the stack hard limit (~64M?)
	logical, parameter :: output_matrix = .false.
	logical, parameter :: output_results = .false.
	integer :: i,j,isize,info,nrhs,ldb
	integer :: sys_time_begin, sys_time_end, clock_rate	!Timing variables
	integer, dimension(:), pointer :: itime
	double precision, dimension(:), pointer :: rhs,sol,sub,dia,sup,wk1,wk2,resid,&
		dia_lapack, rhs_lapack, sol_lapack, resid_lapack
	double precision, dimension(:), pointer :: du, dl
	
	write(*,*) 'Testing Tridiagonal solvers!'
	write(*,*)
	
	allocate(rhs(nz), sol(nz), sub(nz), dia(nz), sup(nz), wk1(nz), wk2(nz), resid(nz), &
		dia_lapack(nz), rhs_lapack(nz), sol_lapack(nz), resid_lapack(nz), &
		du(nz-1), dl(nz-1))
	
	call random_seed(size=isize)
	allocate(itime(isize))
	itime = time()
	!write(*,*) itime
	call random_seed(put=itime)
	deallocate(itime)
	
	rhs(:) = 1d0
	call random_number(rhs)
	call gen_tridiag_matrix(nz,sub,dia,sup)
	
	if(output_matrix) then
		write(*,*) "Here's our matrix:"
		do i=1,nz
			do j=1,i-2
				write(*,'(a12)',advance='no') ''
			end do
			if((i-1).ge.1) write(*,'(ES12.4)',advance='no') sub(i-1)
			write(*,'(ES12.4)',advance='no') dia(i)
			if((i+1).le.nz) write(*,'(ES12.4)',advance='no') sup(i+1)
			write(*,*)
		end do
	endif
	
	!Call the solver from Gary's book:
	call system_clock(sys_time_begin,clock_rate)
	call tridi(nz,rhs,sol,sub,dia,sup,wk1,wk2)
	call system_clock(sys_time_end,clock_rate)
	write(*,*) 'tridi computation time:', (sys_time_end-sys_time_begin)*1d0/clock_rate
	write(*,*)
	
	!Call LAPACK's tridiagonal solver:
	nrhs = 1
	dl = sub(2:nz)
	du = sup(1:nz-1)
	ldb = nz
	dia_lapack = dia
	rhs_lapack = rhs
	call system_clock(sys_time_begin,clock_rate)
	call DGTSV(nz,nrhs,dl,dia_lapack,du,rhs_lapack,ldb,info)
	call system_clock(sys_time_end,clock_rate)
	write(*,*) 'LAPACK computation time:', (sys_time_end-sys_time_begin)*1d0/clock_rate
	write(*,*)
	sol_lapack = rhs_lapack
	
	if(output_results) then
		write(*,'(4a15)') 'Solution', 'Resid', 'LAPACK Resid', 'RHS'
		do i=1,nz
			do j=1,nz
				if(i == 1) then
					resid(i) = rhs(i) - (dia(i)*sol(i) + sup(i)*sol(i+1))
					resid_lapack(i) = rhs(i) - (dia(i)*sol_lapack(i) + sup(i)*sol_lapack(i+1))
				else if(i == nz) then
					resid(i) = rhs(i) - (sub(i)*sol(i-1) + dia(i)*sol(i))
					resid_lapack(i) = rhs(i) - (sub(i)*sol_lapack(i-1) + dia(i)*sol_lapack(i))
				else
					resid(i) = rhs(i) - (sub(i)*sol(i-1) + dia(i)*sol(i) + sup(i)*sol(i+1))
					resid_lapack(i) = rhs(i) - (sub(i)*sol_lapack(i-1) + &
						dia(i)*sol_lapack(i) + sup(i)*sol_lapack(i+1))
				endif
			end do
			write(*,'(4ES15.4)') sol(i), resid(i), resid_lapack(i), rhs(i)
		end do
	endif
	write(*,*)
	
	deallocate(rhs, sol, sub, dia, sup, wk1, wk2, resid, dia_lapack, rhs_lapack, &
		sol_lapack, resid_lapack, du, dl)
end program test_tridiag