subroutine TIM_Central()

! Purpose: Central Difference Method.
	use linear_operators
	use module_parameter
	use module_data
	
	implicit none
	integer i
	real :: c0, c1, c2, c3, t, M_eff(1:ndof,1:ndof), f_eff(1:ndof), &
            f_t(1:ndof), d1(1:ndof), d2(1:ndof), d3(1:ndof), v2(1:ndof), &
            a2(1:ndof), invM(ndof,ndof)
	real PropForce

! Compute integral constants.
	c0=1/dt/dt; c1=0.5/dt; c2=2*c0; c3=1/c2;

! Compute displacement at -dt.
	d1=d0_vector - dt*v0_vector +c3*a0_vector
! Form the equivalent mass matrix
	M_eff = c0*M_matrix + c1*C_matrix  
! Compute the inverse of the equivalent mass matrix
	invM = .i. M_eff  

	d2=d0_vector
	t=0.0
	do i=1,nstep
		f_t = PropForce(t)*NodalForceId  ! force at time t
! Equivalent force at time t
		f_eff = f_t - ((K_matrix - c2*M_matrix) .x. d2)    &
					- ((c0*M_matrix - c1*C_matrix) .x. d1)    
		d3 = invM .x. f_eff  ! Solve displacement at time t+dt

		if(i .gt. 1) then
			a2 = c0*(d1 - 2*d2 +d3) ! acceleration at time t.
			v2 = c1*(d3-d1)         ! velocity at time t.
			call WriteToFile(t,d2,'d')
			call WriteToFile(t,v2,'v')
			call WriteToFile(t,a2,'a')
		endif
		d1=d2
		d2=d3
		t=t+dt
	end do

! output dispalcement at the last time node.
	call WriteToFile(t+dt,d2,'d')  

end subroutine TIM_Central