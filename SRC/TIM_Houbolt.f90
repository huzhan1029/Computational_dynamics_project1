subroutine TIM_Houbolt()

! Purpose: Houbolt Method.
	use linear_operators
	use module_parameter
	use module_data
	use module_ioport
	
	implicit none
	integer i
	real :: c0, c1, c2, c3, c4, c5, c6, c7, cc0,cc1,cc2,cc3,t 
	real :: K_eff(ndof,ndof), invK(ndof,ndof),M_eff(ndof,ndof), &
            invM(ndof,ndof), f_eff(1:ndof), f_t(1:ndof), &
	        d1(1:ndof), d2(1:ndof), d3(1:ndof), d4(1:ndof),   &
			v2(1:ndof), a2(1:ndof)
	real PropForce
	
	c0=2.0/dt/dt; c1=11.0/6.0/dt; c2=5.0/dt/dt; c3= 3.0/dt
	c4=-2.0*c0;	c5=-c3/2.0; c6=c0/2.0; c7=c3/9.0

	write(oport_runmsg, 301) c0, c1, c2, c3, c4, c5, c6, c7

	K_eff=K_matrix + c0*M_matrix + c1*C_matrix
	invK= .i. K_eff

! Start up by using Central Difference Method.
	cc0=1.0/dt/dt; cc1=0.5/dt; cc2=2.0*cc0; cc3=1.0/cc2
	d1=d0_vector - dt*v0_vector +cc3*a0_vector ! displacement at -dt.
	M_eff = cc0*M_matrix + cc1*C_matrix  ! Form the equivalent mass matrix
	invM = .i. M_eff  ! Compute the inverse of the equivalent mass matrix
	d2=d0_vector  ! displacement at 0
	t=0.0
	do i=1,2
		f_t = PropForce(t)*NodalForceId  ! nodal force at time 0
! equivalent force at time t
		f_eff = f_t - ((K_matrix - cc2*M_matrix) .x. d2)       &
					- ((cc0*M_matrix - cc1*C_matrix) .x. d1)   
! Solve displacement at time t+dt
		d3 = invM .x. f_eff  
		select case (i)
		case(1) 
			d1=d2
			d2=d3
		case(2)
			a2= cc0*(d1 - 2*d2 +d3) ! acceleration at time t.
			v2= cc1*(d3-d1)         ! velocity at time t.
			call WriteToFile(t, d2,'d') ! write d(dt) to file
			call WriteToFile(t, v2,'v') ! write v(dt) to file
			call WriteToFile(t, a2,'a') ! write a(dt) to file
		end select
		t=t+dt
	end do
! now, d(dt), v(dt), a(dt), d(2*dt) are known.
!      d(dt), v(dt), a(dt) have been written to files.

! Using Houbolt Method.
	do i=3,nstep
		t=t+dt
		f_t = PropForce(t)*NodalForceId
		f_eff = f_t + (M_matrix .x. (c2*d3 + c4*d2 + c6*d1))  &
				    + (C_matrix .x. (c3*d3 + c5*d2 + c7*d1))
		d4 = invK .x. f_eff

		if (i .eq. 3) then
			a2= cc0*(d1 - 2*d2 +d3) ! acceleration at time 2*dt.
			v2= cc1*(d3-d1)         ! velocity at time 2*t.
			call WriteToFile(2*dt, d3,'d') ! write d(2*dt) to file
			call WriteToFile(2*dt, v2,'v') ! write v(2*dt) to file
			call WriteToFile(2*dt, a2,'a') ! write a(2*dt) to file
		endif
		a2 = c0*d4 - c2*d3 - c4*d2 - c6*d1
		v2 = c1*d4 - c3*d3 - c5*d2 - c7*d1
		call WriteToFile(t, d4,'d') ! write d(t) to file
		call WriteToFile(t, v2,'v') ! write v(t) to file
		call WriteToFile(t, a2,'a') ! write a(t) to file
		d1=d2
		d2=d3
		d3=d4
	end do
	
301 format(4X,'Integration constants :', /, 6X, 'c0=', f10.4, 2X,  &
              ';  c1=', f10.4, 2X, ';  c2=', f10.4, 2X, ';  c3=',  &
               f10.4, /, 6X, 'c4=', f10.4, 2X, ';  c5=',           &
  	           f10.4,2X,';  c6=',f10.4,2X,';  c7=',f10.4)				
end subroutine TIM_Houbolt