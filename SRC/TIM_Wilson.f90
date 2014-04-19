subroutine TIM_Wilson()

! Purpose: Wilson-theta Method.
	use linear_operators
	use module_parameter
	use module_data
	use module_ioport
	
	implicit none
	integer i
	real :: c0, c1, c2, c3, c4, c5, c6, c7,c8,c9,c10, t                
	real ::	K_eff(1:ndof,1:ndof), invK(1:ndof,1:ndof),            &
            f_eff(1:ndof), f_t(1:ndof), f_l(1:ndof), d1(1:ndof),  &
            v1(1:ndof), a1(1:ndof), d2(1:ndof), v2(1:ndof), a2(1:ndof)
	real :: theta
	real PropForce

	if(TIM_para1 .eq. 0.0) then
		theta = 1.4  ! default value
	else
	    theta = TIM_para1
	end if
	write(oport_runmsg,301) theta
	c0=6.0/(theta**2*dt**2); c1=3.0/theta/dt; c2=2.0*c1; 
    c3=2.0; c4=2.0; c5=theta*dt/2.0
	c6=c0/theta; c7=-c2/theta; c8=1.0-3.0/theta; 
    c9=dt/2.0; c10=dt*dt/6.0

	K_eff = K_matrix + c0*M_matrix + c1*C_matrix
	invK = .i. K_eff

	d1=d0_vector; v1=v0_vector; a1=a0_vector
	t=0.0
	do i =1,nstep
		t=t+dt 
		f_t=PropForce(t)*NodalForceId
		f_l=PropForce(t-dt)*NodalForceId
		f_eff=f_l + theta*(f_t - f_l)                       & 
		          + (M_matrix .x. (c0*d1 + c2*v1 + c3*a1))  &
				  + (C_matrix .x. (c1*d1 + c4*v1 + c5*a1))
		d2=invK .x. f_eff  ! solve displacement at 't+theta*dt'
		a2=c6*(d2-d1) + c7*v1 + c8*a1
		v2=v1 + c9*(a2 +a1)
		d2=d1 + dt*v1 +c10*(a2 +2.0*a1)
		call WriteToFile(t,d2,'d')
		call WriteToFile(t,v2,'v')
		call WriteToFile(t,a2,'a')
		d1=d2; v1=v2; a1=a2
	end do
301 format(4X,'theta = ',f7.4,/,4X,  &
          '* Note : (1) theta < 1, method is unstable;',/,  &
      12X,' (2) 1 =< theta <1.3667, method is conditionally stable;',/, &
	  12X,' (3) theta >= 1.3667, method is unconditionally stable.')
end subroutine TIM_Wilson