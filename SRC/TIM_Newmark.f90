subroutine TIM_Newmark()

! Purpose: Newmark Method.
	use linear_operators
	use module_parameter
	use module_data
	use module_ioport
	
	implicit none
	integer i
	real :: c0, c1, c2, c3, c4, c5, c6, c7,t                
	real ::	K_eff(1:ndof,1:ndof), invK(1:ndof,1:ndof),            &
            f_eff(1:ndof), f_t(1:ndof), d1(1:ndof), v1(1:ndof),   &
            a1(1:ndof),d2(1:ndof), v2(1:ndof), a2(1:ndof)
	real :: gamma,beta
	real PropForce

	if(TIM_para1 .eq. 0.0) then
		gamma = 0.5
	else
		gamma=TIM_para1
	endif
	if(TIM_para2 .eq. 0.0) then
		beta = 0.25
	else
		beta=TIM_para2
	endif
	write(oport_runmsg,301) gamma, beta

	c0=1.0/(beta*dt*dt); c1=gamma/beta/dt; 
    c2=1.0/beta/dt; c3=0.5/beta-1.0
	c4=gamma/beta-1.0; c5=dt*(gamma/2.0/beta-1.0); 
    c6=dt*(1-gamma); c7=gamma*dt

	K_eff = K_matrix + c0*M_matrix + c1*C_matrix
	invK = .i. K_eff

	d1=d0_vector; v1=v0_vector; a1=a0_vector
	t=0.0
	do i=1,nstep
		t=t+dt
		f_t=PropForce(t)*NodalForceId
		f_eff = f_t + (M_matrix .x. (c0*d1 + c2*v1 +c3*a1))   &
		            + (C_matrix .x. (c1*d1 + c4*v1 +c5*a1))
		d2=invK .x. f_eff
		a2=c0*(d2-d1) - c2*v1 - c3*a1
		v2=v1 + c6*a1 +c7*a2
		call WriteToFile(t,d2,'d')
		call WriteToFile(t,v2,'v')
		call WriteToFile(t,a2,'a')
		d1=d2; v1=v2; a1=a2
	end do

301 format(4X,'gamma = ',f7.4,8X,'beta = ',f7.4)

end subroutine TIM_Newmark