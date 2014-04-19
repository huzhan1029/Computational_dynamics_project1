subroutine TIM_GW4()

! Purpose: GW3 or GW3ucs method.
	use linear_operators
	use module_parameter
	use module_data
	use module_ioport
	
	implicit none
	integer i,j,ii
	real :: Lk(4,4), Lc(4,4), Lm(4,4)
	real :: Ke(ndof,ndof,4,4), fe1(ndof), fe2(ndof),      & 
            fe3(ndof), fe4(ndof)
	real :: S1(ndof,ndof),S2(ndof,ndof),S3(ndof,ndof),    &
            S4(ndof,ndof),S5(ndof,ndof),S6(ndof,ndof),    &
            S7(ndof,ndof),S8(ndof,ndof),T1(ndof,ndof),    &
            T2(ndof,ndof),T3(ndof,ndof),T4(ndof,ndof),    &
			T5(ndof,ndof),T6(ndof,ndof),T7(ndof,ndof),    &
            T8(ndof,ndof),DD1(ndof,ndof),DD2(ndof,ndof),  &
            DD3(ndof,ndof),DD4(ndof,ndof)
	real :: t, f1(ndof), f2(ndof), fm1(ndof), fm2(ndof) 
	real :: d1(ndof), d2(ndof), v1(ndof), v2(ndof), a2(ndof)
	real :: invKe22(ndof,ndof), invKe33(ndof,ndof),   &
            invDD2(ndof,ndof), invM(ndof,ndof)
	real PropForce
	logical pcomp

! f1 - nodal force at last time (t_n); 
! f2 - nodal force at current time (t_n+1);
! fm1 - nodal force at time (t_n+1/3);
! fm2 - nodal force at time (t_n+2/3);
! d1 - displacement at last time (t_n);
! d2 - displacement at current time (t_n+1);
! v1 - velocity at last time (t_n);
! v2 - velocity at current time (t_n+1);
! a2 - acceleration at current time (t_n+1)

! conditionally stable method - GW4
	if(pcomp(TIM_type, 'GW4 ', 4)) then 
		Lk = reshape(source=(/128, 99, -36, 19,   &
							  99, 648, -81, -36,  &
							  -36, -81, 648, 99,  &
							  19, -36, 99, 128/), &
							  shape=(/4,4/), order=(/2,1/))
		Lk = Lk/840.0
		write(oport_runmsg,101) 
! unconditionally stable method - GW4ucs
	elseif(pcomp(TIM_type, 'GW4u', 4)) then  
		Lk = reshape(source=(/331, 387, -207, 89,   &
							  387, 1539, 81, -207,  &
							  -207, 81, 1539, 387,  &
							  89, -207, 387, 331/), &
							  shape=(/4,4/), order=(/2,1/))
		Lk = Lk/2400.0
		write(oport_runmsg,102)
	else
		write(oport_runmsg,104) TIM_type
		stop '***Error : type of TIM is wrong !'
	endif
	Lc = reshape(source=(/-40, 57, -24, 7,   &
						  -57, 0, 81, -24,   &
						  24, -81, 0, 57,    &
						  -7, 24, -57, 40/), &
						  shape=(/4,4/), order=(/2,1/))
	Lc = Lc/80.0
	Lm = reshape(source=(/148, -189, 54, -13,   &
						  -189, 432, -297, 54,  &
						  54, -297, 432, -189,  &
						  -13, 54, -189, 148/), &
						  shape=(/4,4/), order=(/2,1/))
	Lm = Lm/80.0
	write(oport_runmsg,103) 'Lk', ((Lk(i,j),j=1,4),i=1,4)
	write(oport_runmsg,103) 'Lc', ((Lc(i,j),j=1,4),i=1,4)
	write(oport_runmsg,103) 'Lm', ((Lm(i,j),j=1,4),i=1,4)
	do j=1,4
		do i=1,4
			Ke(:,:,i,j) = Lk(i,j)*K_matrix +            &
                          (2.0/dt*Lc(i,j))*C_matrix -   &
						  (4.0/dt/dt*Lm(i,j))*M_matrix
		end do
	end do
	
	invKe22 = .i. Ke(:,:,2,2)
	S1 = Ke(:,:,3,2) .x. invKe22
	S2 = .i. (Ke(:,:,3,3) - (S1 .x. Ke(:,:,2,3)))
	S3 = Ke(:,:,3,1) - (S1 .x. Ke(:,:,2,1))
	S4 = Ke(:,:,3,4) - (S1 .x. Ke(:,:,2,4))
	invKe33 = .i. Ke(:,:,3,3)
	S5 = Ke(:,:,2,3) .x. invKe33
	S6 = .i. (Ke(:,:,2,2) - (S5 .x. Ke(:,:,3,2)))
	S7 = Ke(:,:,2,1) - (S5 .x. Ke(:,:,3,1))
	S8 = Ke(:,:,2,4) - (S5 .x. Ke(:,:,3,4))

	T1 = Ke(:,:,1,2) .x. S6
	T2 = Ke(:,:,1,3) .x. S2
	T3 = T1 .x. S5
	T4 = T2 .x. S1
	DD1 = Ke(:,:,1,1) - (T1 .x. S7) - (T2 .x. S3)
	DD2 = Ke(:,:,1,4) - (T1 .x. S8) - (T2 .x. S4)
	T5 = Ke(:,:,4,2) .x. S6
	T6 = Ke(:,:,4,3) .x. S2
	T7 = T5 .x. S5
	T8 = T6 .x. S1
	DD3 = Ke(:,:,4,1) - (T5 .x. S7) - (T6 .x. S3)
	DD4 = Ke(:,:,4,4) - (T5 .x. S8) - (T6 .x. S4)

	invDD2 = .i. DD2
	invM = .i. M_matrix 
	d1 = d0_vector; v1 = v0_vector
	t=0.0
	f1 = NodalForceId*PropForce(t)  ! nodal force at time 0
	do ii=1,nstep
		t =t+dt  ! current time
		fm1 = NodalForceId*PropForce(t-dt*2.0/3.0)
		fm2 = NodalForceId*PropForce(t-dt/3.0)
		f2 = NodalForceId*PropForce(t)
		fe1 = Lk(1,1)*f1 + Lk(1,2)*fm1 + Lk(1,3)*fm2 + Lk(1,4)*f2
		fe2 = Lk(2,1)*f1 + Lk(2,2)*fm1 + Lk(2,3)*fm2 + Lk(2,4)*f2
		fe3 = Lk(3,1)*f1 + Lk(3,2)*fm1 + Lk(3,3)*fm2 + Lk(3,4)*f2
		fe4 = Lk(4,1)*f1 + Lk(4,2)*fm1 + Lk(4,3)*fm2 + Lk(4,4)*f2

		d2 = invDD2 .x. (fe1 + ((T4-T1).x.fe2) + ((T3-T2).x.fe3) + &
				2.0/dt*(M_matrix .x. v1) - (DD1 .x. d1))
		v2 = invM .x. (fe4 + ((T8-T5).x.fe2) + ((T7-T6).x.fe3) -   &
		        (DD3.x.d1) - (DD4.x.d2))
		v2 = v2*dt/2.0
		a2 = invM .x. (f2 - (K_matrix.x.d2) - (C_matrix.x.v2))

		call WriteToFile(t,d2,'d')
		call WriteToFile(t,v2,'v')
		call WriteToFile(t,a2,'a')
		d1=d2; v1=v2; f1=f2
	end do

101 format(4X,'Conditionally stable method : GW4')
102 format(4X,'Un-Conditionally stable method : GW4ucs')
103 format(4X,a2, ' is : --------------------------- ', /,(6X, 4f12.4))
104 format(4X,'***Error : type of TIM is wrong :', a)
end subroutine TIM_GW4