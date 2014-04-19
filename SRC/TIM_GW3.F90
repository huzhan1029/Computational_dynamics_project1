subroutine TIM_GW3()

! Purpose: GW3 or GW3ucs method.
	use linear_operators
	use module_parameter
	use module_data
	use module_ioport
	
	implicit none
	integer i,j,ii
	real :: Lk(3,3), Lc(3,3), Lm(3,3)
	real :: Ke(ndof,ndof,3,3), fe1(ndof), fe2(ndof), fe3(ndof)
	real :: Q1(ndof,ndof),Q2(ndof,ndof),R1(ndof,ndof),  &
            R2(ndof,ndof),R3(ndof,ndof),R4(ndof,ndof)
	real :: t, f1(ndof), f2(ndof), fm(ndof) 
	real :: d1(ndof), d2(ndof), v1(ndof), v2(ndof), a2(ndof)
	real :: invKe22(ndof,ndof), invR1(ndof,ndof), invM(ndof,ndof)
	real PropForce
	logical pcomp

! f1 - nodal force at last time (t_n); 
! f2 - nodal force at current time (t_n+1);
! fm - nodal force at middle time (t_n+1/2)
! d1 - displacement at last time (t_n);
! d2 - displacement at current time (t_n+1);
! v1 - velocity at last time (t_n);
! v2 - velocity at current time (t_n+1);
! a2 - acceleration at current time (t_n+1)

! conditionally stable method - GW3
	if(pcomp(TIM_type, 'GW3 ', 4)) then 
		Lk = reshape(source=(/4,2,-1,2,16,2,-1,2,4/),   &
                     shape=(/3,3/), order=(/2,1/))
		Lk = Lk/15.0
		write(oport_runmsg,101) 
! unconditionally stable method - GW3ucs
	elseif(pcomp(TIM_type, 'GW3u', 4)) then  
		Lk = reshape(source=(/2,2,-1,2,8,2,-1,2,2/),   &
                     shape=(/3,3/), order=(/2,1/))
		Lk = Lk/9.0
		write(oport_runmsg,102)
	else
		write(oport_runmsg,104) TIM_type
		stop '***Error : type of TIM is wrong !'
	endif
	Lc = reshape(source=(/-3,4,-1,-4,0,4,1,-4,3/),   &
                 shape=(/3,3/), order=(/2,1/))
	Lc = Lc/6.0
	Lm = reshape(source=(/7,-8,1,-8,16,-8,1,-8,7/),  &
                 shape=(/3,3/), order=(/2,1/))
	Lm = Lm/6.0
	write(oport_runmsg,103) 'Lk', ((Lk(i,j),j=1,3),i=1,3)
	write(oport_runmsg,103) 'Lc', ((Lc(i,j),j=1,3),i=1,3)
	write(oport_runmsg,103) 'Lm', ((Lm(i,j),j=1,3),i=1,3)
	do j=1,3
		do i=1,3
			Ke(:,:,i,j) = Lk(i,j)*K_matrix +             &
                          (2.0/dt*Lc(i,j))*C_matrix -    &
						  (4.0/dt/dt*Lm(i,j))*M_matrix
		end do
	end do
	
	invKe22 = .i. Ke(:,:,2,2)
	Q1 = Ke(:,:,1,2) .x. invKe22
	R1 = Ke(:,:,1,3) - (Q1 .x. Ke(:,:,2,3))
	R2 = Ke(:,:,1,1) - (Q1 .x. Ke(:,:,2,1))
	Q2 = Ke(:,:,3,2) .x. invKe22
	R3 = Ke(:,:,3,1) - (Q2 .x. Ke(:,:,2,1))
	R4 = Ke(:,:,3,3) - (Q2 .x. Ke(:,:,2,3))
	invR1 = .i. R1
	invM = .i. M_matrix 
	d1 = d0_vector; v1 = v0_vector
	t=0.0
	f1 = NodalForceId*PropForce(t)  ! nodal force at time 0
	do ii=1,nstep
		t =t+dt  ! current time
		fm = NodalForceId*PropForce(t-dt/2)
		f2 = NodalForceId*PropForce(t)
		fe1 = Lk(1,1)*f1 + Lk(1,2)*fm + Lk(1,3)*f2
		fe2 = Lk(2,1)*f1 + Lk(2,2)*fm + Lk(2,3)*f2
		fe3 = Lk(3,1)*f1 + Lk(3,2)*fm + Lk(3,3)*f2

		d2 = invR1 .x. (fe1 - (Q1.x.fe2) + 2.0/dt*(M_matrix.x.v1) -  &
             (R2.x.d1))
		v2 = invM .x. (fe3 - (Q2.x.fe2) - (R3.x.d1) - (R4.x.d2))
		v2 = v2*dt/2.0
		a2 = invM .x. (f2 - (K_matrix.x.d2) - (C_matrix.x.v2))

		call WriteToFile(t,d2,'d')
		call WriteToFile(t,v2,'v')
		call WriteToFile(t,a2,'a')
		d1=d2; v1=v2; f1=f2
	end do

101 format(4X,'Conditionally stable method : GW3')
102 format(4X,'Un-Conditionally stable method : GW3ucs')
103 format(4X,a2, ' is : --------------------------- ', /,(6X, 3f12.4))
104 format(4X,'***Error : type of TIM is wrong :', a)
end subroutine TIM_GW3