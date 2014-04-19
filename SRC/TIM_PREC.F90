subroutine TIM_Precise()

! Purpose: Precise integration scheme
	use linear_operators
	use module_parameter
	use module_data
	use module_ioport
	
	implicit none
	integer i, ndof1, ndof2
	real :: t                
	real :: d2(1:ndof), v2(1:ndof), a2(1:ndof)
	real ::	H_matrix(1:2*ndof, 1:2*ndof)	
! matrix H = [A, D; B, G], also used as H^{-1}
	real :: Htau(1:2*ndof, 1:2*ndof), Htau2(1:2*ndof, 1:2*ndof)
! H * tau, (H * tau) ^ 2, these two matrices can be released
! as soon as the matrix T will have been calculated
	real :: f_t(1:ndof), f_tp1(1:ndof)
! force at time level t and time level t+1
	real :: T_matrix(1:2*ndof, 1:2*ndof)
! matrix T = exp(H * dt)
	real :: TmI_Hinv(1:2*ndof, 1:2*ndof)
! matrix (T - I) * H^{-1}, also used as T_a
	real :: r0(1:2*ndof), r1(1:2*ndof)
	real :: Minv(1:ndof, 1:ndof)
! M^{-1}
	real :: u_t(1:2*ndof), u_tp1(1:2*ndof)
! u = {d, M * v + 0.5 * C * d}
	real :: f2(1:2*ndof), f3(1:2*ndof)
! temporary storage used in integration
	real PropForce
	integer :: m, n
	real :: tau

	ndof1 = ndof + 1
	ndof2 = ndof * 2
	if(TIM_para1 .eq. 0.0) then
		m = 20
	else
		m = nint(TIM_para1)
		if (m < 0) then
			write(*,*) 'The parameter ', m, &
                       ' is not valid for precise integration!'
			stop
		end if
	endif	
	n = 2 ** m
	
	write(oport_runmsg, 301) m

! construct H
	Minv = .i. M_Matrix
	H_matrix(1:ndof, 1:ndof) = -0.5 * Minv .x. C_matrix
	H_matrix(1:ndof, ndof1:ndof2) = Minv
	H_matrix(ndof1:ndof2, 1:ndof) = ((0.25 * C_matrix .x. Minv)    &
          .x. C_matrix) - K_matrix
	H_matrix(ndof1:ndof2, ndof1:ndof2) = -0.5 * C_matrix .x. Minv

! H * tau and its square
	tau = dt / dble(n)
	Htau = H_matrix * tau
	Htau2 = Htau .x. Htau

! calculate T_a
	TmI_Hinv = 1.0/3.0 * Htau + 1.0/12.0 * Htau2
	do i = 1, ndof2
		TmI_Hinv(i,i) = TmI_Hinv(i,i) + 1.0
	end do
	TmI_Hinv = Htau + (Htau2 .x. TmI_Hinv * 0.5)
	do i = 1, m
		TmI_Hinv = 2.0 * TmI_Hinv + (TmI_Hinv .x. TmI_Hinv)
	end do

! calculate T
	T_matrix = TmI_Hinv
	do i = 1, ndof2
		T_matrix(i,i) = T_matrix(i,i) + 1.0
	end do

! inversion of H
	H_matrix = .i. H_matrix

! calculate (T - I) * H^{-1}
	TmI_Hinv = TmI_Hinv .x. H_matrix

! initial conditions
	t = 0.0
	f_t = PropForce(t) * NodalForceId
	u_t(1:ndof) = d0_vector
	u_t(ndof+1 : ndof*2) = (M_matrix .x. v0_vector) +    &
         (0.5 * C_matrix .x. d0_vector)
	r0 = 0.0
	r1 = 0.0

! All is ready except for east wind!
	do i = 1, nstep

		t = t+dt
		f_tp1 = PropForce(t) * NodalForceId
		r0(ndof1 : ndof2) = f_t
		r1(ndof1 : ndof2) = (f_tp1 - f_t) / dt	

		f3 = H_matrix .x. r1
		f2 = f3
		f3 = f3 * dt
		f2 = f2 + r0
		f2 = TmI_Hinv .x. f2

		u_tp1 = (T_matrix .x. u_t) + f2 - f3

! calculate displacement, velocity and acceleration
		d2 = u_tp1(1:ndof)
		v2 = (Minv .x. (u_tp1(ndof+1 : ndof*2)) -      &
             (0.5 * C_matrix .x. d2))  
		a2 = (Minv .x. (f_tp1 - (C_matrix .x. v2)) -   &
             (K_matrix .x. d2))
		call WriteToFile(t,d2,'d')
		call WriteToFile(t,v2,'v')
		call WriteToFile(t,a2,'a')

! update variables
		f_t = f_tp1
		u_t = u_tp1

	end do

301 format(4X,'m = ', I5)

end subroutine TIM_Precise