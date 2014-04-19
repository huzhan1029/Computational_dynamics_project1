subroutine TimeIntegration ()


	use linear_operators
	use module_parameter
	use module_data
	use module_ioport
	implicit none
	real :: f0(1:ndof)
	real PropForce
	logical pcomp
	allocate (a0_vector(ndof))

	call OpenOutFiles

	write(oport_runmsg, 302) 

! calculate initial acceleration
	f0 = PropForce(0.0)*NodalForceId
	a0_vector = M_matrix .ix. (f0 - (K_matrix .x. d0_vector)  - (C_matrix .x. v0_vector))
	call WriteToFile(0.0, d0_vector,'d') ! write d(0) to file
	call WriteToFile(0.0, v0_vector,'v') ! write v(0) to file
	call WriteToFile(0.0, a0_vector,'a') ! write a(0) to file
	nstep= int(TotalTime/dt) + 1  ! number of time steps
	write(oport_runmsg, 301) nstep

	if(pcomp(TIM_type, 'cent',4))  then ! Central Difference Method
		write(oport_runmsg,300) 'Central Difference'
		call TIM_Central
	elseif(pcomp(TIM_type, 'houb',4)) then ! Houbolt Method
		write(oport_runmsg,300) 'Houbolt'
		call TIM_Houbolt
	elseif(pcomp(TIM_type, 'newm',4)) then ! Newmark Method
		write(oport_runmsg,300) 'Newmark'
		call TIM_Newmark
	elseif(pcomp(TIM_type, 'wils',4)) then ! Wilson-theta Method
		write(oport_runmsg,300) 'Wilson-theta'
		call TIM_Wilson
	elseif(pcomp(TIM_type, 'GW3',3)) then ! GW3 or GW3ucs
		write(oport_runmsg,300) 'GW3 or GW3ucs'
		call TIM_GW3
	elseif(pcomp(TIM_type, 'GW4',3)) then ! GW4 or GW4ucs
		write(oport_runmsg,300) 'GW4 or GW4ucs'
		call TIM_GW4
! ----- 1120 -----
	elseif(pcomp(TIM_type, 'prec',4)) then ! PRECise integration
		write(oport_runmsg,300) 'Precise Integration'
		call TIM_Precise
! ----- end 1120 -----
	else
		write(oport_runmsg, 304)
		stop '***Error : type fo Time Integration Method is wrong !'
	end if

	call CloseOutFiles
	write(oport_runmsg, 303) 

300 format(/,'Run ', a ,' method for time integration...')
301 format(/,'Total number of time intervals is : ', i)
302 format(/,'====== Begin to compute ... ======')
303 format(/,'------ End of computation ------')
304 format(/,'***Error : type fo Time Integration Method is wrong !')
end subroutine TimeIntegration