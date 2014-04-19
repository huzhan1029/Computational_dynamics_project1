subroutine ReadData()
!-------------------------------------------------
! Purpose: Read data from files specified by K_filename, M_filename, IniD_filename, etc.
!-------------------------------------------------

	use module_parameter
	use module_data
	use module_ioport
	implicit none
	logical :: InquireFile, pcomp
	integer :: i,j
	allocate( K_matrix(ndof,ndof), M_matrix(ndof,ndof), C_matrix(ndof,ndof))
	allocate( d0_vector(ndof), v0_vector(ndof), NodalForceId(ndof))
	
	write(oport_runmsg, 200)
200 format(/, '====== Begin to read data file ... ======') 
! -------Read Stiff Matrix from K_filename.---------------------- 
	if(InquireFile(K_filename)) then ! K_filename exists
		open(unit=iport_data,file=K_filename, status='old')
		read(iport_data,*) ((K_matrix(i,j),j=1,ndof),i=1,ndof)
		close(iport_data)
		write(oport_runmsg, 210)  
		write(oport_runmsg,211) ((K_matrix(i,j),j=1,ndof),i=1,ndof)
		
	endif
210 format(/,'Read data of stiff matrix...')
211	format(4X,'Stiff matrix is :', /,(4X,<ndof>F12.4))
! -------Read Mass Matrix from M_filename.---------------------
	M_matrix=0.0
	if(InquireFile(M_filename)) then ! K_filename exists
		open(unit=iport_data,file=M_filename, status='old')
		if(pcomp(M_type, 'cons',4)) then
			read(iport_data,*) ((M_matrix(i,j),j=1,ndof),i=1,ndof)
			write(oport_runmsg, 220)
			write(oport_runmsg,221) ((M_matrix(i,j),j=1,ndof),i=1,ndof)
		elseif(pcomp(M_type, 'lump',4)) then
			read(iport_data,*) (M_matrix(i,i),i=1,ndof)
			write(oport_runmsg, 220)
			write(oport_runmsg, 222) (M_matrix(i,i),i=1,ndof)
		else
			write(oport_runmsg, 223)
			stop '***Error : The type of Mass Matrix is wrong! Check input file!'
		end if
		close(iport_data)
	endif
220 format(/,'Read data of mass matrix...')
221 format(4X,'Consistent mass matrix is :', /,(4X,<ndof>F12.4))
222 format(4X, 'Lumped mass matrix is :', /, 4X, <ndof>f12.4)
223 format(/,'***Error : The type of Mass Matrix is wrong! Check input file!')
! -------Construct Damping Matrix.-----------------------------
	! If C_type is 'file', read Damping Matrix from C_filename;
	! If C_type is 'rayl', compute Damping Matrix from K_matrix and M_matrix.
	
	if(C_exist) then
		if(pcomp(C_type, 'file',4)) then
			if(InquireFile(C_filename)) then
				open(unit=iport_data,file=C_filename, status='old')
				read(iport_data,*) ((C_matrix(i,j),j=1,ndof),i=1,ndof)
				close(iport_data)
				write(oport_runmsg, 230)
				
			endif
		elseif(pcomp(C_type, 'rayl', 4)) then
			C_matrix = K_matrix*rayl_coef(1) + M_matrix*rayl_coef(2)
		else
			write(oport_runmsg, 232)
			stop '***Error : The type of Damping Matrix is wrong! Check input file!'
		endif
	else
		C_matrix = 0.0
	end if
	write(oport_runmsg,231) ((C_matrix(i,j),j=1,ndof),i=1,ndof)
230 format(/,'Read data of damping matrix...')
231 format(4X,'Damping matrix is :', /,(4X,<ndof>F12.4))
232 format(/,'***Error : The type of Damping Matrix is wrong! Check input file!')
! ------Read Initial Displacement.-----------------------------
	if(pcomp(IniD_type, 'file',4)) then
		if(InquireFile(IniD_filename)) then
			open(unit=iport_data,file=IniD_filename, status='old')
			read(iport_data,*) (d0_vector(i),i=1,ndof)
			close(iport_data)
			write(oport_runmsg, 240) 
			write(oport_runmsg, 241) (d0_vector(i),i=1,ndof)
		endif
	elseif(pcomp(IniD_type, 'zero',4) .or. pcomp(IniV_type, '    ',4)) then
		d0_vector = 0.0
	else
		stop '***Error : The type of u0 is wrong! Check input file!'
	endif
240 format(/,'Read data of initial displacement ...')
241 format(4X, 'Initial displacement is :',/, 4X, <ndof>f12.4)
! ------Read Initial Velocity.-----------------------------
	if(pcomp(IniV_type, 'file',4)) then
		if(InquireFile(IniV_filename)) then
			open(unit=iport_data,file=IniV_filename, status='old')
			read(iport_data,*) (v0_vector(i),i=1,ndof)
			close(iport_data)
			write(oport_runmsg, 250) 
			write(oport_runmsg, 251) (v0_vector(i),i=1,ndof)
		endif
	elseif(pcomp(IniV_type, 'zero',4) .or. pcomp(IniV_type, '    ',4)) then
		v0_vector = 0.0
	else
		stop '***Error : The type of v0 is wrong! Check input file!'
	endif
250 format(/,'Read data of initial velocity ...')
251 format(4X, 'Initial velocity is :',/, 4X, <ndof>f12.4)
! ------Read nodal force ID.-----------------------------

	if(pcomp(F_type, 'file',4)) then
		if(InquireFile(F_filename)) then
			open(unit=iport_data,file=F_filename, status='old')
			read(iport_data,*) (NodalForceId(i),i=1,ndof)
			close(iport_data)
			write(oport_runmsg, 260) 
			write(oport_runmsg, 261) (NodalForceId(i),i=1,ndof)
		endif
	elseif(pcomp(F_type, 'sing',4)) then
		NodalForceId=0.0
		NodalForceId(F_dof)=1.0
	else
		NodalForceId=0.0
	endif
260 format(/,'Read data of nodal force ID ...')
261 format(4X, 'Nodal force ID is :',/, 4X, <ndof>i4)
! --------------------------------------------------------------
	write(oport_runmsg, 201)
201 format(/, '------ End of reading data file ------') 
end subroutine ReadData