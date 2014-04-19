subroutine OpenOutFiles()
! Purpose: open files for outputs
	use module_ioport
	use module_parameter
	implicit none
	integer i

	if(Disp_flag) then
		do i=1,n_Disp_dof
			open(unit=oport_disp(i), file=Disp_filename(i), status='unknown')
		end do
	endif
		
	if(Velo_flag)   then
		do i=1,n_Velo_dof
			open(unit=oport_velo(i), file=Velo_filename(i), status='unknown')
		end do
	endif

	if(Acce_flag)	then
		do i=1,n_Acce_dof
			open(unit=oport_acce(i), file=Acce_filename(i), status='unknown')
		end do
	endif

end subroutine OpenOutFiles