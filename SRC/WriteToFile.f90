subroutine WriteToFile(a,b,c)
!-----------------------------------------------------------------------
! Purpose : write time and value to files
! Inputs:  a - time
!          b - vector of displacement, velocity or acceleration
!          c - 'd' means b is displacement; 'v' means b is velocity; 'a' means acceleration.
!-----------------------------------------------------------------------

	use module_parameter
	use module_ioport
	implicit none
	real :: a,b(*)
	integer i
	character :: c
	logical pcomp

! Write value of displacement to Disp_filename.
	if(Disp_flag .and. pcomp(c,'d',1)) then
		do i=1,n_Disp_dof
			write(oport_disp(i), 201) a,b(Disp_dof(i))
		end do
	end if 

! Write value of displacement to Velo_filename.		
	if(Velo_flag .and. pcomp(c,'v',1)) then
		do i=1,n_Velo_dof
			write(oport_velo(i), 201) a,b(Velo_dof(i))
		end do
	end if 

! Write value of displacement to Disp_filename.
	if(Acce_flag .and. pcomp(c,'a',1)) then
		do i=1,n_Acce_dof
			write(oport_acce(i), 201) a,b(Acce_dof(i))
		end do
	end if 
201 format(f12.6, 4X, f12.6)

end subroutine WriteToFile