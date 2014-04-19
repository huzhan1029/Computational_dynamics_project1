subroutine ReadPara(infilename)
!------------------------------------------------------------------
!	purpose: Read macro commands and parameters for each command.
!            All parameters are stored in module_parameter
!	Input:
!            infilename - name of the marco commands input file		 
!------------------------------------------------------------------
	use module_parameter
	use module_ioport
	implicit none
    character :: infilename*(*)
	character :: macro*4, string*255, string_stripped*255, c*1
	logical :: exist_infile, end_of_file
    integer :: i,j,k,ls
	logical cksep, pcomp
	integer MacroNumber 
!   ls - length of string after stripped

    write(oport_runmsg, 100)	  
    open(unit=iport_macro, file=infilename, status='old')
	end_of_file= .false.
	do while(.not. end_of_file)
!       Read one line from the marco input file.
		read(iport_macro, '(a)') string  
	    macro(1:4) = ' '

!       Strip horizontal tab character (Ctrl-I = ASCII Character 9)
        do i = 1,255
            if(ichar(string(i:i)).eq. 9) string(i:i) = ' '
        end do

!       Strip superfluous blanks and comments
        call pstrip(string_stripped,string, 1,ls)

		if(ls .le. 0) cycle ! Comment or blank line

!       Read macro command from stripped string. 
		i = 1
		do while (.not. cksep(string_stripped(i:i)) .and. i .le. ls)
            i = i+1
		end do
			  
	    macro(1:i-1)=string_stripped(1:i-1)
		write(*,101) macro
		write(oport_runmsg,101) macro

		select case ( MacroNumber(macro) )

			case (1) ! 'end '
				end_of_file = .true.

			case (2) ! 'k    '
				call ReadStiffPara

			case (3) ! 'm    '
			    call ReadMassPara
	        
			case (4) ! 'c    '
			    call ReadDampPara

			case (5) ! 'inid' - initial displacement
				call ReadInitialPara(IniD_type, IniD_filename)

			case (6) ! 'iniv' - initial velocity
				call ReadInitialPara(IniV_type, IniV_filename)

			case (7) ! 'dt  ' 
				read(string_stripped(i+1:ls), *) dt
				write(oport_runmsg, *) '    dt = ', dt

			case (8) ! 'time'
				read(string_stripped(i+1:ls), *) TotalTime
				write(oport_runmsg, *) '    time = ', TotalTime

			case (9) ! 'disp' - parameters of displacement outputs
				call ReadOutputPara(Disp_flag,n_Disp_dof,Disp_dof,Disp_filename)

			case (10) ! 'velo' - parameters of velocity outputs
				call ReadOutputPara(Velo_flag,n_Velo_dof,Velo_dof,Velo_filename)

			case (11) ! 'acce' - parameters of acceleration outputs
				call ReadOutputPara(Acce_flag,n_Acce_dof,Acce_dof,Acce_filename)

			case (12) ! 'ndof' - total number of degree of freedom (in space)
				read(string_stripped(i+1:ls), '(i)') ndof
				write(oport_runmsg,*) '    ndof = ', ndof
 
 			case (13) ! 'meth' - Time Integration method
				call ReadMethodPara

			case (14) ! 'fpro' - parameters of propotional force
				call ReadPropPara

			case (15) ! 'forc' - Read nodal force information 
				call ReadNodalForce      

			case default
				write(*,102) macro
				write(oport_runmsg,102) macro
		
		end select

    end do ! while(.not. exist_infile)
    close(iport_macro)
	write(oport_runmsg,103)	  
100   format(/, '======Begin to read marco commands and parameters======')
101   format(/,'Current macro command is "', a4,'"')
102   format(/,'***Error : "',a4,'" is a wrong macro command ! Check input file!')
103   format(/, '------End of reading marco commands and parameters------')
CONTAINS

! ---------------------------------------------------------------------
	subroutine ReadStiffPara()
		K_filename(1:ls-i) = string_stripped(i+1:ls)
		write(oport_runmsg,131) K_filename
131		format(4X,'Stiff matrix in file --', a)
	end subroutine ReadStiffPara
! ---------------------------------------------------------------------
	subroutine ReadMassPara()
		logical pcomp
		M_type = string_stripped(i+1:i+4)
		write(oport_runmsg,132) M_type
		do j= i+1, ls
		    if( cksep(string_stripped(j:j)) ) then
				M_filename(1:ls-j)=string_stripped(j+1:ls)
				write(oport_runmsg,133) M_filename(1:ls-j)
				return
			endif
		end do
132		format(4X,'Type of mass matrix is :', a4)
133		format(4X,'Mass matrix in file --', a)
	end subroutine ReadMassPara
! ---------------------------------------------------------------------
	subroutine ReadDampPara()
		logical pcomp
		C_exist=.true. ! Damp exists.
		C_type = string_stripped(i+1:i+4)
		write(oport_runmsg,107) C_type
		if(pcomp(C_type,'rayl',4)) then
		    do j=i+1,ls
				if( cksep(string_stripped(j:j)) ) then
					do k=j+1,ls
						if( cksep(string_stripped(k:k)) ) then
		                      read(string_stripped(j+1:k-1), *) rayl_coef(1)
							  read(string_stripped(k+1:ls), *) rayl_coef(2)
					          write(oport_runmsg,108) rayl_coef(1)
					          write(oport_runmsg,109) rayl_coef(2)
					          return
					    end if
					end do
				end if
			end do
		elseif(pcomp(C_type, 'file',4)) then
		    do j= i+1, ls
		        if( cksep(string_stripped(j:j)) ) then
			        C_filename(1:ls-j)=string_stripped(j+1:ls)
					write(oport_runmsg,110) C_filename(1:ls-j)
				    return
				endif
			end do
		endif
107   format(4X, 'Type of damp matrix is :', a4)
108   format(4X, 'The 1st rayleigh coefficient is ', f15.4)
109   format(4X, 'The 2nd rayleigh coefficient is ', f15.4)
110	  format(4X, 'Damping matrix in file --', a)
	end subroutine ReadDampPara
! ---------------------------------------------------------------------
	subroutine ReadInitialPara(IniType, IniFile)
	! Purpose: read parameters for initial conditions
	! Outputs:
	!          IniType - type of initial conditions ( 'zero' or 'file')
	!          IniFile - name of the file where non-zero initial conditions are stored.
		
		character :: IniType*4, IniFile*12
		IniType = string_stripped(i+1:i+4)
		write(oport_runmsg,112) IniType
		do j= i+1, ls
		    if( cksep(string_stripped(j:j)) ) then
				IniFile(1:ls-j)=string_stripped(j+1:ls)
				write(oport_runmsg,111) IniFile(1:ls-j)

				return
			endif
		end do
112   format(4X, 'Method of specifying initial condition is :',a4)
111	  format(4X, 'Data in file --', a)
	end subroutine ReadInitialPara
! ---------------------------------------------------------------------

	subroutine ReadOutputPara(output_flag, n_output_dof, output_dof, output_file)
	! Purpose: read parameters for output macro command, such as 'disp','velo' and 'acce'
	! Outputs: 
	!           output_flag - .true. means results of this type will be outputted. 
	!		  	output_dof  - which degree of freedom do you want ouput?
	!			output_file - name of the file to which the results will be wrote.
		logical output_flag
		integer output_dof(max_output_dof), n_output_dof, k2
		character(12) output_file(max_output_dof)
		output_flag = .true.
		n_output_dof=0
		do j=i+1,ls
			if( cksep(string_stripped(j:j)) ) n_output_dof = n_output_dof + 1
		end do
		n_output_dof = int((n_output_dof +1)/2)  ! total number of dofs to be outputted.
		n_output_dof = min(n_output_dof, max_output_dof)
		write(oport_runmsg,121) n_output_dof
		k2=1
		j=i+1
		do while(k2 .le. n_output_dof)
			j=j+1
			if(cksep(string_stripped(j:j))) then  ! find the separator after the interger No.k2
				k=j+1
				do while (.not. cksep(string_stripped(k:k)))
					k=k+1
				end do
				read(string_stripped(i+1:j-1), '(i)') output_dof(k2)
				output_file(k2)=string_stripped(j+1:k-1)
				write(oport_runmsg, 122) k2,output_dof(k2),output_file(k2)
				k2=k2+1
				i=k
				j=i+1
			endif
		end do
121     format(4X,'Total number of outputted DOF is:', i2, &
				/,8X, 'No.',4X,'Number of DOF',4X,'Filename')
122     format(8X,i2,10X,i4,12X,a12)

	end subroutine ReadOutputPara
! ---------------------------------------------------------------------

	subroutine ReadMethodPara()
		logical pcomp
		integer k2
		TIM_type = string_stripped(i+1:i+4)
		write(*,201) TIM_type
		write(oport_runmsg,201) TIM_type
		
		do j=i+1,ls  ! search for the 2nd separator, after which is the TIM_para1

			if( cksep(string_stripped(j:j)) ) then  ! find the 2nd separator.

				do k=j+1,ls  ! search for the 3rd separator, after which is the TIM_para2

					if( cksep(string_stripped(k:k)) ) then ! find the 3rd separator.
						read(string_stripped(j+1:k-1), '(f)') TIM_para1
						do k2=k+1,ls  ! search for the 4th separator, after which is the TIM_para3
							if( cksep(string_stripped(k2:k2)) ) then ! find the 4th separator.
								read(string_stripped(k+1:k2-1), '(f)') TIM_para2
								read(string_stripped(k2+1:ls), '(f)') TIM_para3
								return
							endif
						end do
						! cann't find the 4th separator, so there are only 2 parameters.
		                read(string_stripped(k+1:ls), '(f)') TIM_para2
						return  
				    end if

				end do
				! cann't find the 3rd separator, so there is only one parameter in this line.
				read(string_stripped(j+1:ls),'(f)') TIM_para1
				return
			end if
		end do
		
201   format(4X,'Selected Time Integration Method (TIM) is :', a4)

	end subroutine ReadMethodPara
! ---------------------------------------------------------------------
	subroutine ReadPropPara()
		implicit none
		integer iprop
		read(string_stripped(i+1:ls), *) (prop_para(iprop),iprop=1,5)
		write(oport_runmsg,211) (prop_para(iprop),iprop=1,5)
211     format(4X, 'Parameters of proptational force :',/,6X,5f12.4)
	end subroutine ReadPropPara
! ---------------------------------------------------------------------
	subroutine ReadNodalForce()
		implicit none
		logical pcomp
		F_type=string_stripped(i+1:i+4)
		write(oport_runmsg,141) F_type
		if(pcomp(F_type, 'file', 4)) then
			do j= i+1, ls
				if( cksep(string_stripped(j:j)) ) then
					F_filename(1:ls-j)=string_stripped(j+1:ls)
					write(oport_runmsg, 142) F_filename(1:ls-j)
					return
				endif
			end do
		elseif(pcomp(F_type, 'sing', 4)) then
			do j= i+1, ls
				if( cksep(string_stripped(j:j)) ) then
					read(string_stripped(j+1:ls), '(i)') F_dof
					write(oport_runmsg,143) F_dof
					return
				endif
			end do
		elseif(pcomp(F_type, 'zero', 4) .or. pcomp(F_type, '    ', 4)) then
			write(oport_runmsg,144)
		else
			write(oport_runmsg,*) '***Error : the type of nodal force is wrong !!!'
			write(oport_runmsg,*) '*** The type of nodal force is : ', F_type
			stop '***Error : the type of nodal force is wrong !!!'
			
			
		endif
141   format(4X, 'Type of nodal force ID is :', a4)
142   format(4X, 'Nodal force ID in file --', a)
143   format(4X, 'Nodal force is imposed on DOF ', i)
144   format(4X, 'All nodal forces are zero')
	end subroutine ReadNodalForce
! ---------------------------------------------------------------------
end subroutine ReadPara