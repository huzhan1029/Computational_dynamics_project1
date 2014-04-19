integer function MacroNumber(MacroName) 

! -------------------------------------------------------------------------
!    Purpose: Calculate the sequence number of a macro in the macro list.
!    Input:
!      MacroName - current macro's name
!    Output:
!      MacroNumber - current macro's sequence number.
!--------------------------------------------------------------------------
	implicit none	
	integer :: i
	integer,parameter:: nmacro = 15
	logical pcomp
	character :: MacroName*4
	character(4),parameter:: macrolist(nmacro) = (/
       'end ', 'k   ', 'm   ', 'c   ', 'inid', &
	   'iniv', 'dt  ', 'time', 'disp', 'velo', &
	   'acce', 'ndof', 'meth', 'fpro', 'forc'
    /)
       		
	do i=1, nmacro
		if (pcomp(MacroName, macrolist(i),4)) then
			MacroNumber = i
			return
		end if
	end do

end function MacroNumber