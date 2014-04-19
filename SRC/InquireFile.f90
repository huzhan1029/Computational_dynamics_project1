logical function InquireFile(filename)
! Purpose: inquire if the file exists in current directory.
	implicit none
	character filename*12
	logical exist_file
	inquire(FILE=filename, EXIST=exist_file)
	if(.not. exist_file) then ! file does not exist
		InquireFile = .false.
		write(*,200) filename
		stop 'Error occurs !!!'
	else ! file exists
		InquireFile = .true.      
    endif
200 format(' *ERROR* : Specified file (',a12,') does not exist ! ')
end function InquireFile

