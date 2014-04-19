logical function cksep(x1)
!----------------------------------------------------------------------------------
!      Purpose: Check for existence of separator characters in data.
!      Inputs:
!         x1  -  Character to check
!      Outputs:
!         cksep - True of character is a valid separator; else false.
!----------------------------------------------------------------------------------
      implicit  none

      character x1*1

!     Input character separators are blank or comma

      cksep = (x1.eq.' ') .or. (x1.eq.',')

end function cksep