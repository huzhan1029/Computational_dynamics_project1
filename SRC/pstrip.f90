subroutine pstrip(xxx,yyy,i,j)
!----------------------------------------------------------------------------------
!      Purpose: Strip off comments (begin with !) & leading blanks.
!      Inputs:
!         yyy(*)    - Input string
!         i         - First character to look for ! on comments
!                     N.B. Command language uses !x and !! for
!                          re-execution of command.
!      Outputs:
!         xxx(*)    - Output string after strips
!         j         - length of string after strips
!----------------------------------------------------------------------------------
      implicit  none

      logical   cksep
      character xxx*255,yyy*255
      integer   i,n,nn,j

!     Strip comments

      do n = i,255
        if(ichar(yyy(n:n)).eq.13) then
          yyy(n:n) = ' '
        elseif(yyy(n:n).eq.'!') then
          yyy(n:255) = ' '
          go to 100
        end if
      end do

!     Strip leading blanks

100   xxx = ' '
      do n = 1,255
        if(yyy(n:n).ne.' ') then
          xxx(1:256-n) = yyy(n:255)
		  go to 200
        end if
      end do
	  j=0
      return

!     Find last character

200   do nn = 255,1,-1
        if(xxx(nn:nn).ne.' ') go to 300
      end do
      nn = 2  

!     Remove extra blanks

300   n = 1
! remove blank before comma and blank
301   if(xxx(n:n).eq.' ' .and. cksep(xxx(n+1:n+1))) then  
        xxx(n:nn-1) = xxx(n+1:nn)
        xxx(nn:nn) = ' '
        nn = nn - 1
        go to 301
      endif
      n = n + 1
      if(n.lt.nn) go to 301
      j=nn

! remove blank after comma
      do n = 1,nn-2
        if(xxx(n:n).eq.',' .and. xxx(n+1:n+1).eq.' ' ) then
          xxx(n+1:nn-1) = xxx(n+2:nn)
          xxx(nn:nn) = ' '
		  j=j-1
        end if
      end do

end subroutine pstrip