logical function pcomp(a,b,n)
!---------------------------------------------------------------------------------
!      Purpose: Compare character strings for match
!               Ignores upper/lower case differences.
!      Inputs:
!         a(*)   - Character string 1
!         b(*)   - Character string 2
!         n      - Number of characters to compare
!      Outputs:
!         pcomp  - Flag, true if a = b
!---------------------------------------------------------------------------------
      implicit  none

      integer   n, inc, i, ia,ib
      character a*(*),b*(*)

!     Compute increment between an upper and lower case letter

      inc = ichar('A') - ichar('a')

!     Compare for match

      pcomp = .false.
      do i = 1,n

        ia = ichar(a(i:i))
        ib = ichar(b(i:i))

!       Test all permutations of characters for match

        if(ia.ne.ib .and. ia+inc.ne.ib .and. ia.ne.ib+inc ) return
      end do ! i

      pcomp = .true.

end function pcomp
