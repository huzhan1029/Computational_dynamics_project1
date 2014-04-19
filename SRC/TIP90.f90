program TIP90

!-----[--.----+----.----+----.-----------------------------------------]

!      * * TIP90 * * A Time Integration Program in FORTRAN 90
!                      -    -           -                  --
!....  Copyright (c) 2005-2005: 
!                Computational Dynamics Group, Tsinghua University
!                All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]

!     The following time integration methods has been included in TIP90
!         Central difference method     - TIM_Central.f90                 
!         Houbolt method                - TIM_Houbolt.f90                   
!         Newmark method                - TIM_Newmark.f90
!         Wilson \theta method          - TIM_Wilson.f90
!         Galerkin Weak form (3 points) - TIM_GW3.f90
!         Galerkin Weak form (4 points) - TIM_GW4.f90
!         Precise integration           - TIM_Precise.f90

!     Programmed by:

!         ZHANG Xiong
!         School of Aerospace
!         Tsinghua University
!         Beijing, 100084
!     E-mail:
!         xzhang@tsinghua.edu.cn
!     Web address:
!         www.dynamics.tsinghua.edu.cn/xzhang/tip90

!-----[--.----+----.----+----.-----------------------------------------]

	use module_ioport
	implicit none

    character :: infilename*12, c
    logical :: exist_infile

	open(unit=oport_runmsg, file='msg.txt', status='unknown')
    exist_infile = .false.

! Read the filename from screen and test existence.
    do while (.not. exist_infile)  
        write(*,100)
		read(*,'(a)') infilename
		inquire(FILE=infilename, EXIST=exist_infile)
		if(.not. exist_infile) then ! Specified file does not exist.
			write(*,105)
        endif
    end do
	write(oport_runmsg,101) infilename

	call ReadPara(infilename)

	call ReadData

	call TimeIntegration

	close(oport_runmsg)

100 format('Specify the Input Filename:')
101 format('The macro command input file is "', a,'"')
105 format(/' *ERROR* FILENAME: Specified input file does not exist, ', &
            'reinput name.'/)

end program TIP90
