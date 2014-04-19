module module_ioport

integer :: iport_macro=1000, iport_data=1001
integer :: oport_runmsg=2000
integer :: oport_disp(1:5)=(/2011,2012,2013,2014,2015/), &
		   oport_velo(1:5)=(/2021,2022,2023,2024,2025/), &
		   oport_acce(1:5)=(/2031,2032,2033,2034,2035/)

end module module_ioport