real function PropForce(t)
!-------------------------------------------------------------
!  Purpose: calculate proportional force at t
!  Input:   t - time
!  Output:  PropForce - value of force at t.
!-------------------------------------------------------------
	use module_parameter, only: a=>prop_para
	implicit none
	real t
	PropForce = a(1) + a(2)*t + a(3)*sin(a(4)*t+a(5))

end function PropForce
	