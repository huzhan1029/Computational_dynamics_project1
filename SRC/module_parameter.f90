module module_parameter

implicit none
character(12) :: K_filename, M_filename, C_filename, IniD_filename,   &
                 IniV_filename, F_filename
character(4) :: M_type, C_type, IniD_type, IniV_type, TIM_type, F_type

integer :: max_output_dof=5, n_Disp_dof, n_Velo_dof, n_Acce_dof,   &
           ndof, F_dof, nstep
integer, dimension (5) :: Disp_dof, Velo_dof, Acce_dof
character(12), dimension (5) :: Disp_filename, Velo_filename,   &
               Acce_filename
logical :: Disp_flag=.false., Velo_flag=.false., Acce_flag=.false.,   &
           C_exist=.false.
real :: rayl_coef(1:2), dt, TotalTime, prop_para(1:5)=0.0
real :: TIM_para1=0.0, TIM_para2=0.0, TIM_para3=0.0

! In next content, DOF means Degree of Freedom.
! K_filename - name the file where stiff matrix K is stored.
! M_filename - name the file where mass matrix M is stored.
! C_filename - name the file where damping matrix C is stored.
! IniD_filename - name the file where initial displacement is stored.
! IniV_filename - name the file where initial velocity is stored.
! F_filename - name the file where nodal force ID is stored.
! M_type - type of the mass matrix. 'lump' for lumped, or 'cons' for consistent.
! C_type - type of the damping matrix. 'file' means C will be imported from file.
!          'rayl' means Rayleigh Damping will be adopted.
! IniD_type, IniV_type - 'file' means initial conditions will be imported from files.
!                        'zero' means initial conditions are all zero.
! TIM_type - Time Integration Method to be used.
! F_type - type of nodal force ID. 'file' means reading ID from file.
!          'sing' (single) means only one node has force on it, number of DOF of this
!          node must be declared. 'zero' means no external force.
! max_output_dof - perhaps you want to output the results of some DOFs, the max number is
!                  defined by this varibale.
! n_Disp_dof - the number of the DOFs whose displacement will be outputted.
! n_Velo_dof, n_Acce_dof - similar to n_Disp_dof
! n_dof - total number of the DOFs in space.
! F_dof - if F_type = 'sing', F_dof is the number of the DOF on which nodal force is imposed
! nstep - total number of time intervals.
! Disp_dof - an array. If you want to output the displacement of DOF 1, 3 and 4, then
!            Disp_dof is /1,3,4,0,0/.
! Velo_dof, Acce_dof - similar to Disp_dof
! Disp_filename - a character-type array. to store the filenames where you want to output
!                 the displacement.
! Velo_filename,Acce_filename - similar to Disp_filename
! Disp_flag - if .true., displacement will be outputted to file.
! Velo_flag - if .true., velocity will be outputted to file.
! Acce_flag - if .true., acceleration  will be outputted to file.
! C_exist - if .true., damping exists in this problem; .false. means no damping.
! rayl_coef - coefficients of Rayleigh damping. C = rayl_coef(1)*M + rayl_coef(2)*K
! dt - length of time step
! TotalTime - total time to be analyzed.
! prop_para - parameters for proportional force. let a(i)=prop_para(i),i=1,5
!             propforce(t) = a(1) + a(2)*t + a(3)*sin( a(4)*t + a(5) )
! TIM_para1,TIM_para2,TIM_para3 - parameters for Time Integration Method
end module module_parameter
