module module_data

implicit none

real,allocatable:: K_matrix(:,:), M_matrix(:,:), C_matrix(:,:), &
				   d0_vector(:), v0_vector(:), a0_vector(:),    &
                   NodalForceId(:)

end module module_data