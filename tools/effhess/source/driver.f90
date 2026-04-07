module hess_driver

use M_strings           ! strings.f90
use get_constants       ! constants.f90
use hess_data           ! readinput.f90
use GAMESS_read_module  ! rdGAMESS.f90
use Molpro_read_module  ! rdMolpro.f90
use QChem_read_module   ! rdQchem.f90
use ORCA_read_module    ! rdORCA.f90
use vibrational         ! vib.f90
use normal_mode         ! normal_mode.f90
use write_data          ! write.f90

implicit none

!---------------------------------------------------------------------------------------------------------
contains
subroutine driver(input_file,rc_file,nast_file,imag_file)
implicit none

character(len=*), intent(in) :: input_file,rc_file,nast_file,imag_file ! from module hess_data (readinput.f90)

!---------------------------------------------------------------------------------------------------------
call read_input_file()  ! Reads the content of the input file. This subroutine is in module hess_data (readinput.f90) 

select case(m) ! The value of 'm' is set in the 'call read_input_file()' subroutine

case(1) ! GAMESS

 open(unit = 10, status = "old", action = "read", file = low_spin_dat)
 open(unit = 11, status = "old", action = "read", file = high_spin_dat)
 open(unit = 12, status = "old", action = "read", file = input_geometry)

case(2) ! Molpro
 
 open(unit = 10, status = "old", action = "read", file = low_spin_grad)
 open(unit = 11, status = "old", action = "read", file = high_spin_grad)
 open(unit = 12, status = "old", action = "read", file = low_spin_hess)
 open(unit = 14, status = "old", action = "read", file = high_spin_hess)

case(3) ! QChem
 open(unit = 10, status = "old", action = "read", file = geometry_in)
 open(unit = 11, status = "old", action = "read", file = low_spin_grad_fchk)
 open(unit = 12, status = "old", action = "read", file = high_spin_grad_fchk)
 open(unit = 14, status = "old", action = "read", file = low_spin_hess_fchk)
 open(unit = 15, status = "old", action = "read", file = high_spin_hess_fchk)

case(4) ! ORCA
 open(unit = 10, status = "old", action = "read", file = low_spin_grad_orca)
 open(unit = 11, status = "old", action = "read", file = high_spin_grad_orca)
 open(unit = 12, status = "old", action = "read", file = low_spin_hess_orca)
 open(unit = 13, status = "old", action = "read", file = high_spin_hess_orca)

end select 

!============= Read input data ====================!

!~~~~~~~ Reading the input depends on the source of
!~~~~~~~ data (GAMESS, Molpro, QChem, or ORCA) since
!~~~~~~~ the programs have different formats for
!~~~~~~~ storing geometry, gradients and Hessian
!~~~~~~~ matrices. For that, there are four
!~~~~~~~ reading-input blocks, executed
!~~~~~~~ via the case statement.

select case(m)

case(1)

!==================== GAMESS ======================!
!                                                  !
!                  rdGAMESS.f90                    !
!                                                  !
!==================================================!

!=========== Determine number of atoms ============!
!================================================= !
! The below subroutine reads the number of atoms   !
! in the molecule from one of the input files.     !
!==================================================!

call GAMESS_read_number_of_atoms(12)

!=============== Read atomic data =================!
!==================================================!
! The below subroutine reads atomic symbols,       !
! charges and atomic coordinates from one of the   !
! input files.                                     !
!==================================================!

call GAMESS_read_atomic_data(12)

!============= Read gradient and Hessian ==========!
!==================================================!
! The below subroutine reads gradient and Hessian  !
! for each of the spin states.                     !
!==================================================!

call GAMESS_read_grad_hess(10,grad1,Hess1)
call GAMESS_read_grad_hess(11,grad2,Hess2)

case(2)

!==================== Molpro ======================!
!                                                  !
!                  rdMolpro.f90                    !
!                                                  !
!==================================================!

call Molpro_read_number_of_atoms(12)

call Molpro_read_atomic_data(12)

call Molpro_read_grad_hess(10,12,grad1,Hess1)
call Molpro_read_grad_hess(11,14,grad2,Hess2)

case(3)

!==================== QChem =======================!
!                                                  !
!                  rdQChem.f90                     !
!                                                  !
!==================================================!

call QChem_read_number_of_atoms(11)

call Qchem_read_atomic_data(10,11)

call QChem_read_grad_hess(11,14,grad1,Hess1)
call QChem_read_grad_hess(12,15,grad2,Hess2)

case(4)

!===================== ORCA =======================!
!                                                  !
!                    rdORCA.f90                    !
!                                                  !
!==================================================!

call ORCA_read_atomic_data(12)

call ORCA_read_grad_hess(10,12,grad1,Hess1)
call ORCA_read_grad_hess(11,13,grad2,Hess2)

end select

!============ End reading input data ==============!

!==== Build and project the Effective Hessian =====!
!                                                  !
!                    vib.f90                       !
!                                                  !
!==================================================!

!~~~~~ Prepare for projection of the Effective Hessian

! 1. Build the Effective Hessian

call construct_eff_hess()

! 2. Calculate center of mass of the molecule

call calc_center_of_mass()

! 3. Mass-weight the Effective Hessian
 
call calc_eff_hess_mw()

! 4. Calculate normalized mass-weighted |g1 - g2|

call calc_norm_deltaG_mw()

! 5. Construct moment of intertia tensor 

call calc_rot()

! 6. Calculate the inverse of the moment of inertia tensor

call calc_inverse()

! 7. Construct projectors to project out translation,
!    rotation, and the reaction coordinate (rc) degrees of freedom 

call proj_construct()

! 8. Project the mass-weighted Effective Hessian

call project_eff_hess()

! 9. Calculate vibrational frequencies

call get_freq(imag_file)

! 10. Calculate reduced masses along normal modes of vibration 

call get_red_mass()

! 11. Calculate the reduced mass along rc

call get_red_mass_along_rc()

!========= End of Effective Hessian block =========!


!= Generate Cartesian displacements along the rc ==! 
!                                                  !
!                 normal_mode.f90                  !
!                                                  !
!==================================================!

open(unit = 33, action = "write", file = rc_file)

call get_cart_disp(33,number_of_atoms,symbol,  &
                      charge,coord,freq_rc,rc_eig_v_amu)


!============ Write the output file ===============!
!                                                  !
!                    write.f90                     !
!                                                  !
!==================================================!

call write_output(input_file)
call write_nast_template(nast_file) 


!~~~~ Formatting

33 format(A1,'    ',F4.1,'   ',F12.9,'   ',F12.9,'    ',F12.9)

end subroutine driver
end module hess_driver
