module hess_driver

use M_strings
use get_constants
use read_geometry ! Added 07/14/2025
use hess_data
use GAMESS_read_module
use Molpro_read_module
use QChem_read_module
use ORCA_read_module
use vibrational
use normal_mode
use write_data

implicit none

!---------------------------------------------------------------------------------------------------------
contains
subroutine driver(input_file,fgr_file,imag_file)
implicit none

character(len=*), intent(in) :: input_file,fgr_file,imag_file

!================ Read input file =================!
!==================================================!
! The subroutine below reads the effhess input     !
! file and determines from which electronic        !
! structure package comes the Hessian matrix.      ! 
!==================================================!

call read_input_file()

select case(m)

case(1) ! GAMESS

 open(unit = 10, status = "old", action = "read", file = hess_gamess)

case(2) ! Molpro
 
 open(unit = 12, status = "old", action = "read", file = hess_molpro)

case(3) ! QChem
 open(unit = 14, status = "old", action = "read", file = hess_qchem)

case(4) ! ORCA
 open(unit = 12, status = "old", action = "read", file = hess_orca)

end select 

!=============== Read atomic data =================!
!==================================================!
! The subroutine below reads the number of atoms,  !
! atomic symbols, charges, and atomic coordinates  !
! from the file "initial_geometry.xyz" and atomic  !            
! coordinates for the final geometry from the file !
! "final_geometry.xyz".                            !
!==================================================!

open(unit = 15, status = "old", action = "read", file = initial_geometry)
open(unit = 16, status = "old", action = "read", file = final_geometry)
call read_atomic_data(15, 16)

!================ Read Hessian ====================!
!==================================================!
! Reading the input Hessian depends on the source  !
! of data (GAMESS, Molpro, QChem, or ORCA) since   !
! the programs have different formats for storing  !
! Hessian matrices, and in general different       !
! formats. For that, there are four subroutines    !
! below for reading input blocks, executed via     ! 
! the case statement.                              !
!==================================================!

select case(m)

case(1)

!==================== GAMESS ======================!
!                                                  !
!                  rdGAMESS.f90                    !
!                                                  !
!==================================================!

!============= Read gradient and Hessian ==========!
!==================================================!
! The below subroutine reads Hessian for the 
! inital state.                                    !
!==================================================!

call GAMESS_read_hess(10,Hess)

case(2)

!==================== Molpro ======================!
!                                                  !
!                  rdMolpro.f90                    !
!                                                  !
!==================================================!

call Molpro_read_hess(12,Hess)

case(3)

!==================== QChem =======================!
!                                                  !
!                  rdQChem.f90                     !
!                                                  !
!==================================================!

call QChem_read_hess(14,Hess)

case(4)

!===================== ORCA =======================!
!                                                  !
!                    rdORCA.f90                    !
!                                                  !
!==================================================!

call ORCA_read_hess(12,Hess)

end select

!============ End reading input data ==============!

!========= Build and project the Hessian ==========!
!                                                  !
!                    vib.f90                       !
!                                                  !
!==================================================!

!~~~~~ Prepare for projection of the Hessian

! 2. Calculate center of mass of the molecule

call calc_center_of_mass()

! 3. Mass-weight the Hessian
 
call calc_hess_mw()

! 5. Construct moment of intertia tensor 

call calc_rot()

! 6. Calculate the inverse of the moment of inertia tensor

call calc_inverse()

! 7. Construct projectors to project out translation,
!    rotation, and the reaction coordinate (rc) degrees of freedom 

call proj_construct()

! 8. Project the mass-weighted Hessian

call project_hess()

! 9. Calculate vibrational frequencies

call get_freq(imag_file)

! 10. Calculate Huang-Rhys factors ! Added 07/08/2025

call Huang_Rhys()

!========= End of Effective Hessian block =========!


!============ Write the output file ===============!
!                                                  !
!                    write.f90                     !
!                                                  !
!==================================================!

call write_output(input_file)
call write_mlj_template(fgr_file) 


end subroutine driver
end module hess_driver
