module hess_driver

use M_strings
use get_constants
use hess_data
use GAMESS_read_module
use Molpro_read_module
use vibrational
use normal_mode
use write_data

implicit none

!---------------------------------------------------------------------------------------------------------
contains
subroutine driver(input_file,rc_file,nast_file,imag_file)
implicit none

character(len=*), intent(in) :: input_file,rc_file,nast_file,imag_file

!---------------------------------------------------------------------------------------------------------
call read_input_file()

select case(m)

case(1)

 open(unit = 10, status = "old", action = "read", file = low_spin_dat)
 open(unit = 11, status = "old", action = "read", file = high_spin_dat)
 open(unit = 12, status = "old", action = "read", file = input_geometry)

case(2)
 
 open(unit = 10, status = "old", action = "read", file = low_spin_grad)
 open(unit = 11, status = "old", action = "read", file = high_spin_grad)
 open(unit = 12, status = "old", action = "read", file = low_spin_hess)
 open(unit = 14, status = "old", action = "read", file = high_spin_hess)

end select 

!~~~~~~~~~~~~~~ Read input data ~~~~~~~~~~~~~~~~~~~~!

!~~~~~~~ Reading the input depends on the source of
!~~~~~~~ data (GAMESS, Molpro) since  two programs
!~~~~~~~ have different formats for storing geometry,
!~~~~~~~ gradients and Hessian matrices. For that,
!~~~~~~~ there ar two reading-input blocks, executed
!~~~~~~~ via the case statement.

select case(m)

case(1)

!~~~~~~~~~~~~~~~~~~~~ GAMESS ~~~~~~~~~~~~~~~~~~~~~~!

!~~~~~~~~~~~ Determine number of atoms ~~~~~~~~~~~~!
!                                                  !
! The 'GAMESS_read_number_of_atoms(file_id,numbr)' !
! subroutine reads the input file of 'file_id'     !
! to determine the number of atoms 'numbr'         !
! in the system of interest                        !
!==================================================!

call GAMESS_read_number_of_atoms(12)

!=============== READ ATOMIC DATA =================!
!==================================================!
! The 'GAMESS_read_atomic(file_id,numbr,symb,      !
!      chrg,coord)' subroutine reads 'file_id1'    !
! to determine atomic symbols, charges and         !
!              coordinates of atoms                !
!==================================================!

call GAMESS_read_atomic_data(12)

!============= READ GRADIENT AND HESSIAN ==========!
!==================================================!
! The 'GAMESS_read_grad_hess(file_id,numbr,        !
!     grad,HESSIAN)' subroutine reads gradient     !
! and Hessian from input file 'file_id'.           !
! For low-spin state, 'file_id' = 10; for the      !
! high-spin = 11. Numbers are assigned above       !
!==================================================!

call GAMESS_read_grad_hess(10,grad1,Hess1)
call GAMESS_read_grad_hess(11,grad2,Hess2)

case(2)

!==================== MOLPRO ======================!

!=========== DETERMINE NUMBER OF ATOMS ============!
!==================================================!
! The 'Molpro_read_number_of_atoms(file_id,numbr)' !
! subroutine reads the input file of 'file_id'     !
! to determine the number of atoms 'numbr'         !
! in the system of interest                        !
!==================================================!

call Molpro_read_number_of_atoms(12)

!=============== READ ATOMIC DATA =================!
!==================================================!
! The 'Molpro_read_atomic(file_id,numbr,symb,      !
!      chrg,coord)' subroutine reads 'file_id1'    !
! to determine atomic symbols, charges and         !
!              coordinates of atoms                !
!==================================================!
call Molpro_read_atomic_data(12)
!============= READ GRADIENT AND HESSIAN ==========!
!==================================================!
! The 'Molpro_read_grad_hess(file_id1,file_id2,    !
! numbr,grad,HESSIAN)' subroutine reads gradient   !
! from input file 'file_id1' and Hessian from      !
! input file 'file_id2'. For low-spin state,       !
! 'file_id1' = 10 and 'file_id2' = 12; for the     !
! high-spin = 11 and 14 correspondingly.           !
!           Numbers are assigned above             !
!==================================================!

call Molpro_read_grad_hess(10,12,grad1,Hess1)
call Molpro_read_grad_hess(11,14,grad2,Hess2)

end select

!=============== END READ INPUT DATA ====================!

!============= CONSTRUCT EFFECTIVE HESSIAN ================!
!==========================================================!
! See vib.f90 for 'construct_eff_hess' subroutine          !
! List of parameters:                                      !
! - numbr: number_of_atoms - required to set size          !
! - grad1,grad2: read above gradients of two spin states   !
! - HESSIAN1, HESSIAN2: HESS1 and HESS2 - state specific   !
!                       Hessian matries                    !
! - dlG: deltaG = grad2 - grad1 - gradient parallel to RC  !
! - n_(dlG,grad1,grad2): norm_(deltaG,grad1,grad2) -       !
!                        Eucledian norms of three gradients!
! - lam: lambda - Lagrangian multiplier                    !
! - grad_mn: gradmean - geometric mean of two gradients    !
! - HESSIAN: HESS - Effective Hessian                      !
! - slope: gradslope - type of intersection (quantity)     !
! - int_type: intersect_type - type of intersection        !
!                                      (quality)           !
!==========================================================! 

call construct_eff_hess ()

!============= END CONSTRUCT EFFECTIVE HESSIAN ==========!

!====== PREPARE FOR PROJECTION OF EFFECTIVE HESSIAN =====!

!===== CALCULATE CENTER OF MASS OF THE MOLECULE ===!
!========================================================!
! The 'calc_com(numbr,coord,chrg,coord_shft,             !
! coord_shft_mw,mass,total_mass,com)' subroutine         !
! calculates the com of the molecule, shifts 'coord' to  !
! the com and then calculates mass-weighted coordinates  !
!========================================================! 

call calc_center_of_mass()
 
!======= CALCULATE MASS-WEIGHTED EFFECTIVE HESSIAN=======!
!========================================================!
! The 'calc_eff_hess_mw(numbr,HESSIAN,mass,HESSIAN_mw)'  !
! subroutine returns the mass-weighted Hessian           !
! 'HESSIAN_mw' for the subsequent vibrational analysis   !
!========================================================!

call calc_eff_hess_mw()

!======= CALCULATE NORMALIZED MASS-WEIGHTED DELTA_G =====! 
!========================================================!
! The 'calc_norm_dlG_mw(numbr,mass,dlG,m_dlG,m_norm_dlG)'!
! subroutine computes the normalized mass-weighted       !
!               parallel gradient                        !
!========================================================!

call calc_norm_deltaG_mw()


!========= CALCULATE MOMENT OF INERTIA TENSOR ===========!
!============= IN MASS-WEIGHTED COORDINATES =============!
!========================================================!
! The calc_rot(numbr,coord_shft_mw,rot,mom_eig,          !
! rot_to_inv) subroutine calculates the moment of        !
! inertia tensor, 'rot' and makes a copy of 'rot' to     !
! 'rot_to_inv' for reasons given in vib.f90. Then,       !
! the MKL subroutine dsyev is called to diagonalize      !
! the 'rot' matrix and get principal moments of inertia, !
!         'moments_eig' as 'rot' eigenvalues             !
!========================================================!

call calc_rot()


!===== CALCULATE INVERSE OF MOMENT OF INERTIA TENSOR ====!
!========================================================!
! The 'dgetrf' subroutine of LAPACK call the LU          !
! decomposition for 'rot_to_inv to calculate the desired !
! inverse of 'rot_to_inv'. The 'dgetri' subroutine       !
! returns the inverse of original 'rot_to_inv'. Thus, the!
! final result is 'rot_to_inv' = inv('rot_to_inv') =     !
!                     inv('rot')                         !
!========================================================!

call calc_inverse()

!==== CONSTRUCT RC, TRANS, ROT AND TOTAL PROJECTORS =====!
!========================================================!
! The 'proj_construct' subroutine constructs the total   !
! projector P to project out six (rot+trans) + one (RC)  !
! degrees of freedom (DOF). Rot,trans and RC components  !
! of P are written in (rot_,trans_,rc_)projector         !
! for debugging. Besides, rc_projector is further used   !
! to determine reduced mass associated with the RC       !
!========================================================!
call proj_construct()

!=== PROJECTION OF THE MASS-WEIGHTED EFFECTIVE HESSIAN ==!
!========================================================!
! The 'project_eff_hess(numbr,HESSIAN_mw,                !
! P,proj,HESSIAN_proj,e)' subroutine returns projected   !
! Effective Hessian 'HESSIAN_proj'. In 'HESSIAN_proj',   !
! seven DOF are washed away. The subroutine then calls   !
! LAPACK diagonalization for 'HESSIAN_proj' and stores   !
! the matrix eigenvalues in 'e'. In addition, the        !
! subroutine returns 'proj = (I - P)', where I is the    !
!                    Identity Matrix                     !
!========================================================!

call project_eff_hess()

!========== CALCULATE VIBRATONAL FREQUENCIES ============!
!========================================================!
! The 'get_freq(numbr,e,freq_array)' subroutine          !
! returns vibrational frequencies associated with the    !
! eigenvalues 'e' of the mass-weighted Effective Hessian !
!========================================================!

call get_freq(imag_file)

!================ CALCULATE REDUCED MASSES ==============!
!========================================================!
! The 'get_red_mass(numbr,mass,HESSIAN_proj,             !
! eig_vec_to_amu,red_m)' subroutine calculates           !
! the reduced mass 'red_m' as an inverse of the dot      !
! product of a mass-weighted eigenvector of HESSIAN_proj !
! with itself: red_m = 1/(vec*vec), where                !
!                vec = eig_vec_to_amu                    !
!========================================================!

call get_red_mass()

!============= DETERMINE THE REDUCED MASS ===============!
!======= ASSOCIATED WITH THE REACTION COORDINATE ========!
!========================================================!
! The 'get_red_mass_along_rc(numbr,mass,rc_proj,         !
! HESSIAN_mw,rc_red_m,HESSIAN_rc)' subroutine calculates !
! the reduced mass along the RC. For that, the mass-     !
! weighted Hessian 'HESS_mw' gets projected onto the RC  !
! projector 'P_RC = rc_projector':                       !
! 'HESS_RC = P_RC*HESS_mw*P_RC'. Then, only the RC DOF is!
! left in HESS_RC with the non-zero eigenvalue.          !
! The corresponding eigenvector is used to calculate the !
! reduced mass                                           !

call get_red_mass_along_rc()

!===== GENERATE CARTESIAN DISPLACEMENTS ALONG THE RC ====! 

open(unit = 33, action = "write", file = rc_file)

call get_cart_disp(33,number_of_atoms,symbol,  &
                      charge,coord,freq_rc,rc_eig_v_amu)

!================ WRITE THE OUTPUT FILE =================!

call write_output(input_file)

!=============== WRITE NAST TEMPLATE FILE ===============!

call write_nast_template(nast_file)


!**************** FORMATTING OPTIONS TO WRITE
!**************** DIFFERENT GAMESS CARDS DATA ************!


33 FORMAT(A1,'    ',F4.1,'   ',F12.9,'   ',F12.9,'    ',F12.9)
21 FORMAT(1X,'$force',1X,'rdhess=.t.',1X,'projct=.t.',1X,'$end',/,                 &
          1X,'$contrl',1X,'runtyp=hessian',1X,'maxit=200',/,                       &
          1X,'icharg=',I2,1X,'mult=',I1,1X'scftyp=uhf',1X,'$end',/,                &
          1X,'$system',1X,'mwords=500',1X,'$end',/,                                &
          1X,'$basis',1X,'gbasis=n31',1X,'ngauss=6',1X,'$END')
!************************* END ***************************!

!**************** DEALLOCATE ALL VARIABLES ***************!
!deallocate(coord,coord_shifted, coord_shifted_mw)
!deallocate(HESS1, HESS2, HESS, HESS_mw, projector)
!deallocate(HESS_proj, m_deltaG, m_norm_deltaG)
!deallocate(rc_projector, eye)
!deallocate(symbol, grad1, grad2, charge, deltaG)
!************************** END **************************!
end subroutine driver
end module hess_driver
