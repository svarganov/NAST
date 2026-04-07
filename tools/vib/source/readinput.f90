module hess_data

implicit none
logical :: GAMESS, Molpro, QChem, ORCA
integer :: m  
character(len=80) :: hess_gamess, hess_molpro, hess_qchem, hess_orca, initial_geometry, final_geometry
!--------------------------------------------------------------------------------------------------------
character(len=80), allocatable :: array(:)
character(len=6) :: intersect_type
character(len=15), dimension(:),  allocatable :: symbol                                       ! 1-D text arrays
double precision, dimension(:), allocatable   :: charge, mass, freq, e, coord_shifted_mw_vec, RM, deployed, & ! deployed added for QChem   
                                                 all_eig_v_amu ! Added on 07/08/2025 to pass for Huang_Rhys subroutine in vib.f90 
!---------------------------------------------------------------------------------------------------------
integer :: number_of_atoms
double precision :: lambda, D_HESS, total_mass
!---------------------------------------------------------------------------------------------------------
double precision :: center_of_mass, rot, moments_eig, rot_to_inv
dimension           center_of_mass(1,3), rot(3,3), moments_eig(3), rot_to_inv(3,3)
!---------------------------------------------------------------------------------------------------------
double precision, dimension(:,:), allocatable :: coord, coord_final, coord_diff, coord_diff_mw, coord_shifted, & ! coord_mw added on 07/08/2025 for HRF
                                                 coord_shifted_mw, Hess, Hess_mw, projector, Hess_projected, eye, &  ! coord_final added on 07/14/2025 for HRF 
                                                 trans_projector, rot_projector, P,      &                  ! coord_diff added on 07/15/2025
                                                 Hess_projected_eig_vec_to_amu
!---------------------------------------------------------------------------------------------------------
double precision, dimension(:), allocatable :: HRF, reorganization_energy ! added on 01/27/2026
contains
!-------------------------------------------------------------------
subroutine read_input_file()

implicit none
namelist /programs/ GAMESS, Molpro, QChem, ORCA
namelist /gamess_input_file/ hess_gamess
namelist /molpro_input_file/ hess_molpro
namelist /qchem_input_file/  hess_qchem
namelist /orca_input_file/   hess_orca
namelist /geometries/ initial_geometry, final_geometry
!-------------------------------------------------------------------
read(11, nml = programs) ! 11 is the input_file

if (GAMESS .and. (.not. Molpro) .and. (.not. QChem) .and. (.not. ORCA)) then
 m = 1
 write(*,*) "Working with GAMESS input files"
 read(11, nml = gamess_input_file)
else if (Molpro .and. (.not. GAMESS) .and. (.not. QChem) .and. (.not. ORCA)) then
 m = 2
 write(*,*) "Working with Molpro input files"
 read(11, nml = molpro_input_file)
else if (QChem .and. (.not. GAMESS) .and. (.not. Molpro) .and. (.not. ORCA)) then
 m = 3
 write(*,*) "Working with QChem input files"
 read(11, nml = qchem_input_file)
else if (ORCA .and. (.not. GAMESS) .and. (.not. Molpro) .and. (.not. QChem)) then
 m = 4
 write(*,*) "Working with ORCA input files"
 read(11, nml = orca_input_file)
else
 write(*,*) "Error! Make a proper selection of a program (GAMESS, Molpro, QChem or ORCA). &
             Only one program can be used at a time"
 call exit
end if 

read(11, nml = geometries)

end subroutine read_input_file
!----------------------------------------------------------------------
end module hess_data
