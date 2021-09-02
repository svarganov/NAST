module hess_data

implicit none
logical :: GAMESS, Molpro
integer :: m   
character(len=80) :: low_spin_dat, high_spin_dat, input_geometry,                       & ! GAMESS files
                     low_spin_grad, high_spin_grad, low_spin_hess, high_spin_hess         ! Molpro files
!--------------------------------------------------------------------------------------------------------
character(len=80), allocatable :: array(:)
character(len=6) :: intersect_type
character(len=15), dimension(:),  allocatable :: symbol                                       ! 1-D text arrays
double precision, dimension(:), allocatable   :: grad1, grad2, charge, deltaG, mass,            &
                                                 freq, e, e_rc, reduced_mass,                   &
                                                 coord_shifted_mw_vec, RM, m_deltaG,            &
                                                 m_norm_deltaG, Eff_Hess_rc_eig_v_amu, rc_eig_v_amu ! number vectors
!---------------------------------------------------------------------------------------------------------
integer :: number_of_atoms
double precision :: norm_deltaG, lambda, gradmean, D_HESS, total_mass, e_rc_max,         &
                    norm_grad1, norm_grad2, gradslope, reduced_mass_rc, freq_rc
!---------------------------------------------------------------------------------------------------------
double precision :: center_of_mass, rot, moments_eig, rot_to_inv
dimension           center_of_mass(1,3), rot(3,3), moments_eig(3), rot_to_inv(3,3)
!---------------------------------------------------------------------------------------------------------
double precision, dimension(:,:), allocatable :: coord, coord_shifted, coord_shifted_mw, &
                                                 Hess1, Hess2, Eff_Hess, Eff_Hess_mw, projector, &
                                                 Eff_Hess_projected, eye, rc_projector,           &
                                                 trans_projector, rot_projector, P,      &
                                                 Eff_Hess_projected_eig_vec_to_amu, Eff_Hess_rc
contains
!-------------------------------------------------------------------
subroutine read_input_file()

implicit none
namelist /programs/ GAMESS, Molpro
namelist /gamess_input_files/ low_spin_dat, high_spin_dat, input_geometry
namelist /molpro_input_files/ low_spin_grad, high_spin_grad,      &
                              low_spin_hess, high_spin_hess
!-------------------------------------------------------------------
read(11, nml = programs)

if (GAMESS .and. (.not. Molpro)) then
 m = 1
 write(*,*) "Working with GAMESS input files"
 read(11, nml = gamess_input_files)
else if (Molpro .and. (.not. GAMESS)) then
 m = 2
 write(*,*) "Working with Molpro input files"
 read(11, nml = molpro_input_files)
else if ((GAMESS .and. Molpro)) then
  write(*,*) "Error! You chose both GAMESS and Molpro to be TRUE. &
             Check your input and return"
  call exit
else if ((.not. GAMESS) .and. (.not. Molpro)) then
  write(*,*) "Error! You chose both GAMESS and Molpro to be FALSE.  &
             Check your input and return"
  call exit
end if 

end subroutine read_input_file
!----------------------------------------------------------------------
end module hess_data
