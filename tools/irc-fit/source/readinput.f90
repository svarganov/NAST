module irc_data

implicit none

logical :: GAMESS, Molpro, ORCA
integer :: m   
character(len=80) :: g_irc_output_reac, g_irc_output_prod, & ! GAMESS files
                     m_irc_output, &                            ! Molpro files
                     g_irc_output_reac_orca, g_irc_output_prod_orca ! ORCA files

integer                       :: n1, n2, numbr
double precision, allocatable :: geom_reac(:,:,:), geom_prod(:,:,:), en_reac(:), en_prod(:), &
                                 irc_reac(:), irc_prod(:), ref(:,:)
double precision              :: l(5)=0.0d0,h(5)=0.0d0

contains
!-------------------------------------------------------------------
subroutine read_input_file(file_inp)

implicit none
character(len=*),intent(in)              :: file_inp
namelist /programs/ GAMESS, Molpro, ORCA
namelist /gamess_input_files/ g_irc_output_reac, g_irc_output_prod
namelist /molpro_input_files/ m_irc_output
namelist /orca_input_files/ g_irc_output_reac_orca, g_irc_output_prod_orca
!-------------------------------------------------------------------
open (unit = 99, file = file_inp)
read(unit = 99, nml = programs)

if (GAMESS .and. (.not. Molpro)) then
 m = 1
 write(*,*) "Working with GAMESS input files"
 read(unit = 99, nml = gamess_input_files)
else if (Molpro .and. (.not. GAMESS)) then
 m = 2
 write(*,*) "Working with Molpro input files"
 read(unit = 99, nml = molpro_input_files)
else if (ORCA .and. (.not. GAMESS) .and. (.not. Molpro)) then ! Added
 m =3                                                         ! Added
 write(*,*) "Working with ORCA input files"                   ! Added
 read(unit = 99, nml = orca_input_files)                      ! Added
else if ((GAMESS .and. Molpro)) then
  write(*,*) "Error! You chose both GAMESS and Molpro to be TRUE. &
             Check your input and return"
  call exit
else if ((.not. GAMESS) .and. (.not. Molpro) .and. (.not. ORCA)) then ! Added
  write(*,*) "Error! You chose GAMESS, Molpro and ORCA all to be FALSE.  &
             Check your input and return" ! Added
  call exit
end if 

end subroutine read_input_file
!----------------------------------------------------------------------
end module irc_data
