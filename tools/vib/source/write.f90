module write_data
implicit none
contains
!--------------------------------------
subroutine write_header()
use M_strings
!--------------------------------------
! The current solution is only for Linux
! as it binds Fortran with Linux shell
! The more appropriate would be connect
! Fortran with C/C++ or Python
!--------------------------------------
implicit none

character(len=120)                      :: line
character(len=32)                       :: username
character(len=120), allocatable         :: array(:)
integer                                 :: date(8)

write(66, 10)
call system('uname -a > _tmp')
open(unit=1,file='_tmp')
read(1,'(A)') line
close(1)                       
call system('rm _tmp')
call split(line,array)
call getlog(username)
write(66, '(1X,A,15X,A)') 'User: ',username
write(66, '(1X,A,15X,A,/,1X,A,13X,A,1X,/,1X,A,2X,A)') 'Host: ',array(2),'System: ',trim(array(15)),&
'Hardware platform: ',array(13)

call date_and_time(values=date)


write(66,'(1X,A,6X,I4,5(A,I2.2))') 'Date and time :', &
date(1), '/', date(2), '/', date(3), ' ', date(5), ':', &
date(6), ':', date(7)

write(66, 12)
write(66, 13)

10 format(/,4X,'*******************************************************', &
           /4X,'~~~~~~~~~~~~ Vib: part of the NAST package ~~~~~~~~~~~~',  &
           /4X,'                   ~~~~~ v. 1.0 ~~~~~               ',  &
           /4X,'*******************************************************',/)

11 format('A')

12 format(/,1X,'   Vib is experimental code to calculate Huang-Rhys factors', &
          /,1X,'   and reorganization energies for Fermi Golden Rule calculations with NAST', &
          /,1X,'   using the Marcus-Levich-Jortner theory. Working equations and the theory', &
          /,1X,'   behind them can be found in the NAST manual, distributed with the package.', &
          /,1X,'   As program distributed free of charge, there is no warranty implied.', &
          /,1X,'   The code distribution is not allowed without written permission from the authors.')

13 format(/,/,1X,'Authors in alphabetical order:',&
          /,1X,'-------------------------------------------------------------', &
          /,1X,'Claudia E. Avalos              New York University, New York',&
          /,1X,'Ilya D. Dergachev,             New York University, New York',&
          /,1X,'Vsevolod D. Dergachev,         University of Nevada, Reno',&
          /,1X,'Yash Patel,                    New York University, New York',&
          /,1X,'Mitra Rooein,                  Uppsala University, Uppsala',&
          /,1X,'Sergey A. Varganov             University of Nevada, Reno',&
          /,1X,'Philip S. Weiss                New York University, New York',&
          /1X,'-------------------------------------------------------------',/)

end subroutine write_header

!-----------------------------------------------------------------------------------------------
subroutine write_output(input)

use hess_data
!-----------------------------------------------------------------------------------------------
implicit none

character(len=*), intent(in)     :: input
integer                          :: i
!------------------------------------------------------------------------------------------------

write(66,12) trim(input)
!write(66,*) ""
write(66,'(/,A)') '...... Performing vib calcualtions'

!write(66,*) ""
write(66,'(/,A)') '....... Vib data calculated'
write(66,*) ""

write(66,13) "Vibrational mode","Frequency, cm-1","Huang-Rhys factor", "Reorganization energy, a.u."
write(66,*) ""

do i = 1, 3*number_of_atoms
 if ( e(i) > 0 .and. freq(i) > 1d-3) then
  write(66,14) i, freq(i), HRF(i), reorganization_energy(i)
 else if ( e(i) < 0 .and. freq(i) > 1d-3) then
  write(66,14) i, -freq(i), HRF(i), reorganization_energy(i)
 else
  write(66,14) i, 0.0d0
 end if
end do

write(66,'(/,a)') 'Vib terminated normally.'

!================= FORMATTING OPTIONS ====================!

12 FORMAT(/,'..... Opening the ',A,' input file')
13 FORMAT(1x,A,2x,A,2x,A,2x,A)
14 FORMAT(I6,f24.2,2es20.3)

end subroutine write_output
!-------------------------------------------------------------------------------
subroutine write_mlj_template(input)

use hess_data

implicit none
character(len=*), intent(in)     :: input
integer                          :: i
!--------------------------------------------
open(unit = 44,access='stream',form='formatted',action = "write", file = input)

write(44,10)
write(44,55)
write(44,12) pack(freq,freq > 1d-3)
write(44,13) pack(HRF(7:),HRF(7:) > 0.00d0)
write(44,14) pack(reorganization_energy(7:),reorganization_energy(7:) > 0.00d0)
write(44,15)

10 format('This is a file automatically genereated by the Vib code.',&
         /,'The file contains Huang-Rhys factors and reorganization energies',&
         ' and can be used as an input file for the MLJ program.',&
         ' Check NAST manual to finish the input file.'/)

55 format('&inputdata',/,'enR = ',/,'enP = ',/,'V = ')

12 format('freR = ',400f7.1)
13 format('S = ',400es10.3)
14 format('lambda_vib = ',400es10.3)
15 format('lambda_ext = ',/,'n_vib = -1',/,'T1 = ',/,'Tstep = ',/,'T2 = ',/,'trsh = 4.796',/,'&end')

end subroutine write_mlj_template
!-----------------------------------
end module
