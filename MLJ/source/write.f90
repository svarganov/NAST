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
!***************************** Now print them to the output file
write(66,'(1X,A,6X,I4,5(A,I2.2))') 'Date and time :', &
date(1), '/', date(2), '/', date(3), ' ', date(5), ':', &
date(6), ':', date(7)

write(66, 12)
write(66, 13)


10 format(/,4X,'*****************************************', &
           /4X,'~~~ MLJ: Marcus-Levich-Jortner theory ~~~',  &
           /4X,'           ~~~~~ v. 1.0 ~~~~~ ',  &
           /4X,'*****************************************',/)

11 format('A')

12 format(/,1X,'   MLJ is experimental code to calculate Fermi Golden Rule-lie rate constants', &
          /,1X,'   in the form of the Marcus-Levich-Jortner theory (MLJ). MLJ theory is an', &
          /,1X,'   extension of the Marcus theory to account for a non-classical conribution', &
          /,1X,'   of high frequency modes for which hv << kT is not true.', &
          /,1X,'   The code distribution is not allowed without written permission from the authors.')

13 format(/,/,1X,'Authors in alphabetical order:',&
          /,1X,'-------------------------------------------------------------', &
          /,1X,'Claudia E. Avalos              New York University, New York',&
          /,1X,'Ilya D. Dergachev,             New York University, New York',&
          /,1X,'Vsevolod D. Dergachev,         University of Nevada, Reno',&
          /,1X,'Yash Patel,                    New York University, New York',&
          /,1X,'Mitra Rooein,                  Uppsala University, Uppsala',&
          /,1X,'Sergey A. Varganov,            University of Nevada, Reno',&
          /,1X,'Philip S. Weiss                New York University, New York',&
          /1X,'-------------------------------------------------------------',/)

end subroutine write_header
!-----------------------------
end module
