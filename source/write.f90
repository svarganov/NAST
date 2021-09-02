module write_nast_data
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
! character(len=120)                      :: array(40)
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


10 format(/,4X,'*******************************************************', &
           /4X,'~~~~~~~ NAST: Nonadiabatic Statistical Theory ~~~~~~~~',  &
           /4X,'                  ~~~~~ v. 2021.1 ~~~~~               ',  &
           /4X,'*******************************************************',/)

11 format('A')

12 format(/,/,1X,'   NAST is experimental code to calculate rate constants', &
          /,1X,'   and transition probabilities for the spin-forbidden chemical reactions', &
          /,1X,'   using equations of nonadiabatic statistical theory.', &
          /,1X,'   As program distributed free of charge, there is no warranty implied.', &
          /,1X,'   The code distribution is not allowed without written permission from the authors.')

13 format(/,/,1X,'Authors in alphabetical order:',&
          /,1X,'-------------------------------------------------------------', &
          /,1X,'Ilya D. Dergachev,             University of Nevada, Reno',&
          /,1X,'Vsevolod D. Dergachev,         University of Nevada, Reno',&
          /,1X,'Aleksandr O. Lykhin,           University of Nevada, Reno',&
          /,1X,'Robert Mauban,                 University of Nevada, Reno',&
          /,1X,'Mitra Rooein,                  University of Nevada, Reno',&
          /,1X,'Sergey A. Varganov             University of Nevada, Reno',&
          /1X,'-------------------------------------------------------------',/)

end subroutine write_header
!-----------------------------
subroutine write_options()
use input_module
implicit none


write(66,'(2X,A)') '---------------------------------------------------------'
write(66,'(8X,A)') 'NAST control parameters and related data'
write(66, '(/,5X,A,I1,5X,A,L,5X,A,L,5X,A,L)') 'zpe = ',zpe,'sp = ',sp,'zn = ',zn, &
                                                     'solution = ',solution
write(66, '(5X,A,L,5X,A,L,5X,A,L,5X,A,L)') 'tst = ',tst,'printmore = ',printmore,'rev = ',rev,'extern = ',extern
write(66,'(2X,A)') '---------------------------------------------------------'

 if (tst) then
    write(66,'(/,1x,a,/)') 'tst is set true. Therefore, calculating &
traditional Transition State rate constant in this run.'
 end if

write(66,'(4xx,a)') 'Polynomials coefficients will be recognized only &
if zn is set true in the input file.'

if (zpe == 1) then
  write(66,'(/,1x,a)') 'zpe = 1: ZPE correction scheme I (eliminates turning points below ZPE).'
else if (zpe == 2) then
   write(66,'(1x,a)') 'zpe = 2: ZPE correction scheme II (accounts for turning points below ZPE).'
else
   write(66,'(1x,a)') 'zpe = 0: ZPE correction is skipped.'
end if

write(66,'(1x,a,i5,a)') 'Electronic barrier from reactant to MECP is ',int(mecp*autocm),' cm-1'
write(66,'(1x,a,i6,a,/,1x,a,i6,a,/,1x,a,i6,a)') 'ZPE of reactants = ',int(zpeR*autocm), ' cm-1', &
'ZPE of MECP = ',int(zpeX*autocm),' cm-1','ZPE-corrected MECP energy bin = ',binX,' cm-1'

if (rev) then
  write(66,'(1x,a,i6,a)') 'ZPE of product is ',int(zpeP*autocm), ' cm-1.'
  write(66,'(1x,a,i6,a)') 'Gap between reactant and product is  ',binGAP,' cm-1.'
end if

end subroutine write_options
!--------------------------------
end module
