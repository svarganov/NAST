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

10 format(/,4X,'*******************************************************', &
           /4X,'~~~~~~~ Effective Hessian: part of NAST package ~~~~~~~',  &
           /4X,'                  ~~~~~ v. 2021.1 ~~~~~               ',  &
           /4X,'*******************************************************',/)

11 format('A')

12 format(/,/,1X,'   Effective Hessian is experimental code to construct effective Hessian matrix', &
          /,1X,'   at the minimum energy crossing point (MECP) of the reaction at study.', &
          /,1X,'   Working equations and the theory behind them can be found in the NAST manual,', &
          /,1X,'   distributed with the package.', &
          /,1X,'     As program distributed free of charge, there is no warranty implied.', &
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
!-----------------------------------------------------------------------------------------------
subroutine write_output(input)

use hess_data
!-----------------------------------------------------------------------------------------------
implicit none

character(len=*), intent(in)     :: input
integer                          :: i
double precision                 :: zpe
!------------------------------------------------------------------------------------------------

write(66,12) trim(input)

!~~~ Write symbols, charges and Cartesian coordinates

write(66,'(/,A,/)') 'Coordinates (bohr):'

do I = 1, number_of_atoms
 write(66,13) symbol(i),charge(i),coord(i,1),      &
              coord(i,2),coord(i,3)
end do

!~~~ Write data

write(66,'(/,A)') '      Effective Hessian data     '

write(66,14) gradmean, norm_deltaG
write(66,15) norm_grad1, norm_grad2, lambda
write(66,16) gradslope, intersect_type
write(66,'(/,A,/)') 'Parallel gradient: '

do i = 1, number_of_atoms
 write(66,17) symbol(i),charge(i),              &
              deltaG(i*3-2),deltaG(i*3-1),      &
              deltaG(i*3)
end do

write(66,'(/,A,/)') 'Principal moments of inertia (amu*bohr**2)'
write(66,18) moments_eig(1), moments_eig(2), moments_eig(3)

write(66,'(/,A,10X,A,/)') 'Frequency (cm-1)', 'Reduced mass (amu)'
zpe = 0.0d0

do i = 1, 3*number_of_atoms
 if ( e(i) > 0 .and. freq(i) > 1d-3) then
  write(66,19) freq(i), reduced_mass(i)
  zpe = zpe + freq(i)
 else if ( e(i) < 0 .and. freq(i) > 1d-3) then
  write(66,19) -freq(i), reduced_mass(i)
 else
  write(66,19) 0.0d0, reduced_mass(i)
 end if
end do

write(66,'(/,a,1x,f10.7)') 'Reduced mass along the RC (amu):',reduced_mass_rc

zpe = zpe/2             ! final ZPE in wavenumbers
zpe = zpe*0.002859144   ! final ZPE in kcal/mol

write(66,'(/,a,1x,f12.3)') 'zpe in kcal mol-1:',zpe

write(66,'(/,a)') 'Effective Hessian terminated normally.'
!================= FORMATTING OPTIONS ====================!

12 FORMAT(/,'Opening the ',A,' input file')
13 FORMAT(A10,4X,F4.1,3X,F12.9,3X,F12.9,4X,F12.9)
14 FORMAT(/,'Gradmean = ',ES16.9,/,'Norm of parallel gradient = ',ES16.9)
15 FORMAT('Norm of low-spin gradient = ',ES16.9,/,'Norm of high-spin gradient = ',ES16.9, &
         /,'Lambda = ',E13.6)
16 FORMAT('Gradients dot product = ',ES16.9,' indicates ',A,' intersection')
17 FORMAT(A8,'    ',F3.0,'  ',ES17.10,'   ',ES17.10,'    ',ES17.10)
18 FORMAT(2X,F14.9,1X,F14.9,1X,F14.9)
19 FORMAT(2X,F7.2,20X,F10.7)

end subroutine write_output
!-------------------------------------------------------------------------------
subroutine write_nast_template(input)

use hess_data

implicit none
character(len=*), intent(in)     :: input
integer                          :: i
!--------------------------------------------
open(unit = 44,access='stream',form='formatted',action = "write", file = input)

! Write heading to the output file
write(44,10)
write(44,11)
write(44,55) pack(int(freq),freq > 1d-3)
write(44,'(A,3f16.5)') 'inertX=',moments_eig
write(44,12) reduced_mass_rc, norm_deltaG, gradmean

55 format('freR = ',/,'freX = ',400i5)
10 format('!This is an NAST template input file automatically ',&
         'genereated by the Effective Hessian code.',&
         /,'!The file contains the MECP part of NAST input.',&
         ' Check NAST manual to finish the input file.')

11 format(/,"&keys",1x,"zpe=1",1x,"&end",/,&
         /,"&inputdata")


12 format("inertR = ",/,"enX = ",1x,"enR = ",/,&
         "redmass = ",f10.5,1x,"soc = ",/,"grad = ",f11.6,1x,&
         "gradmean = ",f11.6,/,"maxn = 10000",/,"&end") 

end subroutine write_nast_template
!-----------------------------------
end module
