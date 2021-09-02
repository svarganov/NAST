module write_irc_data
implicit none
contains
!--------------------------------------
subroutine write_header(n)
use M_strings
!--------------------------------------
! The current solution is only for Linux
! as it binds Fortran with Linux shell
! The more appropriate would be connect
! Fortran with C/C++ or Python
!--------------------------------------
implicit none

integer, intent(in)                     :: n
character(len=120)                      :: line
character(len=32)                       :: username
character(len=120), allocatable         :: array(:)
integer                                 :: date(8)

write(n, 10)
call system('uname -a > _tmp')
open(unit=1,file='_tmp')
read(1,'(A)') line
close(1)
call system('rm _tmp')
call split(line,array)
call getlog(username)
write(n, '(1X,A,15X,A)') 'User: ',username
write(n, '(1X,A,15X,A,/,1X,A,13X,A,1X,/,1X,A,2X,A)') 'Host: ',array(2),'System: ',trim(array(15)),&
'Hardware platform: ',array(13)

call date_and_time(values=date)
!***************************** Now print them to the output file
write(n,'(1X,A,6X,I4,5(A,I2.2))') 'Date and time :', &
date(1), '/', date(2), '/', date(3), ' ', date(5), ':', &
date(6), ':', date(7)

write(n, 12)
write(n, 13)


10 format(/,4X,'*******************************************************', &
           /4X,'~~~~~~~~~ Intrinsic Reaction Coordinate fit ~~~~~~~~~~~',  &
           /4X,'                  ~~~~~ v. 2021.1 ~~~~~               ',  &
           /4X,'*******************************************************',/)

11 format('A')

12 format(/,/,1X,'   IRC is experimental code to fit explicit 1-D crossing potentials', &
          /,1X,'              for the spin-forbidden reaction at hand.', &
          /,1X,'   The program is distributed free of charge and there is no implied warranty.', &
          /,1X,'   The code distribution is not allowed without written permission from the authors.')

13 format(/,/,1X,'Authors in alphabetical order:',&
          /,1X,'-------------------------------------------------------------', &
          /,1X,'Ilya D. Dergachev,             University of Nevada, Reno',&
          /,1X,'Vsevolod D. Dergachev          University of Nevada, Reno',&
          /,1X,'Aleksandr O. Lykhin,           University of Nevada, Reno',&
          /,1X,'Robert Mauban,                 University of Nevada, Reno',&
          /,1X,'Mitra Rooein,                  University of Nevada, Reno',&
          /,1X,'Sergey A. Varganov             University of Nevada, Reno',&
          /1X,'-------------------------------------------------------------',/)

end subroutine write_header
!--------------------------------
subroutine write_output(input_file,n_prod)

use irc_data

implicit none

character(len=*), intent(in)    :: input_file
integer,          intent(in)    :: n_prod
integer                         :: i
double precision                :: en_reac_fit(n1), en_prod_fit(n2)

en_reac_fit = 0.0d0
en_prod_fit = 0.0d0

write(66,12) trim(input_file)

! => fitting of reactan irc points.

if ( (l(1) .eq. 0.0d0) .and. (l(2) .eq. 0.0d0) .and. (l(3) .eq. 0.0d0) ) then

  write(66,13) l(4),l(5),irc_reac(n1),irc_reac(1)
  en_reac_fit(1:n1) = l(4)*irc_reac(1:n1) + l(5)

else if ( (l(1) .eq. 0.0d0) .and. (l(2) .eq. 0.0d0)  ) then
  
  write(66,14) l(3),l(4),l(5),irc_reac(n1),irc_reac(1)
  en_reac_fit(1:n1) = l(3)*irc_reac(1:n1)**2 + l(4)*irc_reac(1:n1) + l(5)

else if ( (l(1) .eq. 0.0d0 ) ) then
 
  write(66,15) l(2),l(3),l(4),l(5),irc_reac(n1),irc_reac(1)
  en_reac_fit(1:n1) = l(2)*irc_reac(1:n1)**3 + l(3)*irc_reac(1:n1)**2 + &
                      l(4)*irc_reac(1:n1) + l(5)

else

  write(66,16) l(1),l(2),l(3),l(4),l(5),irc_reac(n1),irc_reac(1)
  en_reac_fit(1:n1) = l(1)*irc_reac(1:n1)**4 + l(2)*irc_reac(1:n1)**3 + &
                      l(3)*irc_reac(1:n1)**2 + l(4)*irc_reac(1:n1) + l(5) 

end if

! => fitting of product irc points.

if ( (h(1) .eq. 0.0d0) .and. (h(2) .eq. 0.0d0) .and. (h(3) .eq. 0.0d0) ) then

  write(66,17) h(4),h(5),irc_prod(n_prod),irc_prod(1)
  en_prod_fit(1:n_prod) = h(4)*irc_prod(1:n_prod) + h(5)

else if ( (h(1) .eq. 0.0d0) .and. (h(2) .eq. 0.0d0)) then
  
  write(66,18) h(3),h(4),h(5),irc_prod(n1),irc_prod(1)
  en_prod_fit(1:n_prod) = h(3)*irc_prod(1:n_prod)**2 + h(4)*irc_prod(1:n_prod) + h(5)

else if ( (h(1) .eq. 0.0d0 ) ) then
 
  write(66,19) h(2),h(3),h(4),h(5),irc_prod(n_prod),irc_prod(1)
  en_prod_fit(1:n_prod) = h(2)*irc_prod(1:n_prod)**3 + h(3)*irc_prod(1:n_prod)**2 + &
                          h(4)*irc_prod(1:n_prod) + h(5)

else

  write(66,20) h(1),h(2),h(3),h(4),h(5),irc_prod(n_prod),irc_prod(1)
  en_prod_fit(1:n_prod) = h(1)*irc_prod(1:n_prod)**4 + h(2)*irc_prod(1:n_prod)**3 + &
                          h(3)*irc_prod(1:n_prod)**2 + h(4)*irc_prod(1:n_prod) + h(5) 

end if

! => Write raw and fitted data to the output file for comparison

write(66,21)

do i = 1, n1
   write(66,22) i,irc_reac(i),en_reac(i),en_reac_fit(i)
end do

write(66,23)

do i = 1, n_prod
   write(66,22) i,irc_prod(i),en_prod(i),en_prod_fit(i)
end do

12 format(/,'Opening the ',A,' input file')

13 format(/,/,1X,'MECP-reactant minimum energy path has been fit to a line:', &
          /,/,10X,'f(x) = a1*x + a2, where',/,/,2X, &
          'a1 = ',F10.5,' a2 = ',F10.5, &
          /,/,10X,'The polynomial has been fit in the range of [',F5.3,',',F5.3,'] bohr.',&
          /,3x,'Warning! Fitting reactant to a line might inidcate that there are not enough IRC points.',&
          /,3x,'Use fitting results with caution. Consider re-running IRC with smaller search step.')

14 format(/,/,1X,'MECP-reactant minimum energy path has been fit to a parabola:', &
          /,/,10X,'f(x) = a1*x^2 + a2*x + a3, where',/,/,2X, &
          'a1 = ',F10.5,' a2 = ',F10.5,' a3 = ',F10.5, &
          /,/,10X,'The polynomial has been fit in the range of [',F5.3,',',F5.3,'] bohr.',&
          /,3x,'Warning! Fitting reactant to parabola might indicate there are not enough IRC points.',&
          /,3x,'Use fitting results with caution. Consider re-running IRC with smaller search step.')

15 format(/,/,1X,'MECP-reactant minimum energy path has been fit to a cubic polynomial:', &
          /,/,10X,'f(x) = a1*x^3 + a2*x^2 + a3*x + a4, where',/,/,2X, &
          'a1 = ',F10.5,' a2 = ',F10.5,' a3 = ',F10.5,' a4 = ',F10.5, &
          /,/,10X,'The polynomial has been fit in the range of [',F5.3,',',F5.3,'] bohr.')

16 format(/,/,1X,'MECP-reactant minimum energy path has been fit to a quartic polynomial:', &
          /,/,10X,'f(x) = a1*x^4 + a2*x^3 + a3*x^2 + a4*x + a5, where',/,/,2X, &
          'a1 = ',F10.5,' a2 = ',F10.5,' a3 = ',F10.5,' a4 = ',F10.5,' a5 = ',F10.5, &
          /,/,10X,'The polynomial has been fit in the range of [',F5.3,',',F5.3,'] bohr.')

17 format(/,/,1X,'MECP-product minimum energy path has been fit to a line:', &
          /,/,10X,'g(x) = b1*x + b2, where',/,/,5X, &
          'b1 = ',F10.5,' b2 = ',F10.5, &
          /,/,10X,'The polynomial has been fit in the range of [',F5.3,',',F5.3,'] bohr.',&
          /,/,3x,'Warning! Fitting to a line means there are only two IRC points on the product minimum energy path',&
          /,3x,'in the region from MECP to point where product energy is equal to energy at reactant minimum.',&
          /,3x,'It might be alright, but please consider re-running IRC with smaller step size.')

18 format(/,/,1X,'MECP-product minimum energy path has been fit to a parabola:', &
          /,/,10X,'g(x) = b1*x^2 + b2*x + b3, where',/,/,5X, &
          'b1 = ',F10.5,' b2 = ',F10.5,' b3 = ',F10.5, &
          /,/,10X,'The polynomial has been fit in the range of [',F5.3,',',F5.3,'] bohr.',&
          /,/,3x,'Warning! Fitting to parabola might indicate there are not enough IRC points',&
          /,3x,'in the region from MECP to point where product energy is equal to energy at reactant minimum.',&
          /,3x,'It might be alright, but please consider re-running IRC with smaller step size.')

19 format(/,/,1X,'MECP-product minimum energy path has been fit to a cubic polynomial:', &
          /,/,10X,'g(x) = b1*x^3 + b2*x^2 + b3*x + b4, where',/,/,5X, &
          'b1 = ',F10.5,' b2 = ',F10.5,' b3 = ',F10.5,' b4 = ',F10.5, &
          /,/,10X,'The polynomial has been fit in the range of [',F5.3,',',F5.3,'] bohr.')

20 format(/,/,1X,'MECP-product minimum energy path has been fit to a quartic polynomial:', &
          /,/,10X,'g(x) = b1*x^4 + b2*x^3 + b3*x^2 + b4*x + b5, where',/,/,5X, &
          'b1 = ',F10.5,' b2 = ',F10.5,' b3 = ',F10.5,' b4 = ',F10.5,' b5 = ',F10.5, &
          /,/,10X,'The polynomial has been fit in the range of [',F5.3,',',F5.3,'] bohr.')

21 format(/,/,1x,'MECP-reactant irc points, energies and energies fit from f(x)', &
          /,/,6x,'n',9x,'x, Bohr',8x,'E, a.u.',6x, 'E_fit, a.u.',/)

22 format(5x,i3,5x,f10.5,5x,f10.5,5x,f10.5)

23 format(/,/,1x,'MECP-product irc points, energies and energies fit from g(x)', &
          /,/,6x,'n',9x,'x, Bohr',8x,'E, a.u.',6x, 'E_fit, a.u.',/)

end subroutine write_output
!--------------------------------
subroutine write_nast_template(input)

use irc_data
use get_constants

implicit none
character(len=*), intent(in)     :: input
character(len=19)                :: intersect_type
double precision                 :: grad1, grad2, limitL, limitR, E_limitR
!--------------------------------------------
open(unit = 44,access='stream',form='formatted',action = "write", file = input)

! Write heading to the output file

grad1 = 0.0d0
grad2 = 0.0d0
limitL = 0.0d0
limitR = 0.0d0

write(44,10)
call write_header(44)
write(44,'(/,A,f10.5)') 'Value of the reaction coordinate at the MECP: x_M = ', irc_reac(1)

if ( (l(1) .eq. 0.0d0) .and. (l(2) .eq. 0.0d0) .and. (l(3) .eq. 0.0d0) ) then

  write(44,'(/,A,f10.5)') "Derivative of the reactant polynomial, df(x)/dx =  ", l(4)
  
  grad1 = l(4)
  write(44,'(/,A,f10.5)') "df(x)/dx | x_M = ", grad1

else if ( (l(1) .eq. 0.0d0) .and. (l(2) .eq. 0.0d0)  ) then
  
  write(44,'(/,A,f10.5,A,f10.5)') "Derivative of the reactant polynomial, df(x)/dx =  ", &
                          2.0d0*l(3),'x + ',l(4)

  grad1 = 2.0d0*l(3)*irc_reac(1) + l(4)
  write(44,'(/,A,f10.5)') "df(x)/dx | x_M = ", grad1

else if ( (l(1) .eq. 0.0d0 ) ) then
 
  write(44,'(/,A,f10.5,A,f10.5,A,f10.5)') "Derivative of the reactant polynomial, df(x)/dx =  ", &
                       3.0d0*l(2),'x^2 + ',2.0d0*l(3),'x + ',l(4)
  
  grad1 = 3.0d0*l(2)*(irc_reac(1)**2) + 2.0d0*l(3)*irc_reac(1) + l(4)
  write(44,'(/,A,f10.5)') "df(x)/dx | x_M = ", grad1

else

  write(44,'(/,A,f10.5,A,f10.5,A,f10.5,A,f10.5)') "Derivative of the reactant polynomial, df(x)/dx =  ", &
                       4.0d0*l(1),'x^3 + ',3.0d0*l(2),'x^2 + ',2.0d0*l(3),'x + ',l(4)
  
  grad1 = 4.0d0*l(1)*(irc_reac(1)**3) + 3.0d0*l(2)*(irc_reac(1)**2) + 2.0d0*l(3)*irc_reac(1) + l(4)
  write(44,'(/,A,f10.5)') "df(x)/dx | x_M = ", grad1

end if

if ( (h(1) .eq. 0.0d0) .and. (h(2) .eq. 0.0d0) .and. (h(3) .eq. 0.0d0) ) then

  write(44,'(/,A,f10.5)') "Derivative of the product polynomial, dg(x)/dx =  ", h(4)
  
  grad2 = h(4)
  write(44,'(/,A,f10.5)') "dg(x)/dx | x_M = ", grad2

else if ( (h(1) .eq. 0.0d0) .and. (h(2) .eq. 0.0d0)  ) then
  
  write(44,'(/,A,f10.5,A,f10.5)') "Derivative of the product polynomial, dg(x)/dx =  ", &
                          2.0d0*h(3),'x + ',h(4)
  
  grad2 = 2.0d0*h(3)*irc_reac(1) + h(4)
  write(44,'(/,A,f10.5)') "dg(x)/dx | x_M = ", grad2

else if ( (h(1) .eq. 0.0d0 ) ) then
 
  write(44,'(/,A,f10.5,A,f10.5,A,f10.5)') "Derivative of the product polynomial, dg(x)/dx =  ", &
                       3.0d0*h(2),'x^2 + ',2.0d0*h(3),'x + ',h(4)
  
  grad2 = 3.0d0*h(2)*(irc_reac(1)**2) + 2.0d0*h(3)*irc_reac(1) + h(4)
  write(44,'(/,A,f10.5)') "dg(x)/dx | x_M = ", grad2

else

  write(44,'(/,A,f10.5,A,f10.5,A,f10.5,A,f10.5)') "Derivative of the product polynomial, dg(x)/dx =  ", &
                       4.0d0*h(1),'x^3 + ',3.0d0*h(2),'x^2 + ',2.0d0*h(3),'x + ',h(4)
  
  grad2 = 4.0d0*h(1)*(irc_reac(1)**3) + 3.0d0*h(2)*(irc_reac(1)**2) + 2.0d0*h(3)*irc_reac(1) + h(4)
  write(44,'(/,A,f10.5)') "dg(x)/dx | x_M = ", grad2

end if

write(44,'(/,A,f10.5)') "Derivatives product at x_M: ", grad1*grad2

if (grad1*grad2 < 0) then
   intersect_type = 'peaked intersection'
else
   intersect_type = 'sloped intersection'
end if

write(44,'(/,A,A)') "Derivatives product indicates ", intersect_type

write(44,11)

if (intersect_type == 'peaked intersection') then
   limitR = irc_prod(n2)
else
   limitR = 2.0d0*irc_reac(1)
   E_limitR = l(1)*limitR**4 + l(2)*limitR**3 + &
              l(3)*limitR**2 + l(4)*limitR + l(5)
end if


write(44,12) l(1),l(2),l(3),l(4),l(5), &
             h(1),h(2),h(3),h(4),h(5), &
             limitL, limitR 

write(44,13)

if (intersect_type == 'peaked intersection') then
   write(44,14) irc_prod(n2)
else
   write (44,15) irc_reac(1),en_reac(1),en_reac(1)*au2wavenumber, &
                 limitR,E_limitR,E_limitR*au2wavenumber
end if

10 format('=======',/,'=======',/,'=======', &
          /,/,'This is an NAST template input file automatically ',&
         'genereated by the irc fit code.',&
         /,'The file contains the polynomial coefficients',&
         ' needed to run Zhu-Nakamura approach.', &
         /,'Check NAST manual for more details.')

11 format (/,/,'======================================================== ', &
           /,/,"     NAST &polynomial input group for Zhu-Nakamura", &
           /,/,'========================================================')

12 format(/,"&polynomials",/,"hs4 = ",f10.7,/,"hs3 = ",f10.7,/,&
         "hs2 = ",f10.7,/,"hs1 = ",f10.7,/, "hs0 = ",f10.7,/, &
         "ls4 = ",f10.7,/,"ls3 = ",f10.7,/, "ls2 = ",f10.7,/, &
         "ls1 = ",f10.7,/,"ls0 = ",f10.7,/,"limitL = ",f5.2,/, &
         "limitR = ",f5.2,/,"&end")

13 format (/,'=============================== ', &
           /,/,"        limitL & limitR", &
           /,/,'================================', &
           /,/,'limitL and limitR are the left and right limits on the x-axis,', &
           /,'where x is the 1-D reaction coordinate.', &
           /,/,"limitL and limitR given above in the &polynomial group are authors' suggestion.", &
           /,'Feel free to change these defaults at your own risk and understanding.', &
           /,'An explanation of our choice is given below.', &
           /,/,'limitR must be greater than limitL. limitL should be chosen at the reactant minimum.',&
           /,'Therefore, limitL is effectively zero: limitL = 0.0.', &
           /,/,'The choice of limitR depends on the intersection type.')

14 format  (/,'For the peaked intersection, limitR should be chosen as the product minimum,' &
           /,'which is the final point on the MECP -> product path.', &
           /,'In this example, the product minimum is at',f5.2,'.',/,/,'Fly with us again.')

15 format  (/,'For the sloped intersection, limitR should be chosen so that', &
           /,'E1 = f(limitR) is well above MECP to ensure convergence.', &
           /,/,'For that, our default choice is limitR = 2.0 * x_M,', &
           /,'where x_M is the coordinate of MECP.', &
           /,/,'In this example, x_M was ',f5.2,', E1(x_M) = ',f10.7,' a.u. (',f11.2,' cm-1 ).', &
           /,'Then, suggested limitR = ',f5.2,', E1(limitR) = ',f10.7,' a.u. (',f11.2,' cm-1 ).',/,/,'Fly with us again.')

end subroutine write_nast_template
!-----------------------------------
end module
