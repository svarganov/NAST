module irc_driver

use fit
use get_constants
use irc_data
use write_irc_data
use GAMESS_read_module

implicit none
contains
!---------------------------------------------
subroutine driver_irc(input_file,nast_file)
implicit none

character(len=*), intent(in)    :: input_file, nast_file
!---------------------------------------------
integer                         :: i, j, k, n_prod

call read_input_file(input_file)

select case(m)

case(1)

 open(unit = 11, status = "old", action = "read", file = g_irc_output_reac)
 open(unit = 12, status = "old", action = "read", file = g_irc_output_prod)

case(2)
 
 open(unit = 11, status = "old", action = "read", file = m_irc_output)

end select 

select case(m)

case(1)

call GAMESS_read_atomic_data(11,n1,numbr,geom_reac,en_reac)
call GAMESS_read_atomic_data(12,n2,numbr,geom_prod,en_prod)

case(2)

!call Molpro_read_atomic_data(11)

end select

allocate ( ref (numbr,3) )
allocate ( irc_reac(n1), irc_prod(n2) )

en_prod = en_prod - en_reac(n1)
en_reac = en_reac - en_reac(n1)

geom_reac = geom_reac*ang2bohr
geom_prod = geom_prod*ang2bohr

irc_reac = 0.0d0
irc_prod = 0.0d0

ref(1:numbr,1:3) = geom_reac(n1,1:numbr,1:3)

do i = 1, n1
  geom_reac(i,1:numbr,1:3) = geom_reac(i,1:numbr,1:3) - ref(1:numbr,1:3)
end do

do i = 1, n2
  geom_prod(i,1:numbr,1:3) = geom_prod(i,1:numbr,1:3) - ref(1:numbr,1:3)
end do

do i = 1, n1
 do j = 1, numbr
  do k = 1, 3
   irc_reac(i) = irc_reac(i) + geom_reac(i,j,k)**2
  end do
 end do
end do

do i = 1, n2
 do j = 1, numbr
  do k = 1, 3
   irc_prod(i) = irc_prod(i) + geom_prod(i,j,k)**2
  end do
 end do
end do

irc_reac = sqrt(irc_reac)
irc_prod = sqrt(irc_prod)

if      (n1 .eq. 2) then
  call linear_polynom_fit (n1,irc_reac,en_reac,l(4),l(5))
else if (n1 .eq. 3) then
  call quadr_polynom_fit  (n1,irc_reac,en_reac,l(3),l(4),l(5))
else if (n1 .eq. 4) then
  call cubic_polynom_fit  (n1,irc_reac,en_reac,l(2),l(3),l(4),l(5))
else
  call quart_polynom_fit  (n1,irc_reac,en_reac,l(1),l(2),l(3),l(4),l(5))
end if

n_prod = count(en_prod .gt. 0)

! n_prod => We need only that points on the product curve, for which the corresponding energies
! are higher or equal to the energy minimum on the reactant curve, which, in turn, is shiftet to zero.

if (n_prod .eq. 1) then
  n_prod = n_prod + 1   ! To fit a line, we need a minimum of two points.
end if

if      (n_prod .eq. 2) then
  call linear_polynom_fit (2,irc_prod(1:2),en_prod(1:2),h(4),h(5))
else if (n_prod .eq. 3) then
  call quadr_polynom_fit  (3,irc_prod(1:3),en_prod(1:3),h(3),h(4),h(5))
else if (n_prod .eq. 4) then
  call cubic_polynom_fit  (4,irc_prod(1:4),en_prod(1:4),h(2),h(3),h(4),h(5))
else
  call quart_polynom_fit  (n_prod,irc_prod(1:n_prod),en_prod(1:n_prod),&
       h(1),h(2),h(3),h(4),h(5))
end if

call write_output(input_file,n_prod)
call write_nast_template(nast_file)

end subroutine driver_irc
!----------------------------------------------
end module
