module QChem_read_module

use M_strings
use hess_data

implicit none

contains
!----------------------------------------------------------------------------------
subroutine QChem_read_hess(file_id,Hess)

implicit none
integer, intent(in)                                           :: file_id
double precision, dimension(:,:), allocatable, intent(out)    :: Hess
!----------------------------------------------------------------------------------
character(len=80)                                             :: line
character(len=80), allocatable                                :: array(:)
double precision, dimension(:), allocatable                   :: deployed
integer                                                       :: i, j, m, N, last_number, column
!----------------------------------------------------------------------------------
allocate  (Hess (3*number_of_atoms,3*number_of_atoms) )
!----------------------------------------------------------------------------------
! Read Hessian

Hess = 0.0d0

N = int(3*real(number_of_atoms)*((3*real(number_of_atoms) - 1)/2 + 1))

allocate  (deployed (N) )
deployed = 0.0d0

do while (index(line, 'Cartesian Force Constants') == 0)
 read(file_id,'(A)') line
end do                                             ! Now on the Force Constants line

do i = 1, int(ceiling(real(N)/5))

  if (i*5 <= N) then
    column = 5
  else
    column = MOD(N,5)              
  end if
  read(file_id,'(A)') line
  call split(line, array)
  if (column == 5) then
    read(array(1),*) deployed((i-1)*5 + 1)
    read(array(2),*) deployed((i-1)*5 + 2)
    read(array(3),*) deployed((i-1)*5 + 3)
    read(array(4),*) deployed((i-1)*5 + 4)
    read(array(5),*) deployed((i-1)*5 + 5)
  else if (column == 1) then
    read(array(1),*) deployed((i-1)*5 + 1)
  else if (column == 2) then
    read(array(1),*) deployed((i-1)*5 + 1)
    read(array(2),*) deployed((i-1)*5 + 2)
  else if (column == 3) then
    read(array(1),*) deployed((i-1)*5 + 1)
    read(array(2),*) deployed((i-1)*5 + 2)
    read(array(3),*) deployed((i-1)*5 + 3)
  else if (column == 4) then
    read(array(1),*) deployed((i-1)*5 + 1)
    read(array(2),*) deployed((i-1)*5 + 2)
    read(array(3),*) deployed((i-1)*5 + 3)
    read(array(4),*) deployed((i-1)*5 + 4)
  end if
end do

rewind file_id

m = 0
do i=1, 3*number_of_atoms
  do j=1, i
    Hess(j,i)=deployed(j+m)
    Hess(i,j)=Hess(j,i)
  end do
  m=m+i
end do
end subroutine QChem_read_hess
!---------------------------------------------
end module QChem_read_module
