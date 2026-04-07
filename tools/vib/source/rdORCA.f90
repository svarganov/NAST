module ORCA_read_module

use M_strings
use hess_data

implicit none

contains
!--------------------------------------------------------------------------
subroutine ORCA_read_hess(file_id,Hess)

implicit none
integer, intent(in)                                             :: file_id
double precision, dimension(:,:), allocatable, intent(out)      :: Hess
!----------------------------------------------------------------------------------
character(len=132)                                               :: line
character(len=132), allocatable                                  :: array(:)
integer :: i,j,k,n, column, number_of_five_column_blocks
!----------------------------------------------------------------------------------
allocate  (Hess    (3*number_of_atoms,3*number_of_atoms) )
!----------------------------------------------------------------------------------

rewind file_id

! ORCA's .hess coordinates are already in bohr units, no need for conversion.

! Read Hessian
Hess = 0.0d0

do while (index(line, '$hessian') == 0)
 read(file_id,'(A)') line
end do

read(file_id,'(A)') line
read(file_id,'(A)') line

number_of_five_column_blocks = int(ceiling(3*real(number_of_atoms)/5))

do i = 1, number_of_five_column_blocks

 if (i*5 <= 3*number_of_atoms) then
   column = 5
 else
   column = mod(3*number_of_atoms,5) ! 	mod(x,y)=x-y*INT(x/y)
 end if

 if (column == 5) then
   do j = 1, 3*number_of_atoms
     read(file_id,'(A)') line
     call split(line,array)
     read(array(2),*) Hess(j,(i-1)*5 + 1)
     read(array(3),*) Hess(j,(i-1)*5 + 2)
     read(array(4),*) Hess(j,(i-1)*5 + 3)
     read(array(5),*) Hess(j,(i-1)*5 + 4)
     read(array(6),*) Hess(j,(i-1)*5 + 5)
   end do
 else if (column == 1) then
   do j = 1, 3*number_of_atoms
     read(file_id,'(A)') line
     call split(line,array)
     read(array(2),*) Hess(j,(i-1)*5 + 1)
   end do
 else if (column == 2) then
   do j = 1, 3*number_of_atoms
     read(file_id,'(A)') line
     call split(line,array)
     read(array(2),*) Hess(j,(i-1)*5 + 1)
     read(array(3),*) Hess(j,(i-1)*5 + 2)
   end do
 else if (column == 3) then
   do j = 1, 3*number_of_atoms
     read(file_id,'(A)') line
     call split(line,array)
     read(array(2),*) Hess(j,(i-1)*5 + 1)
     read(array(3),*) Hess(j,(i-1)*5 + 2)
     read(array(4),*) Hess(j,(i-1)*5 + 3)
   end do
 else if (column == 4) then
   do j = 1, 3*number_of_atoms
     read(file_id,'(A)') line
     call split(line,array)
     read(array(2),*) Hess(j,(i-1)*5 + 1)
     read(array(3),*) Hess(j,(i-1)*5 + 2)
     read(array(4),*) Hess(j,(i-1)*5 + 3)
     read(array(5),*) Hess(j,(i-1)*5 + 4)
   end do
 end if

 read(file_id,'(A)') line

end do

end subroutine ORCA_read_hess
!--------------------------------------------------
end module ORCA_read_module
