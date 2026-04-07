module Molpro_read_module

use M_strings
use hess_data

implicit none

contains
!----------------------------------------
subroutine Molpro_read_hess(file_id,Hess)

implicit none
integer, intent(in)                                           :: file_id
double precision, dimension(:,:), allocatable, intent(out)    :: Hess
!----------------------------------------------------------------------------------
character(len=80)                                             :: line
character(len=80), allocatable                                :: array(:)
integer                                                       :: i, j, k, l, last_number, &
                                                                 number_of_blocks, column
!----------------------------------------------------------------------------------
allocate  (Hess (3*number_of_atoms,3*number_of_atoms) )
!----------------------------------------------------------------------------------
! Read Hessian
                                                   ! Need to shift position to the
                                                   ! line containing Force Constants
do while ((index(line, 'Force Constants') == 0) .and. (index(line, 'FORCE CONSTANTS') == 0))
 read(file_id,'(A)') line
end do                                             ! Now on the Force Constants line

read(file_id,'(A)') line                          ! Skip header (blocks separation line)

do i = 1, 3*number_of_atoms
 do J = 1, 3*number_of_atoms
    Hess(i,j) = 0.0d0
 end do
end do

number_of_blocks = int(ceiling(3*real(number_of_atoms)/5))

do i = 1, number_of_blocks
   if (i*5 <= 3*number_of_atoms) then
       column = 5
   else
       column = mod(3*number_of_atoms,5)              
   end if
       
   do j = 1, 3*number_of_atoms - 5*(i-1)
      read(file_id,'(A)') line
      call split(line, array)
      if (column == 5) then
        if (j < 5) then
          if (j == 1) then
            read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
          else if (j == 2) then
            read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
            read(array(3),*) Hess(j + (i-1)*5,(i-1)*5 + 2)
          else if (j == 3) then
            read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
            read(array(3),*) Hess(j + (i-1)*5,(i-1)*5 + 2)
            read(array(4),*) Hess(j + (i-1)*5,(i-1)*5 + 3)
          else if (j == 4) then
            read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
            read(array(3),*) Hess(j + (i-1)*5,(i-1)*5 + 2)
            read(array(4),*) Hess(j + (i-1)*5,(i-1)*5 + 3)
            read(array(5),*) Hess(j + (i-1)*5,(i-1)*5 + 4)
          end if
        else
          read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
          read(array(3),*) Hess(j + (i-1)*5,(i-1)*5 + 2)
          read(array(4),*) Hess(j + (i-1)*5,(i-1)*5 + 3)
          read(array(5),*) Hess(j + (i-1)*5,(i-1)*5 + 4)
          read(array(6),*) Hess(j + (i-1)*5,(i-1)*5 + 5)
        end if  
      else if (column == 1) then
        read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
      else if (column == 2) then
        if (j == 1) then
          read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
        else if (J == 2) then
          read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
          read(array(3),*) Hess(j + (i-1)*5,(i-1)*5 + 2)
        end if
      else if (column == 3) then
        if (j == 1) then
          read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
        else if (j == 2) then
          read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
          read(array(3),*) Hess(j + (i-1)*5,(i-1)*5 + 2)
        else if (j == 3) then
          read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
          read(array(3),*) Hess(j + (i-1)*5,(i-1)*5 + 2)
          read(array(4),*) Hess(j + (i-1)*5,(i-1)*5 + 3)
        end if
      else if (column == 4) then
        if (j == 1) then
          read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
        else if (j == 2) then
          read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
          read(array(3),*) Hess(j + (i-1)*5,(i-1)*5 + 2)
        else if (j == 3) then
          read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
          read(array(3),*) Hess(j + (i-1)*5,(i-1)*5 + 2)
          read(array(4),*) Hess(j + (i-1)*5,(i-1)*5 + 3)
        else if (j == 4) then
          read(array(2),*) Hess(j + (i-1)*5,(i-1)*5 + 1)
          read(array(3),*) Hess(j + (i-1)*5,(i-1)*5 + 2)
          read(array(4),*) Hess(j + (i-1)*5,(i-1)*5 + 3)
          read(array(5),*) Hess(j + (i-1)*5,(i-1)*5 + 4)
        end if
      end if
   end do
   read(file_id,'(A)') line                       ! Skip header (blocks separation line)
end do

do i = 1, 3*number_of_atoms
 do j = 1, 3*number_of_atoms
   if (j /= i) then
     Hess(i,j) = Hess(j,i)
   end if
 end do
end do
print*,Hess
 
end subroutine Molpro_read_hess
!-----------------------------------------------
end module Molpro_read_module
