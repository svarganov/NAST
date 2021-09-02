module GAMESS_read_module

use M_strings
use hess_data
implicit none

contains
!------------------
subroutine GAMESS_read_number_of_atoms(file_id)

implicit none
integer, intent(in)                                    :: file_id
!-----------------------------------------------------------------
character(len=80)                                      :: line
character(len=3)                                       :: point_g
!-----------------------------------------------------------------

                                                  ! Need to shift position to the line containing $DATA

read(file_id,'(A)') line                          ! Check from first line

do while (index(line, '$DATA') == 0)
 read(file_id,'(A)') line
end do                                            ! Now on the $DATA line

                                                  ! Now need to shift position to the line containing $END
read(file_id,'(A)') line                          ! The line next to $DATA: mol. spec.
read(file_id,'(A)') line                          ! The line next to mol. spec.:
                                                  ! the point group of the molecule
point_g = trim(adjustl(line))                     ! Assign the point group of the molecule
if (point_g /= 'C1') then                         ! Check the point group
  read(file_id,'(A)') line                        ! Skip one line if not C1
end if

number_of_atoms = 0

do while (index(line, '$END') == 0)
 read(file_id,'(A)') line
 number_of_atoms = number_of_atoms + 1
end do

number_of_atoms = number_of_atoms - 1             ! Looks like DO WHILE counts the extra $END
                                                  ! line which is not the atom containing line
rewind file_id                                    ! Positions the file associated with the
                                                  ! specified unit to its initial point
end subroutine GAMESS_read_number_of_atoms
!--------------------------------------------------------------------------
subroutine GAMESS_read_atomic_data(file_id)

implicit none
integer, intent(in)                                             :: file_id
!-------------------------------------------------------------------------
character(len=80)                                               :: line
character(len=3)                                                :: point_g
character(len=80), allocatable                                  :: array(:)
integer                                                         :: i
!-------------------------------------------------------------------------
allocate  (symbol (number_of_atoms) )
allocate  (charge (number_of_atoms) )
allocate  (coord (number_of_atoms,3) )
!-------------------------------------------------------------------------
                                                   ! Need to shift position to the
                                                   ! line containing $DATA

do while (index(line, '$DATA') == 0)
 read(file_id,'(A)') line
end do                                             ! Now on the $DATA line

read(file_id,'(A)') line                           ! The line next to $DATA: mol. spec.

read(file_id,'(A)') line                           ! The line next to mol. spec.:
                                                   ! the point group of the molecule
point_g = trim(adjustl(line))                      ! Assign the point group of the molecule

if (point_g /= 'C1') then                          ! Check the point group
 read(file_id,'(A)') line                          ! Skip one line if not C1
end if

do i =1, number_of_atoms
 read(file_id,'(A)') line
 call split(line, array)
 read(array(1),*) symbol(i)
 read(array(2),*) charge(i)
 read(array(3),*) coord(i,1)
 read(array(4),*) coord(i,2)
 read(array(5),*) coord(i,3)
end do

end subroutine GAMESS_read_atomic_data
!---------------------------------------------------------------------------------
subroutine GAMESS_read_grad_hess(file_id,grad,Hess)

implicit none
integer, intent(in)                                             :: file_id
double precision, dimension(:), allocatable, intent(out)        :: grad
double precision, dimension(:,:), allocatable, intent(out)      :: Hess
!----------------------------------------------------------------------------------
character(len=80)                                               :: line
character(len=80), allocatable                                  :: array(:)
integer                                                         :: i, j, k, l, last_number
!----------------------------------------------------------------------------------
allocate  (grad    (3*number_of_atoms) )
allocate  (Hess (3*number_of_atoms,3*number_of_atoms) )
!----------------------------------------------------------------------------------
! Read gradient
                                                   ! Need to shift position to the
                                                   ! line containing $GRAD

do while (index(line, '$GRAD') == 0)
 read(file_id,'(A)') line
end do                                             ! Now on the $GRAD line

read(file_id,'(A)') line                           ! Skip extra line containing
                                                   ! converged energy


do i = 1, number_of_atoms
 read(file_id,'(A)') line
 call split(line, array)
 read(array(3),*) grad(i*3-2)
 read(array(4),*) grad(i*3-1)
 read(array(5),*) grad(i*3)
end do

rewind file_id
!------------------------------------------------------------------------------------
! Read Hessian
                                                   ! Need to shift position to the
                                                   ! line containing $HESS

do while (index(line, '$HESS') == 0)
 read(file_id,'(A)') line
end do                                             ! Now on the $HESS line

read(file_id,'(A)') line                           ! Now on the line of
                                                   ! converged energy

do i = 1, 3*number_of_atoms
 do j = 1, ceiling(3*real(number_of_atoms)/5)
  read(file_id,'(A)') line
  if (j /= ceiling(3*real(number_of_atoms)/5)) then
   last_number = 5
  else
   last_number = 5-5*ceiling(3*real(number_of_atoms)/5) &
                   + 3*number_of_atoms
  end if

  do l = 1, last_number
   k = j*5 + l - 5
   read(line((15*(l-1)+6):(15*l+5)),*) Hess(i,k)
  end do

 end do
end do

end subroutine GAMESS_read_grad_hess
!-------------------------------------------------------------------------------------
end module GAMESS_read_module
