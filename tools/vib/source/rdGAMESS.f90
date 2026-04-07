module GAMESS_read_module

use M_strings
use hess_data
implicit none

contains
!---------------------------------------------------------------------------------
subroutine GAMESS_read_hess(file_id,Hess)

implicit none
integer, intent(in)                                             :: file_id
double precision, dimension(:,:), allocatable, intent(out)      :: Hess
!----------------------------------------------------------------------------------
character(len=80)                                               :: line
character(len=80), allocatable                                  :: array(:)
integer                                                         :: i, j, k, l, last_number
!----------------------------------------------------------------------------------
allocate  (Hess (3*number_of_atoms,3*number_of_atoms) )
!----------------------------------------------------------------------------------
! Read Hessian
                                                   ! Need to shift position to the
                                                   ! line containing $HESS

do while ((index(line, '$HESS') == 0) .and. (index(line, '$hess') == 0))
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

end subroutine GAMESS_read_hess
!-------------------------------------------------------------------------------------
end module GAMESS_read_module
