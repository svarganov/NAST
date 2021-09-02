module diag_module
  use precision_module, only : wp
  use error_module

  implicit none
  public

contains
!------------------------------------------------------------------------

!*****************************************************************************        
!subroutine adiabatization(ediab1,ediab2,h12,h13,h14,eig)
!== Old version of the subroutine. Employs diagonalization of 4x4 matrix
!== with specific matrix elements - h12, h13, h14. Coded to model singlet-triplet
!== Intersystem Crossing.
!*****************************************************************************
!
!  integer           ::  OK, INFO
!  complex(wp)    ::  A(4,4), B(4), DUMMY(1,1), WORK(8)
!  complex(wp), intent(in)     ::  h12, h13, h14
!  real(wp), intent(in)   ::  ediab1, ediab2
!  real(wp), intent(out)  ::  eig(4)
!
!! define elements of matrix A for singlet-triplet intersections
!! define diagonal elements of matrix A
!      A(1,1)=ediab1
!      A(2,2)=ediab2
!      A(3,3)=ediab2
!      A(4,4)=ediab2
!! define upper-right elements of matrix A
!      A(1,2)=h12
!      A(1,3)=h13
!      A(1,4)=h14
!      A(2,3)=(0.0d0, 0.0d0)
!      A(2,4)=(0.0d0, 0.0d0)
!      A(3,4)=(0.0d0, 0.0d0)
!! define upper-right elements of matrix A
!      A(2,1)=DCONJG(A(1,2))
!      A(3,1)=DCONJG(A(1,3))
!      A(3,2)=DCONJG(A(2,3))
!      A(4,1)=DCONJG(A(1,4))
!      A(4,2)=DCONJG(A(2,4))
!      A(4,3)=DCONJG(A(3,4))
!
!! diagonalize matrix A and determine eigenvalues (adiabatic energies)
!      CALL ZGEEV('N', 'N', 4, A, 4, B, DUMMY, &
!                1, DUMMY, 1, WORK, 8, WORK, OK)
!
!      eig=real(B)
!
! sort adiabatic energies in increasing ('I') order
!      CALL DLASRT('I',4,eig,INFO)
!
!      IF (OK .ne. 0) THEN
!          call error_message(3)!error code in parentheses   
!      ENDIF
!
subroutine adiabatization(ediab1,ediab2,eig)

use sort
use input_module, only : soc 

real(wp), intent(in)  :: ediab1, ediab2
real(wp), intent(out) :: eig(2)
real(wp)              :: A(2,2)
integer               :: l, inf
real(wp)              :: work(2*(3+2/2))

A(1,1) = ediab1
A(2,2) = ediab2
A(1,2) = soc
A(2,1) = A(1,2)

l = 2*(3+2/2)

call dsyev('V','U',2,A,2,eig,work,l,inf)

! sort adiabatic energies in increasing ('I') order
call insert_sort(eig)

end subroutine adiabatization

!---------------------------------------------------------------------
end module diag_module
