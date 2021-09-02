module rot_module
 use precision_module, only : wp
 use pes_module,  only : f_inertX
 use constants_module
 use input_module

 implicit none
 public

contains
subroutine asym_top(rotR,rotX,rotTP)
  implicit none
  integer                      :: i, k
  real(wp), dimension(:), intent(inout) :: rotR, rotX, rotTP
 
rotR     = 0.0d0
rotX     = 0.0d0
rotTP    = 0.0d0

write(66,'(A)') '.......rotational.'

! Note that in the book Unimolecular Kinetics by Green N.J.B.
! There factor of 1/pi is missing in partition function (eq. 81).
! As a result density of states in eq. 83 has an extra pi which is wrong.

! Rotational DOS for reactant.

do i=1,maxn
  rotR(i) = 4.0d0*sqrt(2.0d0*dble(i)/autocm*&
            inertR(1)*inertR(2)*inertR(3))/autocm
enddo

! Rotational DOS for turning points.

! Rotational DOS are taken from MECP.
! maxn-1 here means that one 1cm-1 should be already in rxn coord 
! and the lowest rate is k(2) and the highest k(maxn)
! the total number of k is equal to maxn-1.

do i=1,maxn-1
  rotTP(i) = 4.0d0*sqrt(2.0d0*dble(i)/autocm*&
             inertX(1)*inertX(2)*inertX(3))/autocm
enddo

! Rotational DOS at the MECP; max[i]+binX=maxn due to conservation of energy.

do i=1,maxn-binX
  rotX(i+binX) = 4.0d0*sqrt(2.0d0*dble(i)/autocm*&
                 inertX(1)*inertX(2)*inertX(3))/autocm
enddo
 
end subroutine asym_top

subroutine rev_asym_top(rotP)
 implicit none
 integer                      :: i
 real(wp), dimension(:), intent(out) :: rotP
 
 rotP  = 0.0d0

 write(66,'(A)') '...........product rotational.'
 do i=1, maxn+binGAP
   rotP(i) = 4.0d0*sqrt(2.0d0*dble(i)/autocm*&
             inertP(1)*inertP(2)*inertP(3))/autocm
 enddo

end subroutine rev_asym_top
end module rot_module
