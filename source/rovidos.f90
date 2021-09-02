module convolution_module
 use precision_module, only : wp
 use error_module
 use input_module

 implicit none
 public

contains
subroutine rovib_dos(vibR,vibX,vibTP,rotR,rotX,rotTP, & 
                                     dosR,dosX,dosTP)

real(wp), dimension(:), intent(in)    :: vibR, vibX, vibTP, rotR, rotX, rotTP 
real(wp), dimension(:), intent(out)   :: dosR, dosX, dosTP
real(wp)                              :: summ
integer                               :: i,j,k

dosR     = 0.0d0
dosX     = 0.0d0
dosTP    = 0.0d0

! Rovibrational DOS for reactant.

write(66,'(A)') '........rovibrational.'
if (solution) then
! No rotational levels in bulk solution.
  dosR = vibR
else
  do i=1,maxn
    summ = 0.0d0
    do j=1,i
      summ = rotR(j)*vibR(i-j+1) + summ ! Accumulate DOS.
    end do
    dosR(i) = summ 
  end do
end if
    
! Rovibrational DOS for turning points (tunneling regime).

if (solution) then
  dosTP = vibTP
else
  do i=1,maxn
    summ = 0.0d0
    do j=1,i
      summ = rotTP(j)*vibTP(i-j+1) + summ 
    enddo
    dosTP(i) = summ
  enddo
end if

! Rovibrational DOS above MECP.

dosX(1+binX:maxn) = dosTP(1:maxn-binX)

! Expression above gives the same result
! as an explicit convolution below.
! However, it is much more compact.

!do i = 1 + binX, maxn
!  summ = zero
!  do j = 1 + binX, i
!    summ = rotX(j)*vibX(i-j+1+binX) + summ
!  end do
!  dosX(i) = summ
!end do

end subroutine rovib_dos

!subroutine rev_rovib_dos(vibP,rotP,dosP)
subroutine rev_rovib_dos(vibP,rotP,dosP,dosP_full,dosX_full,dosX_TP_rev)
 real(wp), dimension(:), intent(in)  :: vibP, rotP
 real(wp), dimension(:), intent(out) :: dosP, dosP_full, dosX_full, dosX_TP_rev 
 real(wp)                            :: conP, summ
 real(wp), allocatable               :: vibX_full(:),rotX_full(:)
 integer                             :: i,j,k,q

dosP_full   = zero
dosX_full   = zero
dosX_TP_rev = zero
allocate (rotX_full (maxn+binGAP))
allocate (vibX_full (maxn+binGAP))
 
write(66,'(A)') '.............product rovibrational.'

do i=1,maxn+binGAP!-binXrev
! rotX_full(i+binXrev) = 4.0d0*sqrt(2.0d0*dble(i)/autocm*&
  rotX_full(i) = 4.0d0*sqrt(2.0d0*dble(i)/autocm*&
                 inertX(1)*inertX(2)*inertX(3))/autocm
enddo

do j = 1, size(freX)
  k=0
  do while (floor(freX(j)*dble(k)*autocm).le.maxn+binGAP)
     k = k + 1
     q = ceiling(freX(j)*dble(k)*autocm)
     if (q .eq. zero .or. q .gt. maxn+binGAP) cycle
       vibX_full(q) = vibX_full(q) + one
  end do
end do
! Add ZPE level.
  vibX_full(1) = vibX_full(1) + one

! Add external vibrational levels.

if (extern) then
  do i = 1, size(extX)
    q = ceiling(extX(i))
    if (q .eq. zero .or. q .gt. maxn+binGAP) cycle
    vibX_full(q) = vibX_full(q) + 1.0d0
  end do
end if
 
if (solution) then
  dosP_full = vibP
else
  do i=1, maxn+binGAP
    conP = 0.0d0
    do j=1,i
      conP=rotP(j)*vibP(i-j+1)+conP
    end do
    dosP_full(i)=conP ! Calculates rovibrational DOS for product.
  end do
end if
 
dosP(1:maxn)=dosP_full(1+binGAP: maxn+binGAP)

if (solution) then
  dosX_TP_rev = vibX_full
else
  do i = 1, maxn + binGAP
    summ = 0.0d0
    do j = 1, i
      summ = rotX_full(j)*vibX_full(i-j+1) + summ 
    end do
    dosX_TP_rev(i) = summ
  end do
end if

! Rovibrational DOS above MECP.

dosX_full(1+binXrev:maxn+binGAP) = dosX_TP_rev(1:maxn+binGAP-binXrev)

open (55,file='dos_product.out') 

write(55,*) '#Energy (cm-1)   ','Product(rot, vib, rovib) - states/cm-1'
   
do i=1, maxn+binGAP
  write(55,71) i,  rotP(i),  vibP(i),  dosP_full(i)
end do

71 format (I6,F30.15,F30.15,F30.15)   
   
close(55)

end subroutine rev_rovib_dos

!---------------------------------------------------------------------------
subroutine write_dos(vibR,vibX,vibTP,rotR,rotX,rotTP,   dosR,dosX,dosTP)
  
 integer                              :: i
 real(wp), dimension(:), intent(in)   :: vibR, vibX, vibTP, rotR, rotX, rotTP
 real(wp), dimension(:), intent(in)   :: dosR, dosX, dosTP 
   
open (11,file='dos_mecp.out') 
open (12,file='dos_reactant.out') 
open (13,file='dos_tp.out')

write(11,*) '#Energy (cm-1)   ','MECP(rot, vib, rovib) - states/cm-1'
write(12,*) '#Energy (cm-1)   ','Reactant(rot, vib, rovib) -  states/cm-1'
write(13,*) '#Energy (cm-1)   ','MECP TP(rot, vib, rovib) - states/cm-1'
   
do i=1, maxn
  write(11,71) i,  rotX(i),  vibX(i),  dosX(i)
  write(12,71) i,  rotR(i),  vibR(i),  dosR(i)
  write(13,71) i, rotTP(i), vibTP(i), dosTP(i)
end do

71 format (I6,F30.15,F30.15,F30.15)   
   
close(11)
close(12)
close(13)
  
end subroutine write_dos
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
end module convolution_module
