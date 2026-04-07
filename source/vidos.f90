module vib_module
 use precision_module, only : wp
 use constants_module
 use input_module
   
 implicit none
 public

contains
subroutine vib_dos(vibR,vibX,vibTP)
 
 implicit none
 integer                                    :: h, i, j, k, q, ignore
 real(wp), dimension(:), intent(inout) :: vibR, vibX, vibTP
   
 vibR  = zero
 vibX  = zero
 vibTP = zero

 write (66,'(A)') '......vibrational.'

! Vibrational DOS for reactant.

do j = 1, size(freR)
  k = 0
  do while ((freR(j)*dble(k)*(1/Estep)*autocm).le.maxn)
    k = k + 1
    q = ceiling(freR(j)*dble(k)*(1/Estep)*autocm)
    if (q .eq. zero .or. q .gt. maxn) cycle
    vibR(q) = vibR(q) + one
  end do
end do

! Add external vibrational levels.

if (externR) then
  do i = 1, size(extR)
    q = ceiling(extR(i)*(1/Estep))
    if (q .eq. zero .or. q .gt. maxn) cycle
    vibR(q) = vibR(q) + one
  end do
end if

! Add ZPE level
vibR(1) = vibR(1) + one

vibR = vibR/Estep

! Vibrational DOS for turning points. 

do j = 1, size(freX)
  k=0
  do while ((freX(j)*dble(k)*(1/Estep)*autocm).le.maxn)
     k = k + 1
     q = ceiling(freX(j)*dble(k)*(1/Estep)*autocm)
     if (q .eq. zero .or. q .gt. maxn) cycle
       vibTP(q) = vibTP(q) + one
  end do
end do

! Add external vibrational levels.

if (externX) then
  do i = 1, size(extX)
    q = ceiling(extX(i)*(1/Estep))
    if (q .eq. zero .or. q .gt. maxn) cycle
    vibTP(q) = vibTP(q) + one
  end do
end if

! Add ZPE level.
vibTP(1) = vibTP(1) + one 

vibTP = vibTP/Estep

! Vibrational DOS above MECP. 

vibX(1+binX:maxn)=vibTP(1:maxn-binX)
 
end subroutine vib_dos

subroutine rev_vib_dos(vibP)
 
 implicit none
 integer       :: h, i, j, k, q, ignore
 real(wp), dimension(:), intent(inout) :: vibP

! Vibrational DOS for product.
   
vibP  = zero
 
write(66,'(/,a)') 'Calculating product densities of states.'
write(66,'(/,a)') '.........vibrational.'

! Vibrational DOS for product.
 
do j = 1, size(freP)
  k = 0
  do while ((freP(j)*dble(k)*(1/Estep)*autocm).le.maxn+binGAP)
   k = k + 1
   q = ceiling(freP(j)*dble(k)*(1/Estep)*autocm)
   if (q .eq. zero .or. q .gt. maxn+binGAP) cycle
     vibP(q) = vibP(q) + one
  end do
end do

! Add external vibrational levels.

if (externP) then
  do i = 1, size(extP)
    q = ceiling(extP(i)*(1/Estep))
    if (q .eq. zero .or. q .gt. maxn+binGAP) cycle
    vibP(q) = vibP(q) + one
 end do
end if

! Add ZPE level.
vibP(1) = vibP(1) + one

vibP = vibP/Estep

end subroutine rev_vib_dos

end module vib_module
