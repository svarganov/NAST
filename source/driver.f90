module driver_module
  use precision_module, only : wp
  use input_module
  use vib_module
  use rot_module
  use convolution_module
  use prob_module
  use rate_module
  
  implicit none
  public

contains
subroutine run
   
 implicit none
 real(wp), allocatable :: vibR(:), vibX(:), vibTP(:), vibP(:)
 real(wp), allocatable :: rotR(:), rotX(:), rotTP(:), rotP(:)
 real(wp), allocatable :: dosR(:), dosP(:), dosP_full(:), dosX(:), dosX1(:), dosTP(:)
 real(wp), allocatable :: dosX_full(:), dosX_TP_rev(:)
 real(wp), allocatable :: bbDia(:), bbDia_rev(:)
 real(wp)              :: aaDia 
 real(wp), allocatable :: probLZ(:), probWC(:), probZN(:), probEckart(:)
 real(wp), allocatable :: probLZS(:,:), probWCS(:,:), probLZM(:,:)
 real(wp), allocatable :: probLZ_rev(:), probWC_rev(:), probEckart_rev(:)
 real(wp), allocatable :: nosLZ(:), nosWC(:), nosZN(:), nosLZ_full(:)
 real(wp), allocatable :: nosWC_full(:), nosLZS(:,:), nosWCS(:,:), nosLZM(:,:)
 real(wp), allocatable :: nosTST(:), nosTST_full(:), nosTST_Eckart(:), nosTST_full_Eckart(:)
  
allocate (vibR  (maxn))
allocate (vibX  (maxn))
allocate (vibTP (maxn))
allocate (rotR  (maxn))
allocate (rotX  (maxn))
allocate (rotTP (maxn))
allocate (dosR  (maxn))
allocate (dosX  (maxn))
allocate (dosX1 (maxn))
allocate (dosTP (maxn))
allocate (bbDia (maxn))

! ~~~~ Program begins.   
  
! Densities of states.

write(66,'(/,2X,A)') '---------------------------------------------------------'
write(66,'(15X,A)') 'Start NAST calculation'
write(66,'(2X,A)') '---------------------------------------------------------'
write(66,'(/,3X,A,/)') '1. Calculating densities of states'
  
call vib_dos(vibR,vibX,vibTP)
call asym_top(rotR,rotX,rotTP)

call rovib_dos(vibR,vibX,vibTP,rotR,rotX,rotTP, &
               dosR,dosX,dosX1,dosTP)

if (printmore) then
  call write_dos(vibR,vibX,vibTP,rotR,rotX,rotTP, &
                dosR,dosX,dosTP)
end if

! Deallocation.

deallocate (vibR)
deallocate (vibX)
deallocate (vibTP)
deallocate (rotR)
deallocate (rotX)
deallocate (rotTP)


! Microcanonical probabilities.
if (.not. tst) then
  write(66,'(/,3X,A)') '2. Calculating microcanonical probabilities'
  allocate (probLZ (maxn))
  allocate (probWC (maxn))
  probLZ = zero
  probWC = zero
  bbDia  = zero
  call diab_param(aaDia,bbDia)
  call LZ(aaDia,bbDia,probLZ)
  call WC(aaDia,bbDia,probWC)
  if (zn) then
    allocate (probZN(maxn))
    call ZhN(probZN)
  else
    allocate (probZN (1))
    call ZNempty(probZN) 
  endif
  if (sp) then
    allocate (probLZS (maxn,icp))
    allocate (probLZM (maxn,4))
    allocate (probWCS (maxn,icp))
    call LZS(probLZS)
    call LZM(probLZM)
    call WCS(probWCS)
  else
    allocate (probLZS (1,1))
    allocate (probLZM (1,1))
    allocate (probWCS (1,1))
  end if
else
! Do formal allocation of non-used array(s).
    allocate (probLZ (1))
    allocate (probWC (1))
    allocate (probZN (1))
    allocate (probLZS (1,1))
    allocate (probLZM (1,1))
    allocate (probWCS (1,1))
end if

if (TST_tunn) then               
    allocate (probEckart (maxn))
    call EckartProb(probEckart) 
else
    allocate (probEckart (1))
end if


! Microcanonical and canonical rate constants.

if (.not. tst) then
  allocate (nosLZ  (maxn))
  allocate (nosWC  (maxn))
  allocate (nosZN  (maxn))
  nosLZ = zero
! Do formal allocation of non-used array(s).
  allocate (nosTST  (1))
  allocate (nosTST_Eckart  (1))
  write(66,'(/,3X,A)') '3. Calculating microcanonical NAST rate constant'
  call nast_en_rate(dosR,dosX,dosTP,dosP_full,probLZ,probWC,probZN,nosLZ,nosWC,nosZN)
  call nast_can_rate(dosR,dosP_full,nosLZ,nosWC,nosZN)
else
    if (TST_tunn) then
      allocate (nosTST (maxn))
      allocate (nosTST_Eckart (maxn))
    else
      allocate (nosTST (maxn))
      allocate (nosTST_Eckart (1))
    end if
! Do formal allocation of non-used array(s)
  allocate (nosLZ (1))
  allocate (nosWC (1))
  allocate (nosZN (1))
  call tst_en_rate(dosR,dosX,dosTP,probEckart, nosTST,nosTST_Eckart)
  call tst_can_rate(dosR,nosTST,nosTST_Eckart)
end if

! Microcanonical rate constants for transitions between MS components.

if (sp) then
  allocate (nosLZS (maxn,icp))
  nosLZS = zero
  allocate (nosLZM (maxn,4))
  nosLZM = zero
  allocate (nosWCS (maxn,icp))
  nosWCS = zero
  call nast_en_rate_LZ_split(dosR,dosX,dosX1,probLZS,nosLZS,probLZM,nosLZM)
  call nast_en_rate_WC_split(dosR,dosX,dosTP,probWCS,nosWCS)
  call nast_can_rate_split(dosR,nosLZS,nosWCS,nosLZM)
else
  allocate (nosLZS (1,1))
  allocate (nosLZM (1,1))
  allocate (nosWCS (1,1))
end if

! Deallocation.
deallocate (dosR)

deallocate (dosTP)
deallocate (dosX1)

! Velocity-averaged transition probabilities (VATP).
 
 if (.not. tst ) then
    call VATP()
 end if

! Reverse rate constants.

if (rev) then
 allocate (vibP          (maxn + binGAP))
 allocate (rotP          (maxn + binGAP))
 allocate (dosP          (maxn))
 allocate (dosP_full     (maxn + binGAP))
 allocate (dosX_full     (maxn + binGAP))
 allocate (dosX_TP_rev   (maxn + binGAP))
 call rev_vib_dos(vibP)
 call rev_asym_top(rotP)
 call rev_rovib_dos(vibP,rotP,dosP,dosP_full,dosX_full,dosX_TP_rev)

 if (.not. tst) then
   allocate (bbDia_rev   (maxn + binGAP))
   allocate (probLZ_rev  (maxn + binGAP))
   allocate (probWC_rev  (maxn + binGAP))
   allocate (nosLZ_full  (maxn + binGAP))
   allocate (nosWC_full  (maxn + binGAP))
   bbDia_rev  = zero
   probLZ_rev = zero
   probWC_rev = zero
   call LZ_rev (aaDia,bbDia_rev,probLZ_rev)
   call WC_rev (aaDia,bbDia_rev,probWC_rev)
   call nast_en_rate_rev(dosP_full,dosX_full,dosX_TP_rev,dosR,probLZ_rev,probWC_rev,nosLZ_full,nosWC_full)
   call nast_can_rate_rev(dosR,dosP_full,nosLZ_full,nosWC_full)
 else
   allocate (nosTST_full (maxn + binGAP))

   if (TST_tunn) then
     allocate (probEckart_rev (maxn+binGAP))
     call EckartProb_rev(probEckart_rev)
     allocate (nosTST_full_Eckart (maxn+binGAP))
   else
     allocate (probEckart_rev (1))
     allocate (nosTST_full_Eckart (1))
   end if

   call tst_en_rate_rev(dosP_full,dosX_full,dosX_TP_rev,probEckart_rev, nosTST_full,nosTST_full_Eckart)
   call tst_can_rate_rev(dosP_full,nosTST_full,nosTST_full_Eckart)
 end if
else
! Do formal allocation of non-used array(s)
  allocate (vibP       (1))
  allocate (rotP       (1))
  allocate (dosP_full  (1))
  allocate (dosX_full  (1))
  allocate (nosLZ_full (1))
  allocate (dosP       (1))
  allocate (bbDia_rev  (1))
  allocate (probLZ_rev (1))
endif

if (.not. tst) then
  call write_prob(aaDia,bbDia,probLZ,probLZ_rev,probLZS,probLZM,probWC,probWC_rev,probWCS,probZN,bbDia_rev)
end if

if (TST_tunn) then
  call write_prob_TST(probEckart,probEckart_rev)
end if

! Deallocation.

 deallocate (vibP)
 deallocate (rotP)
 deallocate (dosP)
 deallocate (dosX)
 deallocate (dosP_full)
 deallocate (probLZ)
 deallocate (probWC)
 deallocate (probZN)
 deallocate (nosLZ)
 deallocate (nosWC)
 deallocate (nosZN)
 deallocate (nosTST)
 deallocate (nosTST_Eckart)
 deallocate (bbDia)
 deallocate (probLZS)
 deallocate (probLZM)
 deallocate (probWCS)
 deallocate (nosLZS)
 deallocate (nosLZM)
 deallocate (nosWCS)

end subroutine run
end module driver_module
