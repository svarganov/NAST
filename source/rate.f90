module rate_module
   use precision_module, only : wp
   use input_module
   use constants_module
   use error_module
   implicit none
   public
 
contains
!==================================================================================
!====================== MICROCANONICAL RATE CONSTANT k(E) =========================
!==================================================================================
subroutine nast_en_rate(dosR,dosX,dosTP,dosP_full,probLZ,probWC,probZN,nosLZ,nosWC,nosZN)
 
real(wp)                             :: WCsto,LZsto,ZNsto
real(wp), allocatable                :: nosWC_TP(:), nosZN_TP(:)
real(wp), allocatable                :: nosWC_X (:), nosZN_X (:)
real(wp), allocatable                :: rateLZ(:), rateWC(:), rateZN(:)
real(wp), dimension(:), intent(in)   :: dosR, dosX, dosTP, dosP_full
real(wp), dimension(:), intent(in)   :: probLZ, probWC, probZN
real(wp), dimension(:), intent(out)  :: nosLZ, nosWC, nosZN
integer                              :: i, j


allocate (nosWC_X  (maxn))
allocate (nosZN_X  (maxn))
allocate (nosWC_TP (maxn))
allocate (nosZN_TP (maxn))
allocate (rateLZ   (maxn)) 
allocate (rateWC   (maxn)) 
allocate (rateZN   (maxn))

nosWC_X  = zero
nosZN_X  = zero
nosWC_TP = zero
nosZN_TP = zero
rateLZ   = zero 
rateWC   = zero
rateZN   = zero
   
!**********************************************************************
!   OVER-THE-BARRIER CONTRIBUTION TO THE RATE
!      (over-the-barrier number of states) 
!**********************************************************************
if (maxn.gt.binX) then 
  do i=1+binX, maxn !start from MECP bin
    LZsto = zero
    WCsto = zero
    ZNsto = zero
    do j=1+binX, i ! NOS above crossing point at i-th bin
      LZsto = probLZ(j)*dosX(i-j+1+binX)*Estep + LZsto 
      WCsto = probWC(j)*dosX(i-j+1+binX)*Estep + WCsto
      ZNsto = probZN(j)*dosX(i-j+1+binX)*Estep + ZNsto
    enddo
    nosLZ(i)   = LZsto/autocm 
    nosWC_X(i) = WCsto/autocm 
    nosZN_X(i) = ZNsto/autocm 
 enddo
endif

!**********************************************************************
!   TUNNELING CONTRIBUTION TO THE RATE         ^_^
!      (tunneling number of states) 
!**********************************************************************
! Number of states at each turning point is approximated by single point with
! properties taken at MECP.

do i=1, maxn
  WCsto = zero
  ZNsto = zero
  if (i.le.binX) then
    do j=1, i
      WCsto = probWC(j)*dosTP(i-j+1)*Estep + WCsto 
      ZNsto = probZN(j)*dosTP(i-j+1)*Estep + ZNsto 
    enddo
  else
    do j=1, binX
      WCsto = probWC(j)*dosTP(i-j+1)*Estep + WCsto
      ZNsto = probZN(j)*dosTP(i-j+1)*Estep + ZNsto
    end do
  end if
   nosWC_TP(i) = WCsto/autocm
   nosZN_TP(i) = ZNsto/autocm
end do
 
!*****************************************
!        TOTAL NUMBER OF STATES 
!*****************************************
nosWC = nosWC_X + nosWC_TP
nosZN = nosZN_X + nosZN_TP


!*****************************************
!      Microcanonical Rate Constant 
!*****************************************
!binZPE = 0 if either zpe = 0 or zpe = 1 correction is used. 
do i=1, maxn-binZPE
  if (dosR(i).ne.zero) then
    rateLZ(i) = degenR*nosLZ(i+binZPE)/(dosR(i)*planck_hsec) 
    rateWC(i) = degenR*nosWC(i+binZPE)/(dosR(i)*planck_hsec)
    rateZN(i) = degenR*nosZN(i+binZPE)/(dosR(i)*planck_hsec)
  end if

end do


open(1451,file='energy_rates.out')

if (zn) then
 write(1451,'(25x,A,5x,A,5x,A,5x,A)') 'E, cm-1','LZ rate constant, s-1','WC rate constant, s-1','ZN rate constant, s-1'

 do i=1, maxn
   write(1451,'(F30.3,3es24.2)') i*Estep, rateLZ(i), rateWC(i), rateZN(i)
 end do

else
 write(1451,'(25x,A,5x,A,5x,A)') 'E, cm-1','LZ rate constant, s-1','WC rate constant, s-1'
 do i=1, maxn
   write(1451,'(F30.3,2es24.2)') i*Estep, rateLZ(i), rateWC(i)
 end do
endif
   
close(1451)
  
if (printmore) then
  open(1415,file='effective_number_of_states.out')
  write(1415,'(25x,A,15x,A,15x,A)') 'E, cm-1','nosLZ','nosWC'
  do i=1, maxn
    write(1415,'(F30.3,2es24.2)') i*Estep, nosLZ(i), nosWC(i)
  enddo

  close(1415)

end if
   
write(66,'(A,I1,A)') '     (The forward rate constant k(E) is multiplied by &
reaction path degeneracy equal to ',int(degenR),' ).'

deallocate (nosWC_X)
deallocate (nosZN_X)
deallocate (nosWC_TP)
deallocate (nosZN_TP)
deallocate (rateLZ) 
deallocate (rateWC) 
deallocate (rateZN)

end subroutine nast_en_rate
!=====================================================================================
subroutine nast_en_rate_LZ_split(dosR,dosX,dosX1,probLZS,nosLZS,probLZM,nosLZM)
   real(wp), allocatable                 :: LZSsto(:), LZMsto(:)
   real(wp), allocatable                 :: rateLZS(:,:), rateLZM(:,:)
   real(wp), dimension(:), intent(in)    :: dosR, dosX, dosX1
   real(wp), dimension(:,:), intent(in)  :: probLZS, probLZM
   real(wp), dimension(:,:), intent(out) :: nosLZS, nosLZM 
   integer                               :: i, j

   allocate (rateLZS (maxn,icp))
   allocate (LZSsto  (icp))
   allocate (rateLZM (maxn,4))
   allocate (LZMsto  (4))
  
   rateLZS = zero
   rateLZM = zero

!**********************************************************************
!   OVER-THE-BARRIER CONTRIBUTION TO THE RATE
!      (over-the-barrier number of states) 
!**********************************************************************
    if (maxn.gt.binX1) then 
      do i=1+binX1, maxn !start from 1st MECP bin
         LZMsto = zero
           do j=1+binX1, i ! NOS above crossing point at i-th bin
              LZMsto(1) = probLZM(j,1)*dosX1(i-j+1+binX1)*Estep + LZMsto(1)
           enddo
         nosLZM(i,1)   = LZMsto(1)/autocm
      enddo
    endif  
   
    if (maxn.gt.binX) then 
      do i=1+binX, maxn !start from MECP bin
         LZSsto = zero
         LZMsto = zero
           do j=1+binX, i ! NOS above crossing point at i-th bin
              LZSsto(:) = probLZS(j,:)*dosX(i-j+1+binX)*Estep + LZSsto(:)
              LZMsto(2) = probLZM(j,2)*dosX(i-j+1+binX)*Estep + LZMsto(2)
           enddo
         nosLZS(i,:)   = LZSsto(:)/autocm
         nosLZM(i,2)   = LZMsto(2)/autocm
      enddo
    endif

    if (maxn.gt.binX3) then 
      do i=1+binX3, maxn !start from 3rd MECP bin
         LZMsto = zero
           do j=1+binX3, i ! NOS above crossing point at i-th bin
              LZMsto(3) = probLZM(j,3)*dosX(i-j+1+binX3)*Estep + LZMsto(3)
           enddo
         nosLZM(i,3)   = LZMsto(3)/autocm
      enddo
    endif

!*****************************************
!      Microcanonical Rate Constant 
!*****************************************
   !binZPE = 0 if no ZPE correction or ZPE1 correction is used 
  do i=1, maxn-binZPE
   if (dosR(i).ne.zero) then 
      rateLZS(i,:) = degenR*nosLZS(i+binZPE,:)/(dosR(i)*planck_hsec)
      rateLZM(i,1) = degenR*nosLZM(i+binZPE,1)/(dosR(i)*planck_hsec)
      rateLZM(i,2) = degenR*nosLZM(i+binZPE,2)/(dosR(i)*planck_hsec)
      rateLZM(i,3) = degenR*nosLZM(i+binZPE,3)/(dosR(i)*planck_hsec)
      rateLZM(i,4) = rateLZM(i,1) + rateLZM(i,2) + rateLZM(i,3)
   end if
  end do


  open(556,file='split_LZ_energy_rate.out')
  write(556,'(25x,a)', advance="no") "Energy(cm-1)"
  write(556, '(4X,I2,A,I2,A,I2,A,I2)', advance="no") mults(1), "_", 1, &
                                                   "  ->",mults(2), "_", 1
  do i=2, icp
   write(556, '(11X,I2,A,I2,A,I2,A,I2)', advance="no") mults(1), "_", 1+floor((real(i)-1)/3), &
                                            " ->", mults(2), "_", mod(i-1,3)+1+floor((real(i)-1)/3)
  end do
  
  write(556,'(/)')
 
  do i=1, maxn
      write(556,'(F30.3,2es24.2)') i*Estep, rateLZS(i,:)
  enddo
  close(556)
  open(557,file='Bz_LZ_energy_rate.out')
      write(557,'(25x,a,4x,a,4x,a,4x,a,4x,a)') "E, cm-1         ", & 
      "k1(r-,Bz)              ", &
      "k2(r0,Bz)              ", &
      "k3(r+,Bz)              ", &
      "k_overal(Bz)              "
  do i=1, maxn
    write(557,'(F30.2,4es24.4)') i*Estep, rateLZM(i,:)
  enddo
  close(557)
  
deallocate (rateLZS)
deallocate (LZSsto)
deallocate (rateLZM)
deallocate (LZMsto)

end subroutine nast_en_rate_LZ_split
!============================================================================
subroutine nast_en_rate_WC_split(dosR,dosX,dosTP,probWCS,nosWCS)
   real(wp), allocatable                 :: WCSsto(:)
   real(wp), allocatable, dimension(:,:) :: nosWCS_TP
   real(wp), allocatable, dimension(:,:) :: nosWCS_X
   real(wp), allocatable                 :: rateWCS(:,:)
   real(wp), dimension(:), intent(in)    :: dosR, dosX, dosTP
   real(wp), dimension(:,:), intent(in)  :: probWCS
   real(wp), dimension(:,:), intent(out) :: nosWCS
   integer                               :: i, j

   allocate (nosWCS_X  (maxn,icp))
   allocate (nosWCS_TP (maxn,icp))
   allocate (rateWCS   (maxn,icp))
   allocate (WCSsto    (icp))

   rateWCS   = zero
   nosWCS_X  = zero
   nosWCS_TP = zero

!**********************************************************************
!   OVER-THE-BARRIER CONTRIBUTION TO THE RATE
!      (over-the-barrier number of states) 
!**********************************************************************
   if (maxn.gt.binX) then 
      do i=1+binX, maxn !start from MECP bin
         WCSsto = zero
           do j=1+binX, i ! NOS above crossing point at i-th bin
              WCSsto(:) = probWCS(j,:)*dosX(i-j+1+binX)*Estep + WCSsto(:)
           enddo
         nosWCS_X(i,:)   = WCSsto(:)/autocm
      enddo
  endif

!**********************************************************************
!   TUNNELING CONTRIBUTION TO THE RATE         ^_^
!      (tunneling number of states)
!**********************************************************************
! Number of states at each turning point is approximated by single point with
! properties taken at MECP.

do i=1, maxn
  WCSsto = zero
  if (i.le.binX) then
    do j=1, i
      WCSsto(:) = probWCS(j,:)*dosTP(i-j+1)*Estep + WCSsto(:) 
    enddo
  else
    do j=1, binX
      WCSsto(:) = probWCS(j,:)*dosTP(i-j+1)*Estep + WCSsto(:)
    end do
  end if
   nosWCS_TP(i,:) = WCSsto(:)/autocm
end do

!*****************************************
!        TOTAL NUMBER OF STATES 
!*****************************************
nosWCS= nosWCS_X + nosWCS_TP

!*****************************************
!    Microcanonical Split Rate Constant 
!*****************************************
   !binZPE = 0 if no ZPE correction or ZPE1 correction is used 
  do i=1, maxn-binZPE
   if (dosR(i).ne.zero) then
      rateWCS(i,:) = degenR*nosWCS(i+binZPE,:)/(dosR(i)*planck_hsec)
   end if
  end do

  open(557,file='split_WC_energy_rate.out')
  write(557,'(25x,a)', advance="no") "Energy(cm-1)"
  write(557, '(4X,I2,A,I2,A,I2,A,I2)', advance="no") mults(1), "_", 1, &
                                                   "  ->",mults(2), "_", 1
  do i=2, icp
   write(557, '(11X,I2,A,I2,A,I2,A,I2)', advance="no") mults(1), "_", 1+floor((real(i)-1)/3), &
                                            " ->", mults(2), "_", mod(i-1,3)+1+floor((real(i)-1)/3)
  end do
  
  write(557,'(/)')
 
  do i=1, maxn
      write(557,'(F30.3,2es24.2)') i*Estep, rateWCS(i,:)
  enddo
close(557)

deallocate (rateWCS)
deallocate (WCSsto)

end subroutine nast_en_rate_WC_split
!============================================================================
subroutine tst_en_rate(dosR,dosX,dosTP,probEckart,nosTST,nosTST_Eckart)
implicit none

   integer                            :: i, j
   real(wp)                           :: TSTsto
   real(wp), allocatable              :: nosTST_X(:), nosTST_TP(:)
   real(wp), allocatable              :: rateTST(:), rateTST_Eckart(:)
   real(wp), dimension(:), intent(in)   :: dosR, dosX, dosTP
   real(wp), allocatable, intent(in)  :: probEckart(:)
   real(wp), dimension(:), intent(out) :: nosTST, nosTST_Eckart
 
!**********************************************************************
!   CONVENTIONAL TST NUMBER-OF-STATES AND MICROCANONICAL RATE CONSTANT
!**********************************************************************
   allocate (rateTST   (maxn))
   
   write(66,'(/,A)') 'Calculating microcanonical TST rate constant.'

   nosTST     = zero
   rateTST    = zero

   do i = 1+binX, maxn ! Start from the MECP bin
      nosTST(i) = sum(dosX(1+binX:i))*Estep
   end do
   nosTST = nosTST/autocm

   where (dosR .ne. zero) 
       rateTST(1:maxn-binZPE) = degenR*nosTST(1+binZPE:maxn)&
       *autoev/(dosR(1:maxn-binZPE)*planck)
   end where

if (TST_tunn) then
!***************************************************************************
!   QUANTUM-MECHANICAL TST NUMBER-OF-STATES AND MICROCANONICAL RATE CONSTANT
!***************************************************************************
   allocate (nosTST_X  (maxn))
   allocate (nosTST_TP (maxn))
   allocate (rateTST_Eckart   (maxn))
   nosTST_X   = zero
   nosTST_TP  = zero
   nosTST_Eckart = zero
   rateTST_Eckart    = zero
!**********************************************************************
!   OVER-THE-BARRIER CONTRIBUTION TO THE RATE
!      (over-the-barrier number of states) 
!**********************************************************************
 if (maxn.gt.binX) then 
   do i=1+binX, maxn !start from MECP bin
     TSTsto = zero
     do j=1+binX, i ! NOS above crossing point at i-th bin
       TSTsto = probEckart(j)*dosX(i-j+1+binX)*Estep + TSTsto 
     enddo
     nosTST_X(i)   = TSTsto/autocm 
  enddo
 endif

!**********************************************************************
!   TUNNELING CONTRIBUTION TO THE RATE         ^_^
!      (tunneling number of states) 
!**********************************************************************
! Number of states at each turning point is approximated by single point with
! properties taken at MECP.

 do i=1, maxn
   TSTsto = zero
   if (i.le.binX) then
     do j=1, i
       TSTsto = probEckart(j)*dosTP(i-j+1)*Estep + TSTsto 
     enddo
   else
     do j=1, binX
       TSTsto = probEckart(j)*dosTP(i-j+1)*Estep + TSTsto
     end do
   end if
    nosTST_TP(i) = TSTsto/autocm
 end do
 
!*****************************************
!        TOTAL NUMBER OF STATES 
!*****************************************
 nosTST_Eckart = nosTST_X + nosTST_TP

!************************************************************************
!      Microcanonical TST Rate Constant with Eckart tunneling probability
!************************************************************************
!binZPE = 0 if either zpe = 0 or zpe = 1 correction is used. 

 do i=1, maxn-binZPE
   if (dosR(i).ne.zero) then
     rateTST_Eckart(i) = degenR*nosTST_Eckart(i+binZPE)/(dosR(i)*planck_hsec) 
   end if

 end do

endif

if (printmore) then

  open(321,file='effective_number_of_states_TST.out')

  if (TST_tunn) then
      write(321,'(25x,A,7x,A,7x,A)') "E, cm-1", "nos TST", "nos TST Eckart"

      do i=1, maxn
         write(321,200) i*Estep, nosTST(i), nosTST_Eckart(i)
      enddo

  else
      write(321,'(25x,A,7x,A)') "E, cm-1", "nos TST"

      do i=1, maxn
         write(321,201) i*Estep, nosTST(i)
      enddo

  endif

  close(321)

endif

open(322,file='energy_tst_rate.out')

if (TST_tunn) then
    write(322,'(25x,A,7x,A,7x,A)') "E, cm-1", "TST rate constant, s-1", "TST rate constant Eckart, s-1"

    do i=1, maxn
       write(322,200) i*Estep, rateTST(i), rateTST_Eckart(i)
    enddo

else
    write(322,'(25x,A,7x,A)') "E, cm-1", "TST rate constant, s-1"

    do i=1, maxn
       write(322,201) i*Estep, rateTST(i)
    enddo

endif

close(322)

  200 format('',F30.3,2es24.2)
  201 format('',F30.3,es24.2)
 
  write(66,'(A)') '................................... Done.'
  write(66,'(/,A,I1,A)') '     (The forward rate constant k(E) is multiplied by &
 reaction path degeneracy equal to ',int(degenR),' )'

if (.not. TST_tunn) then
 deallocate (rateTST)
else
 deallocate (nosTST_X)
 deallocate (nosTST_TP)
 deallocate (rateTST_Eckart)
endif

end subroutine tst_en_rate
!=============================================================================
!============== CANONICAL (THERMAL) RATE CONSTANT k(T) =======================
!=============================================================================
subroutine nast_can_rate(dosR,dosP_full,nosLZ,nosWC,nosZN)
   implicit none
   character(len=33)                       :: printLZ, printWC, printZN
   integer                                 :: i,j,T_n_points
   real(wp)                                :: qR, kbt, Temp
   real(wp)                                :: trateLZ, trateWC, trateZN
   real(wp), dimension(:), intent(in)      :: dosR, dosP_full, nosLZ, nosWC, nosZN
 
   Temp = zero
   kbt = zero

   if (zn) then
     write(66, 201)
   else
     write(66, 202)
   endif

     do j = 0, Tpoints
       Temp = T1 + dble(j)*Tstep
   
     trateLZ = zero
     trateWC = zero
     trateZN = zero
 
       qR = zero
       kbt = kb_hartK*Temp
 
       do i=1, maxn
         qR=qR + dosR(i)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
       enddo
 
       do i=1, maxn
         trateLZ  = trateLZ + nosLZ(i)*exp(-dble(i)*Estep/(autocm*kbt))*Estep 
         trateWC  = trateWC + nosWC(i)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
         trateZN  = trateZN + nosZN(i)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
       end do
       trateLZ = trateLZ/(planck_hsec*qR) 
       trateWC = trateWC/(planck_hsec*qR) 
       trateZN = trateZN/(planck_hsec*qR)
       
       if (zn) then
         write(printLZ,'(es24.2)') trateLZ
         write(printWC,'(es24.2)') trateWC
         write(printZN,'(es24.2)') trateZN
         write(66, 204) Temp, trim(adjustl(adjustr(printLZ))), &
                             trim(adjustl(adjustr(printWC))), &
                             trim(adjustl(adjustr(printZN)))
       else
         write(printLZ,'(es24.2)') trateLZ
         write(printWC,'(es24.2)') trateWC
         write(66, 205) Temp, trim(adjustl(adjustr(printLZ))), &
                             trim(adjustl(adjustr(printWC)))
       end if
     enddo
  
  201 format(/,/,13x,'Canonical rate constant k(T), s-1.',/,/,&
            5x,'T(K)',3x,'Landau-Zener',3x,'Weak Coupling',3x,'Zhu Nakamura',/)
  202 format(/,/,13x,'Canonical rate constant k(T), s-1',/,/,&
            5x,'T(K)',3X,'Landau-Zener',3x,'Weak Coupling',/)
  203 format(4x,f7.1,4x,a)
  204 format(4x,f7.1,4x,a,8x,a,8x,a)
  205 format(4x,f7.1,4x,a,8x,a)

end subroutine nast_can_rate
!================================================================================
!=========== CANONICAL (THERMAL) Split RATE CONSTANT k(T) =======================
!================================================================================
subroutine nast_can_rate_split(dosR,nosLZS,nosWCS,nosLZM)
   implicit none
   integer                                 :: i,j,T_n_points
   real(wp)                                :: qR, kbt, Temp
   real(wp), allocatable                   :: trateLZS(:), trateWCS(:), trateLZM(:)
   real(wp), dimension(:), intent(in)      :: dosR
   real(wp), dimension(:,:), intent(in)    :: nosLZS, nosWCS, nosLZM
 
   allocate (trateLZS  (icp))
   allocate (trateWCS  (icp))
   allocate (trateLZM  (3))

   Temp = zero
   kbt = zero
  
  open(558,file='split_LZ_tempt_rate.out')
  open(559,file='split_WC_tempt_rate.out')
  open(560,file='Bz_LZ_tempt_rate.out')
  write(558,'(a)') "Canonical rate constants between two Ms spin states using LZ formula, s-1"
  write(559,'(a)') "Canonical rate constants between two Ms spin states using WC formula, s-1"
  write(560,'(a)') "Canonical MFE rate constants using LZ formula, s-1"
  write(558, '(6x,a)', advance="no") "T (K)"
  write(559, '(6x,a)', advance="no") "T (K)"
  write(560, '(6x,a)', advance="no") "T (K)"
  write(558, '(4x,I2,a,I2,a,I2,a,I2)', advance="no") mults(1), " _", 1, &
                                      " ->", mults(2), " _", 1
  write(559, '(4x,I2,a,I2,a,I2,a,I2)', advance="no") mults(1), " _", 1, &
                                      " ->", mults(2), " _", 1
  write(560, '(15x,a,15x,a,15x,a)', advance="no") "k1", "k2", "k3"
  do i=2, icp
       write(558, '(11x,I2,a,I2,a,I2,a,I2)', advance="no") mults(1), " _", 1+floor((real(i)-1)/3), &
                                       " ->", mults(2), " _", mod(i-1,3)+1+floor((real(i)-1)/3)
       
       write(559, '(11x,I2,a,I2,a,I2,a,I2)', advance="no") mults(1), " _", 1+floor((real(i)-1)/3), &  
                                       " ->", mults(2), " _", mod(i-1,3)+1+floor((real(i)-1)/3)
  enddo
  write (558, '(/)')
  write (559, '(/)')
  write (560, '(/)')
  
  do j = 0, Tpoints
      Temp = T1 + dble(j)*Tstep
   
      trateLZS = zero
      trateWCS = zero
      trateLZM = zero
      qR      = zero 
      kbt = kb_hartK*Temp
 
      do i=1, maxn
        qR=qR + dosR(i)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
        trateLZS(:)  = trateLZS(:) + nosLZS(i,:)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
        trateWCS(:)  = trateWCS(:) + nosWCS(i,:)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
        trateLZM(1)  = trateLZM(1) + nosLZM(i,1)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
        trateLZM(2)  = trateLZM(2) + nosLZM(i,2)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
        trateLZM(3)  = trateLZM(3) + nosLZM(i,3)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
      end do
 
      trateLZS(:) = trateLZS(:)/(planck_hsec*qR) 
      trateWCS(:) = trateWCS(:)/(planck_hsec*qR) 
      trateLZM(:) = trateLZM(:)/(planck_hsec*qR)
     
      write(558,*) Temp, trateLZS(:)
      write(559,*) Temp, trateWCS(:)
      write(560,333) Temp, trateLZM(:)

  enddo

  close(558)
  close(559)
  close(560)

  deallocate (trateLZS)
  deallocate (trateWCS)
  deallocate (trateLZM)

333 format(4x,f7.1,3es24.2)
end subroutine nast_can_rate_split
!================================================================
subroutine tst_can_rate(dosR,nosTST,nosTST_Eckart)
   implicit none
   character(len=33)                          :: printTST, printTST_Wigner, printTST_Eckart
   integer                                    :: i,j,T_n_points
   real(wp)                                   :: qR, kbt, Temp, kappa_Wigner
   real(wp)                                   :: trateTST, trateTST_Wigner, trateTST_Eckart
   real(wp), allocatable, intent(in)          :: dosR(:)
   real(wp), dimension(:), intent(in)         :: nosTST, nosTST_Eckart
 
   Temp = zero
   kbt = zero

if (.not. TST_tunn) then
   write(66, 202)
else
   write(66, 201)
end if
   
do j = 0, Tpoints

    Temp = T1 + dble(j)*Tstep
    trateTST = zero

    qR = zero
    kbt = kb_hartK*Temp
 
    do i=1, maxn
      qR=qR + dosR(i)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
    enddo
 
    do i=1, maxn
      trateTST  = trateTST + nosTST(i)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
    end do
 
    trateTST = trateTST/(planck_hsec*qR)
      
    if (TST_tunn) then
       kappa_Wigner = zero
       trateTST_Wigner = zero
       trateTST_Eckart = zero
        
       do i=1, maxn
         trateTST_Eckart  = trateTST_Eckart + nosTST_Eckart(i)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
       end do

       trateTST_Eckart = trateTST_Eckart/(planck_hsec*qR)
   
       kappa_Wigner=one+(((planck_hsec*TST_freq*cmtoHz)/kbt)**2)/24
       trateTST_Wigner = trateTST*kappa_Wigner
       write(printTST,'(es24.2)') trateTST
       write(printTST_Wigner,'(es24.2)') trateTST_Wigner
       write(printTST_Eckart,'(es24.2)') trateTST_Eckart
       write(66, 203) Temp, trim(adjustl(adjustr(printTST))), &
                           trim(adjustl(adjustr(printTST_Wigner))), &
                           trim(adjustl(adjustr(printTST_Eckart)))
    else      
       write(printTST,'(es24.2)') trateTST
       write(66, 204) Temp, trim(adjustl(adjustr(printTST)))
    end if

enddo
  
  201 format(/,/,13x,'Canonical TST rate constant, k(T)',/,/,&
            5x,'T(K)',5x,'k, s-1',4x,'k_Wigner, s-1',4x,'k_Eckart, s-1',/)
  202 format(/,/,13x,'Canonical TST rate constant, k(T)',/,/,&
            5x,'T(K)',5x,'k, s-1',/)
  203 format(2x,f7.1,5x,a,4x,a,4x,a,4x,a)
  204 format(2x,f7.1,5x,a,4x,a)

  205 format('',F30.3,es24.7)

end subroutine tst_can_rate
!==================================================================================
!================ MICROCANONICAL REVERSE RATE CONSTANT k(E) =======================
!==================================================================================
subroutine nast_en_rate_rev(dosP_full,dosX_full,dosX_TP_rev,dosR,probLZ_rev,probWC_rev,nosLZ_full,nosWC_full)
   
 real(wp), dimension(:), intent(in)  :: dosP_full,dosX_full,dosX_TP_rev,dosR,probLZ_rev,probWC_rev
 real(wp), dimension(:), intent(out) :: nosLZ_full,nosWC_full
 real(wp), allocatable               :: REVrateLZ(:),REVrateWC(:),nosWC_X(:),nosWC_TP(:)
 integer                             :: i,j
 real(wp)                            :: LZsto,WCsto,p

 
allocate (REVrateLZ     (maxn+binGAP))
allocate (REVrateWC     (maxn+binGAP))
allocate (nosWC_X       (maxn+binGAP)) 
allocate (nosWC_TP      (maxn+binGAP)) 

nosLZ_full = zero
nosWC_full = zero
nosWC_X    = zero
nosWC_TP   = zero
REVrateLZ  = zero
REVrateWC  = zero

!************************************
!   Over the barrier contribution   *
!   to the resverse rate constant.  *
!   Calculating over-the-barrier    *
!   number of states.               *
!************************************

if (maxn+binGAP.gt.binXrev) then 
  do i=1+binXrev, maxn+binGAP !start from MECP bin.
    LZsto = zero
    WCsto = zero
    do j=1+binXrev, i ! NOS above crossing point at i-th bin
     LZsto = probLZ_rev(j)*dosX_full(i-j+1+binXrev)*Estep + LZsto
     WCsto = probWC_rev(j)*dosX_full(i-j+1+binXrev)*Estep + WCsto
    end do
    nosLZ_full(i)   = LZsto/autocm
    nosWC_X(i)      = WCsto/autocm
  end do
end if

 write(66,'(/,a)') 'Calculating reverse rate constant.'
 write(66,'(/,a,I2)') 'The reverse rate constant k(E) is multiplied by &
           reaction path degeneracy equals to ',int(degenP)

!************************************
!   Tunneling contribution          * 
!   to the resverse rate constant.  *
!   Calculating tunneling           *
!   number of states.               *
!************************************
! Number of states at each turning point is approximated by single point
! with properties taken at MECP.

do i = binGAP+1, maxn + binGAP
  WCsto = zero
  if (i .le. binXrev) then
    do j = binGAP+1, i
     WCsto = probWC_rev(j)*dosX_TP_rev(i-j+binGAP+1)*Estep + WCsto
    end do
  else
    do j  = binGAP+1, binXrev
     WCsto = probWC_rev(j)*dosX_TP_rev(i-j+binGAP+1)*Estep + WCsto
    end do
  end if
  nosWC_TP(i) = WCsto/autocm
end do


!*****************************
!   Total number of states
!*****************************

nosWC_full = nosWC_X + nosWC_TP

!*****************************************
!   Microcanonical Reverse Rate Constant 
!*****************************************

do i = 1, maxn + binGAP
  if (dosP_full(i).ne.0.0d0) then 
    REVrateLZ(i) = degenP*nosLZ_full(i)/(dosP_full(i)*planck_hsec)
    REVrateWC(i) = degenP*nosWC_full(i)/(dosP_full(i)*planck_hsec)
  end if
end do


open(144,file='energy_rate_rev.out')

write(144,'(25x,a,4x,a,4x,a)') 'E, cm-1','Reverse LZ rate constant, s-1','Reverse WC rate constant, s-1'
do i=1, maxn+binGAP
  write(144,201) i*Estep, REVrateLZ(i), REVrateWC(i)
end do

if (printmore) then
  open(1444,file='effective_number_of_states_rev.out')
  write(1444,'(25x,a,4x,a,4x,a)') 'E, cm-1','LZ','WC'
   
  do i = 1, maxn + binGAP
    write(1444,'(F30.3,es24.7,4x,es24.7)') i*Estep,nosLZ_full(i),nosWC_full(i)
  end do

  close(1444)
end if
 
deallocate (REVrateLZ)
deallocate (REVrateWC)

close(144)
200 format(F30.3,3es24.7)
201 format(F30.3,2es24.7)

end subroutine nast_en_rate_rev
!==============================================================================
!============= MACROCANONICAL (THERMAL) REVERSE RATE CONSTNAT k(T) ============
!==============================================================================
subroutine nast_can_rate_rev(dosR,dosP_full,nosLZ_full,nosWC_full)
   
implicit none
real(wp), dimension(:), intent(in)  :: dosP_full,nosLZ_full,nosWC_full,dosR
character(len=33)                   :: printLZ, printWC
integer                             :: i,j
real(wp)                            :: qP_full,kbt, temp 
real(wp)                            :: rev_trateLZ, rev_trateWC
 
temp = zero
kbt = zero

write(66, 202)

do j = 0, Tpoints

   temp = T1 + dble(j)*Tstep

   rev_trateLZ = zero
   rev_trateWC = zero
  
   qP_full = zero
   kbt = kb_hartK*Temp
 
   do i=1, maxn+binGAP
     qP_full = qP_full + dosP_full(i)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
   end do

   do i=1, maxn
     rev_trateLZ = rev_trateLZ + nosLZ_full(i)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
     rev_trateWC = rev_trateWC + nosWC_full(i)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
   end do

   rev_trateLZ = rev_trateLZ/(planck_hsec*qP_full)
   rev_trateWC = rev_trateWC/(planck_hsec*qP_full)
 
   write(printLZ,'(es24.2)') rev_trateLZ
   write(printWC,'(es24.2)') rev_trateWC
   write(66, 204) Temp, trim(adjustl(adjustr(printLZ))), &
                        trim(adjustl(adjustr(printWC)))

end do

202 format(/,/,13X,'Canonical reverse rate constant, k(T)',/,/,&
           5x,'T(K)',5x,'Landau-Zener',5x,'Weak Coupling',/)
203 format(4x,f7.1,4x,a)
204 format(4x,f7.1,4x,a,8x,a)

end subroutine nast_can_rate_rev
!============================================================================
! TST microcanonical reverse rate constant.
!============================================================================
subroutine tst_en_rate_rev(dosP_full,dosX_full,dosX_TP_rev,probEckart_rev,nosTST_full,nosTST_full_Eckart)
implicit none

integer                             :: i, j
real(wp)                            :: TSTsto_rev
real(wp), allocatable               :: nosTST_full_X(:), nosTST_TP_rev(:)
real(wp), dimension(:), intent(in)  :: dosP_full, dosX_full, dosX_TP_rev
real(wp), allocatable, intent(in)   :: probEckart_rev(:)
real(wp), dimension(:), intent(out) :: nosTST_full, nosTST_full_Eckart
real(wp), allocatable               :: rateTST_rev(:), rateTST_rev_Eckart(:)

!*****************************************************************************
!   CONVENTIONAL TST REVERSE NUMBER-OF-STATES AND MICROCANONICAL RATE CONSTANT
!*****************************************************************************
allocate (rateTST_rev  (maxn+binGAP))

write(66,'(/,A)') 'Calculating microcanonical TST reverse rate constant.'

nosTST_full = zero
rateTST_rev = zero

do i = 1+binXrev, maxn+binGAP
  nosTST_full(i) = sum(dosX_full(1+binXrev:i))*Estep
end do

nosTST_full = nosTST_full/autocm

where (dosP_full .ne. zero)
     rateTST_rev(1:maxn+binGAP-binZPE) = degenP*nosTST_full(1+binZPE:maxn+binGAP)&
     *autoev/(dosP_full(1:maxn+binGAP-binZPE)*planck)
end where

if (TST_tunn) then
!***********************************************************************************
!   QUANTUM-MECHANICAL TST REVERSE NUMBER-OF-STATES AND MICROCANONICAL RATE CONSTANT
!***********************************************************************************
   allocate (nosTST_full_X  (maxn+binGAP))
   allocate (nosTST_TP_rev (maxn+binGAP))
   allocate (rateTST_rev_Eckart (maxn+binGAP))
   nosTST_full_X       = zero
   nosTST_TP_rev       = zero
   nosTST_full_Eckart  = zero
   rateTST_rev_Eckart  = zero
!**********************************************************************
!   OVER-THE-BARRIER CONTRIBUTION TO THE RATE
!      (over-the-barrier number of states) 
!**********************************************************************
if (maxn+binGAP.gt.binXrev) then 
  do i=1+binXrev, maxn+binGAP !start from MECP bin.
    TSTsto_rev = zero
    do j=1+binXrev, i ! NOS above crossing point at i-th bin
     TSTsto_rev = probEckart_rev(j)*dosX_full(i-j+1+binXrev)*Estep + TSTsto_rev
    end do
    nosTST_full_X(i)   = TSTsto_rev/autocm
  end do
end if

 write(66,'(/,a)') 'Calculating reverse rate constant.'
 write(66,'(/,a,I2)') 'The reverse rate constant k(E) is multiplied by &
           reaction path degeneracy equals to ',int(degenP)

!**********************************************************************
!   TUNNELING CONTRIBUTION TO THE RATE         ^_^
!      (tunneling number of states) 
!**********************************************************************
! Number of states at each turning point is approximated by single point with
! properties taken at MECP.

do i = binGAP+1, maxn + binGAP
  TSTsto_rev = zero
  if (i .le. binXrev) then
    do j = binGAP+1, i
     TSTsto_rev = probEckart_rev(j)*dosX_TP_rev(i-j+binGAP+1)*Estep + TSTsto_rev
    end do
  else
    do j  = binGAP+1, binXrev
     TSTsto_rev = probEckart_rev(j)*dosX_TP_rev(i-j+binGAP+1)*Estep + TSTsto_rev
    end do
  end if
  nosTST_TP_rev(i) = TSTsto_rev/autocm
end do

!*****************************
!   Total number of states
!*****************************

nosTST_full_Eckart = nosTST_full_X + nosTST_TP_rev

!************************************************************************
!  Microcanonical TST Reverse Rate Constant with Eckart tunneling probability
!************************************************************************

do i = 1, maxn + binGAP
  if (dosP_full(i).ne.0.0d0) then 
    rateTST_rev_Eckart(i) = degenP*nosTST_full_Eckart(i)/(dosP_full(i)*planck_hsec)
  end if
end do

endif

if (printmore) then

  open(324,file='effective_number_of_states_TST_rev.out')

  if (TST_tunn) then
      write(324,'(25x,A,7x,A,7x,A)') "E, cm-1", "nos TST reverse", "nos TST Eckart reverse"

      do i=1, maxn
         write(324,200) i*Estep, nosTST_full(i), nosTST_full_Eckart(i)
      enddo

  else
      write(324,'(25x,A,7x,A)') "E, cm-1", "nos TST reverse"

      do i=1, maxn
         write(324,201) i*Estep, nosTST_full(i)
      enddo

  endif

  close(324)

end if

open(325,file='energy_tst_rate_rev.out')

if (TST_tunn) then
    write(325,'(25x,A,7x,A,7x,A)') "E, cm-1", "Reverse TST rate constant, s-1", "Reverse TST Eckart rate constant, s-1"

    do i=1, maxn
       write(325,200) i*Estep, rateTST_rev(i), rateTST_rev_Eckart(i)
    enddo

else
    write(325,'(25x,A,7x,A)') "E, cm-1", "Reverse TST rate constant, s-1"

    do i=1, maxn
       write(325,201) i*Estep, rateTST_rev(i)
    enddo

endif

close(325)

200 format('',F30.3,2es24.7)
201 format('',F30.3,es24.7)

write(66,'(A)') '................................... Done.'
write(66,'(/,A,I1,A)') '     (The reverse TST rate constant k(E) is multiplied by &
reaction path degeneracy equal to ',int(degenP),' ).'

end subroutine tst_en_rate_rev
!=================================================
! Traditional TST canonical reverse rate constant.
!=================================================
subroutine tst_can_rate_rev(dosP_full,nosTST_full,nosTST_full_Eckart)
implicit none

real(wp), dimension(:), intent(in)  :: dosP_full, nosTST_full, nosTST_full_Eckart
character(len=33)                   :: printTST_rev, printTST_rev_Wigner, printTST_rev_Eckart
real(wp)                            :: qP_full, temp, kbt, kappa_Wigner
real(wp)                            :: trateTST_rev, trateTST_rev_Wigner, trateTST_rev_Eckart
integer                             :: i,j

temp = zero
kbt = zero

if (TST_tunn) then
  write (66, 201)
else
  write(66, 202)
end if

do j = 0, Tpoints

  temp = T1 + dble(j)*Tstep
  trateTST_rev = zero
  qP_full = zero
  kbt = kb_hartK*Temp
 
  do i=1, maxn
    qP_full = qP_full + dosP_full(i)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
    trateTST_rev  = trateTST_rev + nosTST_full(i)*exp(-dble(i)*Estep/(autocm*kbt))*Estep 
  end do
 
  trateTST_rev = trateTST_rev/(planck_hsec*qP_full)

  if (TST_tunn) then
     kappa_Wigner = zero
     trateTST_rev_Wigner = zero
     trateTST_rev_Eckart = zero

     do i=binGAP+1, maxn+binGAP
       trateTST_rev_Eckart = trateTST_rev_Eckart + nosTST_full_Eckart(i)*exp(-dble(i)*Estep/(autocm*kbt))*Estep
     end do

     trateTST_rev_Eckart = trateTST_rev_Eckart/(planck_hsec*qP_full)

     kappa_Wigner=one+(((planck_hsec*TST_freq*cmtoHz)/kbt)**2)/24
     trateTST_rev_Wigner = trateTST_rev*kappa_Wigner
     write(printTST_rev,'(es24.2)') trateTST_rev
     write(printTST_rev_Wigner,'(es24.2)') trateTST_rev_Wigner
     write(printTST_rev_Eckart,'(es24.2)') trateTST_rev_Eckart
     write(66, 203) Temp, trim(adjustl(adjustr(printTST_rev))), &
                       trim(adjustl(adjustr(printTST_rev_Wigner))), &
                       trim(adjustl(adjustr(printTST_rev_Eckart)))
  else
   write(printTST_rev,'(es24.2)') trateTST_rev
   write(66, 204) temp, trim(adjustl(adjustr(printTST_rev)))
  end if
end do
 
201 format(/,/,13x,'Canonical TST reverse rate constant, k(T)',/,/,&
           5x,'T(K)',5x,'k, s-1',4x,'k_Wigner, s-1',4x,'k_Eckart, s-1',/)
202 format(/,/,13x,'Canonical TST reverse rate constant, k(T)',/,/,&
           5x,'T(K)',5x,'k, s-1',/)
203 format(2x,f7.1,5x,a,4x,a,4x,a,4x,a)
204 format(2x,f7.1,5x,a,4x,a)

end subroutine tst_can_rate_rev

end module rate_module
