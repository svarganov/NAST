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
subroutine nast_en_rate(dosR,dosX,dosTP,probLZ,probWC,probZN,nosLZ,nosWC,nosZN)
 
real(wp)                             :: WCsto,LZsto,ZNsto
real(wp), allocatable                :: nosWC_TP(:), nosZN_TP(:)
real(wp), allocatable                :: nosWC_X (:), nosZN_X (:)
real(wp), allocatable                :: rateLZ(:),  rateWC(:),  rateZN(:)
real(wp), dimension(:), intent(in)   :: dosR, dosX, dosTP
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
      LZsto = probLZ(j)*dosX(i-j+1+binX) + LZsto 
      WCsto = probWC(j)*dosX(i-j+1+binX) + WCsto
      ZNsto = probZN(j)*dosX(i-j+1+binX) + ZNsto
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
      WCsto = probWC(j)*dosTP(i-j+1) + WCsto 
      ZNsto = probZN(j)*dosTP(i-j+1) + ZNsto 
    enddo
  else
    do j=1, binX
      WCsto = probWC(j)*dosTP(i-j+1) + WCsto
      ZNsto = probZN(j)*dosTP(i-j+1) + ZNsto
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
  write(1451,'(A,5x,A,5x,A,5x,A)') 'E, cm-1','rate LZ','rate WC','rate ZN'

 do i=1, maxn
   write(1451,'(I5,3es24.7)') i, rateLZ(i), rateWC(i), rateZN(i)
 end do
else
 write(1451,'(A,5x,A,5x,A)') 'E, cm-1','rate LZ','rate WC'
 do i=1, maxn
   write(1451,'(I5,2es24.7)') i, rateLZ(i), rateWC(i)
 end do
end if
   
close(1451)
  
if (printmore) then
  open(1415,file='effective_number_of_states.out')
  write(1415,'(A,5x,A,5x,A)') 'E, cm-1','nosLZ','nosWC'
  do i=1, maxn
    write(1415,'(I5,2es24.7)') i, nosLZ(i), nosWC(i)
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
subroutine nast_en_rate_split(dosR,dosX,probLZS,nosLZS)
   real(wp), allocatable                 :: LZSsto(:)
   real(wp), allocatable                 :: rateLZS(:,:)
   real(wp), dimension(:), intent(in)    :: dosR, dosX
   real(wp), dimension(:,:), intent(in)  :: probLZS
   real(wp), dimension(:,:), intent(out) :: nosLZS
   integer                               :: i, j

   allocate (rateLZS  (maxn,icp))
   allocate (LZSsto   (icp))

   rateLZS = zero

   !**********************************************************************
!   OVER-THE-BARRIER CONTRIBUTION TO THE RATE
!      (over-the-barrier number of states) 
!**********************************************************************
   if (maxn.gt.binX) then 
      do i=1+binX, maxn !start from MECP bin
         LZSsto = zero
           do j=1+binX, i ! NOS above crossing point at i-th bin
              LZSsto(:) = probLZS(j,:)*dosX(i-j+1+binX) + LZSsto(:)
           enddo
         nosLZS(i,:)   = LZSsto(:)/autocm
      enddo
  endif

     !*****************************************
   !      Microcanonical Rate Constant 
   !*****************************************
   !binZPE = 0 if no ZPE correction or ZPE1 correction is used 
  do i=1, maxn-binZPE
   if (dosR(i).ne.zero) then !Eq. (1) ref_1
      rateLZS(i,:) = degenR*nosLZS(i+binZPE,:)/(dosR(i)*planck_hsec)
   end if
  end do
  
  open(556,file='split_energy_rate.out')
  write(556,'(a)', advance="no") "Energy(cm-1)"
  write(556, '(4X,I2,A,I2,A,I2,A,I2)', advance="no") mults(1), "_", 1, &
                                                   "  ->",mults(2), "_", 1
  do i=2, icp
   write(556, '(11X,I2,A,I2,A,I2,A,I2)', advance="no") mults(1), "_", 1+floor((real(i)-1)/3), &
                                            " ->", mults(2), "_", mod(i-1,3)+1+floor((real(i)-1)/3)
  end do
  
  write(556,'(/)')
 
  do i=1, maxn
      write(556,*) i, rateLZS(i,:)
  enddo
close(556)
  
deallocate (rateLZS)
deallocate (LZSsto)

end subroutine nast_en_rate_split
!============================================================================
subroutine tst_en_rate(dosR,dosX,nosTST)
implicit none

   integer                            :: i
   real(wp), allocatable              :: rateTST(:)
   real(wp), allocatable, intent(in)  :: dosR(:), dosX(:)
   !real(wp), allocatable, intent(out) :: nosTST(:)
   real(wp), dimension(:), intent(out) :: nosTST
   
   allocate (rateTST  (maxn))
   
   write(66,'(/,A)') 'Calculating microcanonical TST rate constant.'

   nosTST = zero
   rateTST = zero
   
   do i = 1+binX, maxn ! Start from the MECP bin
      nosTST(i) = sum(dosX(1+binX:i))
   end do

   nosTST = nosTST/autocm
   
   !do i = 1, maxn-binZPE
   !   if (dosR(i) .ne. 0.0d0) then
   !      rateTST(i) = degen*nosTST(i+binZPE)*autoev/(dosR(i)*planck)
   !   end if
   !enddo
   
   where (dosR .ne. zero)
       rateTST(1:maxn-binZPE) = degenR*nosTST(1+binZPE:maxn)&
       *autoev/(dosR(1:maxn-binZPE)*planck)
   end where

   open(322,file='energy_tst_rate.out')
      write(322,'(a)') 'E, cm-1     ', &
                    'TST rate'
   do i=1, maxn
      write(322,200) i, rateTST(i)
   enddo
   close(322)
  200 format('',i5,es24.7)
 
  write(66,'(A)') '................................... Done.'
  write(66,'(/,A,I1,A)') '     (The forward rate constant k(E) is multiplied by &
 reaction path degeneracy equal to ',int(degenR),' )'

end subroutine tst_en_rate
!==================================================================================
!============== MACROCANONICAL (THERMAL) RATE CONSTANT k(T) =======================
!==================================================================================
subroutine nast_can_rate(dosR,nosLZ,nosWC,nosZN)
   implicit none
   character(len=33)                       :: printLZ, printWC, printZN
   integer                                 :: i,j,T_n_points
   real(wp)                                :: qR, kbt, Temp 
   real(wp)                                :: trateLZ, trateWC, trateZN
   real(wp), dimension(:), intent(in)      :: dosR, nosLZ, nosWC, nosZN
 
   Temp = zero
   kbt = zero

   if (zn) then
       write(66, 200)
   else
       write(66, 201)
   end if

    do j = 0, Tpoints
      Temp = T1 + dble(j)*Tstep
   
      trateLZ = zero
      trateWC = zero
      trateZN = zero
 
      qR      = zero !Partition function for reactants
      kbt = kb_hartK*Temp !kb is Boltzmann constant
 
      do i=1, maxn
        qR=qR + dosR(i)*exp(-dble(i)/(autocm*kbt))
      enddo
 
      do i=1, maxn
        trateLZ  = trateLZ + nosLZ(i)*exp(-dble(i)/(autocm*kbt)) 
        trateWC  = trateWC + nosWC(i)*exp(-dble(i)/(autocm*kbt)) 
        trateZN  = trateZN + nosZN(i)*exp(-dble(i)/(autocm*kbt)) 
      end do
 
      trateLZ = trateLZ/(planck_hsec*qR) 
      trateWC = trateWC/(planck_hsec*qR) 
      trateZN = trateZN/(planck_hsec*qR) 

      if (zn) then
         write(printLZ,'(es24.2)') trateLZ
         write(printWC,'(es24.2)') trateWC
         write(printZN,'(es24.2)') trateZN
         write(66, 202) Temp, trim(adjustl(adjustr(printLZ))), &
                              trim(adjustl(adjustr(printWC))), &
                              trim(adjustl(adjustr(printZN)))
      else
         write(printLZ,'(es24.2)') trateLZ
         write(printWC,'(es24.2)') trateWC
         write(66, 203) Temp, trim(adjustl(adjustr(printLZ))), &
                              trim(adjustl(adjustr(printWC)))
      end if
    
    enddo
  
  200 format(/,/,13x,'Canonical rate constant, k(T).',/,/,&
            5x,'T(K)',3x,'Landau-Zener',3x,'Weak Coupling',3x,'Zhu Nakamura',/)
  201 format(/,/,13x,'Canonical rate constant, k(T)',/,/,&
            5x,'T(K)',3X,'Landau-Zener',3x,'Weak Coupling',/)

  202 format(4x,f7.1,4x,a,8x,a,8x,a)
  203 format(4x,f7.1,4x,a,8x,a)

end subroutine nast_can_rate
!================================================================
subroutine tst_can_rate(dosR,nosTST)
   implicit none
   character(len=33)                          :: printTST
   integer                                    :: i,j,T_n_points
   real(wp)                                   :: qR, kbt, Temp 
   real(wp)                                   :: trateTST
   real(wp), allocatable, intent(in)          :: dosR(:)
   real(wp), allocatable, intent(in)          :: nosTST(:)
 
   Temp = zero
   kbt = zero

   write(66, 201)
   
   do j = 0, Tpoints
      Temp = T1 + dble(j)*Tstep
   
      trateTST = zero
 
      qR      = zero !Partition function for reactants
      kbt = kb_hartK*Temp !kb is Boltzmann constant
 
      do i=1, maxn
        qR=qR + dosR(i)*exp(-dble(i)/(autocm*kbt))
      enddo
 
      do i=1, maxn
        trateTST  = trateTST + nosTST(i)*exp(-dble(i)/(autocm*kbt)) 
      end do
 
      trateTST = trateTST/(planck_hsec*qR) 
      write(printTST,'(es24.2)') trateTST
      write(66, 203) Temp, trim(adjustl(adjustr(printTST)))

   enddo
  
  201 format(/,/,13x,'Canonical rate constant, k(T)',/,/,&
            5x,'T(K)',5x,'TST rate constant',/)

  203 format(2x,f7.1,5x,a)
end subroutine tst_can_rate
!==================================================================================
!================ MICROCANONICAL REVERSE RATE CONSTANT k(E) =======================
!==================================================================================
subroutine nast_en_rate_rev(dosP_full,dosX_full,dosX_TP_rev,probLZ_rev,probWC_rev,nosLZ_full,nosWC_full)
   
 real(wp), dimension(:), intent(in)  :: dosP_full,dosX_full,dosX_TP_rev,probLZ_rev,probWC_rev
 real(wp), dimension(:), intent(out) :: nosLZ_full,nosWC_full
 real(wp), allocatable               :: REVrateLZ(:), REVrateWC(:),nosWC_X(:),nosWC_TP(:)
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
     LZsto = probLZ_rev(j)*dosX_full(i-j+1+binXrev) + LZsto
     WCsto = probWC_rev(j)*dosX_full(i-j+1+binXrev) + WCsto
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
     WCsto = probWC_rev(j)*dosX_TP_rev(i-j+binGAP+1) + WCsto
    end do
  else
    do j  = binGAP+1, binXrev
     WCsto = probWC_rev(j)*dosX_TP_rev(i-j+binGAP+1) + WCsto
    end do
  end if
  nosWC_TP(i) = WCsto/autocm
end do


!*****************************
!   Total number of states
!*****************************

nosWC_full = nosWC_X + nosWC_TP

! Microcanonical reverse rate constant.

do i = 1, maxn + binGAP
  if (dosP_full(i).ne.0.0d0) then 
    REVrateLZ(i) = degenP*nosLZ_full(i)/(dosP_full(i)*planck_hsec)
    REVrateWC(i) = degenP*nosWC_full(i)/(dosP_full(i)*planck_hsec)
  end if
end do

open(144,file='energy_rate_reverse.out')
write(144,'(a,4x,a,4x,a)') 'E, cm-1','Reverse rate LZ','Reverse rate WC'
   
do i = 1, maxn + binGAP
  write(144,200) i, REVrateLZ(i), REVrateWC(i)
end do

if (printmore) then
  open(1444,file='effective_number_of_states_rev.out')
  write(1444,'(a,4x,a,4x,a)') 'E, cm-1','LZ','WC'
   
  do i = 1, maxn + binGAP
    write(1444,'(i5,es24.7,4x,es24.7)') i,nosLZ_full(i),nosWC_full(i)
  end do

close(1444)
end if
 
deallocate (REVrateLZ)  
deallocate (REVrateWC)  
  
close(144)
200 format(i5,2es24.15)

end subroutine nast_en_rate_rev
!==============================================================================
!============= MACROCANONICAL (THERMAL) REVERSE RATE CONSTNAT k(T) ============
!==============================================================================
subroutine nast_can_rate_rev(dosP_full,nosLZ_full,nosWC_full)
   
implicit none
real(wp), dimension(:), intent(in)  :: dosP_full,nosLZ_full,nosWC_full
character(len=33)                   :: printLZ, printWC
integer                             :: i,j
real(wp)                            :: qP_full,kbt, temp 
real(wp)                            :: rev_trateLZ, rev_trateWC
 
temp = zero
kbt = zero

write(66, 201)

do j = 0, Tpoints

  temp = T1 + dble(j)*Tstep
  rev_trateLZ = zero
  rev_trateWC = zero
  qP_full = zero
  kbt = kb_hartK*Temp !kb is Boltzmann constant
 
  do i=1, maxn+binGAP
    qP_full = qP_full + dosP_full(i)*exp(-dble(i)/(autocm*kbt))
    rev_trateLZ = rev_trateLZ + nosLZ_full(i)*exp(-dble(i)/(autocm*kbt))
    rev_trateWC = rev_trateWC + nosWC_full(i)*exp(-dble(i)/(autocm*kbt)) 
  end do
 
  rev_trateLZ = rev_trateLZ/(planck_hsec*qP_full) 
  rev_trateWC = rev_trateWC/(planck_hsec*qP_full) 

  write(printLZ,'(es24.2)') rev_trateLZ
  write(printWC,'(es24.2)') rev_trateWC
  write(66, 203) Temp, trim(adjustl(adjustr(printLZ))), &
                         trim(adjustl(adjustr(printWC)))

end do
  
201 format(/,/,13X,'Canonical reverse rate constant, k(T)',/,/,&
           5x,'T(K)',5x,'Landau-Zener',5x,'Weak Coupling',/)
203 format(4x,f7.1,4x,a,8x,a)

end subroutine nast_can_rate_rev

! Traditional TST microcanonical reverse rate constant.

subroutine tst_en_rate_rev(dosP_full,dosX_full,nosTST_full)
implicit none

real(wp), dimension(:), intent(in)  :: dosP_full, dosX_full
real(wp), dimension(:), intent(out) :: nosTST_full
real(wp), allocatable               :: rateTST_rev(:)
integer                             :: i

allocate (rateTST_rev  (maxn+binGAP))
nosTST_full = zero
rateTST_rev = zero
write(66,'(/,A)') 'Calculating microcanonical TST reverse rate constant.'

do i = 1+binXrev, maxn+binGAP
  nosTST_full(i) = sum(dosX_full(1+binXrev:i))
end do

nosTST_full = nosTST_full/autocm

where (dosP_full .ne. zero)
     rateTST_rev(1:maxn+binGAP-binZPE) = degenP*nosTST_full(1+binZPE:maxn+binGAP)&
     *autoev/(dosP_full(1:maxn+binGAP-binZPE)*planck)
end where

open (307, file = 'energy_tst_rate_rev.out')
write(307,'(a,4x,a)') 'E, cm-1','TST reverse rate'
do i = 1, maxn + binGAP
  write(307,199) i, rateTST_rev(i)
enddo
close(307)
199 format('',i5,es24.7)
 
write(66,'(A)') '................................... Done.'
write(66,'(/,A,I1,A)') '     (The reverse TST rate constant k(E) is multiplied by &
reaction path degeneracy equal to ',int(degenP),' ).'

end subroutine tst_en_rate_rev

! Traditional TST canonical reverse rate constant.

subroutine tst_can_rate_rev(dosP_full,nosTST_full)
implicit none

real(wp), dimension(:), intent(in)  :: dosP_full, nosTST_full
character(len=33)                   :: printTST_rev
real(wp)                            :: qP_full, trateTST_rev, temp, kbt
integer                             :: i,j

temp = zero
kbt = zero
write (66, 201)

do j = 0, Tpoints
  temp = T1 + dble(j)*Tstep
  trateTST_rev = zero
  qP_full = zero
  kbt = kb_hartK*Temp !kb is Boltzmann constant
 
  do i=1, maxn+binGAP
    qP_full = qP_full + dosP_full(i)*exp(-dble(i)/(autocm*kbt))
    trateTST_rev  = trateTST_rev + nosTST_full(i)*exp(-dble(i)/(autocm*kbt)) 
  end do
 
  trateTST_rev = trateTST_rev/(planck_hsec*qP_full) 
  write(printTST_rev,'(es24.2)') trateTST_rev
  write(66, 203) temp, trim(adjustl(adjustr(printTST_rev)))

end do
 
201 format(/,/,13x,'Canonical TST reverse rate constant, k(T).',/,/&
            5x,'T(K)',3x,'TST reverse',/)
203 format(4x,f7.1,4x,a)

end subroutine tst_can_rate_rev

end module rate_module
