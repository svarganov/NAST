module input_module
 use precision_module, only : wp
 use error_module
 use constants_module
 implicit none
 public
   
 integer                :: zpe = 0, Tpoints
 real(wp), allocatable  :: freX(:), freR(:), freP(:)
 real(wp), allocatable  :: extX(:), extR(:), extP(:),tmp(:)
 real(wp)               :: enX, enR, enP, degenR, degenP
 real(wp)               :: TST_freq = 0
 real(wp)               :: inertX(3), inertR(3), inertP(3)
 complex(wp), allocatable  :: socmat(:), temp(:)
 integer                :: mults(2), i, icp
 real(wp)               :: delta(2), Bz, g
 real(wp)               :: mecp, gap, redmass, soc, grad, gradmean, Estep
 real(wp)               :: hs4,hs3,hs2,hs1,hs0, ls4,ls3,ls2,ls1,ls0
 real(wp)               :: limitL, limitR, T1, T2, Tstep
 complex(wp)            :: h12, h13, h14
 logical                :: zn=.false., rev=.false., upper=.true.
 logical                :: tst=.false., printmore=.false., sp=.false.
 logical                :: externR=.false.,  externX=.false.,  externP=.false.
 logical                :: solution=.false., TST_tunn=.false.
 real(wp) :: hbar, beta, lmda, tda
 integer          :: istat
 integer          :: symX=1, symR=1, symP=1, chir=1
 real(wp) :: zpeX, zpeR, zpeP
 integer          :: maxn, binX, binX1, binX3, binGAP, binZPE, binZPE_P, binXrev

 real(wp) :: cxx(3),cyy(3),czz(3)

 contains
 
 subroutine read_file()
  
namelist /keys/ zn,zpe,rev,tst,sp,solution,printmore,externR, & 
                externX, externP, TST_tunn
namelist /external/ extX, extR, extP 
namelist /inputdata/   symR, symX, chir, &
                       freX, freR, inertX, inertR, enX, enR, &
                       maxn, lmda, tda, T1, T2, Tstep, TST_freq, Estep
namelist /probability/ redmass, soc, grad, gradmean
namelist /reverse/ symP, freP, inertP, enP
namelist /polynomials/ hs4,hs3,hs2,hs1,hs0,ls4,ls3,ls2,ls1,ls0, &
      limitL,limitR,upper
namelist /split/ mults, icp, socmat, delta, Bz, g

! Begin reading input data.
  
read (unit=11, nml=keys)

! Allocate larger arrays to make sure you read
! all frequencies provided by user in an input file.

allocate (freX(400))
allocate (freR(400))

! Initialize arrays to zero. The zero elements will be needed later.

freR = zero
freX = zero

read (unit=11, nml=inputdata)

! Redifine maxn

if (Estep == zero) then
  Estep = one
end if

maxn = nint(maxn/Estep)

! Below is automated search for number of non-zero elements in
! freR, freX and freP. These arrays will be re-allocated to contain
! only non-zero elements.

! For reactant.

! size(pack(freR, freR /= 0)) - number of non-zero elements
! in freR read from the 'inputdata' namelist.

allocate (tmp(size(pack(freR, freR /= 0))))
tmp = pack(freR, freR /= 0)
deallocate (freR)

! Allocate freR again, but this time the size of array is equal to
! the number of non-zero elements in array freR read from namelist.

allocate (freR(size(tmp)))
freR = tmp
deallocate (tmp)

! For MECP.

allocate (tmp(size(pack(freX, freX /= 0))))
tmp = pack(freX, freX /= 0)
deallocate (freX)
allocate (freX(size(tmp)))
freX = tmp
deallocate (tmp)

! ------------------------------------------------------------------
! One line command that does the same as above. No need to use tmp.

!            freR=pack(freR, freR /= 0)
!            freX=pack(freX, freX /= 0)

! However, does not work with (old) Intel compilers.
! Feel free to try it instead of allocate(tmp1( .. deallocate (tmp1)
! and allocate(tmp2( .. deallocate (tmp2).
!------------------------------------------------------------------

if (T1 == zero) then
  T1  = 290_wp
end if
if (T2 == zero) then
   T2 = 300_wp
end if
if (T2 < T1) then
   T2 = T1
end if
if (Tstep == zero) then
  Tstep = one
end if
Tpoints = int((T2 - T1)/Tstep)

if (externR.or.externX.or.externP) then
  if (externR) then
    allocate (extR(400)) 
  end if
 
  if (externX) then
    allocate (extX(400)) 
  end if

  if (externP) then
    allocate (extP(400)) 
  end if

  read (unit=11, nml=external)
  
end if

if (externR) then
  allocate (tmp(size(pack(extR, extR /= 0))))
  tmp = pack(extR, extR /= 0)
  deallocate (extR)
  allocate (extR(size(tmp)))
  extR = tmp
  deallocate (tmp)
else
  allocate (extR(1))
  extR = zero
end if

if (externX) then
  allocate (tmp(size(pack(extX, extX /= 0))))
  tmp = pack(extX, extX /= 0)
  deallocate (extX)
  allocate (extX(size(tmp)))
  extX = tmp
  deallocate (tmp)
else
  allocate (extX(1))
  extX = zero
end if

if (externP) then
  allocate (tmp(size(pack(extP, extP /= 0))))
  tmp = pack(extP, extP /= 0)
  deallocate (extP)
  allocate (extP(size(tmp)))
  extP = tmp
  deallocate (tmp)
else
  allocate (extP(1))
  extP = zero
end if

if (.not. tst) then
  read (unit=11, nml=probability)
end if

if (rev) then
  allocate (freP(400))
  freP = zero
  read (unit=11, nml=reverse)
  allocate (tmp(size(pack(freP, freP /= 0))))
  tmp = pack(freP, freP /= 0)
  deallocate (freP)
  allocate (freP(size(tmp)))
  freP = tmp
  deallocate (tmp)
  if (istat/=0) then
    call error_message(10) 
  end if
else
  allocate(freP(1))
end if

if (sp) then
  allocate(socmat(30))
  read (unit=11, nml=split)
  if (istat/=0) then
    call error_message(16)
  end if
  allocate(temp(icp))
  do i=1, icp
    temp(i) = socmat(i)
  enddo
  deallocate(socmat)
  allocate(socmat(icp))
  socmat = temp
  deallocate(temp)
else
  allocate(socmat(1))
end if

if (zn) then
  read (unit=11, nml=polynomials, iostat=istat)
  if (istat/=0) then
    call error_message(1)
  end if
  if (limitL.ge.limitR) then
    call error_message(5)
  endif
endif

if (TST_tunn.and.TST_freq.eq.0) then
  call error_message(18)
endif

if ((TST_tunn).and.(.not.rev)) then
  call error_message(19)
endif

call atomic_units()

end subroutine read_file
!-------------------------------------------------------------------

!-------------------------------------------------------------------
subroutine atomic_units() 
 implicit none

!zpe corrrection
if (zpe /=0 .and. zpe /=1 .and. zpe /=2) then !only three values are accepted
        call error_message(13)   
endif
!zpe = 0 means no ZPE correction is applied
!zpe = 1 means that no shift will be made for nunmber of states,
! but probability will be shifted by zpeR-zpeX
if (zpe == 1) then
   zpeX = sum(freX)/two
   zpeR = sum(freR)/two
   zpeP = sum(freP)/two
   binZPE = zero
   binZPE_P = zero
!zpe = 2 means that probability will be calculted everywhere along the 
!reaction path and later shift will be made in rate.f90 for nunmber of 
!states using binZPE
elseif (zpe == 2) then
   zpeX = zero
   zpeR = zero
   zpeP = zero
   binZPE = floor((sum(freR)-sum(freX))/two)
   binZPE_P = floor((sum(freP)-sum(freX))/two)
! Warning. Use of zpe == 2 together with extern = .true.
! can cause problems. If extern = .true., user had better
! deleted fundamental frequencies from, say, freR, which are
! replaced by extR. However, since zpeR is calculated from freR,
! and due to frequencies of excluded modes had been removed,
! zpeR can result to be smaller than zpeX. Therefore,
! binZPE above will be negative! It will cause negative
! array indices later.
! One solution is: if binZPE < 0, binZPE = 0.

  if (binZPE < 0) then
    binZPE = 0
  end if

!No ZPE correction is applied
else
   zpeX = zero
   zpeR = zero
   zpeP = zero
   binZPE = zero
   binZPE_P = zero
endif

!Reverse rate constant data
if (rev) then
   gap    = enR-enP !in hartree
   binGAP = floor(gap*autocm+zpeR-zpeP)
   binXrev = ceiling((enX-enP)*autocm+zpeX-zpeP)
       
       if (Estep.ne.one) then
           binXrev = ceiling(binXrev/Estep)
           binGAP = ceiling(binGAP/Estep)
       endif
   
   if (binGAP.lt.zero) then
     call error_message(20)
   end if
else
   gap    = two
   binGAP = 0
   binXrev = 0
endif

mecp = enX-enR
binX = ceiling(mecp*autocm+zpeX-zpeR)

if (binX.lt.zero) then
  call error_message(21)
end if

if (Estep.ne.one) then
    binX = ceiling(binX/Estep)
endif

if (sp) then
   socmat     = socmat/autocm
   delta      = delta/autocm
   binX1      = ceiling(mecp*autocm-g*mu*Bz/2.0+zpeX-zpeR) !number of bins below (MECP-gmu.Bz/2) in cm-1 for p1
   binX3      = ceiling(mecp*autocm+g*mu*Bz/2.0+zpeX-zpeR) !number of bins below (MECP+gmu.Bz/2) in cm-1 for p3
endif

! Degeneracy of the reaction path.

degenR = dble(symR*chir/symX)
degenP = dble(symP*chir/symX)

! Set all properties to atomic units.

inertX   = inertX*amutoau
inertR   = inertR*amutoau 
inertP   = inertP*amutoau
freX     = freX/autocm
freR     = freR/autocm
freP     = freP/autocm
redmass  = redmass*amutoau
soc      = soc/autocm
zpeX     = zpeX/autocm
zpeR     = zpeR/autocm
zpeP     = zpeP/autocm
 
end subroutine atomic_units

end module input_module
