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
 real(wp)               :: inertX(3), inertR(3), inertP(3)
 complex(wp), allocatable  :: socmat(:), temp(:)
 integer                :: mults(2), i, icp
 real(wp)               :: mecp, gap, redmass, soc, grad, gradmean
 real(wp)               :: hs4,hs3,hs2,hs1,hs0, ls4,ls3,ls2,ls1,ls0
 real(wp)               :: limitL, limitR, T1, T2, Tstep
 complex(wp)            :: h12, h13, h14
 logical                :: zn=.false., rev=.false., upper=.true.
 logical                :: tst=.false., printmore=.false., sp=.false.
 logical                :: solution=.false., extern=.false.
   
 real(wp) :: hbar, beta, lmda, tda
 integer          :: istat
 integer          :: symX=1, symR=1, symP=1, chir=1
 real(wp) :: zpeX, zpeR, zpeP
 integer          :: maxn, binX, binGAP, binZPE, binXrev

 real(wp) :: cxx(3),cyy(3),czz(3)

 contains
 
 subroutine read_file()
  
namelist /keys/ zn,zpe,rev,tst,sp,solution,printmore,extern
namelist /external/ extX, extR, extP 
namelist /inputdata/   symR, symX, chir, &
                       freX, freR, inertX, inertR, enX, enR, &
                       maxn, lmda, tda, T1, T2, Tstep
namelist /probability/ redmass, soc, grad, gradmean
namelist /reverse/ symP, freP, inertP, enP
namelist /polynomials/ hs4,hs3,hs2,hs1,hs0,ls4,ls3,ls2,ls1,ls0, &
      limitL,limitR,upper
namelist /split/ mults, icp, socmat

! Begin reading input data.
  
read (unit=11, nml=keys)

! Allocate larger arrays to make sure you read
! all frequencies provided by user in an input file.

allocate (freX(300))
allocate (freR(300))

! Initialize arrays to zero. The zero elements will be needed later.

freR = zero
freX = zero

read (unit=11, nml=inputdata)

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

if (extern) then
  
  allocate (extR(300))
  allocate (extX(300))
  allocate (extP(300))

  read (unit=11, nml=external)

! For reactant
  allocate (tmp(size(pack(extR, extR /= 0))))
  tmp = pack(extR, extR /= 0)
  deallocate (extR)
  allocate (extR(size(tmp)))
  extR = tmp
  deallocate (tmp)

! For MECP
  allocate (tmp(size(pack(extX, extX /= 0))))
  tmp = pack(extX, extX /= 0)
  deallocate (extX)
  allocate (extX(size(tmp)))
  extX = tmp
  deallocate (tmp)

! For product
  allocate (tmp(size(pack(extP, extP /= 0))))
  tmp = pack(extP, extP /= 0)
  deallocate (extP)
  allocate (extP(size(tmp)))
  extP = tmp
  deallocate (tmp)
  
end if

if (.not. tst) then
  read (unit=11, nml=probability)
end if

if (rev) then
  allocate (freP(300))
  freP = zero
  read (unit=11, nml=reverse)
  allocate (tmp(size(pack(freP, freP /= 0))))
  tmp = pack(freP, freP /= 0)
  deallocate (freP)
  allocate (freP(size(tmp)))
  freP = tmp
  deallocate (tmp)
  if (istat/=0) then
    call error_message(10)!error code in parentheses   
  end if
else
  allocate(freP(1))
 end if

if (sp) then
  allocate(socmat(30))
  read (unit=11, nml=split)
  if (istat/=0) then
    call error_message(16)!error code in parentheses
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
    call error_message(1)!error code in parentheses   
  end if
  if (limitL.ge.limitR) then
    call error_message(5)!error code in parentheses   
  endif
endif

call atomic_units()

end subroutine read_file
!-------------------------------------------------------------------

!-------------------------------------------------------------------
subroutine atomic_units() 
 implicit none

!zpe corrrection
if (zpe /=0 .and. zpe /=1 .and. zpe /=2) then !only three values are accepted
        call error_message(13)!   
endif
!zpe = 0 means no ZPE correction is applied
!zpe = 1 means that no shift will be made for nunmber of states,
! but probability will be shifted by zpeR-zpeX
if (zpe == 1) then
   zpeX = sum(freX)/two
   zpeR = sum(freR)/two
   zpeP = sum(freP)/two
   binZPE = 0 
!zpe = 2 means that probability will be calculted everywhere along the 
!reaction path and later shift will be made in rate.f90 for nunmber of 
!states using binZPE
elseif (zpe == 2) then
   zpeX = zero
   zpeR = zero
   zpeP = zero
   binZPE = floor((sum(freR)-sum(freX))/two)

! Warning. Use of zpe == 2 together with extern = .true.
! can cause problems. If extern = .true., user has better
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
   binZPE = 0
endif

!Reverse rate constant data
if (rev) then
   gap    = enR-enP !in hartree
   binGAP = floor(gap*autocm+zpeR-zpeP)  !number of bins for ZPE corr gap in cm-1
   binXrev = ceiling((enX-enP)*autocm+zpeX-zpeP)
else
   gap    = two
   binGAP = 0
   binXrev = 0
endif

mecp = enX-enR
binX = ceiling(mecp*autocm+zpeX-zpeR)

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
socmat     = socmat/autocm
 
end subroutine atomic_units

end module input_module
