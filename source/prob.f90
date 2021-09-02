module prob_module
 use precision_module, only : wp
 use input_module
 use constants_module
 use zn_adiabatic_module, &
   only : zn_path_sloped, zn_path_peaked, set_slope, set_origin
 use zn_param_module,  only : zn_param_sloped, zn_param_peaked
 use zn_phase_module,  only : zn_zero_phase, zn_phase_sloped, zn_phase_peaked
 use zn_prob_module,   only : zn_prob_sloped, zn_prob_peaked

 implicit none
 public

 real(wp) :: kbt

contains

! Spin-diabatic parameters aaDia & bbDia.

subroutine diab_param(aaDia,bbDia)
 implicit none
 integer                             :: i
 real(wp)                            :: a
 real(wp), intent(out)               :: aaDia
 real(wp), dimension(:), intent(out) :: bbDia
 
aaDia=grad*gradmean/(16.0*redmass*soc**3.0)
a=sqrt(aaDia)

if (a .le. one) then
  write(66, '(/,1X,A,1X,f6.2,A)') 'WC formula is valid when a >> 1. Calculated "a" is', a,'.'   
  write(66,'(A)')'WC formula IS NOT VALID for this system.'   
endif

do i=1, maxn
  bbDia(i)=((dble(i)/autocm-mecp+zpeR-zpeX)*grad)/(2.0*soc*gradmean)
enddo

do i=1, maxn
  if (bbDia(i) .gt. one .and. bbDia(i).gt.sqrt(aaDia/bbDia(i))) then
     write(66,'(/,6X,A,I5,A)') 'LZ formula is valid at energies much greater than', i, ' cm-1.'
     exit !stop cycle  LZ valid when bb>>1 and bb>>a/b 
  endif
enddo

end subroutine diab_param

! Landau-Zener double passage formula.

subroutine LZ(aaDia,bbDia,  probLZ)

integer                             :: i
real(wp)                            :: p
real(wp),               intent(in)  :: aaDia
real(wp), dimension(:), intent(in)  :: bbDia
real(wp), dimension(:), intent(out) :: probLZ

probLZ = zero

write(66,'(/,A)') "........Landau-Zener"

where (bbDia .gt. zero) probLZ = one-exp(-pi/(2.0*sqrt(bbDia*aaDia)))
  
end subroutine LZ

! Landau-Zener double passage for reverse rate constant.
! Now the energy is spanned up to maxn+binGAP, so need
! to calculate probability at new points. It is fast enough to
! re-calculate probability at all points (1:maxn+binGAP) so no need
! to copy existing values in 1:maxn.

subroutine LZ_rev(aaDia,bbDia_rev,probLZ_rev)
real(wp),               intent(in)  :: aaDia
real(wp), dimension(:), intent(out) :: bbDia_rev,probLZ_rev
integer                             :: i

do i=1, maxn+binGAP
  bbDia_rev(i)=((dble(i)/autocm-(enX-enP)+zpeP-zpeX)*grad)/(2.0*soc*gradmean)
enddo

do i=1, maxn+binGAP
  if (bbDia_rev(i) .gt. one .and. bbDia_rev(i).gt.sqrt(aaDia/bbDia_rev(i))) then
     write(66,'(/,6X,A,I5,A)') 'Reverse LZ probability is valid at energies much greater than', i, ' cm-1.'
     exit 
  endif
enddo

where (bbDia_rev .gt. zero) probLZ_rev = one-exp(-pi/(2.0*sqrt(bbDia_rev*aaDia)))

end subroutine LZ_rev

!=====================================================================
!======================= LANDAU-ZENER SPIN SPLIT =====================
!=====================================================================
subroutine LZS(probLZS)
!calculates diabatic transition between state specific spin states using double-passage LZ
!Returns 2d array of transition probabilities.
!Requires sp = .true.
integer               :: i
real(wp)              :: const, mecp_zpe, lz_coeff
real(wp), allocatable :: sqnorm(:), p(:)
real(wp), dimension(:,:), intent(out) :: probLZS

allocate (sqnorm(icp))
allocate (p(icp))

mecp_zpe = mecp - zpeR + zpeX
const = (-2.0 * pi / grad)
probLZS = zero
sqnorm = socmat * conjg(socmat)      !sqnorm is squared norm of spin-orbit matrix elements.

!Starts at the MECP, probabilities below MECP are zero.
do i = binX, maxn
   lz_coeff = redmass / (2.0 * (dble(i)/autocm - mecp_zpe))
   lz_coeff = const * sqrt(lz_coeff)
   p = sqnorm * lz_coeff
   p = exp(p)
   p = 1.0 - p*p
   probLZS(i, :) = p(:)
end do

end subroutine LZS

! Weak Coupling double passage formula.

subroutine WC(aaDia,bbDia,  probWC)
 integer                             :: i, NZ, IERR
 complex(wp)                         :: airyvalue
 real(wp)                            :: AII, airytemp
 real(wp),              intent(in)   :: aaDia
 real(wp), dimension(:), intent(in)  :: bbDia
 real(wp), dimension(:), intent(out) :: probWC


probWC    = zero
airytemp  = zero
airyvalue = (zero,zero)
 
write(66,'(/,a,/,a)') "Warning! Rate calculations using transition probabilities that account",&
                      "for quantum tunneling should be run with caution."
write(66,'(a)') "Make sure that in a reaction studied, reactant lies higher in energy than product."
write(66,'(a,/,a,/,a,/,a)') "Otherwise, there will be a region of unphysical tunneling, which should be excluded.",&
                      "If this is the case, run a reverse (rev = .true.) calculation. This will allow NAST to prevent",&
                      "unphysical tunneling. Details on reverse rate calculations can be found in the NAST package manual." 
write(66,'(/,a)') "..........Weak Coupling"

do i=1, maxn            
  airytemp = -bbDia(i)*(aaDia**(-1.0/3.0))
  call ZAIRY(airytemp, zero, 0, 1, airyvalue, AII, NZ, IERR)
  if (NZ .eq. 0 .and. IERR .eq. 0) then
    probWC(i) = pi**2.0*aaDia**(-2.0/3.0)*(abs(airyvalue))**2.0
    if (probWC(i) .gt. one) then
      write(66,'(A,I5,A)') 'Unphysical WC probability greater than 1.0 at ',i,' cm-1'
      stop   
    endif
  end if
end do

end subroutine WC

! Weak Coupling double passage for reverse rate constant.
! Now the energy is spanned up to maxn+binGAP, so need
! to calculate probability at new points. It is fast enough to
! re-calculate probability at all points (1:maxn+binGAP) so no need
! to copy existing values in 1:maxn.

subroutine WC_rev(aaDia,bbDia_rev,  probWC_rev)
implicit none

integer                             :: i, NZ, IERR
complex(wp)                         :: airyvalue
real(wp)                            :: AII, airytemp
real(wp),               intent(in)  :: aaDia
real(wp), dimension(:), intent(in)  :: bbDia_rev
real(wp), dimension(:), intent(out) :: probWC_rev

probWC_rev    = zero
airytemp      = zero
airyvalue     = (zero,zero)
 
write(66,'(A)') "..........Weak Coupling for reverse rate constant."

do i = 1, maxn+binGAP            
  airytemp = -bbDia_rev(i)*(aaDia**(-1.0/3.0))
  call ZAIRY(airytemp, zero, 0, 1, airyvalue, AII, NZ, IERR)
  if (NZ .eq. 0 .and. IERR .eq. 0) then
    probWC_rev(i) = pi**2.0*aaDia**(-2.0/3.0)*(abs(airyvalue))**2.0
    if (probWC_rev(i) .gt. one) then
      write(66,'(A,I5,A)') 'Unphysical WC probability greater than 1.0 at ',i,' cm-1 in reverse direction.'
      stop   
    endif
  end if
end do

probWC_rev(1:binGAP) = zero

end subroutine WC_rev

! Zhu-Nakamura double passage formula

subroutine ZhN(probZN)
 character(6)          :: crosstype
 integer               :: binLower, binUpper
 real(wp)              :: aa, dd, up_origin, down_origin, r0
 real(wp)              :: origin, xcross, mecp
 real(wp), allocatable :: bb(:), arrayT1(:), arrayT2(:)
 real(wp), allocatable :: xpoint(:), up(:), down(:)
 real(wp), allocatable :: sigma0(:), delta0(:)
 real(wp), allocatable :: sigma(:), delta(:)
 real(wp), allocatable :: arT1tun(:), arT2tun(:)
 real(wp), allocatable :: arT1ref(:), arT2ref(:)
 real(wp), dimension(:), intent(inout) :: probZN

probZN = zero

allocate (bb      (maxn))
allocate (xpoint  (maxn))
allocate (up      (maxn))
allocate (down    (maxn))
allocate (sigma   (maxn))
allocate (delta   (maxn))

write(66,'(A)') "...........Zhu-Nakamura"

call set_slope(hs0,hs1,hs2,hs3,hs4,ls0,ls1,ls2,ls3,ls4,limitL,limitR, &
               crosstype,xcross,mecp)

call set_origin(limitL,limitR, origin)
select case (crosstype)
  case ("sloped")
   allocate (arrayT1 (maxn))
   allocate (arrayT2 (maxn))
   allocate (sigma0  (maxn))
   allocate (delta0  (maxn))
! Do formal allocation of non-used array(s).
   allocate (arT1tun (1))
   allocate (arT2tun (1))
   allocate (arT1ref (1))
   allocate (arT2ref (1))

   call zn_path_sloped(origin,  xpoint,up,down,up_origin,arrayT1,arrayT2)
   
   call zn_param_sloped(up_origin,xcross,xpoint,up,down,  &
                        aa,dd,bb,r0,binLower,binUpper)
   
   call zn_zero_phase(aa,dd,bb,  sigma0,delta0)
   
   call zn_phase_sloped(binLower,binUpper,up_origin,r0,   &
                        arrayT1,arrayT2,sigma0,delta0, &
                        sigma,delta)
   
   call zn_prob_sloped(aa,bb,sigma,delta,probZN)
    
  
  case ("peaked")
   call zn_param_peaked(origin,xcross,            &
                        aa,bb,binLower,binUpper,down_origin)

   allocate (arT1tun (binLower-1))
   allocate (arT2tun (binLower-1))
   allocate (arT1ref (maxn-binUpper))
   allocate (arT2ref (maxn-binUpper))
! Do formal allocation of non-used array(s).
   allocate (arrayT1 (1))
   allocate (arrayT2 (1))
   allocate (sigma0  (1))
   allocate (delta0  (1))
 
   call zn_path_peaked(binLower,binUpper,down_origin,      & 
                       arT1tun,arT2tun,arT1ref,arT2ref)


   call zn_phase_peaked(binLower,binUpper,down_origin,aa,bb, &
                        arT1tun,arT2tun,arT1ref,arT2ref,  &
                        sigma,delta)

   call zn_prob_peaked(binLower,binUpper,aa,bb,sigma,delta,probZN) 
end select

! Deallocation.

deallocate (bb     )
deallocate (arrayT1)
deallocate (arrayT2)
deallocate (up     )
deallocate (down   )
deallocate (xpoint )
deallocate (sigma0 )
deallocate (delta0 )
deallocate (sigma  )
deallocate (delta  )
deallocate (arT1tun)
deallocate (arT2tun)
deallocate (arT1ref)
deallocate (arT2ref)

end subroutine ZhN

subroutine ZNempty(prob)
 real(wp), allocatable, intent(out) :: prob(:)

allocate (prob(maxn))
   
prob = 0.0
! Generates empty probabilities to avoid conflicts with printing data.

end subroutine ZNempty

subroutine LZempty(prob)
 real(wp), allocatable, intent(out) :: prob(:,:)

allocate (prob(maxn,4))
mults(1) = 1
mults(2) = 3

prob = 0.0
! Generates empty probabilities to avoid conflicts with printing data.
end subroutine LZempty

! Print probabilities of transition.

subroutine write_prob(aaDia,bbDia,probLZ,probWC,probZN,probLZS,bbDia_rev,probLZ_rev,probWC_rev)
implicit none
 
real(wp),                 intent(in) :: aaDia
real(wp), dimension(:),   intent(in) :: bbDia, probLZ, probWC, probZN
real(wp), dimension(:,:), intent(in) :: probLZS
real(wp), dimension(:),   intent(in) :: bbDia_rev, probLZ_rev, probWC_rev
integer                              :: i

open(15,file='energy_probabilities.out')
write(15,*) "Double Passage Transition Probabilities."
   
if (.not. zn) then
  write(15,'(A,7x,A,7x,A)') "Energy(cm-1)","Landau-Zener","Weak Coupling" 
  do i=1,maxn
    write (15,'(I6,2es24.5)') i, probLZ(i), probWC(i)
  end do
else
  write(15,'(A,7x,A,7x,A,7x,A)') "Energy(cm-1)","Landau-Zener","Weak Coupling","Zhu-Nakamura"
  do i=1,maxn
    write (15,'(I6,3es24.5)') i, probLZ(i), probWC(i),probZN(i)
  end do
end if
 
close(15)

if (sp) then
  open(16,file='split_probabilities.out')
  write(16,'(a)') "Double passage transition probabilities between two Ms spin states"
  write(16, '(a)', advance="no") "Energy (cm-1)"
  write(16, '(4x,I2,a,I2,a,I2,a,I2)', advance="no") mults(1), " _", 1, &
                                    " ->", mults(2), " _", 1
  do i=2, icp
    write(16, '(11x,I2,a,I2,a,I2,a,I2)', advance="no") mults(1), " _", 1+floor((real(i)-1)/3), &
                                    " ->", mults(2), " _", mod(i-1,3)+1+floor((real(i)-1)/3)
  end do
  write(16,'(/)')

  do i=1,maxn
    write (16,*) i, probLZS(i,:)
  end do

close(16)

end if

if (rev) then
  open(17,file='energy_probabilities_rev.out')
  write(17,'(a)') "Double passage LZ and WC transition probabilities in reverse direction."
  write(17,'(a,7x,a,7x,a,7x,a)') "Energy(cm-1)","b^2 diabatic","Landau-Zener","Weak Coupling"
  do i = 1, maxn+binGAP
    write (17,'(i6,3es24.5)') i, bbDia_rev(i), probLZ_rev(i), probWC_rev(i)
  end do
end if

end subroutine write_prob


! Velocity-averaged transition probabilities (VATP).

! Calculate LZ and WC velocity-averaged probabilities
! using Maxwell-Boltzmann and Kuki normalized distributions.

subroutine VATP()

 integer        :: i
 real(wp)       :: Temp, ksi
 real(wp)       :: probLZ_temp_MB
 real(wp)       :: probWC_temp_MB
 real(wp)       :: probLZ_temp_Kk
 real(wp)       :: probWC_temp_Kk

 real ( kind = 4 ), parameter   :: a = 0.0E+00
 real ( kind = 4 )                 abserr
 real ( kind = 4 ), parameter   :: b = 1.0E+00
 real ( kind = 4 )                 resabs
 real ( kind = 4 )                 resasc
 real ( kind = 4 )                 result_1
 real ( kind = 4 )                 result_2
 real ( kind = 4 )                 result_3
 real ( kind = 4 )                 result_4
 
Temp = zero
kbt = zero

probLZ_temp_MB = 0.0
probLZ_temp_Kk = 0.0
probWC_temp_MB = 0.0
probWC_temp_Kk = 0.0
ksi = 0.0

ksi = ((4*pi**2)*soc**2)*((2*redmass/(gradmean*grad))**(2.0/3.0)) 
   
write(66, 9)
write(66, 10)
  
do i = 0, Tpoints
  Temp = T1 + dble(i)*Tstep
  kbt = kb_hartK*dble(Temp)    
! Call # 1 : LZ_MB
  call qk61 (f05, a, b, result_1, abserr, resabs, resasc)
  probLZ_temp_MB = 1 - result_1*sqrt(2/(pi*kbt))
  probLZ_temp_MB = 2*probLZ_temp_MB - probLZ_temp_MB**2
! Call # 2 : LZ_Kk
  call qk61 (f06, a, b, result_2, abserr, resabs, resasc)
  probLZ_temp_Kk = 1 - result_2/kbt
  probLZ_temp_Kk = 2*probLZ_temp_Kk- probLZ_temp_Kk**2
! Call # 3 : WC_MB
  call qk61 (f07, a, b, result_3, abserr, resabs, resasc)
  probWC_temp_MB = (ksi*sqrt(2/(pi*kbt)))*result_3
! Call # 4 : WC_Kk
  call qk61 (f08, a, b, result_4, abserr, resabs, resasc)
  probWC_temp_Kk = (ksi/kbt)*result_4
     
  write(66,11) Temp, probLZ_temp_MB, probLZ_temp_Kk, &
               probWC_temp_MB, probWC_temp_Kk

11 format(3x,f7.1,3x,f6.4,4x,f6.4,4x,f6.4,4x,f6.4)
end do

9  format(/,/,11x,'Landau-Zener (LZ) and Weak Coupling (WC) double passage velocity-averaged', &
          /,27x,'probabilities of transition, p(T)', &
          /,11x,'Calculated using Maxwell-Boltzmann (MB) and Kuki (K)', &
          /,27x,'energy distributions')
10 format (/,6x,'T(K)',4x,'LZ MB',5x,'LZ K',6x,'WC MB',5x,'WC K',/)

end subroutine VATP

function f05 (x)
 real ( kind = 4 ) f05
 real ( kind = 4 ) x

! Note that gradmean is not used here. It is included as argument to
! make the qk61 call unified for LZ and WC formulas

f05 = exp((-2*pi*sqrt(redmass)*soc**2)/(dble(x)*grad)) &
                       *exp(-(dble(x)**2)/(2*kbt))

return
end function

function f06 (x)
 real (kind = 4 ) f06
 real (kind = 4 ) x

! Note that gradmean is not used here. It is included as argument to
! make the qk61 call unified for LZ and WC formulas

f06 = x*exp((-2*pi*sqrt(redmass)*soc**2)/(dble(x)*grad)) &
                         *exp(-(dble(x)**2)/(2*kbt))

return
end function

function f07 (x)
 real (kind = 4) f07
 real (kind = 4) x
 integer        :: NZ, IERR
 complex(wp)    :: airyvalue
 real(wp)       :: AII, airytemp

airytemp  = 0.0
airyvalue = (0.0,0.0) 
airytemp = -0.5*(x**2)*redmass*(((2*redmass*grad**2)/gradmean**4)**(1.0/3.0))
call ZAIRY(airytemp, 0.0, 0, 1, airyvalue, AII, NZ, IERR)
f07 = ((abs(airyvalue))**2.0)*exp(-(dble(x)**2)/(2*kbt))   
return

end function

function f08 (x)
 real (kind = 4) f08
 real (kind = 4) x
 integer        :: NZ, IERR
 complex(wp)    :: airyvalue
 real(wp)       :: AII, airytemp

airytemp  = 0.0
airyvalue = (0.0,0.0)
airytemp = -0.5*(x**2)*redmass*(((2*redmass*grad**2)/gradmean**4)**(1.0/3.0))
call ZAIRY(airytemp, 0.0, 0, 1, airyvalue, AII, NZ, IERR)
f08 = ((abs(airyvalue))**2.0)*(x**2)*exp(-(dble(x)**2)/(2*kbt))   
return

end function
end module prob_module
