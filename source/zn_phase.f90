module zn_phase_module
 use precision_module, only : wp
 use diag_module
 use input_module, only : maxn, redmass, zpeR, zpeX
 use pes_module 
 use zn_param_module, only : target_func 

 implicit none
 public

contains

subroutine zn_zero_phase(aa,dd,bb, sigma0,delta0)
 real(wp), intent(in) :: aa, dd, bb(:)
 real(wp), dimension(:), intent(out) :: sigma0, delta0
 integer          :: i
 real(wp) :: fp, fm, fpc, fmc
 real(wp) :: gam1, gam2, gam2_old, prefactor
 real(wp) :: temp, temm, vx, bbnew, bx
  
gam1=0.9d0*sqrt(dd-1.0d0)
gam2=7.0d0*sqrt(dd)/16.0d0
prefactor=sqrt(2.0d0/aa)*pi/4.0d0
  
do i=1, maxn
  temp = (bb(i)+gam1)**2.0d0+gam2
  temp = sqrt(temp)
  temm = (bb(i)-gam1)**2.0d0+gam2
  temm = sqrt(temm)
  fp = sqrt(temp+(bb(i)+gam1)) + sqrt(temm+(bb(i)-gam1))
  fm = sqrt(temp-(bb(i)+gam1)) + sqrt(temm-(bb(i)-gam1))
  bx = bb(i)-0.9553d0
  gam2_old = gam2
  
! fmc = fm with redefined gam2.
  gam2 = 0.45d0*sqrt(dd)/(1.0d0+1.5d0*dexp(2.2d0*bx*abs(bx)**0.57d0))
  temp = (bb(i)+gam1)*(bb(i)+gam1)+gam2
  temp = sqrt(temp)
  temm = (bb(i)-gam1)*(bb(i)-gam1)+gam2
  temm = sqrt(temm)
  fmc = sqrt(temp-(bb(i)+gam1)) + sqrt(temm-(bb(i)-gam1))
  gam2 = gam2_old
  
! fpc = fp with redefined bb.
  bbnew = bb(i)-0.16d0*bx/sqrt(bb(i)**2.0d0+1.0d0)
  temp = (bbnew+gam1)*(bbnew+gam1)+gam2
  temp = sqrt(temp)
  temm = (bbnew-gam1)*(bbnew-gam1)+gam2
  temm = sqrt(temm)
  
  fpc = sqrt(temp+(bbnew+gam1)) + sqrt(temm+(bbnew-gam1))
  temp = fp**2.0d0+fm**2.0d0
  sigma0(i) = prefactor*fmc/temp
  delta0(i) = prefactor*fpc/temp
enddo

end subroutine zn_zero_phase
!---------------------------------------------------------------------

!---------------------------------------------------------------------
subroutine zn_phase_sloped(binLower,binUpper,shift,r0, &
           arrayT1,arrayT2,sigma0,delta0,  sigma,delta)

 integer              :: i, j, k 
 real(wp)             :: en, t1, t2
 integer,  intent(in) :: binLower, binUpper
 real(wp), intent(in) :: shift, r0
 real(wp), dimension(:), intent(in)  :: arrayT1,arrayT2
 real(wp), dimension(:), intent(in)  :: sigma0, delta0
 real(wp), dimension(:), intent(out) :: sigma, delta

sigma = zero
delta = zero

do i=1, maxn
  t1 = arrayT1(i)
  t2 = arrayT2(i)
  en = dble(i)/autocm+zpeR-zpeX 
  j  = nint(en*autocm)
     
! E < binLower.           
  if (j.lt.dble(binLower)) then
! Lower state.
    k = 1
    delta(i) = -integ(k,en,r0,t1,shift)
! Ipper state.
    k = 4
    delta(i) = delta(i)+integ(k,en,r0,t2,shift)
    delta(i) = delta(i)+delta0(i)
    sigma(i) = sigma0(i)

! binLower <= E < binUpper.           
  else if (j .ge. binLower .and. j .lt. binUpper) then
    k = 1
    sigma(i) = integ(k,en,r0,t1,shift)
    k=4
    delta(i) = integ(k,en,r0,t2,shift)
    delta(i) = delta(i)+delta0(i)
    sigma(i) = sigma(i)+sigma0(i)

! E > binUpper.            
  else if (j .ge. binUpper) then
    k = 1
    sigma(i) = integ(k,en,r0,t1,shift)
    k = 4
    sigma(i) = sigma(i)-integ(k,en,r0,t2,shift)
    delta(i) = delta0(i)
    sigma(i) = sigma(i)+sigma0(i)
  endif
enddo

end subroutine zn_phase_sloped

! Simpson's formula for integration.

function integ(k,en,r0,t,shift)
 integer          :: i, j, nstep, n
 real(wp) :: integ, dx, f, e_state, x, x0, coeff
 real(wp), intent(in) :: en, r0, t, shift
 integer,  intent(in) :: k

integ = zero
nstep = 1000
dx = abs(r0-t)/dble(nstep)

if (r0 .gt. t) then
  x0 = t
else
  x0 = r0
end if

do n=0, nstep
  x = x0 + dx*dble(n)
  e_state = target_func(x,k) - shift
  f = sqrt(abs(two*redmass*(en - e_state)))

! Set coefficient in Simpson's rule.          
  if (n == 0 .or. n == nstep) then
    coeff = one 
  else
    if (mod(n,2) == 1) then
       coeff = four
    else if (mod(n,2) == 0) then
       coeff = two
    end if
  end if
  integ = integ + coeff*f
end do  
integ = integ*dx/three

end function integ

subroutine zn_phase_peaked(binLower,binUpper,shift,aa,bb, &
            arT1tun,arT2tun,arT1ref,arT2ref,  sigma,delta)
 
 integer          :: i, j, k
 real(wp) :: en, t1, t2, temp, temm
 real(wp), intent(in) :: aa, shift
 integer,          intent(in) :: binLower, binUpper
 real(wp), dimension(:), intent(in) :: bb,arT1tun,arT2tun,arT1ref,arT2ref
 real(wp), dimension(:), intent(out) :: sigma, delta
 
sigma = zero
delta = zero

do i=1, maxn
  en = dble(i)/autocm + zpeR - zpeX 
  j = nint(en*autocm)
      
! E < binLower, delta in range from arT1tun to arT2tun.         
  if (j .lt. dble(binLower)) then

! Sigma.    
    temp = sqrt(one - one/bb(i)**two)
    temp = sqrt(6.0d0 + 10.0d0*temp)/(one + temp)
    sigma(i) = pi*temp/(16.0d0*sqrt(abs(bb(i))*aa))

! Delta.   
    t1 = arT1tun(i)
    t2 = arT2tun(i)
    k = 1
    delta(i) = integ(k,en,t1,t2,shift)

! binLower <= E < binUpper, sigma in range from r0 to t1.           
  else if (j .ge. binLower .and. j .le. binUpper) then

    temp = one + bb(i)
    temm = one - bb(i)
    sigma(i) =  temm*sqrt(five + three*bb(i))/sqrt(aa)
    sigma(i) = -sigma(i)*(0.057d0*temp**0.25d0 + one/three)
    delta(i) =  temp*sqrt(five - three*bb(i))/sqrt(aa)
    delta(i) =  delta(i)*(0.057d0*temm**0.25d0 + one/three)


! E > binUpper, sigma in range from arT1ref to arT2ref.           
  else if (j .ge. binUpper) then

! Delta.     
    temp = sqrt(one - one/bb(i)**two)
    temp = sqrt(6.0d0 + 10.0d0*temp)/(one + temp)
    delta(i) = pi*temp/(16.0d0*sqrt(bb(i)*aa))

! Sigma.    
    t1 = arT1ref(i-binUpper)
    t2 = arT2ref(i-binUpper)
    k = 4
    sigma(i) = integ(k,en,t1,t2,shift)
  end if
end do

!print *,'sigma 1',sigma(1)
!print *,'delta 1',delta(1)
!print *,'sigma 2400',sigma(2400)
!print *,'delta 2400',delta(2400)
!print *,'sigma 3100',sigma(3100)
!print *,'delta 3100',delta(3100)
!print *,'bb 1',bb(1)
!print *,'bb 2400',bb(2400)
!print *,'bb 3100',bb(3100)

end subroutine zn_phase_peaked
end module zn_phase_module
