module zn_prob_module
   use precision_module, only : wp
   use constants_module
   use input_module,     only : maxn

   implicit none
   public

contains
!---------------------------------------------------------------------
subroutine zn_prob_sloped(aa,bb,sigma,delta,probZN)
 implicit none
 integer          :: i, k
 real(wp) :: corr, arg, phi, g1, g2, UIm, URe, p, temp
 real(wp) :: x, part_i, part_r, arggamma, phase
 real(wp), intent(in)  :: aa
 real(wp), dimension(:),intent(in)   :: bb, sigma, delta
 real(wp), dimension(:), intent(out) :: probZN
 integer  :: jerr

do i=1,maxn
  if (bb(i) .lt. zero) then
  
    g2 = three*sigma(i)*log(1.2d0 + aa)
    g2 = g2/(delta(i)*pi) - one/aa
    g1 = 1.8d0*(aa**0.23d0)*exp(-delta(i))
    temp = bfunc(sigma(i)/pi)
  
    URe = sqrt(temp)*exp(delta(i)) - g1*sin(sigma(i))**two
    URe = cos(sigma(i))*URe
  
    UIm = -(g1*sin(sigma(i))*cos(sigma(i)))**two
    UIm = UIm*exp(-two*delta(i))/temp
    UIm = UIm + temp*exp(two*delta(i))
    UIm = UIm + two*g1*(cos(sigma(i))**two) - g2
    UIm = sin(sigma(i))*sqrt(UIm)
  
    phase = atan(UIm/URe)
    p = one + temp*exp(two*delta(i))
    p = p - g2*(sin(sigma(i))**two)
    p = one/p
    probZN(i) = four*p*(one - p)*(sin(phase)**two)
  
  else if (bb(i) .ge. zero) then
  
    corr = five*(aa**0.25d0)*(10.0d0**(-sigma(i)))
    corr = corr/(0.8d0+(aa**0.25d0))
    arg = delta(i)*(one + corr)/pi 
    x = zero
    k = 1
    part_r = zero
    part_i = zero
    call cgamma(x,arg,k,part_r,part_i)
    arggamma = atan(part_i/part_r)
    phi = -arg + arg*log(arg) - arggamma - 0.25d0*pi
    corr = one + sqrt(one + (bb(i)**(-two))*(0.4d0*aa + 0.7d0))
    corr = sqrt(two/corr)
    p = exp(-0.25d0*pi*corr/sqrt(bb(i)*aa))
    probZN(i) = four*p*(one - p)
    probZN(i) = probZN(i)*sin(phi + sigma(i))**two
  
  endif
enddo

end subroutine zn_prob_sloped
!------------------------------------------------------------------------
subroutine zn_prob_peaked(binLower,binUpper,aa,bb,sigma,delta,probZN)
implicit none
integer             :: i, k, n
integer, allocatable :: m(:)
real(wp) :: arg, phi, phibar, h1, h2, h3, h4, p, temp, corr
real(wp) :: sigmac, w, t, temp2
real(wp) :: x, part_i, part_r, arggamma, phase
integer,  intent(in)  :: binLower, binUpper
real(wp), intent(in)  :: aa 
real(wp), intent(in)  :: bb(:), sigma(:), delta(:)
real(wp), dimension(:),intent(out) :: probZN
real(wp), allocatable :: grid_small_aa_column(:), integral_matrix(:,:)
  
probZN = 0.0d0

do i=1,maxn
   if(i.lt.binLower) then
  
       sigmac=1.0d0-0.32d0*(10.0d0**(-2.0d0/aa))*exp(-delta(i))
       sigmac=sigmac*sigma(i)
       temp=bfunc(sigmac/pi)*exp(-2.0d0*delta(i))
       probZN(i)=1.0d0+(0.5d0*sqrt(aa)/(1.0d0+sqrt(aa)))*temp
       probZN(i)=temp/(probZN(i)**2.0d0+temp)
 
   elseif (i.ge.binLower .and. i.le.binUpper) then
 
       h4=0.61d0*sqrt(2.0d0+bb(i))
       h3=(sqrt(aa)-3.0d0*bb(i))*sqrt(1.23d0+bb(i))
       h3=h3/(sqrt(aa)+3.0d0)
       h2=(1.0d0+bb(i))**(1.2d0-0.4d0*bb(i))
       h2=1.0d0+0.38d0*h2/aa
       call omega_solver(aa,bb(i),h2,h3,h4,probZN(i))

 
   elseif (i.gt.binUpper) then
 
       arg=delta(i)/pi 
       x=0.0d0
       k=1
       part_r = 0.0d0
       part_i = 0.0d0
       call cgamma(x,arg,k,part_r,part_i)
       arggamma=atan(part_i/part_r)
       phi=-arg+arg*log(arg)-arggamma-0.25d0*pi
       h1=0.23d0*(aa**0.25d0)*(40.0d0**(-sigma(i)))
       h1=h1/(aa**0.25d0+0.75d0)
       phibar=phi+h1
       p=bb(i)**(-2.0d0)*(0.72d0-0.62d0*(aa**(1.43d0/2.0d0)))
       p=1.0d0+sqrt(1.0d0-p)
       p=sqrt(2.0d0/p)
       p=exp(-0.25d0*pi*p/sqrt(bb(i)*aa))
       probZN(i)=4.0d0*(cos(sigma(i)-phibar)**2.0d0)
       probZN(i)=probZN(i)/(probZN(i)+p*p/(1.0d0-p))
 
   endif
 enddo

if ( sqrt(aa) .lt. one ) then 
 allocate ( m((binUpper-1) - (binLower+1)) )
 allocate ( grid_small_aa_column(101), integral_matrix(100,101) )
 open( 12, file="omega_grid_small_aa" )
 read( 12,*) integral_matrix
 close(12)
 n=int( sqrt(aa)/0.01 )
 m=int( 50*bb(binLower+1:binUpper-1) + 1)
 grid_small_aa_column = integral_matrix(n,:)
 deallocate (integral_matrix)
 
 do i = binLower+1, binUpper-1
  corr = 0.0d0
  corr = (h2/(aa**(1.0d0/3.0d0)))*grid_small_aa_column(m(i))
  corr = (corr**2)/((corr**2)+1)
  probZN(i) = probZN(i) + corr
 end do
 
end if

end subroutine zn_prob_peaked
!-------------------------------------------------------------------
subroutine omega_solver(aa,bb,h2,h3,h4,prob)
! ********
! *
! Below is the heuristic algorithm how to solve a highly
! oscillatory integral using the QUADPACK integration library of Fortran.
! Wolfram Mathematica solution of this integral has been used as a reference
! to estimate an error of the below Fortran alogirthm. The error has been
! estimated to be of the order of 10-4.
! *
! The original integral is W = (pre-factor)*INT(0,oo,t,aa,bb)dt
! What we do is divide W in two integrals, I1 + I2: INT(0,1) + INT(1,oo)
! Then, we transform I2: INT(1,oo) -> INT(0,1) by changing the integration variable.
! Function f01 implements I2 integrand, while f02 - I1. 
! *
! ********
implicit none
real(wp), intent(in) :: aa,bb,h2,h3,h4
real(wp), intent(out) :: prob
real ( kind = 4 ), parameter :: a1 = 0.02e+00
real ( kind = 4 ), parameter :: b1 = 0.1e+00
real ( kind = 4 ) abserr, resabs, resasc
integer ( kind = 4 ) ier
integer ( kind = 4 ) neval
real ( kind = 4 ) result, res, tmp
!real ( kind = 4 ), external :: f01, f02
real ( kind = 4 ) dq, q1, q2, dd, n1, n2
integer :: i, j
integer, parameter :: nq = 1000000, nn = 8

res = 0.0e+00
dd = (b1-a1)/nn
! *
! Get I2 ( 0.02 , 0.1 )
! *
do j = 0, nn-1

  n1 = a1 + dble(j)*dd
  n2 = a1 + dble(j+1)*dd
  dq = (n2 - n1)/nq
  tmp = 0.0e+00
  result = 0.0e+00
  
  do i=0, nq-1
     q1 = n1 + dble(i)*dq
     q2 = n1 + dble(i+1)*dq
     call qk61 (f01, q1, q2, result, abserr, resabs, resasc)
     tmp = tmp + result
  end do

  res = res + tmp
end do
! *
! Get I2 ( 0.1 , 0.2 )
! *     
dd = 0.0e+00
tmp = 0.0e+00
result = 0.0e+00

dd = (0.2e+00 - 0.1e+00)/100000
n1 = 0.0e+00
n2 = 0.0e+00

do i = 0, 99999
  n1 = 0.1e+00 + dble(i)*dd
  n2 = 0.1e+00 + dble(i+1)*dd
  call qk61 (f01, n1, n2, result, abserr, resabs, resasc)
  tmp = tmp + result
end do

res = res + tmp
result = 0.0e+00
tmp = 0.0e+00
! *
! Get I2 ( 0.2 , 0.6 ) and I2 ( 0.6 , 1.0 )
! *
call qk61 (f01, 0.2e+00, 0.6e+00, result, abserr, resabs, resasc)
res = res + result
result=0.0e+00
call qk61 (f01, 0.6e+00, 1.0e+00, result, abserr, resabs, resasc)
res = res + result
result = 0.0e+00
! *
! Divide I1 ( 0.0, 1.0 ) in two parts to reduce the error
! *
call qk61 (f02, 0.0e+00, 0.5e+00, result, abserr, resabs, resasc)
tmp = result
result = 0.0e+00
call qk61 (f02, 0.5e+00, 1.0e+00, result, abserr, resabs, resasc)
tmp = tmp + result
res = res + tmp
result = 0.0e+00
tmp = 0.0e+00

prob = (h2/(aa**(1.0d0/3.0d0)))*res
prob = (prob**2)/((prob**2)+1)
contains
!-------------------------------------------------------------------
real ( kind = 4 ) function f01(x)
implicit none
real ( kind = 4 ) x

f01 = cos(0.333333/(x**3) - (bb/(aa**(1.0d0/3.0d0)))/x - &
         (h3/(aa**(1.0d0/3.0d0)))/(h4*x+(sqrt(aa)**(1.0d0/3.0d0))))/(x**2)

end function
!-------------------------------------------------------------------
real ( kind = 4 ) function f02(x)
implicit none
real ( kind = 4 ) x

f02 = cos(0.333333*(x**3) - (bb/(aa**(1.0d0/3.0d0)))*x - &
          (h3/(aa**(1.0d0/3.0d0)))*x/(h4+(sqrt(aa)**(1.0d0/3.0d0))*x))

end function
!-------------------------------------------------------------------
end subroutine omega_solver
!-------------------------------------------------------------------
function bfunc(x)
   implicit none
   real(wp) :: bfunc,x,g,temp
   
   call vdtgamma(1,x,g) !MKL Gamma function for real values only
   temp=2.d0*pi*x**(2.0d0*x)*exp(-2.0d0*x)
   bfunc=temp/(x*g**2.0d0)
end function
!---------------------------------------------------------------------
end module zn_prob_module
