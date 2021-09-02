module pes_module 
   use precision_module, only : wp
   use input_module
   implicit none
   public

contains
!-------------------------------------------------------------------
function poly_high(x)
   real(wp) :: x, poly_high

   poly_high=hs4*x**4.0d0+hs3*x**3.0d0+hs2*x**2.0d0+hs1*x+hs0      
!   poly_high=hs(1)*x**3.0d0+hs(2)*x**2.0d0+hs(3)*x+hs(4)      
end function poly_high
!-------------------------------------------------------------------

function poly_low(x)
   real(wp) :: x, poly_low 

   poly_low=ls4*x**4.0d0+ls3*x**3.0d0+ls2*x**2.0d0+ls1*x+ls0      
!   poly_low=ls(1)*x**3.0d0+ls(2)*x**2.0d0+ls(3)*x+ls(4)      
end function poly_low

!-------------------------------------------------------------------

function f_inertX(coef,k)
   real(wp) ::  f_inertX, coef(3) 
   integer          ::  k 

   f_inertX=coef(1)*dble(k)**2.0d0+coef(2)*dble(k)+coef(3)      
   f_inertX=f_inertX*amutoau      
!PRINT *,'HERE',f_inertX
end function f_inertX

!-------------------------------------------------------------------
end module pes_module
