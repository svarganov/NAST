module root_module
   use precision_module, only : wp
   use error_module
   use constants_module

   implicit none
   public

contains


function root(cf,limitL,limitR)
   implicit none
   real(wp) :: cf(5), limitL, limitR
   real(wp) :: rcf1(2),rcf2(3),rcf3(4), rcf4(5)
   real(wp) :: root

   if     (cf(1).eq.0.0d0 .and. & 
           cf(2).eq.0.0d0 .and. &
           cf(3).eq.0.0d0) then !linear  potentials
       rcf1=cf(4:5) !get rid of first three elements  
       rcf1=rcf1/rcf1(1) !get reduced polynomial 
       root=poly1_root(rcf1,limitL,limitR) !find root of cubic polynomial 

   elseif (cf(1).eq.0.0d0 .and. & 
           cf(2).eq.0.0d0) then !quadratic polynomial
       rcf2=cf(3:5) !get rid of first element  
       rcf2=rcf2/rcf2(1) !get reduced polynomial 
       root=poly2_root(rcf2,limitL,limitR) !find root of cubic polynomial 

   elseif (cf(1).eq.0.0d0) then !cubic polynomial 
       rcf3=cf(2:5) !get rid of the first element  
       rcf3=rcf3/rcf3(1) !get reduced polynomial 
       root=poly3_root(rcf3,limitL,limitR) !find root of cubic polynomial 

   else !quiartic polynomial 
       rcf4=cf/cf(1) ! get reduced polynomial 
       root=poly4_root(rcf4,limitL,limitR) !find root of quartic polynomial
   endif
end function root

!---------------------------------------------------------------------
!---------------------------------------------------------------------
function poly1_root(rcf,limitL,limitR)
   implicit none
   real(wp) :: rcf(2), limitL, limitR
   real(wp) :: x, poly1_root

     x=-rcf(2)
     if (x.gt.limitL .and. x.lt.limitR) then
         poly1_root=x
     else
        call error_message(8)!error code in parentheses   
     endif
end function poly1_root

!---------------------------------------------------------------------
!---------------------------------------------------------------------
function poly2_root(rcf,limitL,limitR)
   implicit none
   real(wp) :: rcf(3), limitL, limitR
   real(wp) :: b,c,dd,x1,x2, poly2_root
   integer          :: i

     b=rcf(2)
     c=rcf(3)
     dd=b**2.0d0-4.0d0*c !calculate determinant of reduced polynomial
         
     i=0 !number of roots
     if (dd.gt.0.0d0) then !two real root 
         x1=(-b+sqrt(dd))/2
         x2=(-b-sqrt(dd))/2
         if (x1.gt.limitL .and. x1.lt.limitR) then
            poly2_root=x1
            i=i+1
         elseif (x2.gt.limitL .and. x2.lt.limitR) then
            poly2_root=x1
            i=i+1
         endif 
     elseif (dd.eq.0.0d0) then
         x1=-b/2
         if (x1.gt.limitL .and. x1.lt.limitR) then
            poly2_root=x1
            i=i+1
         endif
     endif
  
   if (i.gt.1) then !multiple crossings in the region of interest  
        call error_message(9)   
   elseif (i.eq.0) then !no crossings in the region of interest  
        call error_message(8)   
   endif 
end function poly2_root

!---------------------------------------------------------------------
!---------------ROOTS ARE COMPUTED USING CARDANO METHOD---------------
!----see details at http://mathworld.wolfram.com/CubicFormula.html-----
function poly3_root(rcf,limitL,limitR)
   implicit none
   real(wp) :: rcf(4), limitL, limitR
   real(wp) :: poly3_root
   real(wp) :: y1, y2, y3, factor, phi
   real(wp) :: q, r, det, s, t
   complex(wp)   :: y2i, y3i
   integer          :: i
  
!   rcf=cf/cf(1) ! coefficient in front of x^3 is set to one 
   q=(3.0d0*rcf(3)-rcf(2)**2.0d0)/9.0d0
   r=(-2.0d0*rcf(2)**3.0d0+9.0d0*rcf(2)*rcf(3)-27.0d0*rcf(4))/54.0d0
   det=q**3.0d0+r**2.0d0
  
   if (det.gt.0.0d0) then !Only one real root
    s=r+sqrt(det)
     if (s.lt.0.0d0) then
        s=-((ABS(s))**(1.0d0/3.0d0))
     else
        s=s**(1.0d0/3.0d0)
     endif
    t=r-sqrt(det)
     if (t.lt.0.0d0) then
        t=-(ABS(t)**(1.0d0/3.0d0))
     else
        t=t**(1.0d0/3.0d0)
     endif
    y1=-rcf(2)/3.0d0+s+t
  
   !!!!!!!!!!!!!SKIP COMPEX ROOTS!!!!!!!!!!!!!!!!!!!!!!!!!
   !   re=-rcf(2)/3.0d0-(s+t)/2.0d0
   !   im=(sqrt(3.0d0)/2.0d0)*(s-t)
   !  y2i=cmplx(re,im)
   !  y3i=cmplx(re,-im)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
   elseif (det.eq.0.0d0) then ! three real roots (two are equal)
     if (s.lt.0.0d0) then
        s=-((abs(r))**(1.0d0/3.0d0))
     else
        s=r**(1.0d0/3.0d0)
     endif
    y1=-rcf(2)/3.0d0+s*2.0d0
    y2=-rcf(2)/3.0d0-s
    y3=y2

   elseif (det.lt.0.0d0) then ! three real unequal roots 
      phi=r/sqrt(abs(q**3.0d0))
      phi=1.0d0/(cos(phi))
      factor=2.0d0*sqrt(-q)
    y1=factor*cos(phi/3.0d0)-rcf(2)/3.0d0
    y2=factor*cos((phi+2.0d0*pi)/3.0d0)-rcf(2)/3.0d0
    y3=factor*cos((phi+4.0d0*pi)/3.0d0)-rcf(2)/3.0d0
   endif
   
   i=0 !number of roots
   if (det.gt.0.0d0) then !
     if (y1.gt.limitL .and. y1.lt.limitR) then
        poly3_root=y1
        i=i+1
     else
        call error_message(7)!error code in parentheses   
     endif 
   elseif (det.eq.0.0d0) then ! 
     if (y1.gt.limitL .and. y1.lt.limitR) then
        poly3_root=y1
        i=i+1
     elseif (y2.gt.limitL .and. y2.lt.limitR) then
        poly3_root=y2 ! y3 is equal to y2, so there's no need to track y3
        i=i+1
     endif 
   elseif (det.lt.0.0d0) then!  
     if (y1.gt.limitL .and. y1.lt.limitR) then
        poly3_root=y1
        i=i+1
     elseif (y2.gt.limitL .and. y2.lt.limitR) then
        poly3_root=y2
        i=i+1
     elseif (y3.gt.limitL .and. y3.lt.limitR) then
        poly3_root=y3
        i=i+1
     else
        call error_message(8)!error code in parentheses   
     endif 
   endif
  
   if (i.gt.1) then !you have multiple crossings in region of interest  
        call error_message(9)!error code in parentheses   
   endif 
end function poly3_root


!---------------------------------------------------------------------
!---------------ROOTS ARE COMPUTED USING FERRARI'S METHOD-------------
!----see details at https://en.wikipedia.org/wiki/Quartic_function----      
function poly4_root(rcf,limitL,limitR)
   implicit none
   real(wp) :: rcf(5), limitL, limitR
   real(wp) :: poly4_root
   real(wp) :: a,b,c,d,e
   real(wp) :: p,q, dd,dd0,dd1,T,tmp
   complex(wp)   :: S,x(4)
   integer          :: i,j


   a=rcf(1) !simplify nottion for polynomial coefficients
   b=rcf(2)
   c=rcf(3)
   d=rcf(4)
   e=rcf(5)

   dd=256d0*e**3.0d0-192d0*b*d*e**2.0d0-128d0*c**2.0d0*e**2.0d0
   dd=dd+144*c*d**2.0d0*e-27d0*d**4.0d0+144*b**2.0d0*c*e**2.0d0
   dd=dd-6.0d0*b**2.0d0*d**2.d0*e-80.0d0*b*c**2.0d0*d*e
   dd=dd+18.0d0*b*c*d**3.0d0+16.0d0*c**4.0d0*e-4.0d0*c**3.0d0*d**2.0d0
   dd=dd-27.0d0*b**4.0d0*e**2.0d0+18.0d0*b**3.0d0*c*d*e
   dd=dd-4.0d0*b**3.0d0*d**3.0d0-4.0d0*b**2.0d0*c**3.0d0*e
   dd=dd+b**2.0d0*c**2.0d0*d**2.0d0

   dd0=c**2.0d0-3.0d0*b*d+12.0d0*e
   dd1=2.0d0*c**3.0d0-9.0d0*b*c*d+27.0d0*b**2.0d0*e
   dd1=dd1+27.0d0*d**2.0d0-72.0d0*c*e
   p=(8.0d0*c-3.0d0*b**2.0d0)/8
   q=(b**3.0d0+8.0d0*d-4.0d0*b*c)/8
   
!   if (dd.lt.0) then

     T=(sqrt(-27.0d0*dd)+dd1)/2
       if (T.lt.0.0d0) then 
           T=-(abs(T)**(1.0d0/3.0d0))
       else 
           T=T**(1.0d0/3.0d0)
       endif
   
     S=-2.0d0/3.0d0*p+(T+dd0/T)/3.0d0
     S=0.5d0*sqrt(S)

!print *,'a,b,c,d,e',a,b,c,d,e
!print *,'dd,dd1,dd0',dd,dd1,dd0
!print *,'p,q,T,S',p,q,T,S

   x(1)=-4.0d0*S**2.0d0-2.0d0*p+q/S
   x(1)=-b/4.0d0-S+0.5d0*sqrt(x(1))

   x(2)=-4.0d0*S**2.0d0-2.0d0*p+q/S
   x(2)=-b/4.0d0-S-0.5d0*sqrt(x(2))

   x(3)=-4.0d0*S**2.0d0-2.0d0*p-q/S
   x(3)=-b/4.0d0+S+0.5d0*sqrt(x(3))

   x(4)=-4.0d0*S**2.0d0-2.0d0*p-q/S
   x(4)=-b/4.0d0+S-0.5d0*sqrt(x(4))

!print *,'ROOOOOTS',x(1),x(2),x(3),x(4)
!print *,'IMAG PART',aimag(x(3)),aimag(x(4))
   i=0
   do j=1,4
      if (aimag(x(j)).eq.0.0d0) then !real root
          tmp=real(x(j))
!print *,'tmp',tmp
          if (tmp.gt.limitL .and. tmp.lt.limitR) then
              poly4_root=tmp
              i=i+1
          endif
      endif 
   enddo
!print *,'iiiiiiii',i
!print *,'root',poly4_root
          
   if (i.gt.1) then !multiple crossings in the region of interest  
        call error_message(9)   
   elseif (i.eq.0) then !no crossings in the region of interest  
        call error_message(8)   
   endif 
end function poly4_root
!---------------------------------------------------------------------
end module root_module
