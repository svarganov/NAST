module zn_param_module
 use precision_module, only : wp
 use input_module,     only : maxn, zpeR, zpeX
 use diag_module 
 use pes_module 
 
 implicit none
 public
 
contains

subroutine zn_param_sloped(up_origin,xcross,xpoint,up,down, &
                         aa,dd,bb,r0,binLower,binUpper)
 
 real(wp), intent(in)     :: up_origin, xcross 
 real(wp), dimension(:), intent(inout)  :: xpoint, up, down
 real(wp), intent(out)    :: aa, dd, r0 
 integer,          intent(out)    :: binLower, binUpper
 real(wp), dimension(:), intent(out) :: bb 
 integer          :: i, r0_index 
 integer          :: t1_index, t2_index, itarget,k, ispec 
 real(wp) :: Ex, E1t1, E1t2, E2t1, E2t2, E1r0, E2r0
 real(wp) :: bb_tmp
 real(wp) :: gap, x0, E1x0, E2x0, t1, t2
 real(wp) :: xval, Eup, Edown 
 real(wp), allocatable :: diffUpDown(:)
 
allocate(diffUpDown(maxn))
   
! Find approximate r0 value.

! Find gap between two states.
diffUpDown = abs(up-down)
r0_index = minloc(diffUpDown,1)
r0 = xpoint(r0_index)
E2r0 = up(r0_index)
E1r0 = down(r0_index)
gap  = E2r0-E1r0
Ex   = (E2r0+E1r0)/2.0d0  
 
! Find precise r0, E1r0, and E2r0 values.

! call find_extremum(r0,3, x0,E1x0,E2x0) 
! find_extremum subroutine works, but 
! highly depends on an initial guess; minimize_func is used instead.
k = 0
! Means that target function is a gap.
 call minimize_func(k,xcross, x0,E1x0,E2x0) 
 xpoint(r0_index) = x0
! Replace E1r0 value in down array.
 down(r0_index) = E1x0
! Replace E2r0 value in up array. 
 up(r0_index) = E2x0
! Replace r0 value in xpoint array.
 r0 = xpoint(r0_index)
! Update E1r0.
 E1r0 = down(r0_index)
! Update E2r0
 E2r0  =  up(r0_index)
 Ex=(E2r0+E1r0)/2.0d0
 
! Find approximate t1,t2, E1t1 and E2t2 values.

! Initial guess.
t2_index = r0_index
t2 = xpoint(t2_index)
do while (up(t2_index) .gt. Ex)
  t2_index = t2_index - 1
end do

! Initial guess. 
t1_index = r0_index
t1 = xpoint(t1_index)
do while ((down(t1_index) .lt. Ex) .and. (t1_index .lt. size(down)))
  t1_index = t1_index+1
enddo
 
t1   =  xpoint(t1_index)
t2   =  xpoint(t2_index)
E1t1 =    down(t1_index)
E2t1 =      up(t1_index)
E2t2 =      up(t2_index)
E1t2 =    down(t2_index)
 
! Find precise t1, E1t1 and E2t1 values.

k = 1
 call find_x(t1,Ex,k, xval,Eup,Edown) 
! Replace t1 value in xpoint array.
 xpoint(t1_index) = xval
! Replace E1t1 value in down array.
 down(t1_index) = Edown
! Replace t1 value in up array.
 up(t1_index) = Eup
! Update t1.
 t1   = xpoint(t1_index)
! Update E1t1.
 E1t1 =   down(t1_index)
! Update E2t1.
 E2t1 =   up(t1_index)
   
! Find precise t2 and E2t2 values.

k = 4
 call find_x(t2,Ex,k, xval,Eup,Edown) 
! Replace t2 value in xpoint array.
 xpoint(t2_index) = xval
! Replace E1t2 value in down array.
 down(t2_index) = Edown
! Replace E2t2 value in up array.
 up(t2_index) = Eup
! Update t2.
 t2   = xpoint(t2_index)
! Update E1t2.
 E1t2 =   down(t2_index)
! Update E2t2.
 E2t2 =   up(t2_index)
 
! Shift the origin to the reactant.

up   = up - up_origin
down = up - up_origin
E1r0 = E1r0 - up_origin
E2r0 = E2r0 - up_origin
E1t1 = E1t1 - up_origin
E2t1 = E2t1 - up_origin
E1t2 = E1t2 - up_origin
E2t2 = E2t2 - up_origin
 
! Determine energy regimes.

binLower = ceiling(E1r0*autocm)
binUpper = ceiling(E2r0*autocm)
 
! Compute aa and bb adiabatic parameters.

bb = 0.0d0
dd = (E2t1-E1t1)*(E2t2-E1t2)/(E2r0-E1r0)**2.0d0
aa = sqrt(dd-1.0d0)/(redmass*(t2-t1)**2.0d0*(E2r0-E1r0))
do i= 1, maxn
  bb_tmp = sqrt(dd-1.0d0)*(dble(i)/autocm-(E2r0+E1r0)/2.0d0+zpeR-zpeX)
  bb(i)  = 2.0d0*bb_tmp/(E2r0-E1r0)  
enddo
 
! Print critical points & adiabatic parameters.

!   print *,'bb(1) is              ', bb(1)
!   print *,'aa and dd are      ',aa, dd
!   print *,'Ex is              ',Ex
!   print *,'r0 and r0_index are',r0, r0_index
!   print *,'t2 and t2_index are',t2, t2_index
!   print *,'t1 and t1_index are',t1, t1_index
!   print *,'E1t1 and E1t2 are', E1t1, E1t2
!   print *,'E2t1 and E2t2 are', E2t1, E2t2
!   print *,'E2r0 and E1r0 are', E2r0, E1r0
 
! Print adiabatic energies.

open(21,file='adiabatic_energies.out')
write(21,*) "#bin     ","xpoint    ", "up      ","down      "
do i=1, maxn
  write (21,*) i, xpoint(i), up(i), down(i) 
enddo
close(21)

end subroutine zn_param_sloped

subroutine zn_param_peaked(origin,xcross, &
                        aa,bb,binLower,binUpper,down_origin)

 real(wp), intent(in)                :: origin, xcross
 real(wp), intent(out)               :: aa, down_origin 
 integer,  intent(out)               :: binLower, binUpper
 real(wp), dimension(:), intent(out) :: bb 
 integer          :: i, k
 real(wp) :: x0, deltaE, deltaR, der1, der2, d_up, d_down
 real(wp) :: E1x0, E2x0, Et, Eb, Rt, Rb, h, E1h, E2h, gam


! Determine the adiabatic origin.

! Origin is set to lower state by default. 
k = 1
 down_origin = target_func(origin,k)

! Find precise r0, E1r0, and E2r0 values.

! Target function now is an upper state.
k = 4
 call minimize_func(k,xcross, x0,E1x0,E2x0)
 Rb = x0
 Eb = E2x0
! Target function now is a lower state.
k = 1
 call minimize_func(k,xcross, x0,E1x0,E2x0)
 Rt = x0
 Et = E1x0
 h = (Rt + Rb)/two
 k = 1
 E1h =  target_func(h,k) 
 k = 4
 E2h =  target_func(h,k)
 deltaE = Eb-Et
 deltaR = Rb-Rt
 gam    = deltaE/(E2h - E1h)

! Shift the origin to the reactant.

Et = Et - down_origin
Eb = Eb - down_origin

! Determine energy regimes.

binLower = ceiling(Et*autocm)
binUpper = ceiling(Eb*autocm)

! Compute aa and bb adiabatic parameters.

if (abs(gam - one) .lt. 1.0d-8) then
  k = 1
    call derivatives(Rt,k, der1,der2)
    d_down = der2
  k = 4
    call derivatives(Rb,k, der1,der2)
    d_up=der2
    aa = (d_up-d_down)/(4.0d0*redmass*deltaE)
else
    aa = (1.0d0-gam**2.0d0)/(redmass*deltaR**2.0d0*deltaE)
endif

do i=1, maxn
  bb(i) = dble(i)/autocm - (Eb + Et)/two + zpeR - zpeX
  bb(i) = two*bb(i)/deltaE
enddo

! Print critical points and adiabatic parameters.

!   print *,'bb(1) is              ', bb(1)
!   print *,'aa and gamma are      ',aa, gam
!   print *,'deltE  is',deltaE
!   print *,'deltaR is',deltaR
!   print *,'Et and Eb are', Et, Eb
!   print *,'Rt and Rb are', Rt, Rb
!   print *,'down origin is', down_origin

end subroutine zn_param_peaked

!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

!----------------------GOLDEN SECTION-------------------------------
subroutine minimize_func(k,xcross, x0,E1x0,E2x0)
   real(wp) :: a, b, x1, x2 
   real(wp) :: fun1, fun2 
   real(wp) :: gratio, threshold, shift 
   real(wp) :: der1, der2   
   real(wp) :: der1a, der2a
   real(wp) :: der1b, der2b 
   integer          :: kk, i 
   integer,          intent(in)   :: k 
   real(wp), intent(in)   :: xcross 
   real(wp), intent(out)  :: x0, E1x0, E2x0
 
 !FUNCTION is subject to minimization
 !Determine the range where the smallest function value is expected
    shift=1.0d-10
    der1a=1.0d0 !set guess for derivative of FUNCTION at a
    der1b=1.0d0 !set guess for derivative of FUNCTION at b
    do while (der1a*der1b.gt.0) !till dirivatives have oposite signs 
       shift=shift*10.d0
       a = xcross-shift    !set guess for a
             call derivatives(a,k, der1,der2)
             der1a = der1
             der2a = der2
       b = xcross+shift    !set guess for b
             call derivatives(b,k, der1,der2)
             der1b = der1
             der2b = der2
    enddo
 
    gratio=(1.0d0+dsqrt(5.0d0))/2.0d0
    threshold=1.0d-16 

    i=0 !Number of the cycle 
    do while (ABS(a-b).gt.threshold)

      !if solution is not found within 100 cycles
      !threshold will be rised by one order of magnitude
      i=i+1
      if (mod(i,100).eq.0) then
       threshold=threshold*10.d0
!       print *,'Threshold for extremum has been raised to',threshold
      endif    

      !body of golden section cycle 
       x1 = b-(b-a)/gratio
       x2 = a+(b-a)/gratio
         fun1 = target_func(x1,k)
         fun2 = target_func(x2,k)

         !Minimize convex potential
         if (der2a.gt.0.0d0 .and. der2b.gt.0.0d0) then
             if (fun1.gt.fun2) then
                 a=x1
             else 
                 b=x2
             endif

         !Maximize concave potential
         elseif (der2a.lt.0.0d0 .and. der2b.lt.0.0d0) then
             if (fun1.lt.fun2) then
                 b=x1
             else 
                 a=x2
             endif
         endif
    enddo
    
    x0=a
   kk=4 !4 means that target function is an upper state
      E2x0 = target_func(x0,kk)
   kk=1 !1 means that target function is a lower state
      E1x0 = target_func(x0,kk) 
end subroutine minimize_func
!------------------------------------------------------------------

!------------------------------------------------------------------
function target_func(x,k) !returns target function
   integer, intent(in)           :: k 
   real(wp), intent(in)  :: x
   real(wp)              :: target_func
   real(wp)              :: e1, e2, eig(2) !eig(4)


   if     (k.eq.1) then   !lower state is a target function  
             e1=poly_high(x) 
             e2=poly_low(x)
!            call adiabatization(e1,e2,h12,h13,h14,eig)
             call adiabatization(e1,e2,eig)
             target_func=eig(1) 

   elseif (k.eq.4) then   !upper state is a target function
             e1=poly_high(x) 
             e2=poly_low(x)
!            call adiabatization(e1,e2,h12,h13,h14,eig)
             call adiabatization(e1,e2,eig)
!            target_func=eig(4)
             target_func=eig(2)

   elseif (k.eq.0) then   !gap is a target function
             e1=poly_high(x) 
             e2=poly_low(x)
!            call adiabatization(e1,e2,h12,h13,h14,eig)
             call adiabatization(e1,e2,eig)
!            target_func=eig(4)-eig(1)
             target_func=eig(2)-eig(1)
   endif

end function target_func
!------------------------------------------------------------------

!--------------SECANT METHOD---------------------------------------
subroutine find_x(x0,const,k, xval,Eup,Edown)
   integer, intent(in)           :: k 
   integer                       :: kk, i 
   real(wp), intent(in)  :: x0, const
   real(wp), intent(out) :: xval, Eup, Edown
   real(wp) :: delta, threshold, Ex
   real(wp) :: dx, x, old_x, new_x, f, old_f 
   
     dx=0.001d0 !small arbitrary shift to get new x
     threshold=1.0d-14 ! threshold for deviation in x
     delta=1.0d0 ! arbitrary value to initiate while cycle
     i=0 

       Ex=target_func(x0,k)  
       old_f=Ex-const
       old_x=x0
   
       x=x0+dx
       Ex=target_func(x,k)
       f=Ex-const
          do while (delta.gt.threshold)

          !if solution is not found within 100 cycles
          !threshold will be rised by one order of magnitude
          i=i+1
          if (mod(i,100).eq.0) then
           threshold=threshold*10.d0
!           print *,'Threshold for turning point & 
!                              has been raised to',threshold
          endif    

             new_x=old_x-(old_f*(x-old_x))/(f-old_f)
               old_f=f
               Ex=target_func(new_x,k)
               f=Ex-const
               old_x=x
               x=new_x
             delta=ABS(x-old_x)
          enddo 
      xval=x
   kk=4 !4 means that target function is an upper state
      Eup=target_func(xval,kk)
   kk=1 !1 means that target function is a lower state
      Edown=target_func(xval,kk) 
end subroutine find_x


!--------------NEWTON'S METHOD--------------------------------
subroutine find_extremum(xini,k, x0,E1x0,E2x0)
   real(wp), intent(in)  :: xini
   integer,          intent(in)  :: k
   real(wp), intent(out) :: x0, E1x0, E2x0
   real(wp) :: delta, threshold
   real(wp) :: x, old_x, gap 
   real(wp) :: der1, der2 
   integer          :: kk
   
 !    dx=0.001d0 !small arbitrary shift to get new x
     threshold=1.0d-14 ! threshold for deviation in x
     delta=1.0d0 ! arbitrary value to initiate while cycle
       x=xini
! print *,'threshold >>>>> ',threshold
          do while (delta.gt.threshold)
             call derivatives(x,k, der1,der2)
             old_x=x
             !f=f/(f)
             x=x-der1/der2
             gap=target_func(x,k)
 
             delta=ABS(x-old_x)
! print *,'REFINED GAP and r0 and delta are ',gap, x, delta
! print *,'!!! 1st and 2nd derivatives are ',der1, der2
          enddo 
      x0=x
       kk=4 ! 
       E2x0=target_func(x0,kk)
       kk=1 ! 
       E1x0=target_func(x0,kk)
end subroutine find_extremum

subroutine derivatives(x,k, der1,der2)
   real(wp) :: dx, h, f, y0, y1, y2, delta
   real(wp) :: y11, y22, y111, y222
   real(wp), intent(out) :: der1, der2
   real(wp), intent(in)  :: x
   integer,          intent(in)  :: k
           
   dx=0.0001d0
           h=x
               y0   =  target_func(h,k)

           h=x+dx
               y1   =  target_func(h,k)
           h=x-dx
               y2   =  target_func(h,k)
   
           h=x+2.0d0*dx
               y11  =  target_func(h,k)
           h=x-2.0d0*dx
               y22  =  target_func(h,k)
   
           h=x+3.0d0*dx
               y111 =  target_func(h,k)
           h=x-3.0d0*dx
               y222 =  target_func(h,k)
   
   der1=(y111-9.0d0*y11+45.0d0*y1-45.0d0* &
          y2+9.0d0*y22-y222)/(60.0d0*dx)
   der2=(2.0d0*y111-27.0d0*y11+270.0d0*y1-490.0d0*y0+ &
         270.0d0*y2-27.0d0*y22+2.0d0*y222)/(180.d0*dx*dx)

end subroutine derivatives

subroutine der_spin_diab(x,n, der)
   integer          :: i
   real(wp) :: dx, h 
   real(wp) :: p(3), m(3)
   real(wp), intent(out) :: der 
   real(wp), intent(in)  :: x
   integer,          intent(in)  :: n
           
   dx=0.0001d0
   if (n==2) then
        do i=1, 3 
           h = x-dx*dble(i)
           m(i)  =  poly_high(h)
        enddo

        do i=1, 3 
           h = x+dx*dble(i)
           p(i)  =  poly_high(h)
        enddo

   elseif (n==1) then
        do i=1, 3 
           h = x-dx*dble(i)
           m(i)  =  poly_low(h)
        enddo

        do i=1, 3 
           h = x+dx*dble(i)
           p(i)  =  poly_low(h)
        enddo
   endif
   
   der = (p(3)-9.0d0*p(2)+45.0d0*p(1) &
         -m(3)+9.0d0*m(2)-45.0d0*m(1))/(60.0d0*dx)
!           m(1)+9.0d0*m(2)-m(3))/(60.0d0*dx)

end subroutine der_spin_diab
!---------------------------------------------------------------------
end module zn_param_module
