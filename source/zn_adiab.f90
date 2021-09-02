module zn_adiabatic_module
   use precision_module, only : wp
   use error_module
   use input_module
   use pes_module,   only : poly_high, poly_low 
   use root_module,  only : root 
   use diag_module 
   use zn_param_module 

   implicit none
   public
!   include 'constants'
 
contains
!------------------------------------------------------------------------
subroutine zn_path_sloped(origin,xpoint,up,down,up_origin,arrayT1,arrayT2)

 real(wp), intent(in)                :: origin
 real(wp), intent(out)               :: up_origin 
 real(wp), dimension(:), intent(out) :: arrayT1, arrayT2, xpoint, up, down
 real(wp) :: x, e_const, xval
 real(wp) :: Eup, Edown
 integer          :: i, j, k

! Calculate diabatic and adiabatic energies.

k = 4 ! Origin is set to upper state by default.
up_origin = target_func(origin,k)
do i=1, maxn ! Set guess point x 
  if (i .eq. 1) then 
    x = origin
  else
    x = arrayT2(i-1)
  endif
  e_const = up_origin+dble(i)/autocm+zpeR-zpeX
  j=0 ! number of cycles to find 'x' for given energy bin 'i'
  do
   call find_x(x,e_const,k, xval,Eup,Edown) 
   if ((origin.eq.limitL .and. xval.lt.limitL) .or. &
      (origin.eq.limitR .and. xval.gt.limitR)) then
     j = j + 1 
     x = origin + (-1.0d0)**dble(j)*(abs(limitR-limitL)/10.d0)  
     if (j .gt. 100) then
! Too many cycles.
       call error_message(12)   
     endif
   else
! xcval was found within region of interest.
     exit 
   endif 
  enddo 
  up(i) = Eup
  down(i) = Edown
  arrayT2(i) = xval
enddo

! Lower state.      
k = 1
do i=1, maxn
  if (i.eq.1) then 
    x = (limitL+limitR)/2.0d0
  else
    x = arrayT1(i-1)
  endif
  e_const = up_origin+dble(i)/autocm+zpeR-zpeX
  j=0
  do
   call find_x(x,e_const,k, xval,Eup,Edown) 
   if ((origin.eq.limitL .and. xval.lt.limitL) .or. &
      (origin.eq.limitR .and. xval.gt.limitR)) then
     j = j + 1 
     x = (limitL+limitR)/2.0d0 + (-1.0d0)**dble(j)*((limitR-limitL)/10.d0)  
     if (j .gt. 100) then
! Too many cycles.
       call error_message(12)   
     endif
   else
! xcval was found within region of interest.
     exit 
   endif 
 enddo 
 arrayT1(i) = xval
enddo
xpoint = arrayT2
         
! Lower potential approaches to the limit faster.

if ((origin.eq.limitL .and. arrayT1(maxn).gt.limitR) .or. &
   (origin.eq.limitR .and. arrayT1(maxn).lt.limitL)) then
  write(66,*) 'For a given maximum energy of',maxn,'cm-1, the last & 
         turning point',arrayT1(maxn),'is beyond the given limit'   
  write(66,*) 'Make sure that potentials are still reasonable there'  
endif
 
end subroutine zn_path_sloped
!-------------------------------------------------------------------------
        
subroutine zn_path_peaked(binLower,binUpper,down_origin, &
                          arT1tun,arT2tun,arT1ref,arT2ref)

 integer,          intent(in)  :: binLower, binUpper
 real(wp), intent(in)  :: down_origin
 real(wp), dimension(:), intent(out) :: arT1tun,arT2tun,arT1ref,arT2ref
 real(wp), allocatable :: up(:)
 real(wp) :: x, e_const, xval
 real(wp) :: Eup, Edown
 integer  :: i, j, k


! Calculate diabatic and adiabatic energies.
 
! Upper potential.
k = 4 
do i=binUpper+1, maxn
  if (i-binUpper .eq. 1) then 
    x = limitL
  else
    x = arT1ref(i-binUpper-1)
  end if
  e_const = down_origin + dble(i)/autocm + zpeR - zpeX
! Number of cycle to find 'x' for given energy bin 'i'.
  j = 0
  do
   call find_x(x,e_const,k, xval,Eup,Edown) 
   if (xval .lt. limitL .or. xval .gt. limitR) then
     j = j + 1 
     x = limitL + (-one)**dble(j)*((limitR - limitL)/10.d0)  
     if (j .gt. 100) then
! Too many cycles.
       call error_message(12)   
     end if
   else
! xval was found within region of interest.
     exit 
   end if 
  end do 
  arT1ref(i-binUpper) = xval
end do

! Upper state. 
do i=binUpper+1, maxn 
  if (i-binUpper .eq. 1) then 
    x = limitR
  else
    x = arT2ref(i-binUpper-1)
  endif
  e_const = down_origin+dble(i)/autocm + zpeR - zpeX
  j = 0
  do
   call find_x(x,e_const,k, xval,Eup,Edown) 
   if (xval .lt. limitL .or. xval .gt. limitR) then
     j = j + 1 
     x = limitR + (-one)**dble(j)*((limitR - limitL)/10.d0)  
     if (j .gt. 100) then
! Too many cycles.
       call error_message(12)
     end if
   else
! xval was found within region of interest.
     exit !xval was found within region of interest 
   end if 
  end do 
  arT2ref(i-binUpper) = xval
end do
         
! Lower potential.

k = 1 
! Lower state. Left.
do i=1, binLower-1
  if (i .eq. 1) then 
    x = limitL
  else
    x = arT1tun(i-1)
  endif
  e_const = down_origin + dble(i)/autocm + zpeR - zpeX
  j = 0
  do
   call find_x(x,e_const,k, xval,Eup,Edown) 
!  if ((origin.eq.limitL .and. xval.lt.limitL) .or. &
!     (origin.eq.limitR .and. xval.gt.limitR)) then
   if (xval .lt. limitL) then
     j = j + 1 
     x = limitL + (-one)**dble(j)*((limitR - limitL)/10.d0)  
     if (j .gt. 100) then
! Too many cycles
     call error_message(12)   
     end if
   else
! xval was found within region of interest.
     exit
   end if 
  end do 
  arT1tun(i) = xval
end do

! Lower state. Right.
do i=1, binLower-1
  if (i .eq. 1) then 
     x = limitR
  else
     x = arT2tun(i-1)
  end if
  e_const = down_origin + dble(i)/autocm + zpeR - zpeX
  j = 0
  do
   call find_x(x,e_const,k, xval,Eup,Edown) 
   if (xval .gt. limitR) then
     j = j + 1 
     x = limitR + (-one)**dble(j)*((limitR - limitL)/10.d0)  
     if (j .gt. 100) then
! Too many cycles.
       call error_message(12)!too many cycles   
     end if
   else
! xval was found within region of interest.
     exit 
   end if 
  end do 
  arT2tun(i) = xval
end do

open(17,file='arTXtun.out')
write(17,'(A,7x,A,7x,A)') "Energy(cm-1)","arT1tun","arT2tun"
do i=1, binLower-1
  write (17,*) i, arT1tun(i), arT2tun(i)
end do
close(17)

open(18,file='arTXref.out')
write(18,'(A,7x,A,7x,A)') "Energy(cm-1)","arT1ref","arT2ref"
do i=1, maxn-binUpper
  write (18,*) i, arT1ref(i), arT2ref(i)
end do
close(18)

end subroutine zn_path_peaked


!-------------------------------------------------------------------------
subroutine set_slope(hs0,hs1,hs2,hs3,hs4,ls0,ls1,ls2,ls3,ls4,limitL,limitR, &
                                             crosstype,xcross,mecp)
   real(wp), intent(in)  :: limitL, limitR
   real(wp), intent(in)  :: hs4,hs3,hs2,hs1,hs0, ls4,ls3,ls2,ls1,ls0
   character(6),     intent(out) :: crosstype
   real(wp), intent(out) :: xcross, mecp
   real(wp)              :: d1, d2, der
   integer                       :: n
   real(wp)              :: prod, gradmean, grad
   real(wp)              :: cf(5)

   cf(1)=hs4-ls4
   cf(2)=hs3-ls3
   cf(3)=hs2-ls2
   cf(4)=hs1-ls1
   cf(5)=hs0-ls0
   xcross=root(cf,limitL,limitR) !returns root of polynomial 
   mecp=poly_high(xcross) !locate absolute mecp barrier

   n=1 !lower state
    call der_spin_diab(xcross,n, der)
    d1=der
   n=2 !upper state
    call der_spin_diab(xcross,n, der)
    d2=der

   prod=d1*d2
   gradmean=sqrt(abs(prod))
   grad=abs(d1-d2)

     if (prod.gt.0.0d0) then
        write (66,'(/,10x,a)') "Sloped intersection of PESs"
        crosstype='sloped'
     else if (prod.lt.0.0d0) then
        write (66,'(/,10x,a)') "Peaked intersection of PESs"
        crosstype='peaked'
     else
        call error_message(4)   
     endif
end subroutine set_slope

!------------------------------------------------------------------------
subroutine set_origin(limitL,limitR,  origin)
   real(wp)              :: eL_hs, eL_ls, eR_hs, eR_ls
   real(wp)              :: ar(4), EreacDiab
   integer                       :: INFO
   real(wp), intent(in)  :: limitL, limitR
   real(wp), intent(out) :: origin
  
   eL_hs = poly_high (limitL)
   eL_ls = poly_low  (limitL)
   eR_hs = poly_high (limitR)
   eR_ls = poly_low  (limitR)

   ar(1) = eL_hs
   ar(2) = eL_ls
   ar(3) = eR_hs
   ar(4) = eR_ls

   CALL DLASRT('D',4,ar,INFO) !sort in decreasing order
   EreacDiab = ar(3) ! 
   
   if (EreacDiab.eq.eL_hs) then
      origin = limitL 
      write (66,'(10x,a)') 'Origin is High-Spin state at limitL'
   elseif (EreacDiab.eq.eL_ls) then
      origin = limitL 
      write (66,'(10x,a)') 'Origin is Low-Spin state at limitL'
   elseif (EreacDiab.eq.eR_hs) then
      origin = limitR 
      write (66,'(10x,a)') 'Origin is High-Spin state at limitR'
   elseif (EreacDiab.eq.eR_ls) then
      origin = limitR 
      write (66,'(10x,a)') 'Origin is Low-Spin state at limitR'
   endif

    write (66,'(10x,a,f5.2,a)') 'Origin is at ',origin,' bohr'

end subroutine set_origin
!------------------------------------------------------------------------


!------------------------------------------------------------------------
end module zn_adiabatic_module
