module rate_module
  use precision_module, only : wp
  use input_module
  use constants_module
  use error_module
  implicit none
  public
 
contains
!=============================================================================
!============== Marcus-Levich-Jortner rate constant k(T) =======================
!=============================================================================
subroutine MLJ_rate()
  implicit none
  character(len=33)                       :: printMLJ
  integer                                 :: i, j, k, j_fact, n
  real(wp)                                :: Temp, kbt, trateMLJ
  real(wp)                                :: freq_eff, S_eff, delta_freq
  real(wp)                                :: lambda_eff, lambda_int, lambda
  real(wp)                                :: a, b, c, d

Temp = zero
kbt = zero

write(66, 200)
write(67, 201)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!~~~~~~~~~~~~ Begiining of the T cycle ~~~~~~~~~~~!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

do i = 0, Tpoints

   Temp = T1 + dble(i)*Tstep

   kbt = kb_hartK*Temp ! Gives Hartree

   trateMLJ = zero

   ! Calculate internal reorganization energy lambda_int from low-frequency modes

   lambda_int = zero

   do j = 1, size(freR)
     if (freR(j) < trsh*kbT) then
       lambda_int = lambda_int + lambda_vib(j)
     end if
   end do

   lambda = lambda_int + lambda_ext

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !~~ Calcualte effective Huang-Rhys factor S_eff ~~!
   !~~ Calculate effective internal reorganization ~~!
   !~~~ energy lambda_eff and effective frequency ~~~! 
   !~~~~~~~ freq_eff for high-frequncy modes ~~~~~~~~!
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

   freq_eff   = zero
   S_eff      = zero
   lambda_eff = zero

   ! Check hv/kbT for every frequency at every T.
   ! Threshold chosen as hv/kbT < trsh for a classical mode
   ! and hv/kbT > trsh for a quantum mode.
   ! Modify if needed with trsh in the input file. 

   do j = 1, size(freR)
     if (freR(j) > trsh*kbT) then
       S_eff = S_eff + S(j)
       lambda_eff = lambda_eff + lambda_vib(j)
     end if
   end do

   if (S_eff.eq.0) then ! Classical limit, Marcus theory.

       a = 4*(pi**2)*(V**2)*(1/planck_hsec)
       b = 1/sqrt(4*pi*lambda*kbt)
       d = exp(-((gap+lambda)**2)/(4*lambda*kbt))

       trateMLJ = a*b*d

   else ! MLJ theory.

      freq_eff = lambda_eff/S_eff

      ! Calculate the max vibrational quantum number n for the effective vibrational mode.
      if (n_vib .eq. -1) then
        n = int(1/(exp(freq_eff/kbt)-1)) 
      else
        n = n_vib
      end if

      ! Calculate n! (factorial).

      do j = 0, n
     
        if (j .eq. zero) then
           j_fact = 1
        else 
           j_fact = 1
           do k = 1, j
             j_fact = j_fact*k
           end do 
        end if

      ! Calculate the MLJ rate constant at a given temperature.

        a = 4*(pi**2)*(V**2)*(1/planck_hsec)
        b = 1/sqrt(4*pi*lambda*kbt)
        c = exp(-S_eff)*((S_eff**j)/j_fact)
        d = exp(-((gap+lambda+freq_eff*j)**2)/(4*lambda*kbt))

        trateMLJ = trateMLJ + a*b*c*d

      end do

   end if

   write(printMLJ,'(es24.2)') trateMLJ
   write(66, 202) Temp, trim(adjustl(adjustr(printMLJ)))
   write(67, 203) Temp, freq_eff, S_eff, lambda_eff, lambda_int, lambda, n

end do 

write(67, 204) 

200 format(/,/,13x,'Canonical rate constant k(T), s-1.',/,/,&
            5x,'T(K)',3x,'Marcus-Levich-Jortner',/)
201 format(2x, "T, K", 7x, 'v_eff, a.u.', 6x, 'S_eff',&
                       7x, 'l_eff, a.u.', 4x, 'l_int, a.u.', 5x, 'l_tot, a.u.', 3x, 'n')
202 format(4x,f7.1,4x,a)
203 format(f7.1,5es15.3,4x,i2)
204 format(/,/,'v_eff is the frequency of the effective high-frequency vibrational mode.',/ &
               'S_eff is the Huang-Rhys factor for the effective high-frequency mode.',/ &
               'l_eff is the reorganization energy of the effective high-frequency mode.',/ &
               'l_int is the internal reorganization energy due to low-frequency modes.',/ &
               'l_tot is the total reorganization energy.',/ &
               'n is the vibrational quantum number of the effective high-frequency mode populated at a given temperature.',/,/ &

               'v_eff equal to zero means that all the vibrational modes were classical at a given temperature,',/ &
               'and the rate at this temperature was calculated with the classical Marcus theory.')

end subroutine MLJ_rate

end module rate_module
