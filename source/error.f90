module error_module

  implicit none
  public

  integer e

contains
!---------------------------------------------------------------------
subroutine error_message(e)

  integer e

  if (e .eq. 1) then
   write (*,*) "Polynomial coefficient are not recognized. Check the input file"
  
  elseif (e .eq. 2) then
   write (*,*) "SOC matrix elements are not recognized. Check the input file"
  
  elseif (e .eq. 3) then
   write (*,*) "Error in adiabatization"

  elseif (e .eq. 4) then
   write (*,*) "One of the gradients at the crossing point is equal to zero"

  elseif (e .eq. 5) then
   write (*,*) "Left limit is equal to or greater than the right limit"

  elseif (e .eq. 6) then
   write (*,*) "Something wrong with crosstype"

  elseif (e .eq. 7) then
   write (*,*) "Crossing point of two states is outside the giving limits"
   write (*,*) "Only one real crossing point, while two others are complex"

  elseif (e .eq. 8) then
   write (*,*) "Crossing point(s) are outside the given limits"

  elseif (e .eq. 9) then
   write (*,*) "You have multiple crossings within the given limits"

  elseif (e .eq. 10) then
   write (*,*) "Reverse rate input is not recognized. Check the input file"

  elseif (e .eq. 11) then
   write (*,*) "WC probability is greater than 100%"

  elseif (e .eq. 12) then
   write (*,*) "Couldn't determine reaction coordinate for a bin"

  elseif (e .eq. 13) then
   write (*,*) "Only ZPE values equal to 0, 1 or 2 are accepted (ZPE=0 or ZPE=1 or ZPE=2)"

  elseif (e .eq. 14) then
   write (*,*) "Neither 'classic' or 'quantum' rotational model is chosen"

  elseif (e .eq. 15) then
   write (*,*) "TP data is not recognized. Check the input file"

  elseif (e .eq. 16) then
   write (*,*) "HSO components are not recognized. Check the input file"

  elseif (e .eq. 17) then
   write (*,*) "The T_final value in the input file is set negative"

  endif

  STOP

end subroutine error_message
        
!---------------------------------------------------------------------
end module error_module
