module error_module

  implicit none
  public

  integer e

contains
!---------------------------------------------------------------------
subroutine error_message(e)

  integer e

  if (e .eq. 1) then
   write (*,*) "freR, lambda_vib, and S must be arrays of the same size."
  endif

  STOP

end subroutine error_message
        
!---------------------------------------------------------------------
end module error_module
