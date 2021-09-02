program effective_hessian

!**************************************************!
!       This code implements calculation           !
!       of the Effective Hessian out of two        !
!       spin states Hessian matrices               !
!                       By                         !
!   Aleksandr O. Lykhin and Vsevolod D. Dergachev  !
!                                                  !
!         University of Nevada, Reno, 2020         !
!                                                  !
!**************************************************!

use M_strings
use write_data
use hess_driver

implicit none
integer :: point_pos
character(len=80) :: input_file, output_file, rc_file, &
                     nast_file, imag_file
double precision :: t1, t2

call get_command_argument(1,  input_file)
if (len_trim(input_file) == 0) then
 write (*,*) 'No input file'
 call exit
end if

point_pos = scan(trim(input_file),".", back = .true.)
if ( point_pos > 0 ) then
  output_file = input_file(1:point_pos)//'out'
  nast_file = input_file(1:point_pos)//'nast'
  imag_file = input_file(1:point_pos)//'imag'
  rc_file = input_file(1:point_pos)//'rc'
else
  output_file = input_file//'.out'
  nast_file = input_file//'.nast'
  imag_file = input_file//'.imag'
  rc_file = input_file//'.rc'
end if

open (unit = 11, file = input_file)
open (unit = 66, file = output_file)

call cpu_time(t1)

call write_header()
call driver(input_file,rc_file,nast_file,imag_file)
call cpu_time(t2)

write(66,'(/,1x,a,f11.2)') 'Total CPU time = ',t2-t1
write(66,'(1x,a)') 'Effective Hessian terminated now.'

close(66)
close(99)

end program effective_hessian
