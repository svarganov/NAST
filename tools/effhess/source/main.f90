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
!         Last updated: December 2022              !
!                                                  !
!**************************************************!

use M_strings   ! strings.f90
use write_data  ! write.f90 
use hess_driver ! driver.f90                                                      |
!                                                                                 |
implicit none   !                                                                 V
integer :: point_pos ! Position of . in the suffix of the like ../checkcomp/exam01.com
character(len=80) :: input_file, output_file, rc_file, &
                     nast_file, imag_file
double precision :: t1, t2

call get_command_argument(1,  input_file) ! saves the name of the input file in 'input_file'. Findes the name of the input file from ./effhess.x input_file.com
if (len_trim(input_file) == 0) then
 write (*,*) 'No input file' ! Try './effhess.x' one folder above (effhess)
 call exit
end if

point_pos = scan(trim(input_file),".", back = .true.) ! scans for "." in the name of input file
if ( point_pos > 0 ) then
  output_file = input_file(1:point_pos)//'out' ! If "." is found "out" is added. The output_file "input_file.out" is created
  nast_file = input_file(1:point_pos)//'nast' ! The nast_file "input_file.nast" is created
  imag_file = input_file(1:point_pos)//'imag' ! The imag_file "input_file.imag" is created
  rc_file = input_file(1:point_pos)//'rc' ! The rc_file "input_file.rc" is created
else
  output_file = input_file//'.out' ! If for some reason extension is missing (missing ".") it will be fixed
  nast_file = input_file//'.nast'
  imag_file = input_file//'.imag'
  rc_file = input_file//'.rc'
end if

open (unit = 11, file = input_file)
open (unit = 66, file = output_file)

call cpu_time(t1)

call write_header()
call driver(input_file,rc_file,nast_file,imag_file) ! The body of the effhess.x program is executed
call cpu_time(t2)

write(66,'(/,1x,a,f11.2)') 'Total CPU time = ',t2-t1
write(66,'(1x,a)') 'Effective Hessian terminated now.'

close(66)
close(99)

end program effective_hessian

! Summary of the module
! 1. Gets the name of the input file
! 2. Creates the input_file.out, input_file.nast, input_file.imag, input_file.rc files
! 3. Writes the header (call write_header())
! 4. Performs the main job (call driver)
! 5. Calculates the job time 
! 6. Job is finished
