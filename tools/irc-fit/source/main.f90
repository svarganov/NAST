program irc_fit

use irc_driver
use M_strings
use write_irc_data
use irc_data
use GAMESS_read_module

implicit none
!---------------------------------------------------------------------------------------------------------
character(len=80)             :: input_file, output_file, nast_file
integer                       :: point_pos
double precision              :: t1, t2
!---------------------------------------------------------------------------------------------------------
! Start the program

call cpu_time(t1)

call get_command_argument(1,  input_file)
if (len_trim(input_file) == 0) then
 write (66,'(/,A)') 'No input file'
 call exit
end if

point_pos = scan(trim(input_file),".", back = .true.)
if ( point_pos > 0 ) then
  output_file = input_file(1:point_pos)//'out'
  nast_file = input_file(1:point_pos)//'nast'
else
  output_file = input_file//'.out'
  nast_file = input_file//'.nast'
end if

open (unit=66, file = output_file)

call write_header(66)

open (unit = 99, file = input_file)

call driver_irc(input_file,nast_file)

call cpu_time(t2)

write(66,'(/,1X,A,F6.2)') 'Total CPU time = ',t2-t1
write(66,'(/,a)') 'Fly with us again.'
close(11)
close(12)
close(66)

end program irc_fit
