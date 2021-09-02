program nast
  use input_module, only : read_file
  use driver_module
  use write_nast_data 
  implicit none

  double precision  :: t_1, t_2
  character(len=80) :: input_file
  
  open(unit=66, file = 'nast.out')
  call get_command_argument(1,  input_file)
  if (len_trim(input_file) == 0) then
     write (66,*) 'No input file'
     call exit
  end if
 
  write(66,'(A,A)') 'Input file: ', input_file
  open (unit=11, file=input_file )

  call cpu_time(t_1)

  call write_header()
  call read_file()
  call write_options()
  call run()

  call cpu_time(t_2)
  write(66,'(/,1X,A,F7.2)') 'Total CPU time = ',t_2-t_1
  write(66,'(1X,A)') 'NAST terminated now.'

  close (11)
  close (66)

end program nast
