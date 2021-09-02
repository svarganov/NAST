module GAMESS_read_module

use M_strings
use get_constants
implicit none

contains
!------------------
subroutine GAMESS_read_atomic_data(file_id,npoint,numbr,geom,irc_eng)

implicit none
integer, intent(in)                                          :: file_id
integer, intent(out)                                         :: npoint, numbr
double precision, allocatable, intent(out)                   :: geom(:,:,:), irc_eng(:)
!-----------------------------------------------------------------
character(len=99)                                            :: line
character(len=3)                                             :: point_g
character(len=99), allocatable                               :: array(:)
integer                                                      :: i, end_of_file, num
double precision, allocatable                                :: geom_tmp(:,:,:), eng_tmp(:)
!-----------------------------------------------------------------

read(file_id,'(A)') line                          ! Check from first line

do while (index(line, 'TOTAL NUMBER OF ATOMS') == 0)
 read(file_id,'(A)') line
end do                                            ! Now on the '...TOTAL NUMBER...' line

call split(line, array)
read(array(6),*) numbr
rewind(file_id)

allocate ( geom_tmp (300,numbr,3)  )                   ! 150 is to excede real number of IRC points (npoint)
allocate ( eng_tmp (300) )

npoint = 0

do
 read(file_id,'(A)',iostat=end_of_file) line
 if (end_of_file > 0) then
   print*, 'Something went wrong when reading from file'
   exit
 else if (end_of_file < 0) then
   exit
 else
   if (index(line, ' ON THE REACTION PATH') .ne. 0) then
     npoint = npoint + 1
     call split (line, array)
     read (array(3),*) num
     if ( num == 0 ) then
       read(file_id,'(A)') line
       read(file_id,'(A)') line
       read(file_id,'(A)') line
       read(file_id,'(A)') line
       read(file_id,'(A)') line
       read(file_id,'(A)') line
       do i = 1, numbr
         read(file_id,'(A)') line
         call split (line, array)
         read(array(3),*) geom_tmp (npoint,i,1)
         read(array(4),*) geom_tmp (npoint,i,2)
         read(array(5),*) geom_tmp (npoint,i,3)
       end do
       do while (index(line, 'ENERGY IS ') == 0)
         read(file_id,'(A)') line
       end do
        call split (line, array)
        read (array(5),*) eng_tmp(npoint)
     else
       read(file_id,'(A)') line
       read(file_id,'(A)') line
       read(file_id,'(A)') line
       call split (line, array)
       read(array(4),*) eng_tmp(npoint)
       read(file_id,'(A)') line
       read(file_id,'(A)') line
       read(file_id,'(A)') line
       read(file_id,'(A)') line
       read(file_id,'(A)') line
       do i = 1, numbr
         read(file_id,'(A)') line
         call split (line, array)
         read(array(3),*) geom_tmp (npoint,i,1)
         read(array(4),*) geom_tmp (npoint,i,2)
         read(array(5),*) geom_tmp (npoint,i,3)
       end do
     end if 
   end if
 end if
end do

allocate ( irc_eng(npoint), geom(npoint,numbr,3) )

irc_eng = eng_tmp(1:npoint)

do i = 1, npoint
  geom(i,1:numbr,1:3) = geom_tmp(i,1:numbr,1:3)
end do

deallocate ( eng_tmp, geom_tmp )

end subroutine GAMESS_read_atomic_data
!-------------------------------------------------------------------------------------
end module GAMESS_read_module
