module ORCA_read_module

use M_strings
use hess_data

implicit none

contains
!-----------------------------------------------
subroutine ORCA_read_atomic_data(file_id)

use get_constants

implicit none
integer, intent(in)                                    :: file_id
!-----------------------------------------------------------------
character(len=132)                                      :: line
character(len=132), allocatable                            :: array(:)
integer                                                    :: i
!-----------------------------------------------------------------

! For ORCA, the atomic data is written in the following format in the .hess file:
!
!     $atoms
!     N                               ('numbr' of atoms)
!     atom1   mass1   x1  y1  z1
!     atom2   mass2   x2  y2  z2
! ...
!     atomN   chargeN   xN  yN  zN
!
! The code below reads the atomic data from the ORCA .hess file.
! Note that usually programs print nuclear charge next to atom, like 'atom1  charge1 ...'
! instead of 'atom1  mass1 ...' as ORCA does. For other programs like GAMESS, Molpro, and
! QChem, we read atomic charges to search atom's mass using it's charge as index in the ma
! ss array stored in ams.f90.
! Therefore, for ORCA, we can read the masses from the .hess file and use those instead
! of the ams.f90 array.

do while (index(line, '$atoms') == 0) ! By ORCA's format, the coordinates are in bohr
 read(file_id,'(A)') line
end do

read(file_id,*) number_of_atoms 

allocate  (symbol (number_of_atoms) )
allocate  (charge (number_of_atoms) )
allocate  (mass (number_of_atoms) )
allocate  (coord (number_of_atoms,3) )

mass = 0.0d0
coord = 0.0d0
charge = 0.0d0

do i = 1, number_of_atoms
  read(file_id,'(A)') line
  call split(line, array)
  read(array(1),*) symbol(i)
  read(array(2),*) mass(i)
  read(array(3),*) coord(i,1)
  read(array(4),*) coord(i,2)
  read(array(5),*) coord(i,3)
end do

mass = mass*amu2au

end subroutine ORCA_read_atomic_data
!--------------------------------------------------------------------------
subroutine ORCA_read_grad_hess(file_id1,file_id2,grad,Hess)

implicit none
integer, intent(in)                                             :: file_id1, file_id2
double precision, dimension(:), allocatable, intent(out)        :: grad
double precision, dimension(:,:), allocatable, intent(out)      :: Hess
!----------------------------------------------------------------------------------
character(len=132)                                               :: line
character(len=132), allocatable                                  :: array(:)
integer :: i,j,k,n, column, number_of_five_column_blocks
!----------------------------------------------------------------------------------
allocate  (grad    (3*number_of_atoms) )
allocate  (Hess    (3*number_of_atoms,3*number_of_atoms) )
!----------------------------------------------------------------------------------
! Read gradient

rewind file_id2

! ORCA's .hess coordinates are already in bohr units, no need for conversion.
! Reading the cartesian gradient next. Needed to construct the Effective Hessian.

do while (index(line, 'CARTESIAN GRADIENT') == 0)
 read(file_id1,'(A)') line
end do

read(file_id1,'(A)') line
read(file_id1,'(A)') line

do i = 1, number_of_atoms
  read(file_id1,'(A)') line
  call split(line, array)
  read(array(4),*) grad(i*3-2)
  read(array(5),*) grad(i*3-1)
  read(array(6),*) grad(i*3)
end do

! Read Hessian
Hess = 0.0d0

do while (index(line, '$hessian') == 0)
 read(file_id2,'(A)') line
end do

read(file_id2,'(A)') line
read(file_id2,'(A)') line

number_of_five_column_blocks = int(ceiling(3*real(number_of_atoms)/5))

do i = 1, number_of_five_column_blocks

 if (i*5 <= 3*number_of_atoms) then
   column = 5
 else
   column = mod(3*number_of_atoms,5) ! 	mod(x,y)=x-y*INT(x/y)
 end if

 if (column == 5) then
   do j = 1, 3*number_of_atoms
     read(file_id2,'(A)') line
     call split(line,array)
     read(array(2),*) Hess(j,(i-1)*5 + 1)
     read(array(3),*) Hess(j,(i-1)*5 + 2)
     read(array(4),*) Hess(j,(i-1)*5 + 3)
     read(array(5),*) Hess(j,(i-1)*5 + 4)
     read(array(6),*) Hess(j,(i-1)*5 + 5)
   end do
 else if (column == 1) then
   do j = 1, 3*number_of_atoms
     read(file_id2,'(A)') line
     call split(line,array)
     read(array(2),*) Hess(j,(i-1)*5 + 1)
   end do
 else if (column == 2) then
   do j = 1, 3*number_of_atoms
     read(file_id2,'(A)') line
     call split(line,array)
     read(array(2),*) Hess(j,(i-1)*5 + 1)
     read(array(3),*) Hess(j,(i-1)*5 + 2)
   end do
 else if (column == 3) then
   do j = 1, 3*number_of_atoms
     read(file_id2,'(A)') line
     call split(line,array)
     read(array(2),*) Hess(j,(i-1)*5 + 1)
     read(array(3),*) Hess(j,(i-1)*5 + 2)
     read(array(4),*) Hess(j,(i-1)*5 + 3)
   end do
 else if (column == 4) then
   do j = 1, 3*number_of_atoms
     read(file_id2,'(A)') line
     call split(line,array)
     read(array(2),*) Hess(j,(i-1)*5 + 1)
     read(array(3),*) Hess(j,(i-1)*5 + 2)
     read(array(4),*) Hess(j,(i-1)*5 + 3)
     read(array(5),*) Hess(j,(i-1)*5 + 4)
   end do
 end if

 read(file_id2,'(A)') line

end do

end subroutine ORCA_read_grad_hess
!--------------------------------------------------
end module ORCA_read_module
