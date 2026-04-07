module read_geometry ! Created 07/14/2025, added coord_diff on 07/15/2025

use M_strings
use hess_data

implicit none

contains
!-----------------------------------------------
subroutine read_atomic_data(file_id1, file_id2)

use get_constants
use get_atomic_masses

implicit none
integer, intent(in)                                    :: file_id1, file_id2 ! Commented 07/14/2025
!-----------------------------------------------------------------
character(len=132)                                     :: line
character(len=132), allocatable                        :: array(:)
integer                                                :: i
!-----------------------------------------------------------------

! The format of the 'initial_geometry' file is as follows:
!
!     N                               ('number' of atoms)
!     atom1   charge1   x1  y1  z1    (x, y, z coordinates in Angstroms)
!     atom2   charge2   x2  y2  z2
! ...
!     atomN   chargeN   xN  yN  zN
!
! The code below reads in the atomic data from the input geometry file in this format.
! Reading the first line gives the number of atoms in a molecule. The each from the rest of the
! lines is split to five tokens. The second token "chargeN" is used to find the corresponding
! atomic mass for an element using the ams.f90 array.

! Read full atomic data for the initial geometry

read(file_id1,*) number_of_atoms 

allocate  (symbol (number_of_atoms) )
allocate  (charge (number_of_atoms) )
allocate  (mass (number_of_atoms) )
allocate  (coord (number_of_atoms,3) )

mass = 0.0d0
coord = 0.0d0
charge = 0.0d0

do i = 1, number_of_atoms
  read(file_id1,'(A)') line
  call split(line, array)
  read(array(1),*) symbol(i)
  read(array(2),*) charge(i)
  read(array(3),*) coord(i,1)
  read(array(4),*) coord(i,2)
  read(array(5),*) coord(i,3)
end do

do i = 1, number_of_atoms
  mass(i) = ams(int(charge(i)))
end do

! Read coordinates of the final geometry

allocate  (coord_final (number_of_atoms,3) )

coord_final = 0.0d0

do i = 1, number_of_atoms
  read(file_id2,'(A)') line
  call split(line, array)
  read(array(2),*) coord_final(i,1)
  read(array(3),*) coord_final(i,2)
  read(array(4),*) coord_final(i,3)
end do

! Find the difference of the two geometries 

allocate  (coord_diff (number_of_atoms,3) )

coord_diff = 0.0d0

coord_diff = coord - coord_final

! Unit conversion

coord = coord*ang2bohr
coord_final = coord_final*ang2bohr
coord_diff = coord_diff*ang2bohr
mass = mass*amu2au

end subroutine read_atomic_data
!--------------------------------------------------------------------------
end module read_geometry
