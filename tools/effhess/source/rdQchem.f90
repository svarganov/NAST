module QChem_read_module

use M_strings
use hess_data

implicit none

contains
!-------------------------------------
subroutine QChem_read_number_of_atoms(file_id)
implicit none
integer, intent(in)                           :: file_id
character(len=80)                             :: line

do while (index(line, 'Atomic numbers') == 0)
 read(file_id,'(A)') line
end do                                             ! Now on the Force Constants line
call split(line,array)
read(array(5),*) number_of_atoms
rewind file_id

end subroutine QChem_read_number_of_atoms
!----------------------------------------
subroutine QChem_read_atomic_data(file_id1,file_id2)

use get_constants
use get_atomic_masses

implicit none
integer, intent(in)                                             :: file_id1, file_id2
!-------------------------------------------------
character(len=80)                                               :: line
character(len=80), allocatable                                  :: array(:)
integer                                                         :: i, j, m
!-------------------------------------------------
allocate  (symbol (number_of_atoms)   )
allocate  (charge (number_of_atoms)   )
allocate  (coord  (number_of_atoms,3) )
allocate  (mass   (number_of_atoms)   )
!-------------------------------------------------

charge = 0.0d0
coord = 0.0d0

do while (index(line, 'molecule') == 0)
 read(file_id1,'(A)') line
end do

read(file_id1,'(A)') line

do i = 1, number_of_atoms
 read(file_id1,'(A)') line
 call split(line, array)
 read(array(1),*) symbol(i)
 read(array(2),*) coord(i,1)
 read(array(3),*) coord(i,2)
 read(array(4),*) coord(i,3)
end do

rewind file_id1

do while (index(line, 'Atomic numbers') == 0)
  read(file_id2,'(A)') line
end do                              

m=0
do i=1, int(ceiling(real(number_of_atoms)/6))
  read(file_id2,'(A)') line
  call split(line,array)
  if (i*6 <= number_of_atoms) then
    do j=1,6
     read(array(j),*) charge(j+m)
    end do
    m=m+6
  else 
    do j=1, mod(number_of_atoms,6)
     read(array(j),*) charge(j+m)
    end do
  end if
end do
 
rewind file_id2

do i = 1, number_of_atoms
  mass(i) = ams(int(charge(i)))
end do

coord = coord*ang2bohr
mass = mass*amu2au

end subroutine QChem_read_atomic_data
!--------------------------
subroutine QChem_read_grad_hess(file_id1,file_id2,grad,Hess)

implicit none
integer, intent(in)                                           :: file_id1, file_id2
double precision, dimension(:), allocatable, intent(out)      :: grad
double precision, dimension(:,:), allocatable, intent(out)    :: Hess
!----------------------------------------------------------------------------------
character(len=80)                                             :: line
character(len=80), allocatable                                :: array(:)
double precision, dimension(:), allocatable                   :: deployed
integer                                                       :: i, j, m, N, last_number, column
!----------------------------------------------------------------------------------
allocate  (grad    (3*number_of_atoms) )
allocate  (Hess (3*number_of_atoms,3*number_of_atoms) )
!----------------------------------------------------------------------------------
! Read gradient

grad= 0.d0

do while (index(line, 'Cartesian Forces') == 0)
 read(file_id1,'(A)') line
end do

m=0
do i=1, int(ceiling(real(3*number_of_atoms)/5))
  read(file_id1,'(A)') line
  call split(line,array)
  if (i*5 <= 3*number_of_atoms) then
    do j=1,5
     read(array(j),*) grad(j+m)
    end do
    m=m+5
  else 
    do j=1, mod(3*number_of_atoms,5)
     read(array(j),*) grad(j+m)
    end do
  end if
end do 
rewind file_id1

!------------------------------------------------
! Read Hessian

Hess = 0.0d0

N = int(3*real(number_of_atoms)*((3*real(number_of_atoms) - 1)/2 + 1))

allocate  (deployed (N) )
deployed = 0.0d0

do while (index(line, 'Cartesian Force Constants') == 0)
 read(file_id2,'(A)') line
end do                                             ! Now on the Force Constants line

do i = 1, int(ceiling(real(N)/5))

  if (i*5 <= N) then
    column = 5
  else
    column = MOD(N,5)              
  end if
  read(file_id2,'(A)') line
  call split(line, array)
  if (column == 5) then
    read(array(1),*) deployed((i-1)*5 + 1)
    read(array(2),*) deployed((i-1)*5 + 2)
    read(array(3),*) deployed((i-1)*5 + 3)
    read(array(4),*) deployed((i-1)*5 + 4)
    read(array(5),*) deployed((i-1)*5 + 5)
  else if (column == 1) then
    read(array(1),*) deployed((i-1)*5 + 1)
  else if (column == 2) then
    read(array(1),*) deployed((i-1)*5 + 1)
    read(array(2),*) deployed((i-1)*5 + 2)
  else if (column == 3) then
    read(array(1),*) deployed((i-1)*5 + 1)
    read(array(2),*) deployed((i-1)*5 + 2)
    read(array(3),*) deployed((i-1)*5 + 3)
  else if (column == 4) then
    read(array(1),*) deployed((i-1)*5 + 1)
    read(array(2),*) deployed((i-1)*5 + 2)
    read(array(3),*) deployed((i-1)*5 + 3)
    read(array(4),*) deployed((i-1)*5 + 4)
  end if
end do

rewind file_id2

m = 0
do i=1, 3*number_of_atoms
  do j=1, i
    Hess(j,i)=deployed(j+m)
    Hess(i,j)=Hess(j,i)
  end do
  m=m+i
end do
end subroutine QChem_read_grad_hess
end module QChem_read_module
