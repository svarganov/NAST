module vibrational

use hess_data

implicit none
contains
!------------------
subroutine calc_center_of_mass()

implicit none
integer       :: i, j

allocate  (coord_shifted    (number_of_atoms,3))
allocate  (coord_shifted_mw (number_of_atoms,3))

total_mass = 0.0d0
center_of_mass = 0.0d0
do i = 1, number_of_atoms
 center_of_mass(1,:) = center_of_mass(1,:) + coord(i,:)*mass(i)
 total_mass = total_mass + mass(i)
end do

center_of_mass = center_of_mass/total_mass

do i = 1, number_of_atoms
 coord_shifted(i,:)=coord(i,:)-center_of_mass(1,:)
 coord_shifted_mw(i,:)=coord_shifted(i,:)*sqrt(mass(i))
end do

end subroutine calc_center_of_mass
!-----------------------------------------------------------------------------
subroutine calc_hess_mw()

implicit none
integer      :: i, j
!-----------------------------------------------------------------------------
allocate  (Hess_mw (3*number_of_atoms,3*number_of_atoms))
!-----------------------------------------------------------------------------

do i = 1, 3*number_of_atoms
 do j = 1, i
  Hess_mw(i,j) = Hess(i,j)/sqrt(mass(ceiling(real(i)/3))*mass(ceiling(real(j)/3)))
  if (i /= j) then
    Hess_mw(j,i) = Hess_mw(i,j)
  end if 
 end do
end do

end subroutine calc_hess_mw
!-----------------------------------------------------------------------------
subroutine calc_rot()

use get_constants

implicit none
integer             :: l, inf
double precision    :: work(3*(3+3/2))
integer      :: i

l = 3*(3+3/2)

do i = 1, number_of_atoms
  rot(1,1) = rot(1,1) + coord_shifted_mw(i,2)**2 + coord_shifted_mw(i,3)**2
  rot(1,2) = rot(1,2) - coord_shifted_mw(i,1)*coord_shifted_mw(i,2)
  rot(1,3) = rot(1,3) - coord_shifted_mw(i,1)*coord_shifted_mw(i,3)
  rot(2,2) = rot(2,2) + coord_shifted_mw(i,1)**2 + coord_shifted_mw(i,3)**2
  rot(2,3) = rot(2,3) - coord_shifted_mw(i,2)*coord_shifted_mw(i,3)
  rot(3,3) = rot(3,3) + coord_shifted_mw(i,1)**2 + coord_shifted_mw(i,2)**2
end do
rot(2,1) = rot(1,2)
rot(3,1) = rot(1,3)
rot(3,2) = rot(2,3)

rot_to_inv = rot
call dsyev('V','U',3,rot,3,moments_eig,work,l,inf)
moments_eig = moments_eig/amu2au
end subroutine calc_rot
!---------------------------------------------------------------------------------
subroutine calc_inverse()

implicit none
integer                :: info,lwork
integer                :: ipiv(3)
double precision       :: work(3)

call dgetrf(3,3,rot_to_inv,3,ipiv,info)
lwork = size(work)
call dgetri(3,rot_to_inv,3,ipiv,work,lwork,info)

end subroutine calc_inverse
!---------------------------------------------------------------------------------
subroutine proj_construct()

implicit none
!---------------------------------------------------------------------------------
integer :: i, j, k
integer :: ip,jp,ia,ib,ic,ii,ja,jb,jc,jj,kk,ll,indx,jndx,kndx,lndx,jend
double precision :: tens(3,3,3),coord_shifted_mw_vec(3*number_of_atoms),summ
!---------------------------------------------------------------------------------
allocate (trans_projector  (3*number_of_atoms,3*number_of_atoms))
allocate (rot_projector    (3*number_of_atoms,3*number_of_atoms))
allocate (RM               (3*number_of_atoms))
allocate (P                (3*number_of_atoms,3*number_of_atoms))

do i = 1, number_of_atoms
    coord_shifted_mw_vec(3*i-2) = coord_shifted_mw(i,1) 
    coord_shifted_mw_vec(3*i-1) = coord_shifted_mw(i,2)
    coord_shifted_mw_vec(3*i)   = coord_shifted_mw(i,3)
end do

do i = 1, number_of_atoms
  RM(3*i-2) = sqrt(mass(i))
  RM(3*i-1) = sqrt(mass(i))
  RM(3*i)   = sqrt(mass(i))
end do

tens = 0.0d0
tens(1,2,3) = 1.0d0
tens(2,3,1) = 1.0d0
tens(3,1,2) = 1.0d0
tens(3,2,1) = -1.0d0
tens(1,3,2) = -1.0d0
tens(2,1,3) = -1.0d0

do ip = 1, number_of_atoms
 indx = 3*(ip-1)
 kndx = max(3*(ip-1),6*(ip-1)-3*number_of_atoms)
 do jp = 1, ip
   jndx = 3*(jp-1)
   lndx = max(3*(jp-1),6*(jp-1)-3*number_of_atoms)
   do ic = 1, 3
     ii = indx + ic
     kk = kndx + ic
     jend = 3
     if (jp == ip) then
       jend = ic
     end if
     do jc = 1, jend
       jj = jndx + jc
       ll = lndx + jc
       summ = 0.0d0
       do ia = 1, 3
         do ib =  1, 3
           if (tens(ia,ib,ic) == 0) then
            cycle
           end if
           do ja = 1, 3
             do jb = 1, 3
               if (tens(ja,jb,jc) == 0) then
                cycle
               end if
             summ = summ + tens(ia,ib,ic)*tens(ja,jb,jc)*rot_to_inv(ia,ja) &
                    *coord_shifted_mw_vec(indx+ib)*coord_shifted_mw_vec(jndx+jb)
            end do
           end do
         end do
       end do
     P(kk,ll) = summ
     rot_projector(kk,ll) = P(kk,ll)
     if (ic == jc) then
       trans_projector(kk,ll) = (RM(ii)*RM(jj))/total_mass
       P(kk,ll) = P(kk,ll) + trans_projector(kk,ll)
     end if
     end do
   end do 
 end do
end do

do i = 1, 3*number_of_atoms
 do j = 1, 3*number_of_atoms
  if (j /= i) then
    trans_projector(i,j) = trans_projector(j,i)
    rot_projector(i,j)   = rot_projector(j,i)
    P(i,j) = P(j,i)
  end if
 end do
end do

end subroutine proj_construct
!--------------------------------------------------------------------------------
subroutine project_hess()

implicit none

integer            :: i, l, inf
double precision   :: eye(3*number_of_atoms,3*number_of_atoms)
double precision   :: work(3*number_of_atoms*(3+number_of_atoms/2))
!---------------------------------------------------------------------------------
allocate  (projector          (3*number_of_atoms, 3*number_of_atoms))
allocate  (Hess_projected     (3*number_of_atoms, 3*number_of_atoms))
allocate  (e                  (3*number_of_atoms))
!---------------------------------------------------------------------------------
eye = 0.0d0
l = 3*number_of_atoms*(3+number_of_atoms/2)

forall (i = 1:3*number_of_atoms)  eye(i,i) = 1.0d0

projector = eye - P
Hess_projected = matmul(projector,matmul(Hess_mw,projector))
call dsyev('V','U',3*number_of_atoms,Hess_projected,3*number_of_atoms,e,work,l,inf)

end subroutine project_hess
!--------------------------------------------------------------------------------
subroutine get_freq(input)

use get_constants
use normal_mode

implicit none
character(len=*),intent(in)               :: input
!--------------------------------------------------------------------------------
integer                                   :: n_imag,i,j,point_pos
double precision                          :: imag_freq, imag_eig_v_amu(3*number_of_atoms)
!--------------------------------------------------------------------------------
allocate  (freq (3*number_of_atoms))
allocate  (all_eig_v_amu(3*number_of_atoms))
!--------------------------------------------------------------------------------
n_imag = 0

do i = 1, 3*number_of_atoms
 imag_freq = 0.0d0
 imag_eig_v_amu = 0.0d0
 if (e(i) < 0 .and. abs(e(i)) > 1d-15) then
  n_imag = n_imag + 1
  if (n_imag == 1) then
   open (unit = 22, action = "write", file = input)
  end if

  imag_freq = sqrt(abs(e(i)))/freq_conv
  do j = 1, 3*number_of_atoms
   imag_eig_v_amu(j) = Hess_projected(j,i) &
                       /sqrt(mass(ceiling(real(j)/3))/amu2au)
  end do
  call get_cart_disp(22,number_of_atoms,symbol,charge,coord, &
                     imag_freq,imag_eig_v_amu,i)
 end if
end do

if (n_imag > 0) then
 
 write(66,12)

end if

12 format(/,4X,'    Warning! There is (are) imaginary frequency(ies) found!        ',&
          /,4X,'            Vibrational analysis is not valid.                     ',&
          /,4X,'        (Unless you have reaaally large molecule and               ',&
          /,4X,'       imaginary frequency is very small (< 5-10 cm-1)             ',&
          /,4X,'       Cartesian coordinates corresponding to displacements        ',&
          /,4X,'       along each imaginary frequency mode will be generated       ',&
          /,4X,'                     in the .imag file .                           ',&
          /,4X,'       You can use those coordinates as a starting guess           ',&
          /,4X,'                     for the re-optimization.                      ',/)

freq = sqrt(abs(e))/freq_conv
end subroutine get_freq
!----------------------------------------------------------------------------------
subroutine Huang_Rhys()

use get_constants
use read_geometry

implicit none
!--------------------------------------------------------------------------------
integer                                   :: i, j 
double precision                          :: dq_coord_v (3*number_of_atoms), dQ (3*number_of_atoms)
!--------------------------------------------------------------------------------

allocate  (coord_diff_mw         (number_of_atoms,3))
allocate  (HRF (3*number_of_atoms))
allocate  (reorganization_energy (3*number_of_atoms))

do i = 1, number_of_atoms
  coord_diff_mw(i,:) = coord_diff(i,:)*sqrt(mass(i))
end do

dq_coord_v = 0.d0
dq_coord_v = reshape(transpose(coord_diff_mw),[3*number_of_atoms])

dQ = matmul(transpose(Hess_projected),dq_coord_v)
HRF = 0.0d0
HRF = 0.5*sqrt(e)*dQ**2
reorganization_energy = 0.5*e*dQ**2

end subroutine Huang_Rhys
!----------------------------------------------------------------------------------
end module vibrational
