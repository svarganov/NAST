module vibrational

use hess_data

implicit none
contains
!------------------
subroutine construct_eff_hess()    

implicit none
!----------------------------------------------------------------------------------
allocate  (deltaG (3*number_of_atoms), Eff_Hess(3*number_of_atoms,3*number_of_atoms))

deltaG = grad1 - grad2
norm_deltaG = norm2(deltaG)
norm_grad1 = norm2(grad1)
norm_grad2 = norm2(grad2)
lambda = dot_product(deltaG,grad1)/(norm_deltaG**2)
gradmean = sqrt(norm_grad1*norm_grad2)
Eff_Hess = Hess1 - lambda*(Hess1 - Hess2)
!Eff_Hess = Hess2 +(dot_product(deltaG,grad2)/(norm_deltaG**2))*(Hess2-Hess1)
!Eff_Hess = norm_grad1*Hess2-norm_grad2*Hess1
!Eff_Hess = Eff_Hess/norm_deltaG
gradslope = dot_product(grad1,grad2)

if (gradslope < 0) then
  intersect_type = 'peaked'
else
  intersect_type = 'sloped'
end if

end subroutine construct_eff_hess
!-------------------------------------------
subroutine calc_center_of_mass()

use get_constants
use get_atomic_masses

implicit none
integer       :: i, j

allocate  (coord_shifted    (number_of_atoms,3))
allocate  (coord_shifted_mw (number_of_atoms,3))
allocate  (mass             (number_of_atoms))

coord = coord*ang2bohr                                 ! Convert atomic coordinates to bohr 

do i = 1, number_of_atoms
 mass(i) = ams(int(charge(i)))                         ! ams is array of atomic masses. Comes from ams.f90
end do

mass = mass*amu2au                                     ! Convert mass from atomic mass units to atomic units

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
subroutine calc_eff_hess_mw()

implicit none
integer      :: i, j
!-----------------------------------------------------------------------------
allocate  (Eff_Hess_mw (3*number_of_atoms,3*number_of_atoms))
!-----------------------------------------------------------------------------

do i = 1, 3*number_of_atoms
 do j = 1, i
  Eff_Hess_mw(i,j) = Eff_Hess(i,j)/sqrt(mass(ceiling(real(i)/3))*mass(ceiling(real(j)/3)))
  if (i /= j) then
    Eff_Hess_mw(j,i) = Eff_Hess_mw(i,j)
  end if 
 end do
end do

end subroutine calc_eff_hess_mw
!-----------------------------------------------------------------------------
subroutine calc_norm_deltaG_mw()

implicit none
integer      :: i
!-----------------------------------------------------------------------------
allocate  (m_deltaG      (3*number_of_atoms))
allocate  (m_norm_deltaG (3*number_of_atoms))
!-----------------------------------------------------------------------------

do i = 1, number_of_atoms
 m_deltaG(i*3 - 2) = deltaG(i*3 - 2)/sqrt(mass(i))
 m_deltaG(i*3 - 1) = deltaG(i*3 - 1)/sqrt(mass(i))
 m_deltaG(i*3)     = deltaG(I*3)/sqrt(mass(i))
end do

m_norm_deltaG = m_deltaG/norm2(m_deltaG)

end subroutine calc_norm_deltaG_mw
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

rot_to_inv = rot                                    ! Need to make a copy of 'rot'
                                                    ! because the following diagonalization of 'rot'
                                                    ! will overwrite 'rot' with 'rot' eigenvectors,
                                                    ! but we still need 'rot' to calculate
                                                    ! its inverse. So we'll use 'rot_to_inv' to make inverse
                                                    
                                                    ! Diagonalize the Moment of Inertia Tensor (Matrix)
                                                    ! Eigenalues (principal moments of inertia) are stored in moments_eig 
call dsyev('V','U',3,rot,3,moments_eig,work,l,inf)
moments_eig = moments_eig/amu2au                    ! Convert moments from au*bohr^2 to amu*bohr^2
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
allocate (rc_projector     (3*number_of_atoms,3*number_of_atoms))
allocate (trans_projector  (3*number_of_atoms,3*number_of_atoms))
allocate (rot_projector    (3*number_of_atoms,3*number_of_atoms))
allocate (RM               (3*number_of_atoms))
allocate (P                (3*number_of_atoms,3*number_of_atoms))

!~~~~~ Reshape coord_shifted_mw  (number_of_atoms,3) matrix
!~~~~~ to a coord_shifted_mw_vec (3*number_of_atoms) vector

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
     rc_projector(kk,ll) = m_norm_deltaG(kk)*m_norm_deltaG(ll)
     P(kk,ll) = P(kk,ll) + rc_projector(kk,ll)
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
    rc_projector(i,j)    = rc_projector(j,i)       ! RC projector
    trans_projector(i,j) = trans_projector(j,i)    ! Translational projector
    rot_projector(i,j)   = rot_projector(j,i)      ! Rotational projector
    P(i,j) = P(j,i)                                ! Total projector
  end if
 end do
end do

end subroutine proj_construct
!--------------------------------------------------------------------------------
subroutine project_eff_hess()

implicit none

integer            :: i, l, inf
double precision   :: eye(3*number_of_atoms,3*number_of_atoms)
double precision   :: work(3*number_of_atoms*(3+number_of_atoms/2))
!---------------------------------------------------------------------------------
allocate  (projector          (3*number_of_atoms, 3*number_of_atoms))
allocate  (Eff_Hess_projected (3*number_of_atoms, 3*number_of_atoms))
allocate  (e                  (3*number_of_atoms))
!---------------------------------------------------------------------------------

eye = 0.0d0
l = 3*number_of_atoms*(3+number_of_atoms/2)

forall (i = 1:3*number_of_atoms)  eye(i,i) = 1.0d0       ! Construct Identity matrix(3N, 3N)

projector = eye - P
Eff_Hess_projected = matmul(projector,matmul(Eff_Hess_mw,projector))
call dsyev('V','U',3*number_of_atoms,Eff_Hess_projected,3*number_of_atoms,e,work,l,inf)

end subroutine project_eff_hess
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
!--------------------------------------------------------------------------------
n_imag = 0

do i = 1, 3*number_of_atoms
 imag_freq = 0.0d0
 imag_eig_v_amu = 0.0d0
 if (e(i) < 0 .and. abs(e(i)) > 1d-15) then
! > 1d-15 means it is not an eigenvalue of projected
! translational, rotational or rc mode, because those 
! projected eigenvalues should be zero
! (in practice, ~ 1d-20 - 1d-22 
! Therefore, one more imaginary frequency is found
  n_imag = n_imag + 1
  if (n_imag == 1) then
! Create file for imaginary frequency geometries just once
! We don't want to re-create the file each time a new imaginary
! frequency is found
   open (unit = 22, action = "write", file = input)
  end if
! freq_conv converts eigenvalues of the projected
! Hessian matrix to vibrational frequencies (cm-1)

  imag_freq = sqrt(abs(e(i)))/freq_conv
  do j = 1, 3*number_of_atoms
   imag_eig_v_amu(j) = Eff_Hess_projected(j,i) &
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
          /,4X,'                                                                   ',&
          /,4X,'       Cartesian coordinates corresponding to displacements        ',&
          /,4X,'       along each imaginary frequency mode will be generated       ',&
          /,4X,'                     in the .imag file .                           ',&
          /,4X,'       You can use those coordinates as a starting guess           ',&
          /,4X,'                     for the MECP re-search.                       ')

freq = sqrt(abs(e))/freq_conv
end subroutine get_freq
!----------------------------------------------------------------------------------
subroutine get_red_mass()

use get_constants
use get_atomic_masses

implicit none
!----------------------------------------------------------------------------------
integer                             :: i,j
!----------------------------------------------------------------------------------
allocate  (reduced_mass                      (3*number_of_atoms))
allocate  (Eff_Hess_projected_eig_vec_to_amu (3*number_of_atoms, 3*number_of_atoms))
!----------------------------------------------------------------------------------

do j = 1, 3*number_of_atoms             ! j is the eigenvector number
 do i = 1, 3*number_of_atoms            ! i is the index of element in eigenvector j
  Eff_Hess_projected_eig_vec_to_amu(i,j) = &
  Eff_Hess_projected(i,j)/sqrt(mass(ceiling(real(i)/3))/amu2au) 
 end do
 reduced_mass(j) = &
 1/dot_product(Eff_Hess_projected_eig_vec_to_amu(:,j),&
 Eff_Hess_projected_eig_vec_to_amu(:,j)) 
end do

end subroutine get_red_mass
!----------------------------------------------------------------------------------
subroutine get_red_mass_along_rc()

use get_constants
use get_atomic_masses

implicit none
!-----------------------------------------------------------------------------------
integer             :: i, e_rc_max_indx, l, inf
double precision    :: work(3*number_of_atoms*(3+number_of_atoms/2))
double precision    :: e_rc(3*number_of_atoms)
!-----------------------------------------------------------------------------------
allocate  (Eff_Hess_rc     (3*number_of_atoms,3*number_of_atoms))
allocate  (rc_eig_v_amu    (3*number_of_atoms))
!-----------------------------------------------------------------------------------

l = 3*number_of_atoms*(3+number_of_atoms/2)

Eff_Hess_rc = matmul(rc_projector,matmul(Eff_Hess_mw,rc_projector))
call dsyev('V','U',3*number_of_atoms,Eff_Hess_rc,3*number_of_atoms,e_rc,work,l,inf)

! Fortran sorts eigenvalues in the ascending order. When projecting on the RC,
! we expect only !one! non-zero eigenvalue of Eff_Hess_rc, since all other DOS are
! projected out (note that here our projector is not '1-P', but 'P').
! Thus, if the eigenvalue we seek for is positive, it would have
! the index of '3*number_of_atoms', that is, it is the last (highest) in 'e_rc' array.
! If it is negative, it will be the most negative, and then it will be #1 in array.

e_rc_max_indx = 0

if (abs(e_rc(1)) > 0.0000001) then
 e_rc_max_indx = 1
else
 e_rc_max_indx = 3*number_of_atoms
end if
freq_rc = sqrt(abs(e_rc(e_rc_max_indx)))/freq_conv
do i = 1, 3*number_of_atoms         ! i is the index of an element in an eigenvector
 rc_eig_v_amu(i) = Eff_Hess_rc(i,e_rc_max_indx)/sqrt(mass(ceiling(real(i)/3))/amu2au) 
end do

reduced_mass_rc = 1/dot_product(rc_eig_v_amu,rc_eig_v_amu)

end subroutine get_red_mass_along_rc
!-------------------------------------------------------------------------------------
end module vibrational
