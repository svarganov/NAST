module normal_mode

use get_constants

implicit none
double precision :: hbar, amu2kg, C
!---------------------------------------------------------------------
parameter (hbar = 1.054571725d-14)                      ! (kg*Ang**2/s**2)*s
parameter (amu2kg = 1.660539d-27)                       ! kg
parameter (C = 299792458*100.d0)                        ! speed of light in cm/s
!---------------------------------------------------------------------
contains
!---------------------------------------------------------------------
subroutine get_cart_disp(file_id,numbr,symb,chrg,coord,freq,eig_v_amu,n)

implicit none
integer,intent(in)                 :: file_id,numbr
character(len=15),intent(in)       :: symb(numbr)
double precision,intent(in)        :: coord(numbr,3),         &
                                      chrg(numbr),freq,       &
                                      eig_v_amu(3*numbr)
integer,intent(in),optional        :: n
!---------------------------------------------------------------------
character(len=32)                  :: output
integer                            :: i,j,point_pos
double precision                   :: coord_ang_eq(3*numbr),  &
                                      eig_v_to_kg(3*numbr),   &
                                      qmin,qmax,dq,q,qd,omega,          &
                                      displacement(3*numbr),  &
                                      cartesian_displacement(3*numbr)
!---------------------------------------------------------------------
if (present(n)) write(file_id,12) n

!~~~ Transform coordinates from au to Angstrom
coord_ang_eq = 0.0d0

do i = 1, numbr
 coord_ang_eq(i*3-2) = coord(i,1)/ang2bohr
 coord_ang_eq(i*3-1) = coord(i,2)/ang2bohr
 coord_ang_eq(i*3)   = coord(i,3)/ang2bohr
end do

!========== DIMENSION OF EIG_V_TO_KG ==============!
!== sqrt(1/amu)*sqrt(1/kg/amu)=sqrt(1/kg)
eig_v_to_kg = eig_v_amu*sqrt(1/amu2kg)     

!== Transform vibrational frequency to Hertz (s**(-1))
omega = 2.0d0*pi*C*freq 

!======= DIMENSIONLESS NORMAL MODE LOOP
!======= (LOOP OVER DISPLACEMENTS Q
!======  TO GENERATE CARTESIAN GEOMETRIES) =========! 

!== qmin and qmax are displacement sampling borders
Qmin = -5.0d0
Qmax =  5.0d0

!== dq is a displacement step along the RC
dq   = (qmax - qmin)/20

do i = 0, 20
!== q is the displacement at iteration # I
 q = qmin + dble(i)*dq 
 
!================ DIMENSION OF QD ===================!
!== sqrt(kg*(Ang**2)*(s**-1)/(s**-1)) = sqrt(kg)*Ang

 qd = sqrt(hbar/omega)*q

!===== NORMAL MODES TO CARTESIAN TRANSFORMATION =====!
!============= DIMENSION OF RC_DISP =================! 
!== sqrt(1/kg)*sqrt(kg)*Ang = Ang        
   
  displacement = eig_v_to_kg*qd

!========== GEOMETRY OUTPUT IN XYZ FORMAT ===========!

cartesian_displacement = coord_ang_eq + displacement

!*** Write geometries in file at each loop step ****!

write(file_id,10) i

do j = 1, numbr
 write(file_id,11) symb(j),chrg(j),&
 cartesian_displacement(j*3-2),        &
 cartesian_displacement(j*3-1),        &
 cartesian_displacement(j*3)
end do
!************************ End **********************!

end do

!=================== END LOOP ======================!

!************ Write format specification ***********!
  
10 FORMAT(/,'Geometry # ',I2)
11 FORMAT(A10,4X,F4.1,3X,F17.9,3X,F17.9,4X,F17.9)
12 FORMAT(/,'A new imaginary mode is found. The mode # is Q',I2)

end subroutine get_cart_disp
!---------------------------------------------------
end module normal_mode
