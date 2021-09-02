module fit

implicit none

contains
!----------------------------------------------------
subroutine linear_polynom_fit(n,x,y,p1,p2)

implicit none
integer, intent(in)             :: n
double precision, intent(in)    :: x(n), y(n)
double precision, intent(out)   :: p1, p2
!----------------------------------------------------
double precision                :: sumx, sumy, sumx2, sumxy

sumx = 0.0d0
sumy = 0.0d0
sumx2 = 0.0d0
sumxy = 0.0d0

sumx = sum(x)
sumy = sum(y)
sumx2 = sum(x**2)
sumxy = sum(x*y)

p1 = (n*sumxy - sumx*sumy)/(n*sumx2 - sumx**2)
p2 = (sumy - p1*sumx)/n

end subroutine
!----------------------------------------------------
subroutine quadr_polynom_fit(n,x,y,p1,p2,p3)

implicit none
integer, intent(in)             :: n
double precision, intent(in)    :: x(n), y(n)
double precision, intent(out)   :: p1, p2, p3
!----------------------------------------------------
double precision                :: S1,S2,S3,S4,Sy0,Sy1,Sy2
double precision                :: A(3,3),B(3)
integer                         :: info,lwork,ipiv(3)
double precision                :: work(3)
!----------------------------------------------------

A = 0.0d0
B = 0.0d0

S1 = 0.0d0
S2 = 0.0d0
S3 = 0.0d0
S4 = 0.0d0
Sy0 = 0.0d0
Sy1 = 0.0d0
Sy2 = 0.0d0

S1 = sum(x)
S2 = sum(x*x)
S3 = sum(x*x*x)
S4 = sum(x*x*x*x)

Sy0 = sum(y)
Sy1 = sum(y*x)
Sy2 = sum(y*x*x)

A(1,:) = (/S4, S3, S2/)
A(2,:) = (/S3, S2, S1/)
A(3,:) = (/S2, S1, dble(n)/)

B(:) = (/Sy2, Sy1, Sy0/)

! ==== Calculate inverse of matrix A ==== !

call dgetrf(3,3,A,3,ipiv,info)
lwork = size(work)
call dgetri(3,A,3,ipiv,work,lwork,info)

p1 = dot_product(A(1,:),B)
p2 = dot_product(A(2,:),B)
p3 = dot_product(A(3,:),B)

end subroutine
!----------------------------------------------------
subroutine cubic_polynom_fit(n,x,y,p1,p2,p3,p4)

implicit none
integer, intent(in)             :: n
double precision, intent(in)    :: x(n), y(n)
double precision, intent(out)   :: p1,p2,p3,p4
!----------------------------------------------------
double precision                :: S1,S2,S3,S4,S5,S6
double precision                :: Sy0,Sy1,Sy2,Sy3
double precision                :: A(4,4),B(4)
integer                         :: info,lwork,ipiv(4)
double precision                :: work(4)
!----------------------------------------------------

A = 0.0d0
B = 0.0d0

S1 = sum(x)
S2 = sum(x*x)
S3 = sum(x*x*x)
S4 = sum(x*x*x*x)
S5 = sum(x*x*x*x*x)
S6 = sum(x*x*x*x*x*x)

Sy0 = sum(y)
Sy1 = sum(y*x)
Sy2 = sum(y*x*x)
Sy3 = sum(y*x*x*x)

A(1,:) = (/S6, S5, S4, S3/)
A(2,:) = (/S5, S4, S3, S2/)
A(3,:) = (/S4, S3, S2, S1/)
A(4,:) = (/S3, S2, S1, dble(n)/)

B(:) = (/Sy3, Sy2, Sy1, Sy0/)

! ==== Calculate inverse of matrix A ==== !

call dgetrf(4,4,A,4,ipiv,info)
lwork = size(work)
call dgetri(4,A,4,ipiv,work,lwork,info)

p1 = dot_product(A(1,:),B)
p2 = dot_product(A(2,:),B)
p3 = dot_product(A(3,:),B)
p4 = dot_product(A(4,:),B)
 
end subroutine
!----------------------------------------------------
subroutine quart_polynom_fit(n,x,y,p1,p2,p3,p4,p5)

implicit none
integer, intent(in)             :: n
double precision, intent(in)    :: x(n), y(n)
double precision, intent(out)   :: p1,p2,p3,p4,p5
!----------------------------------------------------
double precision                :: S1,S2,S3,S4,S5,S6
double precision                :: S7,S8,Sy0,Sy1,Sy2
double precision                :: Sy3,Sy4
double precision                :: A(5,5),B(5)
integer                         :: info,lwork,ipiv(5)
double precision                :: work(5)
!----------------------------------------------------

A = 0.0d0
B = 0.0d0

S1 = sum(x)
S2 = sum(x*x)
S3 = sum(x*x*x)
S4 = sum(x*x*x*x)
S5 = sum(x*x*x*x*x)
S6 = sum(x*x*x*x*x*x)
S7 = sum(x*x*x*x*x*x*x)
S8 = sum(x*x*x*x*x*x*x*x)

Sy0 = sum(y)
Sy1 = sum(y*x)
Sy2 = sum(y*x*x)
Sy3 = sum(y*x*x*x)
Sy4 = sum(y*x*x*x*x)

A(1,:) = (/S8, S7, S6, S5, S4/)
A(2,:) = (/S7, S6, S5, S4, S3/)
A(3,:) = (/S6, S5, S4, S3, S2/)
A(4,:) = (/S5, S4, S3, S2, S1/)
A(5,:) = (/S4, S3, S2, S1, dble(n)/)

B(:) = (/Sy4, Sy3, Sy2, Sy1, Sy0/)

! ==== Calculate inverse of matrix A ==== !

call dgetrf(5,5,A,5,ipiv,info)
lwork = size(work)
call dgetri(5,A,5,ipiv,work,lwork,info)

p1 = dot_product(A(1,:),B)
p2 = dot_product(A(2,:),B)
p3 = dot_product(A(3,:),B)
p4 = dot_product(A(4,:),B)
p5 = dot_product(A(5,:),B)
 

end subroutine
!-------------------------
end module
