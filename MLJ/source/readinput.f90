module input_module
 use precision_module, only : wp
 use error_module
 use constants_module
 implicit none
 public
   
 integer                :: Tpoints
 real(wp), allocatable  :: freR(:), S(:), lambda_vib(:)
 real(wp), allocatable  :: tmp(:)
 real(wp)               :: enR, enP, gap
 real(wp)               :: V, lambda_ext
 real(wp)               :: T1, T2, Tstep, trsh = 4.80d0
 integer                :: symR=1, chir=1, n_vib=-1

 contains
 
 subroutine read_file()
  
 namelist /inputdata/ symR, chir, freR, enR, enP, S, lambda_vib, & 
                      lambda_ext, V, T1, T2, Tstep, n_vib, trsh

! Begin reading input data.
  
! Allocate larger arrays to make sure you read
! all frequencies and reorganization energies 
! provided by user in an input file.

allocate (freR(400))
allocate (lambda_vib(400))
allocate (S(400))

freR = zero
lambda_vib = zero

read (unit=11, nml=inputdata)

! Below is automated search for number of non-zero elements in
! freR and lambda_vib. These arrays will be re-allocated to contain only non-zero elements.

allocate (tmp(size(pack(freR, freR /= 0))))
tmp = pack(freR, freR /= 0)
deallocate (freR)

! Allocate freR and lambda_vib again, but this time the size of array is equal to
! the number of non-zero elements in arrays freR and lambda_vib read from namelist.

allocate (freR(size(tmp)))
freR = tmp
deallocate (tmp)

allocate (tmp(size(pack(lambda_vib, lambda_vib /= 0))))
tmp = pack(lambda_vib, lambda_vib /= 0)
deallocate (lambda_vib)

allocate (lambda_vib(size(tmp)))
lambda_vib = tmp
deallocate (tmp)

allocate (tmp(size(pack(S, S /= 0))))
tmp = pack(S, S /= 0)
deallocate (S)

allocate (S(size(tmp)))
S = tmp
deallocate (tmp)

! Check if the freR, lambda_vib, and S are the arrays of the same size
if (size(lambda_vib) .ne. size(freR)) then
  call error_message(1)
end if

if (size(S) .ne. size(freR)) then
  call error_message(1)
end if

if (T1 == zero) then
  T1  = 290_wp
end if
if (T2 == zero) then
   T2 = 300_wp
end if
if (T2 < T1) then
   T2 = T1
end if
if (Tstep == zero) then
  Tstep = one
end if
Tpoints = int((T2 - T1)/Tstep)

gap = enP - enR

! Units conversions

freR       = freR/autocm     ! Converts frequency in cm-1 to hv in hartree
V          = V/autocm

 end subroutine read_file

end module input_module
