module get_constants
implicit none

double precision :: ang2bohr, amu2au, time_au2s, speed_of_light,  &
                    pi, freq_conv
!---------------------------------------------------------------------
parameter (ang2bohr = 1d0/0.529177249238d0)             ! GAMESS conversion factor
parameter (amu2au = 1822.888484995826d0)                ! GAMESS conversion factor
parameter (time_au2s = 2.418884326505d-17)
parameter (speed_of_light = 299792458)
parameter (pi = 3.141592653589793238462643383279502d0)
parameter (freq_conv = 2.0d0*pi*100.d0*speed_of_light*time_au2s)
end module get_constants
