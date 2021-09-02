module get_atomic_masses
implicit none

double precision, dimension(43) :: ams  = (/ 1.007825d0, 4.0026d0, 7.01600d0, 9.01218d0, 11.00931d0, 12.0d0, &
         14.00307d0, 15.99491d0, 18.99840d0, 19.99244d0, 22.9898d0,      &
         23.98504d0, 26.98153d0, 27.97693d0, 30.97376d0, 31.97207d0,     &
         34.96885d0, 39.948d0, 38.96371d0, 39.96259d0, 44.95592d0,       &
         47.90d0, 50.9440d0, 51.9405d0, 54.9381d0, 55.9349d0, 58.9332d0, &
         57.9353d0, 62.9298d0, 63.9291d0, 68.9257d0, 73.9219d0,          &
         74.9216d0, 79.9165d0, 78.9183d0, 83.9115d0, 84.9117d0,          &
         87.9056d0, 89.9054d0, 89.9043d0, 92.9060d0, 97.9055d0, 97.0d0 /)
end module get_atomic_masses

!===========================================!

                ! SOURCE !

!Hydrogen    H    1    1.007825d0
!Helium      He   2    4.0026
!Lithium     Li   3    7.01600
!Beryllium   Be   4    9.01218
!Boron       B    5    11.00931
!Carbon      C    6    12.0
!Nitrogen    N    7    14.00307d0
!Oxygen      O    8    15.99491d0
!Fluorine    F    9    18.99840
!Neon        Ne   10   19.99244
!Sodium      Na   11   22.9898
!Magnesium   Mg   12   23.98504
!Aluminium   Al   13   26.98153
!Silicon     Si   14   27.97693
!Phosphorus  P    15   30.97376
!Sulfur      S    16   31.97207
!Chlorine    Cl   17   34.96885
!Argon       Ar   18   39.948
!Potassium   K    19   38.96371
!Calcium     Ca   20   39.96259
!Scandium    Sc   21   44.95592
!Titanium    Ti   22   47.90
!Vanadium    V    23   50.9440
!Chromium    Cr   24   51.9405
!Manganese   Mn   25   54.9381
!Iron        Fe   26   55.9349
!Cobalt      Co   27   58.9332
!Nickel      Ni   28   57.9353
!Copper      Cu   29   62.9298
!Zinc        Zn   30   63.9291
!Gallium     Ga   31   68.9257
!Germanium   Ge   32   73.9219
!Arsenic     As   33   74.9216
!Selenium    Se   34   79.9165
!Bromine     Br   35   78.9183
!Krypton     Kr   36   83.9115
!Rubidium    Rb   37   84.9117
!Strontium   Sr   38   87.9056
!Yttrium     Y    39   89.9054
!Zirconium   Zr   40   89.9043
!Niobium     Nb   41   92.9060
!Molybdenur  Mo   42   97.9055
!Technetium  Tc   43   97.0
!======================================!
