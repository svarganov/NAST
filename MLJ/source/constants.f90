module constants_module
      use precision_module, only : wp

      real(wp), parameter :: pi=3.1415926535897932d0
      real(wp), parameter :: planck = 4.1356676623401643d-15  ! eV/Hz
      real(wp), parameter :: planck_hsec = 1.519829845d-16    ! Hartree/Hz 
      real(wp), parameter :: kb_hartK = 3.16681053d-6         ! Hartree/K 
      real(wp), parameter :: autocm = 219474.63137172544d0    ! cm-1/Hartree 
      real(wp), parameter :: amutoau = 1822.8884853323708d0   ! Dalton/au
      real(wp), parameter :: autoev = 27.21138602d0           ! eV/Hartree  
      real(wp), parameter :: autokcalpermol = 627.50947377753740d0 ! (kcal/mol)/Hartree
      real(wp), parameter :: cmtohz = 29979245800d0           ! cm-1/Hz
      real(wp), parameter :: mu = 0.46686451877418d0          ! Bohr magneton, cm-1/Tesla
      real(wp), parameter :: zero =  0.0d0
      real(wp), parameter :: one =   1.0d0
      real(wp), parameter :: two =   2.0d0
      real(wp), parameter :: three = 3.0d0
      real(wp), parameter :: four  = 4.0d0
      real(wp), parameter :: five  = 5.0d0
      real(wp), parameter :: half =  0.5d0
end module
