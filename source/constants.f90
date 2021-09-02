module constants_module
      use precision_module, only : wp

      real(wp), parameter :: pi=3.1415926535897932d0
      real(wp), parameter :: planck = 4.1356676623401643d-15
      real(wp), parameter :: planck_hsec = 1.519829845d-16
      real(wp), parameter :: kb_hartK = 3.16681053d-6
      real(wp), parameter :: autocm = 219474.63137172544d0
      real(wp), parameter :: amutoau = 1822.8884853323708d0
      real(wp), parameter :: autoev = 27.21138602d0
      real(wp), parameter :: autokcalpermol = 627.50947377753740d0
      real(wp), parameter :: zero =  0.0d0
      real(wp), parameter :: one =   1.0d0
      real(wp), parameter :: two =   2.0d0
      real(wp), parameter :: three = 3.0d0
      real(wp), parameter :: four  = 4.0d0
      real(wp), parameter :: five  = 5.0d0
      real(wp), parameter :: half =  0.5d0
end module
