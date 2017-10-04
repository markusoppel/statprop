C---------------------------------------------------------
C
C     This file defines all the conversion units to/from
C     atomic units.
C
C---------------------------------------------------------
C
C     convert from atomic units
C
      real*8 TODEBYE
      real*8 TOGIGAV
      real*8 TOEV
      real*8 TOFEMTS

      parameter (TODEBYE = 2.541761D+00)
      parameter (TOGIGAV = 5.142216D+02)
      parameter (TOEV    = 2.721150D+01)
      parameter (TOFEMTS = 2.418895D-02)

C
C    convert to atomic units
C 
      real*8 FROMDEBYE
      real*8 FROMGIGAV
      real*8 FROMEV
      real*8 FROMFEMTS

      parameter (FROMDEBYE = 1.0d0/TODEBYE)
      parameter (FROMGIGAV = 1.0d0/TOGIGAV)
      parameter (FROMEV    = 1.0d0/TOEV   )
      parameter (FROMFEMTS = 1.0d0/TOFEMTS) 


