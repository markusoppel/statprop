C     
C     the sizes   
C
      integer MAXVEC
      parameter (MAXVEC=128)

C
C     the file-handles
C
      integer OPT,LAS,WAVE,DIP,NOUT
      integer EDAT,LDAT,NDAT,PDAT
      integer WOUT
      parameter (OPT  = 10)
      parameter (LAS  = 11)
      parameter (WAVE = 12)
      parameter (DIP  = 13)
      parameter (NOUT = 20)
      parameter (EDAT = 30)
      parameter (NDAT = 31)
      parameter (LDAT = 32)
      parameter (PDAT = 33)
      parameter (WOUT = 40)
C
C     the constants
C

      complex*16 CZERO,CONE,CIM
      parameter (CZERO=(0.0d0,0.0d0),CONE=(1.0d0,0.0d0),
     &           CIM=(0.0d0,1.0d0))

      


