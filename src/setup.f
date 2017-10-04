      subroutine setup
C---------------------------------------------------------------------
C
C     open the outputfile, set up date adn timing
C     
C     Files:
C
C     stateprop.out
C
C---------------------------------------------------------------------    
      implicit none
      include 'param.inc'

C     date and time 
      integer imon,iday,iyear
      integer iarray(3)

C     local variables
      integer ifail

      call idate(imon,iday,iyear)
      call itime(iarray)

      open (NOUT,file='stateprop.out',status='unknown',
     &      IOSTAT=ifail)
      if (ifail.ne.0) then
        write (0,*) 'Error while opening output'
        call exit (1)
      endif

      write (NOUT,10)
      write (NOUT,20) iday,imon,iyear,iarray

C 
C     open outputfiles
C
      open(EDAT,file='energy.dat',status='unknown',
     &     IOSTAT=ifail)
      if (ifail.ne.0) then
	write (0,*)'Error while opening energydat'
	call exit(1)
      endif
      open(LDAT,file='laser.dat',status='unknown',
     &     IOSTAT=ifail)
      if (ifail.ne.0) then
	write (0,*)'Error while opening laser.dat'
	call exit(1)
      endif
      open(NDAT,file='norm.dat',status='unknown',
     &     IOSTAT=ifail)
      if (ifail.ne.0) then
	write (0,*)'Error while opening norm.dat'
	call exit(1)
      endif
      open(PDAT,file='pop.dat',status='unknown',
     &     IOSTAT=ifail)
      if (ifail.ne.0) then
	write (0,*)'Error while opening pop.dat'
	call exit(1)
      endif


      return
      
10    format(/1x,'!! Welcome to Stateprop !!',/
     &        1x,'==========================',//
     &        1x,'(c) M. Oppel, 1997 ')
20    format(/1x,'Started at',t40,': ',
     &        I2,'.',I2,'.',I2,', '     
     &        I2,':',I2,':',I2,/)

      end

