      subroutine laser
C--------------------------------------------------------------
C     read the parameters of the laserfield
C     from file laser.in
C
C     File: laser.in
C
C     ezero1,ezero2,estep	: Amplitude of the laserfield
C                                 (start,stop, step)
C     freq1,freq2,fstep		: Frequency of the laserfield
C                   		  (start,stop,step)
C     tpulse			: Time-duration
C
C-------------------------------------------------------------
      implicit none
      include'param.inc'
      include'units.inc'

C     local variables
      integer ifail

C     the commonblock for the laserfield
      real*8 ezero, freq,tpulse
      common/lasfield/ezero,freq,tpulse

      open (LAS,file='laser.in',STATUS='old',IOSTAT=ifail)
      	if (ifail.ne.0) then
		write (0,*) 'Error while opening laser.in'
		call exit (1)
      	endif
      	read (LAS,*,IOSTAT=ifail) ezero
      	read (LAS,*,IOSTAT=ifail) freq 
      	read (LAS,*,IOSTAT=ifail) tpulse
	if (ifail.ne.0) then
		write (0,*) 'Error while reading laser-paramters'
		call exit (1)
	endif
      close (LAS)

C     write laserparmeters to output

      write (NOUT,10) 
      write (NOUT,20) ezero
      write (NOUT,30) freq
      write (NOUT,40) tpulse
C
C     convert to atomic units
C
      ezero=ezero*FROMGIGAV
      freq=freq*FROMEV
      tpulse=tpulse*FROMFEMTS


      return

10    format (//1x,'Laser Parameters:',/
     &          1x,'=================',/)
20    format (1x,'Amplitude of laserfield',t40,': ',d10.4)
30    format (1x,'Frequency of laserfield',t40,': ',d10.4)
40    format (1x,'Time-duration',t40,': ',d10.4)

      end


      function lasval(time)
C------------------------------------------------------
C
C     this function returns the value of the laserfield
C     at time t, assuming a sin**2-pulse.
C
C     E(t)=E0*cos(wt)*sin**2(pi*t/tpulse)
C     the laser-parameters are stored in the commonblock
C     laserfield.
C     
C     ezero:	Amplitude of the laserfield
C     freq:	Frequency of the laserfield
C     tpulse:	pulse-duration
C 
C------------------------------------------------------------
      implicit none

      include'param.inc'
      real*8 time,lasval
      real*8 pi

      real*8 ezero, freq,tpulse
      common/lasfield/ezero,freq,tpulse

      pi=acos(-1.0d0)
      lasval =  ezero * dcos (freq*time) *
     &          (dsin(pi*time/tpulse))**2
      return
      end
