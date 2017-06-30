      subroutine analyse(ntot,psi,fstate,hamil)
      implicit none

      include 'param.inc'
      include 'units.inc'

      integer ntot
      integer fstate
      complex*16 psi(MAXVEC)
      real*8 hamil(MAXVEC)
      real*8 unit(MAXVEC)
      real*8 norm,energy
      real*8 expval
      integer i
      integer ifail

      external function expval

      do i=1,ntot
	unit(i)=1.0d0
      enddo

      norm=expval(ntot,unit,psi)
      energy=expval(ntot,hamil,psi)

      write (NOUT,30) norm
      write (NOUT,40) energy*TOEV
      write (NOUT,20) ABS(psi(fstate))**2
    
      open(WOUT,file='wave.out',status='unknown',
     &     IOSTAT=ifail)
      if (ifail.ne.0) then
	write (0,*)'Error while opening wave.out'
	call exit(1)
      endif

      do i=1,ntot
	write (WOUT,*) psi(i),zabs(psi(i))
      enddo
    
      close (WOUT)


      close (NOUT)
      close (LDAT)
      close (EDAT)
      close (NDAT)
      close (PDAT)
      return

20    format (1x,'Population of final state',t40,': ',d10.4)
30    format (1x'Norm of final wavefunction',t40,': ',d10.4)
40    format (1x'Final energy',t40,': ',d10.4)
      end 

