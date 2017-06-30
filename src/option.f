      subroutine option(ntot,cart,tottime,timestep,istate,fstate)
C---------------------------------------------------------------------   
C
C     read in options for the propagation of a wavefunction
C     in state-space
C
C     File: option.in
C
C     ncpu		:   Number of Cpu's to be used
C     ntot             	:   Total number of states
C     cart      	:   Cartesian component of dipol-matrix
C     tottime,timestep 	:   Total timestep of proagation (integer),
C			    timestep (real*8) 
C     istate		:   The inital state (=initial wavefunction)
C--------------------------------------------------------------------
      implicit none
      include 'param.inc'
      include 'units.inc'

C     options to be read
      integer ncpu
      integer ntot
      integer istate
      integer fstate
      character*1 cart
      integer tottime
      real*8 timestep

C     local variables
      integer ifail

C     open file and read in the options
      open(OPT,FILE='option.in',STATUS='old',IOSTAT=ifail)
      	if (ifail.ne.0) then
		write (0,*)'Error while opening option file'
		call exit(1)
      	endif
        read(OPT,*,IOSTAT=ifail) ncpu
      	read(OPT,*,IOSTAT=ifail) ntot
      	read(OPT,*,IOSTAT=ifail) cart
      	read(OPT,*,IOSTAT=ifail) tottime,timestep
	read(OPT,*,IOSTAT=ifail) istate
        read(OPT,*,IOSTAT=ifail) fstate
	if (ifail.ne.0) then
		write (0,*) 'Error while reading options'
		call exit (1)
	endif
      close (OPT)


C     check for plausibility
      if (ntot.gt.MAXVEC) then
	write (0,*) 'Error. ntot greater then MAXVEC: ',ntot,MAXVEC
	call exit (1)
      endif

C     write options to  output
      write (NOUT,5) ncpu
      write (NOUT,10) ntot
      write (NOUT,20) cart
      write (NOUT,30) timestep
      write (NOUT,40) tottime
      write (NOUT,50) istate
      write (NOUT,60) fstate
C
C     convert to atomic units
C
      timestep=timestep*FROMFEMTS

C     set up parallel-options
      call mp_set_numthreads(ncpu)

      return

5     format(1x,'Numbers of CPUs to be used ',t40,': ',i10)
10    format(1x,'Total number of states',t40,': ',I10)
20    format(1x,'Cartesian component of dipolmatrix',t40,': ',A10)
30    format(1x,'Timestep',t40,': ',d10.4)
40    format(1x,'Total number of timesteps',t40,': ',I10)
50    format(1x'Initial state (wavefunction)',t40,': ',I10)
60    format(1x,'Final state to analyse',t40,': ',I10)

      end
