      program stateprop
C-----------------------------------------------------
C
C     Propagate a wavefunction on one potential surface
C     under the action of a laserfield in the
C     state-space-representation (i.e. in the basis of
C     the eigenfunctions of that potential) 
C
C     The time-evolution-operator   -i*H*dt    -i*(H0 - E(t)*mue)dt
C                                  e        = e
C     is evaluated using the split-operator-technique
C     (see subroutine propagate)
C
C
C     Markus Oppel, FU Berlin/Uni Wien, 1997 - 2017
C
C     Files 
C     =====
C
C     Input: 
C
C-----------------------------------------------------------------------
C     option.in
C
C     ncpu		:   Number of cpus to be used
C     ntot              :   Total number of states
C     cart              :   Cartesian component of dipol-matrix
C     tottime,timestep  :   Total timestep of propagation (integer),
C                           timestep (real*8)
C     istate		:   the initial state (=initial wavefunction)
C     fstate		:   the final state interested in
C------------------------------------------------------------------------
C     laser.in
C
C     ezero             : Amplitude of the laserfield
C     freq              : Frequency of the laserfield
C     tpulse            : Time-duration
C
C------------------------------------------------------------------------
C     hamilton.in
C
C     file containing the hamiltonmatrix (diagonal) 
C
C     do i=1,ntot
C	write (*,*) hamil(i)
C     enddo
C
C------------------------------------------------------------------------
C     dipol_[x,y,z].in
C
C     file containing the dipolmatrix-elements 
C     for the [x,y,z]-component of the dipol-moment
C
C      do i=1,ntot
C        do j=1,ntot
C        write(*,*) dipol(j,i)   
C        enddo
C      enddo
C
C------------------------------------------------------------------------
C
C      Output:
C
C     stateprop.out
C    
C     general outpufile
C
C-----------------------------------------------------------------------
      implicit none
      include 'param.inc'
    
C     global options
      integer ntot
      integer istate
      integer fstate
      real*8 tottime,timestep
C     the dipol matrix
      character*1 cart
      real*8 dipol(MAXVEC,MAXVEC)
C     the Hamiltonmatrix (diagonal !!) 
      real*8 hamil(MAXVEC)
C     the wavefunction
      complex*16 psi(MAXVEC)
 

      call setup    
      call option(ntot,cart,tottime,timestep,istate,fstate)
      call laser
      call gethamil(ntot,hamil)
      call getdip(ntot,cart,dipol)
      call propagate(ntot,hamil,dipol,tottime,timestep,
     &               istate,psi,fstate)
      call analyse(ntot,psi,fstate,hamil)
      end

