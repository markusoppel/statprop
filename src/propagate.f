      subroutine propagate(ntot,hamil,dipol,ttime,dt,istate,psi,fstate)
C-----------------------------------------------------------------------
C
C     propagate a wavefunction in state-space, using 
C     the split-operator technique
C
C        -iHdt          -i(H0-mue*E(t))
C      e       |psi> = e                |psi>  =
C
C        -iH0*dt/2    i*mue*E(t))*dt    -iH0*dt/2
C      e           * e                 * e          |psi>
C
C----------------------------------------------------------------------
      implicit none
      include 'param.inc'
      include 'units.inc'
C----------------------------------------------------------------------
C
C     the calling options
C
C----------------------------------------------------------------------
C     the hamiltionmatrix 
      integer ntot
      real*8 hamil(MAXVEC)

C     the dipolmatrix
      real*8 dipol(MAXVEC,MAXVEC)
      real*8 pot(MAXVEC,MAXVEC)
C     the unit-matrix
      real*8 unit(MAXVEC)

C     the timing
      real*8 dt
      integer ttime

C     the initial and final state
      integer istate,fstate
C------------------------------------------------------------------------
C
C     local variables
C
C------------------------------------------------------------------------ 

C     the wavefunction and the exponential factors
      complex*16 psi(MAXVEC)       ! diabatic representation
      complex*16 adia(MAXVEC)      ! adiabatic representation
      complex*16 epsi(MAXVEC)      ! exp(-iH0dt)
      complex*16 ehalfpsi(MAXVEC)  ! exp(-iH0dt/2)

C     the dipolmatrix in the adiabatic representation
C     and the exponential factors
      real*8 dipdiag(MAXVEC)   ! in adiabatic representation
      complex*16 potfac(MAXVEC)    ! exp(iE(t)*mue) (adiabatic)

C     the tranformation matrix diabatic <-> adiabatic
      real*8 trans(MAXVEC,MAXVEC)
      real*8 inverse(MAXVEC,MAXVEC)

C     the laserfield
      real*8 efield
      real*8 lasval

C     the workspace for the diagonalization-routine
      integer lwork
      parameter(lwork=3*MAXVEC)
      real*8 work(lwork)

C     the timing
      real*8 time
      integer itime
C
C     timing infomation
      real tarray (2)
      real stime
      real dtime

C     help variables
      integer i,j,ifail,k

C     expectation values
C
      real*8 norm,energy
      real*8 expval

      complex*16 help(MAXVEC,MAXVEC)
C     the external subroutines and functions
C
C     dsyev:  diagonalization of a real symmetric matrix, 
C	      taken from lapack
C
C     zhemv:  hermitian matrix times complex vector,
C	      taken from blas3
C
C     dlincg: calculate the inverse of a real matrix,
C	      taken from imsl
C
C     lasval: function which gives the value of the laserfield
C             at time t (i.e. lasval=lasval(t))
C
C     expval: function to calculate the expectation-value of a 
C             wavefunction with an Operator
C
C     dtime: used for timing information, see man dtime
C
      external function dtime 
      external function lasval
      external function expval
     
C-----------------------------------------------------------------
C
C     Here starts the work
C
C-----------------------------------------------------------------
C
C     Start timing
      stime=dtime(tarray)
C
C     First create the exponential factors, which form the 
C     constant part of the Time-evolution-operator:
C     exp(-iH0dt) and exp(-iH0dt/2)
C     also create the unit-matrix
C
      do i=1,ntot 
	unit(i)=1.0d0
	epsi(i)     = exp(-(0.0d0,1.0d0)*hamil(i)*dt)
	ehalfpsi(i) = exp(-(0.0d0,1.0d0)*hamil(i)*dt/2)
      enddo
C
C     setup the initial wavefunction, 
C     calculate initial norm and energy
C
      do i=1,ntot
        psi(i)=CZERO
      enddo
      psi(istate)=CONE

      norm   = dsqrt(expval(ntot,unit ,psi))
      energy =       expval(ntot,hamil,psi)
      
      write (NOUT,*) 
      write (NOUT,20) norm
      write (NOUT,30) energy*TOEV 
      write (NOUT,*) 
      write (NOUT,40)
C
C     prepare initial wavefunction by multiplying exp(-iH0dt/2)
C  
      call vtv(ntot,ehalfpsi,psi)
C
C       calculate the transformation-matrix
C
	do i=1,ntot
		do j=1,ntot
			pot(i,j)=dipol(i,j)
		enddo
	enddo
        call devcsf(ntot,pot,MAXVEC,dipdiag,trans,MAXVEC)
C
C       calculate the inverse of the transformation matrix
C
        call DLINRG(ntot,trans,MAXVEC,inverse,MAXVEC)

C
C    here starts the actual propagation
C
      do itime=1,ttime-1
C  
C       setup the time and the laserfield
C
	time=real(itime*dt)
        efield=lasval(time)
C
C       calculate the exponential factors exp(i*mue*E(t)) 
C
        call exppot(ntot,efield,dipdiag,potfac,dt)
C
C       change to the adiabatic representation
C       and let the mue*E(T) of the time-evolution-operator act
C       on the wavefuntion
C
	
        call mtv(ntot,inverse,psi,adia,1.0d0)
        call vtv(ntot,potfac,adia)
C
C       switch back to the diabatic picture and let
C       the H0 part of the time-evolution-operator act
C 
        call mtv(ntot,trans,adia,psi,1.0d0)
        call vtv(ntot,ehalfpsi,psi)
C
C       calculate the norm and the energy
        norm=dsqrt(expval(ntot,unit,psi))
        energy=expval(ntot,hamil,psi)
        
        write (PDAT,*) time*TOFEMTS,zabs(psi(fstate))
        write (NDAT,*) time*TOFEMTS,norm
	write (EDAT,*) time*TOFEMTS,energy*TOEV
        write (LDAT,*) time*TOFEMTS,efield*TOGIGAV 
C
C       prepare the next timestep
C
        call vtv(ntot,ehalfpsi,psi)
C
C     end of loop over timesteps
C
      enddo

C     the last timestep
C     we finally have to multiply exp(i*mue*E(t))*exp(-iH0dt/2) 
      itime=ttime 
      time=real(itime*dt)
      efield=lasval(time)
      call exppot(ntot,efield,dipdiag,potfac,dt)
C
C       change to the adiabatic representation
C       and let the mue*E(T) of the time-evolution-operator act
C       on the wavefuntion
C
      call mtv(ntot,inverse,psi,adia,1.0d0)
      call vtv(ntot,potfac,adia)
C
C       switch back to the diabatic picture and let
C       finally exp(-iH0dt/2) act on the wavefunction
C 
      call mtv(ntot,trans,adia,psi,1.0d0)
      call vtv(ntot,ehalfpsi,psi)
C
C     stop timing
      stime=dtime(tarray)
      write (NOUT,100) 
      write (NOUT,10) tarray,stime
      write (NOUT,50)
    
      return

100   format (/1x,'Finished !!')
10    format(/1x,'Usertime',t40,': ',F10.4,' sec.',/
     &        1x,'Sys-Time',t40,': ',F10.4,' sec.',/
     &        1x,'Tot.Time',t40,': ',F10.4,' sec.')
50    format (/1x,'==========================================',/)

20    format (1x,'Norm of initial wavefunction',t40,': ',d10.4)
30    format (1x,'Initial energy',t40,': ',d10.4)
40    format (1x,'===========================================',//
     &        1x,'Starting')
      end


      subroutine exppot(ntot,efield,diag,fac,dt)
C-------------------------------------------------------------
C
C     calculate the exponential factors 
C     exp(i*E(t)*mue*dt)
C
C-------------------------------------------------------------
      implicit none
      include 'param.inc'

C     the dipol-matrix in adiabatic representation
      integer ntot
      real*8 diag(MAXVEC)

C     the exponential factor (return-value of the subroutine)
      complex*16 fac(MAXVEC)

C     the laserfield and the time-step
      real*8 efield, dt

C     local variables
      integer i

      do i=1,ntot
	fac(i)=exp((0.0d0,1.0d0)*efield*diag(i)*dt)
      enddo

      return
      end

      subroutine mtv(ntot,A,V,W,fac)
C-----------------------------------------------------------
C
C     subroutine for (complex) Matrix * Vector Multiplicition
C     w(*) = A(*,*) * v(*)
C
C-----------------------------------------------------------
      implicit none

      include 'param.inc'
      integer ntot
      real*8 A(MAXVEC,MAXVEC)       
      real*8 fac
      complex*16 V(MAXVEC)
      complex*16 W(MAXVEC)

      integer i,j

      do i=1,ntot
	w(i)=0.0d0
	do j=1,ntot
          w(i)=w(i)+dcmplx(a(i,j))*v(j)*fac
	enddo
      enddo

      return
      end	  		


      real*8 function expval(ntot,Operator,w)
C------------------------------------------------
C
C     function to calculate the expectation-value of 
C     a wavefunction with an operator
C
C              ntot   *
C     expval = Sum w(i) * Operator(i) * w(i)
C              i=1 
C
C     assuming the operator to be diagonal 
C
C---------------------------------------------------------
      implicit none
      include 'param.inc'

      integer ntot
      real*8 Operator(MAXVEC)
      complex*16 w(MAXVEC)

      real*8 help 
      integer i

      help=0.0d0
      do i=1,ntot
    	help=help+Operator(i)*zabs(w(i))**2
      enddo

      expval=help

      return
      end
    
      subroutine vtv(ntot,Operator,psi)
C-------------------------------------------------------------
C
C     subroutine to evaluate the action of a (diagonal) operator
C     on an wavefunction
C 
C     psi = Operator * psi 
C
C-------------------------------------------------------------------
      implicit none
      include 'param.inc'

      integer ntot
      complex*16 Operator(MAXVEC)
      complex*16 psi(MAXVEC)

      integer i

      do i=1,ntot
	psi(i)=Operator(i)*psi(i)
      enddo

      return
      end
  
