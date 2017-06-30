      subroutine gethamil(ntot,hamil)
C--------------------------------------------------------------------
C
C     Read in the initial wavefunction
C
C     File:
C
C     psi0.in
C
C     do i=1,ntot
C       write (*,*) psi0(i)
C     enddo
C
C------------------------------------------------------------------------
      implicit none
      include 'param.inc'
      include 'units.inc'

C     the hamiltonmatrix 
      integer ntot
      real*8 hamil(MAXVEC)

C     local variables 
      integer ifail
      integer i

      open (WAVE,FILE='hamil.in',STATUS='old',IOSTAT=ifail)
      	if (ifail.ne.0) then
		write (0,*)'Error while opening hamiltonmatrix'
		call exit (1)
      	endif
      	do i=1,ntot
      		read(WAVE,*) hamil(i)
C               convert to atomic units
		hamil(i)=hamil(i)*FROMEV
      	enddo
      close(WAVE)

      return
      end


      subroutine getdip(ntot,cart,dipol)
C---------------------------------------------------------------------
C
C     Read in the dipol matrix
C
C     Files:
C
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
      implicit none
      include 'param.inc'
      include 'units.inc'

C     the dipol matrix
      integer ntot
      character*1 cart
      real*8 dipol(MAXVEC,MAXVEC)

C     local variables
      integer i,j,ifail

      write (NOUT,10) 'dipol_'//cart//'.in'
      open (DIP,FILE='dipol_'//cart//'.in',STATUS='old',IOSTAT=ifail)
      	if (ifail.ne.0) then
		write (0,*) 'Error while opening dipol-file' 
		call exit (1)
      	endif
      	do i=1,ntot
		do j=1,ntot
    			read (DIP,*) dipol(j,i)    
C                       convert to atomic units
			dipol(j,i)=dipol(j,i)*FROMDEBYE
		enddo
      	enddo
      close(DIP)

10    format(/1x,'Dipol-File',t40,': ',A10)
      return
      end
