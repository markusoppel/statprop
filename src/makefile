<<<<<<< HEAD
#FOR= ifort
FOR= gfortran
#FOROPT=    -fpp -O2  -g -qopenmp -mkl
FOROPT=   -cpp -O2  
LIB=lapack
=======
FOR= gfortran
FOROPT=    -O2 
>>>>>>> bd603ec41ba8817d57c1a84d19e07957e570d4b5

OBJ= stateprop.o propagate.o option.o \
     getfile.o analyse.o      \
     setup.o laser.o



stateprop: $(OBJ) $(LASSIN2) 
	$(FOR) $(FOROPT) -l $(LIB) -o stateprop $(OBJ)   

.f.o:   param.inc units.inc
	$(FOR)  $(FOROPT)   -c $*.f 

clean:
	rm -rf *.o stateprop


