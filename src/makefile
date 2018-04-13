#FOR= ifort
FOR= gfortran
#FOROPT=    -fpp -O2  -g -qopenmp -mkl
FOROPT=   -cpp -O2  
LIB=-llapack -lblas

OBJ= stateprop.o propagate.o option.o \
     getfile.o analyse.o      \
     setup.o laser.o



stateprop: $(OBJ) $(LASSIN2) 
	$(FOR) $(FOROPT) -o stateprop $(OBJ) $(LIB)

.f.o:   param.inc units.inc
	$(FOR)  $(FOROPT)   -c $*.f 

clean:
	rm -rf *.o stateprop


