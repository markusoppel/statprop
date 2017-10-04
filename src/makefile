FOR= ifort
FOROPT=    -fpp -O2  -g -qopenmp -mkl

OBJ= stateprop.o propagate.o option.o \
     getfile.o analyse.o      \
     setup.o laser.o



stateprop: $(OBJ) $(LASSIN2) 
	$(FOR) $(FOROPT) -o stateprop $(OBJ)   

.f.o:   param.inc units.inc
	$(FOR)  $(FOROPT)   -c $*.f 

clean:
	rm -rf *.o stateprop


