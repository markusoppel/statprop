FOR= gfortran
FOROPT=    -O2 

OBJ= stateprop.o propagate.o option.o \
     getfile.o analyse.o      \
     setup.o laser.o



stateprop: $(OBJ) $(LASSIN2) 
	$(FOR) $(FOROPT) -o stateprop $(OBJ)   

.f.o:   param.inc units.inc
	$(FOR)  $(FOROPT)   -c $*.f 

clean:
	rm -rf *.o statesin2 stategaus


