
OBJECTS=main.o \
        init.o \
        eval.o \
        bulb.o \
        bicg.o \
        pois.o \
        pde.o

run: $(OBJECTS)
	gfortran $(OBJECTS) -O2 -mcmodel=medium -fopenmp -o run
	rm -f *.o 
.f.o:
	gfortran -c -O2 -mcmodel=medium -fopenmp $<

