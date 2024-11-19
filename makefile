FC      = gfortran
OLEVEL  = -O2
FOPTS   = -mcmodel=medium -fopenmp -Wall -Wextra -Werror -pedantic -fcheck=all,no-array-temps, -ffpe-trap=invalid,zero,overflow -ffpe-summary=all -g -fbacktrace
FFLAGS	= $(OLEVEL) $(FOPTS)

SOURCE  = module.f90   \
		  main.f90 \
		  init.f90 \
          eval.f90 \
		  output.f90   \
          bulb.f90 \
          pde.f90  \
		  bubble.f90   \


OBJECTS = ${SOURCE:.f90=.o}

EXEC	= run.exe

# Rule to compile .f90 to .o
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

# Default target to link the program together 
run: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS)
#	mkdir -p output
#	rm -f *.o

# Subtargets
clean:
	@del *.o
	@del .\output\*.dat
	@del run.exe
	@echo Completed makefile
#	rm -f *.o run

