FC = gfortran
FFLAGS = -O2

# List of all object files
OBJS = arvis.o boundary2d.o calculate_variables.o cfl.o file_output.o mlw.o model.o main.o

all: main

# Executable depends on all object files
main: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS)

# Specific dependency: model.o depends on boundary2d.o
model.o: model.f90 file_output.o
	$(FC) $(FFLAGS) -c model.f90

# main.o depends on all other object files
main.o: main.f90 arvis.o boundary2d.o calculate_variables.o cfl.o file_output.o mlw.o model.o
	$(FC) $(FFLAGS) -c main.f90

# Generic rule for compiling .f90 files into .o files
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Clean up compiled files
clean:
	rm -f *.o main *.mod *.dat


