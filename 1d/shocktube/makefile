FC = gfortran
FFLAGS = -O2

SRCS = initial_condition.f90 file_output.f90 main.f90

OBJS = initial_condition.o file_output.o main.o

EXEC = calc

all : $(EXEC)

$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

main.o: initial_condition.mod file_output.mod main.f90
	$(FC) $(FFLAGS) -c main.f90

initial_condition.o initial_condition.mod: initial_condition.f90
	$(FC) $(FFLAGS) -c initial_condition.f90

file_output.o file_output.mod: file_output.f90
	$(FC) $(FFLAGS) -c file_output.f90

clean:
	rm -f *.o *.mod $(EXEC)

