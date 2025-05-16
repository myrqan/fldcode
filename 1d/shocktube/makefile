FC = gfortran
FFLAGS = -O2

SRCS = initial_condition.f90 file_output.f90 mlw.f90 boundary.f90 arvis.f90 main.f90

OBJS = initial_condition.o file_output.o mlw.o boundary.o arvis.o main.o

EXEC = a.out

all : $(EXEC)

$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

main.o: initial_condition.mod file_output.mod mlw.mod boundary.mod arvis.mod main.f90
	$(FC) $(FFLAGS) -c main.f90

initial_condition.o initial_condition.mod: initial_condition.f90
	$(FC) $(FFLAGS) -c initial_condition.f90

file_output.o file_output.mod: file_output.f90
	$(FC) $(FFLAGS) -c file_output.f90

mlw.o mlw.mod: mlw.f90
	$(FC) $(FFLAGS) -c mlw.f90

boundary.o boudary.mod: boundary.f90
	$(FC) $(FFLAGS) -c boundary.f90

arvis.o arvis.mod: arvis.f90
	$(FC) $(FFLAGS) -c arvis.f90

clean:
	rm -f *.o *.mod $(EXEC)

