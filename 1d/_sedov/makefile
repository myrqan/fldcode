FC = gfortran
FFLAGS = -O2

SRCS = $(wildcard *.f90)
MODS = $(SRCS:.f90=.o)
OBJS = $(filter-out main.o, $(SRCS:.f90=.o)) main.o

MAIN = main.o

EXEC = a.out

all : $(EXEC)

$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(EXEC)

