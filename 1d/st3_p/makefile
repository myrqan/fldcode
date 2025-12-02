FC     = mpifort
FFLAGS = -O2 -Wall -std=f2018
EXEC   = a.out

VISUAL = visual.py

MPI_NUM = 4
CPU_LIMIT = 50

SRCS = $(wildcard *.f90)
OBJS = $(SRCS:.f90=.o)

MOD_ROOT = vars.o
OBJ_MAIN = main.o
OBJ_SUBS = $(filter-out $(MOD_ROOT) $(OBJ_MAIN), $(OBJS))

.PHONY: all clean run vis

all: $(EXEC)

run: $(TARGET)
	mpirun -np $(MPI_NUM) cpulimit -l $(CPU_LIMIT) ./$(EXEC)

vis: visual.py
	python3 $(VISUAL)

$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $<


$(OBJ_MAIN) $(OBJ_SUBS): $(MOD_ROOT)

$(OBJ_MAIN): $(OBJ_SUBS)

clean:
	rm -f $(EXEC) *.o *.mod dat/*.dat
