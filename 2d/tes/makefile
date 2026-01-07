FC     = mpifort
FFLAGS = -O2 -std=f2018
FFLAGS += -Wall
EXEC   = a.out

VISUAL = visual.py

MPI_NUM = 6
CPU_LIMIT = 50

SRCS = $(wildcard *.f90)
OBJS = $(SRCS:.f90=.o)

MOD_ROOT = vars.o const.o
MOD_ROOT2 = datw.o bnd2d.o
OBJ_MAIN = main.o
OBJ_SUBS = $(filter-out $(MOD_ROOT) $(MOD_ROOT2) $(OBJ_MAIN), $(OBJS))

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


$(OBJ_MAIN) $(OBJ_SUBS): $(MOD_ROOT) $(MOD_ROOT2)

$(MOD_ROOT2): $(MOD_ROOT)

$(OBJ_MAIN): $(OBJ_SUBS)

clean:
	rm -f $(EXEC) *.o *.mod dat/*.dat graph/*.png
	rm -rf __pycache__/
	rm -f param.txt
