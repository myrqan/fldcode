FC = mpifort
FFLAGS = -O2 -Wall -std=f2018
TARGET = main
SRCS= main.f90
OBJS = $(SRCS:.f90=.o)

MPI_NUM = 8
CPU_LIMIT = 50

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<


run: $(TARGET)
	mpirun -np $(MPI_NUM) cpulimit -l $(CPU_LIMIT) ./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean run
