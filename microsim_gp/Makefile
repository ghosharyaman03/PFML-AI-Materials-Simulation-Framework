CC:=mpicc
CFLAGS :=-Ifunctuions -Isolverloop $(shell pkg-config --cflags hdf5)
LIBS	:=$(shell pkg-config --libs hdf5) -lm -lgsl -lgslcblas

all : microsim_gp

microsim_gp : microsim_gp.c
	$(CC) -o microsim_gp microsim_gp.c $(CFLAGS) $(LIBS)

run: microsim_gp
	mpirun -np 16 ./microsim_gp Inp.py Fill.py outpuut

.PHONY : clean

clean :
	rm microsim_gp microsim_gp.o

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/gsl-intel/lib:/opt/hdf5-intel/lib:/opt/intel/oneapi/2025.2/lib
