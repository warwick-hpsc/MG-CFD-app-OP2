MPICC = /mnt/lustre/indy2lfs/sw/openmpi/4.1.4-cuda-11.8/bin/mpicxx

CC = g++

CCOPTS = -I/mnt/lustre/indy2lfs/sw/nvidia/hpcsdk-2211/Linux_x86_64/22.11/cuda/include -O3 -Wall

LOPTS =

CULIBS = -L/mnt/lustre/indy2lfs/sw/nvidia/hpcsdk-2211/Linux_x86_64/22.11/cuda/lib64 -lcudart -L/mnt/lustre/indy2lfs/sw/nvidia/hpcsdk-2211/Linux_x86_64/22.11/REDIST/math_libs/11.8/targets/x86_64-linux/lib -lcufft

NVOPTS = -I${MPI_HOME}/include  -L${MPI_HOME}/lib -lmpi -arch sm_70



# Linkage rules

all : gpupbpic2 lib

lib :   gpupbpic2.o gpulib2.o gpupbpush2.o gpupfft2.o cpbpush2.o \
		 cpplib2.o cgpplib2.o dtimer.o cpx_utils.o
		ar rcs libgpupbpic.a gpupbpic2.o gpulib2.o gpupbpush2.o \
			gpupfft2.o cpbpush2.o cpplib2.o cgpplib2.o dtimer.o\
			cpx_utils.o 

gpupbpic2 : execgpupbpic2.o gpulib2.o gpupbpush2.o gpupfft2.o cpbpush2.o \
             cpplib2.o cgpplib2.o dtimer.o cpx_utils.o
	$(MPICC) $(CCOPTS) $(LOPTS) -o gpupbpic2 execgpupbpic2.o \
    gpulib2.o gpupbpush2.o gpupfft2.o cpbpush2.o cpplib2.o cgpplib2.o \
    dtimer.o cpx_utils.o $(CULIBS) -lm

cpx_utils.o : cpx_utils.cpp
	$(MPICC) $(CCOPTS) -c cpx_utils.cpp

dtimer.o : dtimer.cpp
	$(CC) $(CCOPTS) -c dtimer.cpp

cgpplib2.o : gpplib2.cpp
	$(MPICC) $(CCOPTS) -o cgpplib2.o -c gpplib2.cpp

gpulib2.o : gpulib2.cu
	nvcc $(NVOPTS) -o gpulib2.o -c gpulib2.cu

gpupbpush2.o : gpupbpush2.cu
	nvcc $(NVOPTS) -o gpupbpush2.o -c gpupbpush2.cu

gpupfft2.o : gpupfft2.cu
	nvcc $(NVOPTS) -o gpupfft2.o -c gpupfft2.cu

cpplib2.o : pplib2.cpp
	$(MPICC) $(CCOPTS) -o cpplib2.o -c pplib2.cpp

cpbpush2.o : pbpush2.cpp
	$(MPICC) $(CCOPTS) -o cpbpush2.o -c pbpush2.cpp

execgpupbpic2.o : execgpupbpic2.cpp
	$(MPICC) $(CCOPTS) -o execgpupbpic2.o -c execgpupbpic2.cpp

gpupbpic2.o : gpupbpic2.cpp
	$(MPICC) $(CCOPTS) -o gpupbpic2.o -c gpupbpic2.cpp

clean :
	rm -f *.o *.mod gpupbpic2 *.a
