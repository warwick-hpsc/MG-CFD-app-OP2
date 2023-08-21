MPICC = cc

CC = cc

CCOPTS = -O3 -Wall

LOPTS =
CULIBS = -nocudalib -nocudainc

NVOPTS = -O3 -arch sm_20 -nocudalib -nocudainc



# Linkage rules

all : cgpupbpic2

cgpupbpic2 : cgpupbpic2.o gpulib2.o gpupbpush2.o gpupfft2.o cpbpush2.o \
             cpplib2.o cgpplib2.o dtimer.o
	$(MPICC) $(CCOPTS) $(LOPTS) -o cgpupbpic2 cgpupbpic2.o \
    gpulib2.o gpupbpush2.o gpupfft2.o cpbpush2.o cpplib2.o cgpplib2.o \
    dtimer.o $(CULIBS)


dtimer.o : dtimer.cpp
	$(CC) $(CCOPTS) -c dtimer.cpp

cgpplib2.o : gpplib2.cpp
	$(MPICC) $(CCOPTS) -o cgpplib2.o -c gpplib2.cpp

gpulib2.o : gpulib2.cu
	cc $(NVOPTS) -c gpulib2.cu -I/usr/local/cuda/include

gpupbpush2.o : gpupbpush2.cu
	cc $(NVOPTS) -c gpupbpush2.cu -I/usr/local/cuda/include

gpupfft2.o : gpupfft2.cu
	cc $(NVOPTS) -c gpupfft2.cu -I/usr/local/cuda/include

cpplib2.o : pplib2.cpp
	$(MPICC) $(CCOPTS) -o cpplib2.o -c pplib2.cpp

cpbpush2.o : pbpush2.cpp
	$(MPICC) $(CCOPTS) -o cpbpush2.o -c pbpush2.cpp

cgpupbpic2.o : gpupbpic2.cpp
	$(MPICC) $(CCOPTS) -o cgpupbpic2.o -c gpupbpic2.cpp

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f cgpupbpic2
