MPICC = cc

CC = cc

CCOPTS = -O3 -fPIC
LOPTS = -lm

bin : cpbpic2.o cpbpush2.o cpplib2.o dtimer.o cpx_utils.o
	$(MPICC) $(CCOPTS) $(LOPTS) -o cpbpic2 \
        cpbpic2.o cpbpush2.o cpplib2.o dtimer.o cpx_utils.o

lib : cpbpic2_lib.o cpbpush2.o cpplib2.o dtimer.o cpx_utils.o
	ar rcs libpbpic.a cpbpic2_lib.o cpbpush2.o \
		cpplib2.o dtimer.o cpx_utils.o 


# Compilation rules

cpx_utils.o: cpx_utils.cpp
	$(MPICC) $(CCOPTS) -c cpx_utils.cpp

dtimer.o : dtimer.cpp
	$(CC) $(CCOPTS) -c dtimer.cpp

cpplib2.o : pplib2.cpp
	$(MPICC) $(CCOPTS) -o cpplib2.o -c pplib2.cpp

cpbpush2.o : pbpush2.cpp
	$(MPICC) $(CCOPTS) -o cpbpush2.o -c pbpush2.cpp

cpbpic2.o : pbpic2.cpp
	$(MPICC) $(CCOPTS) -o cpbpic2.o -c pbpic2.cpp

cpbpic2_lib.o : pbpic2_lib.cpp
	$(MPICC) $(CCOPTS) -o cpbpic2_lib.o -c pbpic2_lib.cpp

clean :
	rm -f *.o *.mod libpbpic.a cpbpic2
