#
# set paths for header files and libraries
#
# Some compilers/linkers are happy searching LD_LIBRARY_PATH var.
# For others (cough PGI), need to explicitly state each library
# directory.
#


ifdef OP2_INSTALL_PATH
  OP2_INC = -I$(OP2_INSTALL_PATH)/c/include
  OP2_LIB = -L$(OP2_INSTALL_PATH)/c/lib
endif

ifdef CUDA_INSTALL_PATH
  CUDA_INC = -I$(CUDA_INSTALL_PATH)/include
  CUDA_LIB = -L$(CUDA_INSTALL_PATH)/lib64
endif

ifdef HDF5_INSTALL_PATH
  HDF5_INC = -I$(HDF5_INSTALL_PATH)/include
  HDF5_LIB = -L$(HDF5_INSTALL_PATH)/lib
endif
HDF5_LIB += -lhdf5 -lz

PARMETIS_VER=4
ifdef PARMETIS_INSTALL_PATH
  PARMETIS_INC = -I$(PARMETIS_INSTALL_PATH)/include
  PARMETIS_LIB = -L$(PARMETIS_INSTALL_PATH)/lib
endif
PARMETIS_INC += -DHAVE_PARMETIS
PARMETIS_LIB += -lparmetis -lmetis
ifeq ($(PARMETIS_VER),4)
  PARMETIS_INC += -DPARMETIS_VER_4
endif

ifdef PTSCOTCH_INSTALL_PATH
  PTSCOTCH_INC 	= -I$(PTSCOTCH_INSTALL_PATH)/include
  PTSCOTCH_LIB 	= -L$(PTSCOTCH_INSTALL_PATH)/lib
endif
PTSCOTCH_INC += -DHAVE_PTSCOTCH
PTSCOTCH_LIB += -lptscotch -lscotch -lptscotcherr



ifdef DEBUG
  OPTIMISE := -pg -g -O0
else
  OPTIMISE := -O3
  # OPTIMISE := -O1
  #OPTIMISE := -O0
endif

BIN_DIR = bin
OBJ_DIR = obj
SRC_DIR = src


#
# Locate MPI compilers:
#
ifdef MPI_INSTALL_PATH
  ifneq ("","$(wildcard $(MPI_INSTALL_PATH)/bin/mpicxx)")
    MPICPP := $(MPI_INSTALL_PATH)/bin/mpicxx
  else
  ifneq ("","$(wildcard $(MPI_INSTALL_PATH)/intel64/bin/mpicxx)")
    MPICPP := $(MPI_INSTALL_PATH)/intel64/bin/mpicxx
  else
    MPICPP := mpicxx
  endif
  endif

  ifneq ("","$(wildcard $(MPI_INSTALL_PATH)/bin/mpicc)")
    MPICC := $(MPI_INSTALL_PATH)/bin/mpicc
  else
  ifneq ("","$(wildcard $(MPI_INSTALL_PATH)/intel64/bin/mpicc)")
    MPICC := $(MPI_INSTALL_PATH)/intel64/bin/mpicc
  else
    MPICC := mpicc
  endif
  endif
else
  MPICPP := mpicxx
  MPICC  := mpicc
endif

ifeq ($(OP2_COMPILER),gnu)
  CPP := g++
  CFLAGS	= -fPIC -DUNIX -Wall -Wextra
  CPPFLAGS 	= $(CFLAGS)
  OMPFLAGS 	= -fopenmp
  MPIFLAGS 	= $(CPPFLAGS)
else
ifeq ($(OP2_COMPILER),intel)
  CPP = icpc
  CFLAGS = -DMPICH_IGNORE_CXX_SEEK -restrict -fno-alias -inline-forceinline -parallel -DVECTORIZE -qopt-report=5
  #CFLAGS += -fmax-errors=1
  CPPFLAGS = $(CFLAGS)
  OMPFLAGS = -qopenmp
  OMPOFFLOAD = -qopenmp
  # NVCCFLAGS += -ccbin=$(MPICPP)
  MPIFLAGS	= $(CPPFLAGS)
  ifdef ISET
    OPTIMISE += -x$(ISET)
  else
    OPTIMISE += -xHost
  endif
else
ifeq ($(OP2_COMPILER),xl)
  CPP		 = xlc++
  CFLAGS	 = -qarch=pwr8 -qtune=pwr8 -qhot
  CPPFLAGS 	 = $(CFLAGS)
  OMPFLAGS	 = -qsmp=omp -qthreaded
  OMPOFFLOAD = -qsmp=omp -qoffload -Xptxas -v -g1
  MPIFLAGS	 = $(CPPFLAGS)
else
ifeq ($(OP2_COMPILER),pgi)
  CPP       	= pgc++
  CFLAGS  	=
  CPPFLAGS 	= $(CFLAGS)
  OMPFLAGS 	= -mp
  MPIFLAGS 	= $(CPPFLAGS)
  # NVCCFLAGS	+= -ccbin=$(MPICPP)
  # ACCFLAGS      = -acc -Minfo=acc -ta=tesla:cc35 -DOPENACC
  # ACCFLAGS      = -acc -DOPENACC -Minfo=acc
  ACCFLAGS      = -v -acc -DOPENACC -Minfo=acc
else
ifeq ($(OP2_COMPILER),cray)
  CPP           = CC
  CFLAGS        = -h fp3 -h ipa5
  CPPFLAGS      = $(CFLAGS)
  OMPFLAGS      = -h omp
  MPICPP        = CC
  MPIFLAGS      = $(CPPFLAGS)
else
# print:
# 	@echo "unrecognised value for OP2_COMPILER"
$(error unrecognised value for OP2_COMPILER)
endif
endif
endif
endif
endif

ifdef CPP_WRAPPER
  CPP := $(CPP_WRAPPER)
endif
ifdef MPICPP_WRAPPER
  MPICPP := $(MPICPP_WRAPPER)
endif

#
# set flags for NVCC compilation and linking
#
ifndef NV_ARCH
  MESSAGE=select an NVIDA device to compile in CUDA, e.g. make NV_ARCH=KEPLER
  NV_ARCH=Pascal
endif
ifeq ($(NV_ARCH),Fermi)
  CODE_GEN_CUDA=-arch=sm_20
else
ifeq ($(NV_ARCH),Kepler)
  CODE_GEN_CUDA=-gencode arch=compute_35,code=sm_35
else
ifeq ($(NV_ARCH),Maxwell)
  CODE_GEN_CUDA=-gencode arch=compute_50,code=sm_50
else
ifeq ($(NV_ARCH),Pascal)
  CODE_GEN_CUDA=-gencode arch=compute_60,code=sm_60
endif
endif
endif
endif

NVCCFLAGS += $(CODE_GEN_CUDA) -m64 -Xptxas -dlcm=ca -Xptxas=-v -use_fast_math -O3


MGCFD_INCS := -Isrc -Isrc/Kernels
ifdef PAPI
  MGCFD_INCS += -DPAPI
  MGCFD_LIBS := -lpapi -lpfm
endif


## Enable VERIFY_OP2_TIMING to perform timing measurements external to
## those performed by OP2 internally. Intended to verify whether OP2 timers
## are correct, particularly for MPI sync time.
# MGCFD_INCS += -DVERIFY_OP2_TIMING


OP2_MAIN_SRC = $(SRC_DIR)_op/euler3d_cpu_double_op.cpp


# all: mgcfd_seq mgcfd_openmp mgcfd_openacc mgcfd_openmp4 mgcfd_cuda \
# 	 mgcfd_mpi mgcfd_mpi_openmp mgcfd_mpi_cuda

#all: mgcfd_seq mgcfd_openmp mgcfd_mpi mgcfd_mpi_openmp mgcfd_cuda

all: mgcfd_mpi_vec 

OP2_SEQ_OBJECTS := $(OBJ_DIR)/mgcfd_seq_main.o \
                   $(OBJ_DIR)/mgcfd_seq_kernels.o

OP2_MPI_OBJECTS := $(OBJ_DIR)/mgcfd_mpi_main.o \
                   $(OBJ_DIR)/mgcfd_mpi_kernels.o

OP2_MPI_VEC_OBJECTS := $(OBJ_DIR)/mgcfd_mpi_vec_main.o \
                       $(OBJ_DIR)/mgcfd_mpi_vec_kernels.o

OP2_OMP_OBJECTS := $(OBJ_DIR)/mgcfd_openmp_main.o \
                   $(OBJ_DIR)/mgcfd_openmp_kernels.o

OP2_MPI_OMP_OBJECTS := $(OBJ_DIR)/mgcfd_mpi_openmp_main.o \
                       $(OBJ_DIR)/mgcfd_mpi_openmp_kernels.o

OP2_CUDA_OBJECTS := $(OBJ_DIR)/mgcfd_cuda_main.o \
                    $(OBJ_DIR)/mgcfd_kernels_cu.o

OP2_MPI_CUDA_OBJECTS := $(OBJ_DIR)/mgcfd_mpi_cuda_main.o \
                        $(OBJ_DIR)/mgcfd_mpi_kernels_cu.o

OP2_OMP4_OBJECTS := $(OBJ_DIR)/mgcfd_omp4_main.o \
                    $(OBJ_DIR)/mgcfd_omp4_kernel_funcs.o \
                    $(OBJ_DIR)/mgcfd_omp4_kernels.o

OP2_OPENACC_OBJECTS := $(OBJ_DIR)/mgcfd_openacc_main.o \
                       $(OBJ_DIR)/mgcfd_openacc_kernels.o



## User-friendly wrappers around actual targets:
mgcfd_seq: $(BIN_DIR)/mgcfd_seq
mgcfd_openmp: $(BIN_DIR)/mgcfd_openmp
mgcfd_mpi: $(BIN_DIR)/mgcfd_mpi
mgcfd_mpi_vec: $(BIN_DIR)/mgcfd_mpi_vec
mgcfd_mpi_openmp: $(BIN_DIR)/mgcfd_mpi_openmp
mgcfd_cuda: $(BIN_DIR)/mgcfd_cuda
mgcfd_openmp4: $(BIN_DIR)/mgcfd_openmp4
mgcfd_openacc: $(BIN_DIR)/mgcfd_openacc
mgcfd_mpi_cuda: $(BIN_DIR)/mgcfd_mpi_cuda


## SEQUENTIAL
$(OBJ_DIR)/mgcfd_seq_main.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OPTIMISE) $(MGCFD_INCS) \
	    $(OP2_INC) $(HDF5_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) \
		-c -o $@ $(OP2_MAIN_SRC)
$(OBJ_DIR)/mgcfd_seq_kernels.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OPTIMISE) $(MGCFD_INCS) \
	    $(OP2_INC) $(HDF5_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) \
		-c -o $@ $(SRC_DIR)/../seq/_seqkernels.cpp
$(BIN_DIR)/mgcfd_seq: $(OP2_SEQ_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OPTIMISE) $(MGCFD_LIBS) $^ \
		-lm $(OP2_LIB) -lop2_seq -lop2_hdf5 $(HDF5_LIB) $(PARMETIS_LIB) $(PTSCOTCH_LIB) \
		-o $@


## OPENMP
$(OBJ_DIR)/mgcfd_openmp_main.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_INCS) \
		$(OP2_INC) $(HDF5_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) \
		-c -o $@ $(OP2_MAIN_SRC)
$(OBJ_DIR)/mgcfd_openmp_kernels.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_INCS) \
		$(OP2_INC) $(HDF5_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) \
		-c -o $@ $(SRC_DIR)/../openmp/_kernels.cpp
$(BIN_DIR)/mgcfd_openmp: $(OP2_OMP_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $^ $(OPTIMISE) $(MGCFD_LIBS) \
		-lm $(OP2_LIB) -lop2_openmp -lop2_hdf5 $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) \
		-o $@


## MPI
$(OBJ_DIR)/mgcfd_mpi_kernels.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	     -DMPI_ON \
	    -c -o $@ $(SRC_DIR)/../seq/_seqkernels.cpp
$(OBJ_DIR)/mgcfd_mpi_main.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $^ $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	     -DMPI_ON \
	    -c -o $@ $(OP2_MAIN_SRC)
$(BIN_DIR)/mgcfd_mpi: $(OP2_MPI_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $^ $(OPTIMISE) $(MGCFD_LIBS) \
		-lm $(OP2_LIB) -lop2_mpi $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) \
		-o $@


## MPI_VEC
$(OBJ_DIR)/mgcfd_mpi_vec_kernels.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
        -DMPI_ON \
        -c -o $@ $(SRC_DIR)/../vec/_veckernels.cpp
$(OBJ_DIR)/mgcfd_mpi_vec_main.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $^ $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
        -DMPI_ON \
        -c -o $@ $(OP2_MAIN_SRC)
$(BIN_DIR)/mgcfd_mpi_vec: $(OP2_MPI_VEC_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $^ $(OPTIMISE) $(MGCFD_LIBS) \
        -lm $(OP2_LIB) -lop2_mpi $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) \
        -o $@


## MPI + OPENMP
$(OBJ_DIR)/mgcfd_mpi_openmp_kernels.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -DMPI_ON \
	    -c -o $@ $(SRC_DIR)/../openmp/_kernels.cpp
$(OBJ_DIR)/mgcfd_mpi_openmp_main.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $^ $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -DMPI_ON \
	    -c -o $@ $(OP2_MAIN_SRC)
$(BIN_DIR)/mgcfd_mpi_openmp: $(OP2_MPI_OMP_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $^ $(OPTIMISE) $(MGCFD_LIBS) \
		-lm $(OP2_LIB) -lop2_mpi $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) \
		-o $@

## CUDA
$(OBJ_DIR)/mgcfd_kernels_cu.o:
	mkdir -p $(OBJ_DIR)
	nvcc $(NVCCFLAGS) $(MGCFD_INCS) $(OP2_INC) \
		-c -o $@ $(SRC_DIR)/../cuda/_kernels.cu
$(OBJ_DIR)/mgcfd_cuda_main.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CFLAGS) $^ $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -DCUDA_ON \
	    -c -o $@ $(OP2_MAIN_SRC)
$(BIN_DIR)/mgcfd_cuda: $(OP2_CUDA_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CFLAGS) $^ $(OPTIMISE) $(MGCFD_LIBS) \
	    $(CUDA_LIB) -lcudart $(OP2_LIB) -lop2_cuda $(HDF5_LIB) -lop2_hdf5 \
	    -o $@


## OPENMP4
$(OBJ_DIR)/mgcfd_omp4_kernel_funcs.o:
	mkdir -p $(OBJ_DIR)
	$(CPP) $(CPPFLAGS) $(OMPOFFLOAD) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -Iopenmp4/ -c -o $@ $(SRC_DIR)/../openmp4/_omp4kernel_funcs.cpp
$(OBJ_DIR)/mgcfd_omp4_kernels.o:
	mkdir -p $(OBJ_DIR)
	$(CPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -Iopenmp4/ -c -o $@ $(SRC_DIR)/../openmp4/_omp4kernels.cpp
$(OBJ_DIR)/mgcfd_omp4_main.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $^ $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -Iopenmp4/ -c -o $@ $(OP2_MAIN_SRC)
$(BIN_DIR)/mgcfd_openmp4: $(OP2_OMP4_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPOFFLOAD) $^ $(OPTIMISE) $(MGCFD_LIBS) \
	    $(OP2_LIB) -lop2_openmp4 $(CUDA_LIB) -lcudart \
		$(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) -lop2_hdf5 \
	    -o $@


## OPENACC
# Compile failing: nvlink unable to link _acckernels.o to extern variable definitions
$(OBJ_DIR)/mgcfd_openacc_kernels.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(ACCFLAGS) $(OMPFLAGS) $(OPTIMISE) $^ \
	    $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) -Iopenacc \
	    -c -o $@ $(SRC_DIR)/../openacc/_acckernels.c
$(OBJ_DIR)/mgcfd_openacc_main.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(ACCFLAGS) $(OMPFLAGS) $(OPTIMISE) $^ \
	    $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) -Iopenacc \
	    -c -o $@ $(OP2_MAIN_SRC)
$(BIN_DIR)/mgcfd_openacc: $(OP2_OPENACC_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(ACCFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_LIBS) $^ \
	    $(CUDA_LIB) -lcudart $(OP2_LIB) -lop2_cuda $(HDF5_LIB) -lop2_hdf5 \
	    -o $@


## MPI CUDA
$(OBJ_DIR)/mgcfd_mpi_kernels_cu.o:
	mkdir -p $(OBJ_DIR)
	nvcc $(NVCCFLAGS) $(MGCFD_INCS) $(OP2_INC) -I $(MPI_INSTALL_PATH)/include \
        -c -o $@ $(SRC_DIR)/../cuda/_kernels.cu
$(OBJ_DIR)/mgcfd_mpi_cuda_main.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CFLAGS) $^ $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
        -DCUDA_ON \
        -c -o $@ $(OP2_MAIN_SRC)
$(BIN_DIR)/mgcfd_mpi_cuda: $(OP2_MPI_CUDA_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CFLAGS) $(OPTIMISE) $(MGCFD_LIBS) $^ \
	    $(CUDA_LIB) -lcudart $(OP2_LIB) -lop2_mpi_cuda $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) \
        -o $@


clean:
	rm -f $(BIN_DIR)/* $(OBJ_DIR)/*
