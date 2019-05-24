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



OPTIMISE := -O2
# OPTIMISE := -pg -g -O0

BIN_DIR = bin
OBJ_DIR = obj
SRC_DIR = src


#
# Locate MPI compilers:
#
ifdef MPI_INSTALL_PATH
  ifneq ("","$(wildcard $(MPI_INSTALL_PATH)/bin/mpic++)")
    MPICPP := $(MPI_INSTALL_PATH)/bin/mpic++
  else
  ifneq ("","$(wildcard $(MPI_INSTALL_PATH)/intel64/bin/mpic++)")
    MPICPP := $(MPI_INSTALL_PATH)/intel64/bin/mpic++
  else
    MPICPP := mpic++
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
  MPICPP := mpic++
  MPICC  := mpicc
endif

ifeq ($(COMPILER),gnu)
  CPP := g++
  CFLAGS	= -fPIC -DUNIX -Wall -Wextra
  CPPFLAGS 	= $(CFLAGS)
  OMPFLAGS 	= -fopenmp
  MPIFLAGS 	= $(CPPFLAGS)
else
ifeq ($(COMPILER),intel)
  CPP = icpc
  CFLAGS = -DMPICH_IGNORE_CXX_SEEK -restrict -fno-alias -inline-forceinline -parallel -DVECTORIZE #-parallel #-DCOMM_PERF #-DDEBUG #-qopt-report=5
  CFLAGS += -fmax-errors=1
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
ifeq ($(COMPILER),xl)
  CPP		= xlc++
  CFLAGS	= -O5 -qarch=pwr8 -qtune=pwr8 -qhot
  CPPFLAGS 	= $(CFLAGS)
  OMPFLAGS	= -qsmp=omp -qthreaded
  OMPOFFLOAD    = -qsmp=omp -qoffload -Xptxas -v -g1
  MPIFLAGS	= $(CPPFLAGS)
else
ifeq ($(COMPILER),pgi)
  CPP       	= pgc++
  CFLAGS  	= -O3
  CPPFLAGS 	= $(CFLAGS)
  OMPFLAGS 	= -mp
  MPIFLAGS 	= $(CPPFLAGS)
  # NVCCFLAGS	+= -ccbin=$(MPICPP)
  # ACCFLAGS      = -acc -Minfo=acc -ta=tesla:cc35 -DOPENACC
  # ACCFLAGS      = -acc -DOPENACC -Minfo=acc
  ACCFLAGS      = -v -acc -DOPENACC -Minfo=acc
else
ifeq ($(COMPILER),cray)
  CPP           = CC
  CFLAGS       = -O3 -h fp3 -h ipa5
  CPPFLAGS      = $(CFLAGS)
  OMPFLAGS      = -h omp
  MPICPP        = CC
  MPIFLAGS      = $(CPPFLAGS)
else
# print:
# 	@echo "unrecognised value for COMPILER"
$(error unrecognised value for COMPILER)
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


OP2_MAIN_SRC = $(SRC_DIR)_op/euler3d_cpu_double_op.cpp

# all: mgcfd_seq mgcfd_openmp mgcfd_openacc mgcfd_openmp4 mgcfd_cuda \
# 	 mgcfd_mpi_genseq mgcfd_mpi_openmp mgcfd_mpi_cuda

all: mgcfd_seq mgcfd_openmp mgcfd_mpi_genseq mgcfd_mpi_openmp mgcfd_cuda

OP2_SEQ_SOURCES := $(OP2_MAIN_SRC) \
                   $(SRC_DIR)/../seq/_seqkernels.cpp

OP2_OMP_SOURCES := $(OP2_MAIN_SRC) \
                   $(SRC_DIR)/../openmp/_kernels.cpp

OP2_CUDA_SOURCES := $(OP2_MAIN_SRC) \
                    $(OBJ_DIR)/mgcfd_kernels_cu.o

OP2_OPENACC_SOURCES := $(OP2_MAIN_SRC) \
                       $(SRC_DIR)/../openacc/_acckernels.c

# OP2_OMP4_SOURCES := $(OP2_MAIN_SRC) \
#                     $(OBJ_DIR)/mgcfd_omp4_kernel_funcs.o \
#                     $(OBJ_DIR)/mgcfd_omp4_kernels.o
OP2_OMP4_SOURCES := $(OBJ_DIR)/mgcfd_omp4_main.o \
                    $(OBJ_DIR)/mgcfd_omp4_kernel_funcs.o \
                    $(OBJ_DIR)/mgcfd_omp4_kernels.o

## SEQUENTIAL
mgcfd_seq: $(OP2_SEQ_SOURCES)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $^ $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
		$(PARMETIS_INC) $(PTSCOTCH_INC) \
		-lm $(OP2_LIB) -lop2_seq -lop2_hdf5 \
		$(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB)  -o $(BIN_DIR)/mgcfd_seq

## OPENMP
mgcfd_openmp: $(OP2_OMP_SOURCES)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $^ $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
		$(PARMETIS_INC) $(PTSCOTCH_INC) \
		-lm $(OP2_LIB) -lop2_openmp -lop2_hdf5 \
		$(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB)  -o $(BIN_DIR)/mgcfd_openmp

## MPI
mgcfd_mpi_genseq: $(OP2_SEQ_SOURCES)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $^ $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
		$(PARMETIS_INC) $(PTSCOTCH_INC) \
		-lm $(OP2_LIB) -lop2_mpi -DMPI_ON \
		$(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB)  -o $(BIN_DIR)/mgcfd_mpi_genseq

## MPI + OPENMP
mgcfd_mpi_openmp: $(OP2_OMP_SOURCES)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $^ $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
		$(PARMETIS_INC) $(PTSCOTCH_INC) \
		-lm $(OP2_LIB) -lop2_mpi -DMPI_ON \
		$(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB)  -o $(BIN_DIR)/mgcfd_mpi_openmp

## CUDA
$(OBJ_DIR)/mgcfd_kernels_cu.o:
	mkdir -p $(OBJ_DIR)
	nvcc $(NVCCFLAGS) $(MGCFD_INCS) $(OP2_INC) \
		 -c -o $(OBJ_DIR)/mgcfd_kernels_cu.o $(SRC_DIR)/../cuda/_kernels.cu
mgcfd_cuda: $(OP2_CUDA_SOURCES)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CFLAGS) $^ $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    $(CUDA_LIB) -lcudart \
	    $(OP2_LIB) -lop2_cuda -DCUDA_ON \
	    $(HDF5_LIB) -lop2_hdf5 \
	    -o $(BIN_DIR)/mgcfd_cuda

## OPENMP4
$(OBJ_DIR)/mgcfd_omp4_kernel_funcs.o:
	mkdir -p $(OBJ_DIR)
	$(CPP) $(CPPFLAGS) $(OMPOFFLOAD) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -Iopenmp4/ -c $(SRC_DIR)/../openmp4/_omp4kernel_funcs.cpp -o $@
$(OBJ_DIR)/mgcfd_omp4_kernels.o:
	mkdir -p $(OBJ_DIR)
	$(CPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -Iopenmp4/ -c $(SRC_DIR)/../openmp4/_omp4kernels.cpp -o $@
$(OBJ_DIR)/mgcfd_omp4_main.o:
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $^ $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -Iopenmp4/ -c $(OP2_MAIN_SRC) -o $@
mgcfd_openmp4: $(OP2_OMP4_SOURCES)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPOFFLOAD) $^ $(OPTIMISE) $(MGCFD_INCS) \
	    $(OP2_INC) $(OP2_LIB) -lop2_openmp4 \
		$(PARMETIS_INC) $(PTSCOTCH_INC) $(PARMETIS_LIB) $(PTSCOTCH_LIB) \
		$(HDF5_INC) $(HDF5_LIB) -lop2_hdf5 \
	    $(CUDA_LIB) -lcudart \
	    -o $(BIN_DIR)/mgcfd_openmp4


## OPENACC
# Compile failing: nvlink unable to link _acckernels.o to extern variable definitions
mgcfd_openacc: $(OP2_OPENACC_SOURCES)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(ACCFLAGS) $(OMPFLAGS) $^ $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -Iopenacc $(CUDA_LIB) -lcudart \
	    $(OP2_LIB) -lop2_cuda \
	    $(HDF5_LIB) -lop2_hdf5 \
	    -o $(BIN_DIR)/mgcfd_openacc


## MPI CUDA
$(OBJ_DIR)/mgcfd_mpi_kernels_cu.o:
	mkdir -p $(OBJ_DIR)
	nvcc $(NVCCFLAGS) $(MGCFD_INCS) $(OP2_INC) -I $(MPI_INSTALL_PATH)/include \
                 -c -o $(OBJ_DIR)/mgcfd_mpi_kernels_cu.o $(SRC_DIR)/../cuda/_kernels.cu
mgcfd_mpi_cuda: $(OP2_CUDA_SOURCES)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CFLAGS) $^ $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    $(CUDA_LIB) -lcudart \
            $(OP2_LIB) -lop2_mpi_cuda -DCUDA_ON \
	    $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) \
            -o $(BIN_DIR)/mgcfd_mpi_cuda


clean:
	rm -f $(BIN_DIR)/* $(OBJ_DIR)/*
