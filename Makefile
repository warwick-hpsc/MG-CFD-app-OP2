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

ifdef PETSC_INSTALL_PATH
    PETSC_INC = -I$(PETSC_INSTALL_PATH)/include
    PETSC_LIB = -L$(PETSC_INSTALL_PATH)/lib
endif

ifdef DOLFINX_INSTALL_PATH
    DOLFINX_INC = -I$(DOLFINX_INSTALL_PATH)/include
    DOLFINX_LIB = -L$(DOLFINX_INSTALL_PATH)/lib64
endif

ifdef BOOST_INSTALL_PATH
    BOOST_LIB = -L$(BOOST_INSTALL_PATH)/lib
endif

ifdef SQLITE_INSTALL_PATH
    SQLITE_INC = -I$(SQLITE_INSTALL_PATH)/include
    SQLITE_LIB = -L$(SQLITE_INSTALL_PATH)/lib
endif

ifdef TREETIMER_INSTALL_PATH
    TREETIMER_INC = -I$(TREETIMER_INSTALL_PATH)/include/timing_library/interface
    TREETIMER_LIB = -L$(TREETIMER_INSTALL_PATH) -ltt
endif

ifdef SIMPIC_INSTALL_PATH
    SIMPIC_INC = -I$(SIMPIC_INSTALL_PATH)
    SIMPIC_LIB = -L$(SIMPIC_INSTALL_PATH)/libsimpic.a
endif

ifdef FENICS
    FENICS_DEF = -Ddeffenics
    DOLFINX_LIB += -ldolfinx
    FENICS_LIB = $(BIN_DIR)/libdolfinx_cpx.a
    PETSC_INC += -DHAVE_PETSC
    PETSC_LIB += -lpetsc
    BOOST_LIB += -lboost_filesystem -lboost_program_options -lboost_timer
endif

ifdef SIMPIC
    SIMPIC_DEF = -Ddefsimpic
	SIMPIC_LIB = $(SIMPIC_INSTALL_PATH)/libsimpic.a
endif

ifdef DEBUG
  OPTIMISE := -pg -g -O0
else
  OPTIMISE := -O3
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

ifdef OP2_COMPILER
  ifeq ($(COMPILER),)
    COMPILER=$(OP2_COMPILER)
  endif
endif
ifeq ($(COMPILER),gnu)
  CPP := CC
  CFLAGS	= -fPIC -DUNIX -Wall -Wextra
  CPPFLAGS 	= $(CFLAGS)
  OMPFLAGS 	= -fopenmp
  MPIFLAGS 	= $(CPPFLAGS)
  MPICPP = CC
else
ifeq ($(COMPILER),intel)
  CPP = icpc
  CFLAGS = -DMPICH_IGNORE_CXX_SEEK -inline-forceinline -DVECTORIZE -qopt-report=5
  CFLAGS += -restrict
  # CFLAGS += -parallel ## This flag intoduces a significant slowdown into 'vec' app
  # CFLAGS += -fno-alias ## This flag causes 'vec' app to fail validation, do not enable
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
  CPP		 = xlc++
  CFLAGS	 = -qarch=pwr8 -qtune=pwr8 -qhot
  CPPFLAGS 	 = $(CFLAGS)
  OMPFLAGS	 = -qsmp=omp -qthreaded
  OMPOFFLOAD = -qsmp=omp -qoffload -Xptxas -v -g1
  MPIFLAGS	 = $(CPPFLAGS)
else
ifeq ($(COMPILER),pgi)
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
ifeq ($(COMPILER),cray)
  CPP           = CC
  CFLAGS        = -h fp3 -h ipa5
  CPPFLAGS      = $(CFLAGS)
  OMPFLAGS      = -h omp
  MPICPP        = CC
  MPIFLAGS      = $(CPPFLAGS)
else
  $(error unrecognised value for COMPILER: $(COMPILER))
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
  MESSAGE=select an NVIDA device to compile in CUDA, e.g. make NV_ARCH=Pascal
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
else
ifeq ($(NV_ARCH),Volta)
  CODE_GEN_CUDA=-gencode arch=compute_70,code=sm_70
endif
endif
endif
endif
endif

NVCCFLAGS = 
ifdef NVCC_BIN
	NVCCFLAGS = -ccbin $(NVCC_BIN)
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

## Enable DUMP_EXT_PERF_DATA to write out externally-collected 
## performance data. Included number of loop iterations counts of 
## each kernel, and if VERIFY_OP2_TIMING is enabled then also 
## its compute and sync times.
# MGCFD_INCS += -DDUMP_EXT_PERF_DATA

all: seq openmp mpi mpi_vec mpi_openmp
# all: seq openmp mpi mpi_vec mpi_openmp cuda mpi_cuda
# all: seq openmp mpi mpi_vec mpi_openmp cuda mpi_cuda openacc openmp4

parallel: N = $(shell nproc)
parallel:; @$(MAKE) -j$(N) -l$(N) all

## User-friendly wrappers around actual targets:
seq: $(BIN_DIR)/mgcfd_seq
openmp: $(BIN_DIR)/mgcfd_openmp
mpi: $(BIN_DIR)/mgcfd_mpi
vec: mpi_vec
mpi_vec: $(BIN_DIR)/mgcfd_mpi_vec
mpi_openmp: $(BIN_DIR)/mgcfd_mpi_openmp
cuda: $(BIN_DIR)/mgcfd_cuda
mpi_cuda: $(BIN_DIR)/mgcfd_mpi_cuda
openmp4: $(BIN_DIR)/mgcfd_openmp4
openacc: $(BIN_DIR)/mgcfd_openacc
mpi_cpx: $(BIN_DIR)/mgcfd_cpx.a $(BIN_DIR)/mgcfd_cpx_runtime
mpi_cuda_cpx: $(BIN_DIR)/mgcfd_cpx_cuda.a $(BIN_DIR)/mgcfd_cpx_cuda_runtime


OP2_MAIN_SRC = $(SRC_DIR)_op/euler3d_cpu_double_op.cpp

OP2_SEQ_OBJECTS := $(OBJ_DIR)/mgcfd_seq_main.o \
                   $(OBJ_DIR)/mgcfd_seq_kernels.o

OP2_MPI_OBJECTS := $(OBJ_DIR)/mgcfd_mpi_main.o \
                   $(OBJ_DIR)/mgcfd_mpi_kernels.o

OP2_MPI_CPX_OBJECTS := $(OBJ_DIR)/mgcfd_mpi_cpx_main.o \
                   $(OBJ_DIR)/mgcfd_mpi_cpx_kernels.o

OP2_MPI_CPX_CUDA_OBJECTS := $(OBJ_DIR)/mgcfd_mpi_cuda_cpx_main.o \
                        $(OBJ_DIR)/mgcfd_mpi_cpx_kernels_cu.o

OP2_MPI_CPX_MAIN := $(SRC_DIR)_op/coupler.cpp

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


KERNELS := calc_rms_kernel \
	calculate_cell_volumes \
	calculate_dt_kernel \
	compute_bnd_node_flux_kernel \
	compute_flux_edge_kernel \
	compute_step_factor_kernel \
	copy_double_kernel \
	count_bad_vals \
	count_non_zeros \
	dampen_ewt \
	down_kernel \
	down_v2_kernel_post \
	down_v2_kernel_pre \
	down_v2_kernel \
	get_min_dt_kernel \
	identify_differences \
	initialize_variables_kernel \
	residual_kernel \
	time_step_kernel \
	up_kernel \
	up_post_kernel \
	up_pre_kernel \
	zero_1d_array_kernel \
	zero_5d_array_kernel
SEQ_KERNELS := $(patsubst %, $(SRC_DIR)/../seq/%_seqkernel.cpp, $(KERNELS))
OMP_KERNELS := $(patsubst %, $(SRC_DIR)/../openmp/%_kernel.cpp, $(KERNELS))
CUDA_KERNELS := $(patsubst %, $(SRC_DIR)/../cuda/%_kernel.cu, $(KERNELS))
VEC_KERNELS := $(patsubst %, $(SRC_DIR)/../vec/%_veckernel.cpp, $(KERNELS))
ACC_KERNELS := $(patsubst %, $(SRC_DIR)/../openacc/%_acckernel.c, $(KERNELS))
OMP4_KERNELS := $(patsubst %, $(SRC_DIR)/../openmp4/%_omp4kernel.cpp, $(KERNELS))
OMP4_KERNEL_FUNCS := $(patsubst %, $(SRC_DIR)/../openmp4/%_omp4kernel_func.cpp, $(KERNELS))



## SEQUENTIAL
$(OBJ_DIR)/mgcfd_seq_main.o: $(OP2_MAIN_SRC)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OPTIMISE) $(MGCFD_INCS) \
	    $(OP2_INC) $(HDF5_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) \
		-c -o $@ $^
$(OBJ_DIR)/mgcfd_seq_kernels.o: $(SRC_DIR)/../seq/_seqkernels.cpp $(SEQ_KERNELS)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OPTIMISE) $(MGCFD_INCS) \
	    $(OP2_INC) $(HDF5_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) \
		-c -o $@ $(SRC_DIR)/../seq/_seqkernels.cpp
$(BIN_DIR)/mgcfd_seq: $(OP2_SEQ_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OPTIMISE) $^ $(MGCFD_LIBS) \
		-lm $(OP2_LIB) -lop2_seq -lop2_hdf5 $(HDF5_LIB) $(PARMETIS_LIB) $(PTSCOTCH_LIB) \
		-o $@


## OPENMP
$(OBJ_DIR)/mgcfd_openmp_main.o: $(OP2_MAIN_SRC)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_INCS) \
		$(OP2_INC) $(HDF5_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) \
		-c -o $@ $^
$(OBJ_DIR)/mgcfd_openmp_kernels.o: $(SRC_DIR)/../openmp/_kernels.cpp $(OMP_KERNELS)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_INCS) \
		$(OP2_INC) $(HDF5_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) \
		-c -o $@ $(SRC_DIR)/../openmp/_kernels.cpp
$(BIN_DIR)/mgcfd_openmp: $(OP2_OMP_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $^ $(MGCFD_LIBS) \
		-lm $(OP2_LIB) -lop2_openmp -lop2_hdf5 $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) \
		-o $@


## MPI
$(OBJ_DIR)/mgcfd_mpi_main.o: $(OP2_MAIN_SRC)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	     -DMPI_ON -c -o $@ $^
$(OBJ_DIR)/mgcfd_mpi_kernels.o: $(SRC_DIR)/../seq/_seqkernels.cpp $(SEQ_KERNELS)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	     -DMPI_ON -c -o $@ $(SRC_DIR)/../seq/_seqkernels.cpp
$(BIN_DIR)/mgcfd_mpi: $(OP2_MPI_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OPTIMISE) $^ $(MGCFD_LIBS) \
		-lm $(OP2_LIB) -lop2_mpi $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) \
		-o $@

## MPI_CPX
$(OBJ_DIR)/mgcfd_mpi_cpx_main.o: $(OP2_MAIN_SRC)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	     -DMPI_ON -c -o $@ $^ 

$(OBJ_DIR)/mgcfd_mpi_cpx_kernels.o: $(SRC_DIR)/../seq/_seqkernels.cpp $(SEQ_KERNELS)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	     -DMPI_ON -c -o $@ $(SRC_DIR)/../seq/_seqkernels.cpp
		 
$(BIN_DIR)/mgcfd_cpx.a: $(OP2_MPI_CPX_OBJECTS)
	mkdir -p $(BIN_DIR)
	ar rcs $(BIN_DIR)/mgcfd_cpx.a $(OP2_MPI_CPX_OBJECTS)

$(BIN_DIR)/mgcfd_cpx_runtime: $(OP2_MPI_CPX_MAIN) $(BIN_DIR)/mgcfd_cpx.a $(FENICS_LIB) $(SIMPIC_LIB)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OPTIMISE) $^ $(BIN_DIR)/mgcfd_cpx.a $(SIMPIC_LIB) $(MGCFD_LIBS) \
        -lm $(OP2_LIB) -lop2_mpi $(PARMETIS_LIB) $(DOLFINX_LIB) $(PETSC_LIB) $(BOOST_LIB)\
		$(SQLITE_LIB) $(TREETIMER_INC) $(TREETIMER_LIB) \
        $(PTSCOTCH_LIB) $(HDF5_LIB) $(FENICS_DEF) $(SIMPIC_DEF) -o $@ 


## MPI_CPX CUDA

$(OBJ_DIR)/mgcfd_mpi_cpx_kernels_cu.o: $(SRC_DIR)/../cuda/_kernels.cu $(CUDA_KERNELS)
	mkdir -p $(OBJ_DIR)
	nvcc $(NVCCFLAGS) $(MGCFD_INCS) $(OP2_INC) -I $(MPI_INSTALL_PATH)/include \
        -c -o $@ $(SRC_DIR)/../cuda/_kernels.cu

$(OBJ_DIR)/mgcfd_mpi_cuda_cpx_main.o: $(OP2_MAIN_SRC)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
        -DCUDA_ON -c -o $@ $^

$(BIN_DIR)/mgcfd_cpx_cuda.a: $(OP2_MPI_CPX_CUDA_OBJECTS)
	mkdir -p $(BIN_DIR)
	ar rcs $(BIN_DIR)/mgcfd_cpx_cuda.a $(OP2_MPI_CPX_CUDA_OBJECTS)

$(BIN_DIR)/mgcfd_cpx_cuda_runtime: $(OP2_MPI_CPX_MAIN) $(BIN_DIR)/mgcfd_cpx_cuda.a $(FENICS_LIB) $(SIMPIC_LIB)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OPTIMISE) $^ $(BIN_DIR)/mgcfd_cpx_cuda.a $(SIMPIC_LIB) $(MGCFD_LIBS) \
        $(CUDA_LIB) -lcudart $(OP2_LIB) -lop2_mpi_cuda $(PARMETIS_LIB) $(DOLFINX_LIB) $(PETSC_LIB) \
		$(SQLITE_LIB) $(TREETIMER_INC) $(TREETIMER_LIB) \
        $(PTSCOTCH_LIB) $(HDF5_LIB) $(FENICS_DEF) $(SIMPIC_DEF) -o $@ 

## MPI_VEC
$(OBJ_DIR)/mgcfd_mpi_vec_main.o: $(OP2_MAIN_SRC)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
        -DMPI_ON -c -o $@ $^
$(OBJ_DIR)/mgcfd_mpi_vec_kernels.o: $(SRC_DIR)/../vec/_veckernels.cpp $(VEC_KERNELS)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
        -DMPI_ON -c -o $@ $(SRC_DIR)/../vec/_veckernels.cpp
$(BIN_DIR)/mgcfd_mpi_vec: $(OP2_MPI_VEC_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $^ $(MGCFD_LIBS) \
        -lm $(OP2_LIB) -lop2_mpi $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) \
        -o $@


## MPI + OPENMP
$(OBJ_DIR)/mgcfd_mpi_openmp_main.o: $(OP2_MAIN_SRC)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -DMPI_ON -c -o $@ $^
$(OBJ_DIR)/mgcfd_mpi_openmp_kernels.o: $(SRC_DIR)/../openmp/_kernels.cpp $(OMP_KERNELS)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -DMPI_ON -c -o $@ $(SRC_DIR)/../openmp/_kernels.cpp
$(BIN_DIR)/mgcfd_mpi_openmp: $(OP2_MPI_OMP_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $^ $(MGCFD_LIBS) \
		-lm $(OP2_LIB) -lop2_mpi $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) \
		-o $@


## CUDA
$(OBJ_DIR)/mgcfd_cuda_main.o: $(OP2_MAIN_SRC)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -DCUDA_ON -c -o $@ $^
$(OBJ_DIR)/mgcfd_kernels_cu.o: $(SRC_DIR)/../cuda/_kernels.cu $(CUDA_KERNELS)
	mkdir -p $(OBJ_DIR)
	nvcc $(NVCCFLAGS) $(MGCFD_INCS) $(OP2_INC) \
		-c -o $@ $(SRC_DIR)/../cuda/_kernels.cu
$(BIN_DIR)/mgcfd_cuda: $(OP2_CUDA_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CFLAGS) $(OPTIMISE) $^ $(MGCFD_LIBS) \
	    $(CUDA_LIB) -lcudart $(OP2_LIB) -lop2_cuda $(HDF5_LIB) -lop2_hdf5 \
	    -o $@


## MPI CUDA
$(OBJ_DIR)/mgcfd_mpi_kernels_cu.o: $(SRC_DIR)/../cuda/_kernels.cu $(CUDA_KERNELS)
	mkdir -p $(OBJ_DIR)
	nvcc $(NVCCFLAGS) $(MGCFD_INCS) $(OP2_INC) -I $(MPI_INSTALL_PATH)/include \
        -c -o $@ $(SRC_DIR)/../cuda/_kernels.cu
$(OBJ_DIR)/mgcfd_mpi_cuda_main.o: $(OP2_MAIN_SRC)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
        -DCUDA_ON -c -o $@ $^
$(BIN_DIR)/mgcfd_mpi_cuda: $(OP2_MPI_CUDA_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CFLAGS) $(OPTIMISE) $^ $(MGCFD_LIBS) \
	    $(CUDA_LIB) -lcudart $(OP2_LIB) -lop2_mpi_cuda $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) \
        -o $@


## OPENMP4
$(OBJ_DIR)/mgcfd_omp4_main.o: $(OP2_MAIN_SRC)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -Iopenmp4/ -c -o $@ $^
$(OBJ_DIR)/mgcfd_omp4_kernel_funcs.o: $(SRC_DIR)/../openmp4/_omp4kernel_funcs.cpp $(OMP4_KERNEL_FUNCS)
	mkdir -p $(OBJ_DIR)
	$(CPP) $(CPPFLAGS) $(OMPOFFLOAD) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -Iopenmp4/ -c -o $@ $(SRC_DIR)/../openmp4/_omp4kernel_funcs.cpp
$(OBJ_DIR)/mgcfd_omp4_kernels.o: $(SRC_DIR)/../openmp4/_omp4kernels.cpp $(OMP4_KERNELS)
	mkdir -p $(OBJ_DIR)
	$(CPP) $(CPPFLAGS) $(OMPFLAGS) $(OPTIMISE) $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) \
	    -Iopenmp4/ -c -o $@ $(SRC_DIR)/../openmp4/_omp4kernels.cpp
$(BIN_DIR)/mgcfd_openmp4: $(OP2_OMP4_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(OMPOFFLOAD) $(OPTIMISE) $^ $(MGCFD_LIBS) \
	    $(OP2_LIB) -lop2_openmp4 $(CUDA_LIB) -lcudart \
		$(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) -lop2_hdf5 \
	    -o $@


## OPENACC
$(OBJ_DIR)/mgcfd_openacc_main.o: $(OP2_MAIN_SRC)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(ACCFLAGS) $(OMPFLAGS) $(OPTIMISE) \
	    $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) -Iopenacc \
	    -c -o $@ $^
# Compile failing: nvlink unable to link _acckernels.o to extern variable definitions
$(OBJ_DIR)/mgcfd_openacc_kernels.o: $(SRC_DIR)/../openacc/_acckernels.c $(ACC_KERNELS)
	mkdir -p $(OBJ_DIR)
	$(MPICPP) $(CPPFLAGS) $(ACCFLAGS) $(OMPFLAGS) $(OPTIMISE) \
	    $(MGCFD_INCS) $(OP2_INC) $(HDF5_INC) -Iopenacc \
	    -c -o $@ $(SRC_DIR)/../openacc/_acckernels.c
$(BIN_DIR)/mgcfd_openacc: $(OP2_OPENACC_OBJECTS)
	mkdir -p $(BIN_DIR)
	$(MPICPP) $(CPPFLAGS) $(ACCFLAGS) $(OMPFLAGS) $(OPTIMISE) $^ $(MGCFD_LIBS) \
	    $(CUDA_LIB) -lcudart $(OP2_LIB) -lop2_cuda $(HDF5_LIB) -lop2_hdf5 \
	    -o $@


clean:
	rm -f $(BIN_DIR)/* $(OBJ_DIR)/*
clean_seq:
	rm -f $(BIN_DIR)/mgcfd_seq $(OP2_SEQ_OBJECTS)
clean_mpi:
	rm -f $(BIN_DIR)/mgcfd_mpi $(OP2_MPI_OBJECTS)
clean_mpi_cpx:
	rm -f $(BIN_DIR)/mgcfd_cpx.a $(OP2_MPI_CPX_OBJECTS)
clean_mpi_vec:
	rm -f $(BIN_DIR)/mgcfd_mpi_vec $(OP2_MPI_VEC_OBJECTS)
clean_openmp:
	rm -f $(BIN_DIR)/mgcfd_openmp $(OP2_OMP_OBJECTS)
clean_mpi_openmp:
	rm -f $(BIN_DIR)/mgcfd_mpi_openmp $(OP2_MPI_OMP_OBJECTS)
clean_cuda:
	rm -f $(BIN_DIR)/mgcfd_cuda $(OP2_CUDA_OBJECTS)
clean_mpi_cuda:
	rm -f $(BIN_DIR)/mgcfd_mpi_cuda $(OP2_MPI_CUDA_OBJECTS)
clean_openacc:
	rm -f $(BIN_DIR)/mgcfd_openacc $(OP2_OPENACC_OBJECTS)
clean_openmp4:
	rm -f $(BIN_DIR)/mgcfd_openmp4 $(OP2_OMP4_OBJECTS)


