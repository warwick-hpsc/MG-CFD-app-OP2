include $(OP2_INSTALL_PATH)/../makefiles/common.mk

ifdef HDF5_INSTALL_PATH
  HDF5_SEQ_INSTALL_PATH ?= $(HDF5_INSTALL_PATH)
  HDF5_PAR_INSTALL_PATH ?= $(HDF5_INSTALL_PATH)
endif

APP_BIN_DIR := bin

APP_NAME := mgcfd
APP_INC := -Isrc/ -Isrc/Kernels -I$(HDF5_SEQ_INSTALL_PATH)/include
APP_SRC := src/euler3d_cpu_double.cpp


OP2_LIBS_WITH_HDF5 := true
VARIANT_FILTER_OUT := mpi_%

include $(OP2_INSTALL_PATH)/../makefiles/c_app.mk

APP_NAME := mgcfd_par
APP_INC := -Isrc/ -Isrc/Kernels -I$(HDF5_PAR_INSTALL_PATH)/include
APP_SRC := src/euler3d_cpu_double.cpp

VARIANT_FILTER := mpi_%

include $(OP2_INSTALL_PATH)/../makefiles/c_app.mk