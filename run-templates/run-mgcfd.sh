set -e
set -u

# Compilation variables:
compiler=<COMPILER>
cpp_wrapper="<CPP_WRAPPER>"
mpicpp_wrapper="<MPICPP_WRAPPER>"
mpi=<MPI>
cuda=<CUDA>
openmp=<OPENMP>
openmp4=<OPENMP4>
openacc=<OPENACC>

# File/dir paths:
run_outdir="<RUN_OUTDIR>"
app_dirpath="<APP_DIRPATH>"
data_dirpath="<DATA_DIRPATH>"

# MG-CFD run variables:
nthreads=<NTHREADS>
ntasks=<NTASKS>
mg_cycles=<MG_CYCLES>
partitioner=<PARTITIONER>
validate_solution=<VALIDATE_SOLUTION>


## Exit early if output csv files already exist.
if ls ${run_outdir}/*.csv 1> /dev/null 2>&1; then
  echo "Output CSV files already present, no need to execute."
  exit 0
fi

if $mpi ; then
  if $cuda ; then
    bin_filename=mgcfd_mpi_cuda
  elif $openmp ; then
  	bin_filename=mgcfd_mpi_openmp
  else
  	bin_filename=mgcfd_mpi_genseq
  fi
else
  if $cuda ; then
  	bin_filename=mgcfd_cuda
  elif $openmp ; then
  	bin_filename=mgcfd_mpi_openmp
  elif $openmp4 ; then
  	bin_filename=mgcfd_mpi_openmp4
  elif $openacc ; then
  	bin_filename=mgcfd_mpi_openacc
  else
  	bin_filename=mgcfd_seq
  fi
fi
bin_filepath="${app_dirpath}/bin/${bin_filename}"



# Run:

# if [ ! -f "$bin_filepath" ]; then
  ## Try compiling anyway, source files may have changed
  # if [[ `hostname` == *"login"* ]]; then
    echo "Attempting to compile ..."
    cd "$app_dirpath"
    make_cmd="COMPILER=${compiler} "
    if [ "$cpp_wrapper" != "" ]; then
      make_cmd+="CPP_WRAPPER=$cpp_wrapper "
    fi
    if [ "$mpicpp_wrapper" != "" ]; then
      make_cmd+="MPICPP_WRAPPER=$mpicpp_wrapper "
    fi
    make_cmd+="make -j4 $bin_filename"
    echo "$make_cmd"
    eval "$make_cmd"
  # else
  #   echo "ERROR: Binary does not exist: $bin_filepath"
  #   exit 1
  # fi
# fi

if [[ `hostname` == *"login"* ]]; then
  ## Assume on a login node, do not execute the code.
  echo "Detected presence on login node, aborting before app execution."
  exit 0
fi

if [ ! -z ${RUN_CMD+x} ]; then
  exec_command="$RUN_CMD"
else
  if $mpi ; then
    exec_command="mpirun -n $ntasks"
  else
    exec_command=""
  fi
fi
exec_command+=" $bin_filepath OP_MAPS_BASE_INDEX=1 -i input.dat -p ${run_outdir}/papi.conf -o ${run_outdir}/ -g $mg_cycles -m $partitioner"
if $validate_solution ; then
  exec_command+=" -v"
fi
echo "EXECUTING $bin_filepath"
export OMP_NUM_THREADS=$nthreads
cd "${data_dirpath}"
echo "$exec_command"
eval "$exec_command"
