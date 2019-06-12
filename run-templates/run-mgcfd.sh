set -e
set -u

# Compilation variables:
compiler=<COMPILER>
debug=<DEBUG>
papi=<PAPI>
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
validate_solution=<VALIDATE_SOLUTION>

partitioner=<PARTITIONER>
partitioner_method=<PARTITIONER_METHOD>

touch "${run_outdir}"/job-is-running.txt
if [ -f "${run_outdir}"/job-in-queue.txt ]; then
    rm "${run_outdir}"/job-in-queue.txt
fi

## Exit early if output csv files already exist.
if ls ${run_outdir}/*.csv 1> /dev/null 2>&1; then
  echo "Output CSV files already present, meaning this job has already run."
  rm "${run_outdir}"/job-is-running.txt
  exit 0
fi

# Compile:
if $mpi ; then
  if $cuda ; then
    bin_filename=mgcfd_mpi_cuda
  elif $openmp ; then
  	bin_filename=mgcfd_mpi_openmp
  else
  	bin_filename=mgcfd_mpi
  fi
else
  if $cuda ; then
  	bin_filename=mgcfd_cuda
  elif $openmp ; then
  	bin_filename=mgcfd_openmp
  elif $openmp4 ; then
  	bin_filename=mgcfd_openmp4
  elif $openacc ; then
  	bin_filename=mgcfd_openacc
  else
  	bin_filename=mgcfd_seq
  fi
fi
bin_filepath="${app_dirpath}/bin/${bin_filename}"

# if [ ! -f "$bin_filepath" ]; then
  ## Try compiling anyway, source files may have changed
  if [[ `hostname` == *"login"* ]] || [ `basename "$0"` = run.sh ]; then
    ## On login node, compile
    cd "${app_dirpath}"
    make_cmd="COMPILER=${compiler} "
    if [ "$cpp_wrapper" != "" ]; then
      make_cmd+="CPP_WRAPPER=$cpp_wrapper "
    fi
    if [ "$mpicpp_wrapper" != "" ]; then
      make_cmd+="MPICPP_WRAPPER=$mpicpp_wrapper "
    fi
    if $papi ; then
      make_cmd+="PAPI=1 "
    fi
    if $debug ; then
      make_cmd+="DEBUG=1 "
    fi
    make_cmd+="make -j4 $bin_filename"
    eval "$make_cmd"
    chmod a+x "$bin_filepath"
  elif [ ! -f "$bin_filepath" ]; then
    echo "ERROR: Cannot find binary: $bin_filepath"
    rm "${run_outdir}"/job-is-running.txt
    exit 1
  fi
# fi

if [[ `hostname` == *"login"* ]]; then
  ## Assume on a login node, do not execute the code.
  echo "Detected presence on login node, aborting before app execution."
  rm "${run_outdir}"/job-is-running.txt
  exit 0
fi

# Run:

cd "${data_dirpath}"
input_dat_filename="input-mgcfd.dat"
if [ ! -f "$input_dat_filename" ]; then
  input_dat_filename="input.dat"
  if [ ! -f "$input_dat_filename" ]; then
    echo "ERROR: Cannot find input .dat file"
    exit 1
  fi
fi

exec_command=""
if [ ! -z ${RUN_CMD+x} ]; then
  exec_command+="$RUN_CMD"
else
  if $mpi ; then
    exec_command+="mpirun -n $ntasks"
    if $debug ; then
      # exec_command+=" xterm -e gdb --args"
      exec_command+=" xterm -e gdb -ex run --args"
    fi
  else
    exec_command+=""
    if $debug ; then
      exec_command+=" gdb --args"
    fi
  fi
fi
exec_command+=" $bin_filepath OP_MAPS_BASE_INDEX=1 -i $input_dat_filename -o ${run_outdir}/ -g $mg_cycles -m $partitioner -r $partitioner_method"
if $papi ; then
  exec_command+=" -p ${run_outdir}/papi.conf"
fi
if $validate_solution ; then
  exec_command+=" -v"
fi
echo "$exec_command"
eval "$exec_command"
rm "${run_outdir}"/job-is-running.txt
