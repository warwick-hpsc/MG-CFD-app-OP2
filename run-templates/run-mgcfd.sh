# set -e
# set -u

# Compilation variables:
compiler="<COMPILER>"
debug=<DEBUG>
papi=<PAPI>
cpp_wrapper="<CPP_WRAPPER>"
mpicpp_wrapper="<MPICPP_WRAPPER>"
mpi=<MPI>

# File/dir paths:
run_outdir="<RUN_OUTDIR>"
app_dirpath="<APP_DIRPATH>"
data_dirpath="<DATA_DIRPATH>"

# MG-CFD run variables:
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

################################
# Compile:
################################

bin_filename="<BIN_FILENAME>"
bin_filepath="<BIN_FILEPATH>"
make_target="<MAKE_TARGET>"

recompile_required=false
if [ -f "$bin_filepath" ]; then
  if grep -q "PAPI_NULL" "$bin_filepath" && ! $papi ; then
    echo "Binary compiled with PAPI, recompile required to remove it"
    recompile_required=true
  fi
  if $papi && ! grep -q "PAPI_NULL" "$bin_filepath" ; then
    echo "Binary not compiled with PAPI, recompile required for papi"
    recompile_required=true
  fi
fi

if [[ `hostname` == *"login"* ]] || [ "`head -n1 "$0"`" = "#!/bin/bash" ]; then
  ## On login node, or executing local run script, so compile
  cd "${app_dirpath}"
  if $recompile_required ; then
    make clean_"${make_target}"
  fi
  make_cmd="COMPILER=${compiler} "
  if [ ! -z ${cpp_wrapper+x} ] && [ "$cpp_wrapper" != "" ]; then
    make_cmd+="CPP_WRAPPER=$cpp_wrapper "
  fi
  if [ ! -z ${mpicpp_wrapper+x} ] && [ "$mpicpp_wrapper" != "" ]; then
    make_cmd+="MPICPP_WRAPPER=$mpicpp_wrapper "
  fi
  if $papi ; then
    make_cmd+="PAPI=1 "
  fi
  if $debug ; then
    make_cmd+="DEBUG=1 "
  fi
  make_cmd+="make -j4 $make_target"
  eval "$make_cmd"
  chmod a+x "$bin_filepath"
elif [ ! -f "$bin_filepath" ]; then
  echo "ERROR: Cannot find binary: $bin_filepath"
  rm "${run_outdir}"/job-is-running.txt
  exit 1
fi

if [[ `hostname` == *"login"* ]]; then
  ## Assume on a login node, do not execute the code.
  echo "Detected presence on login node, aborting before app execution."
  rm "${run_outdir}"/job-is-running.txt
  exit 0
fi

################################
# Run:
################################

cd "${data_dirpath}"
input_dat_filename="input-mgcfd.dat"
if [ ! -f "$input_dat_filename" ]; then
  input_dat_filename="mgcfd-input.dat"
  if [ ! -f "$input_dat_filename" ]; then
    input_dat_filename="input.dat"
    if [ ! -f "$input_dat_filename" ]; then
      echo "ERROR: Cannot find input .dat file"
      exit 1
    fi
  fi
fi

exec_command=""
if [ ! -z ${RUN_CMD+x} ]; then
  exec_command+="$RUN_CMD"
else
  if $mpi ; then
    exec_command+="mpirun -n <NTASKS>"
    if $debug ; then
      exec_command+=" xterm -e gdb -ex run --args"
    fi
  else
    exec_command+=""
    if $debug ; then
      exec_command+=" gdb --args"
    fi
  fi
fi
exec_command+=" $bin_filepath -i $input_dat_filename -o ${run_outdir}/ -g $mg_cycles -m $partitioner -r $partitioner_method"
if $papi ; then
  exec_command+=" -p ${run_outdir}/papi.conf"
fi
if $validate_solution ; then
  exec_command+=" -v"
fi
echo "$exec_command"
eval "$exec_command"
rm "${run_outdir}"/job-is-running.txt
