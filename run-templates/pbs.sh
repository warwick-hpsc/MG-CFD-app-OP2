#!/bin/bash

#PBS -N MG-CFD.<RUN ID>
#PBS -o pbs.stdout
#PBS -e pbs.stderr

#PBS -l select=<NODES>
#PBS -l walltime=<HOURS>:<MINUTES>:00
#PBS -A <PROJECT CODE>
#PBS -q <PARTITION>

## Load compiler module(s):
compiler="<COMPILER>"
if [ "$compiler" = "intel" ]; then
  module swap PrgEnv-cray PrgEnv-intel
elif [ "$compiler" = "gnu" ]; then
  module swap PrgEnv-cray PrgEnv-gnu
# elif [ "$compiler" = "cray" ]; then
# elif [ "$compiler" = "clang" ]; then
fi

module load papi
module load cray-hdf5-parallel

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
cd $PBS_O_WORKDIR

# RUN_CMD="aprun -n <NTASKS> -N <TPN> -d <NTHREADS>"
RUN_CMD="aprun -n <NTASKS> -N <TPN>"

## Handle the threading aspect of aprun:
if [ "$compiler" = "intel" ]; then
	tpn=<TPN>
	nt=<NTHREADS>
	NUM_CORES=24
	NUM_CORES_PER_SOCKET=$((NUM_CORES / 2))

	CC_LIST=""
	if $openmp ; then
		if $mpi && [ <TPN> -gt 1 ] ; then
			ntasks_socket_b=$((tpn/2))
			ntasks_socket_a=$((tpn - ntasks_socket_b))
			CC_LIST=""
			core_id=0
			for (( p=0 ; p<${ntasks_socket_a} ; p+=1 )); do
				for (( t=0 ; t<${nt} ; t+=1 )); do
					if [ "$t" = "0" ]; then
						CC_LIST+="$core_id,$((core_id+NUM_CORES))"
					else
						CC_LIST+=",$core_id"
					fi
					core_id=$((core_id + 1))
				done
				CC_LIST+=":"
			done

			core_id=$NUM_CORES_PER_SOCKET
			for (( p=0 ; p<${ntasks_socket_b} ; p+=1 )); do
				for (( t=0 ; t<${nt} ; t+=1 )); do
					if [ "$t" = "0" ]; then
						CC_LIST+="$core_id,$((core_id+NUM_CORES))"
					else
						CC_LIST+=",$core_id"
					fi
					core_id=$((core_id + 1))
				done
				if [ "$p" != "$((ntasks_socket_b-1))" ]; then
					CC_LIST+=":"
				fi
			done
		else
			num_threads_socket_b=$((<NUM_THREADS>/2))
			num_threads_socket_a=$((<NUM_THREADS> - num_threads_socket_b))
			CC_LIST="0,${NUM_CORES}"
			for (( i=1 ; i<${num_threads_socket_a} ; i+=1 )); do
				CC_LIST+=",$i"
			done
			for (( i=0 ; i<${num_threads_socket_b} ; i+=1 )); do
				CC_LIST+=",$((i+${NUM_CORES_PER_SOCKET}))"
			done
		fi
	fi

	if [ "$CC_LIST" != "" ]; then
		# RUN_CMD="aprun -n <NTASKS> -j 2 -d $((2*<NTHREADS>))"
		# RUN_CMD="aprun -n <NTASKS> -N <TPN> -j 2 -d $((2*<NTHREADS>))"
		# RUN_CMD+=" -cc ${CC_LIST}"
		RUN_CMD+=" -j 2 -d $((2*<NTHREADS>)) -cc ${CC_LIST}"
		export KMP_AFFINITY=disabled
	fi
else
	RUN_CMD+=" -d <NTHREADS>"
fi