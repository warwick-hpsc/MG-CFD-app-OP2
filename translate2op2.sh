#!/bin/bash

_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

function _main {
	set -e

	if [ "$OP2_INSTALL_PATH" = "" ]; then
		echo "ERROR: OP2_INSTALL_PATH not set"
		return
	fi

	cd "${_dir}"
	src_dir="./src"

	files=""
	files+="${src_dir}/euler3d_cpu_double.cpp"
	files+=" ${src_dir}/const.h"
	files+=" ${src_dir}/global.h"
	files+=" ${src_dir}"
	files+=" ${src_dir}/Kernels"

	## Default OP2 translator only inserts timers for individual ranks, ignores threads. 
	## Thus cannot detect load imbalance within a team of OpenMP threads. 
	## To add thread timers, set OP_TIME_THREADS environment variable:
	# OP_TIME_THREADS=1 python "${OP2_INSTALL_PATH}"/../translator/c/python/op2.py $files
	if [[ -v OP2_OLD ]]; then
		python3 "${OP2_INSTALL_PATH}"/../translator/c/python/op2.py $files
	else
		python3 "${OP2_INSTALL_PATH}"/../translator/c/op2.py $files
	fi
	sed -i 's/SIMD_VEC 4/SIMD_VEC 8/g' vec/_veckernels.cpp
} 

_main

cd "$_dir"

