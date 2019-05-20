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

	python "${OP2_INSTALL_PATH}"/../translator/c/python/op2.py $files

	## Now apply hacks to get around OP2 compilations issues:
	# 1. Disable op_reductions in cuda kernels:
	# echo "WARNING: Disabling op_reduction() calls in CUDA code."
	# ls "${_dir}"/cuda | while read F ; do
	# 	sed -i "s|op_reduction|// op_reduction|g" cuda/"$F"
	# done
} 

_main

cd "$_dir"

