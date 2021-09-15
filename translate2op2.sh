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

	python3 "${OP2_INSTALL_PATH}"/../translator/c/python/op2.py $files
} 

_main

cd "$_dir"

