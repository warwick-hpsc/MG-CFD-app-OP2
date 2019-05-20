#!/bin/bash

set -e

test_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# input_data_root_dir="$(cd ${test_dir}/../input_data && pwd)"
input_data_root_dir="../input_data"

miniapp_plain_dir="${HOME}/Working.Github/MG-CFD-app-plain"

####################
## Input settings ##
####################

input_orig_data_dir="${input_data_root_dir}/m6wing/rodinia.regenerated"
input_mini_dat=input.dat

input_op2_data_dir="${input_data_root_dir}/m6wing/hdf5.original.OP2-compatible.legacy-numbering"
input_op2_file=input.dat

LEVELS=(0 1 2 3)

####################

###################
## Test settings ##
###################

arrays_to_compare_exact=()
arrays_to_compare_tolerable=()
arrays_to_compare_tolerable+=(variables)

precision="'%.17e'"

# OPTS_master="-DOMP"
OPTS_master="-DLEGACY_ORDERING"

####################

test_dir=`dirname $0`
test_dir=`cd ${test_dir} && pwd`

miniapp_op2_dir=`cd ../../ ; pwd`
miniapp_op2_bin_dir="${miniapp_op2_dir}/bin"
miniapp_op2_src_dir="${miniapp_op2_dir}/src"

miniapp_plain_bin_dir="${miniapp_plain_dir}/bin/`hostname`"
miniapp_plain_src_dir="${miniapp_plain_dir}/src"

output_data_dir="${test_dir}/data"
mkdir -p "${output_data_dir}"

cycles=10
# cycles=1

config_plain="${test_dir}/input-plain.config"
echo "input_file = $input_mini_dat" > "$config_plain"
echo "input_file_directory = ${input_orig_data_dir}" >> "$config_plain"
echo "output_file_prefix = ${output_data_dir}/" >> "$config_plain"
echo "output_variables = Y" >> "$config_plain"
echo "cycles = $cycles" >> "$config_plain"

config_op2="${test_dir}/input-op2.config"
echo "input_file = $input_op2_file" > "$config_op2"
echo "input_file_directory = ${input_op2_data_dir}" >> "$config_op2"
# echo "output_file_prefix = ${output_data_dir}/" >> "$config_op2"
## NOTE: If I use the FULL filepath for 'output_file_prefix', it must 
##       trigger a buffer overflow in HDF5 because it seg faults upon 
##       final cleanup. So must use a shorter (relative) filepath:
echo "output_file_prefix = ./data/" >> "$config_op2"
echo "output_variables = Y" >> "$config_op2"
echo "cycles = $cycles" >> "$config_op2"

source "${test_dir}/../Scripts/fn_verify.sh"

compile() {
	set -e

	cd "${miniapp_plain_dir}"
	export BUILD_FLAGS="${OPTS_master}"
	export CFLAGS="-fp-model precise"
	export OPT_LEVEL="2"
	CC=intel make -j4
	export BUILD_FLAGS=""

	cd "${miniapp_op2_dir}"
	make mgcfd_seq
}

clean() {
	set -e

	rm -f "${output_data_dir}"/*

	cd "${test_dir}/${input_orig_data_dir}"
	for A in ${arrays_to_compare_exact[@]}; do
		rm -f "${A}.L{0,1,2,3,4}.-1"
	done
	for A in ${arrays_to_compare_tolerable[@]}; do
		rm -f "${A}.L{0,1,2,3,4}.-1"
	done
}

grab_output_dataset() {
	set -e

	L=$1
	arr=$2
	suffix=$3

	arr_filepath="${output_data_dir}/${arr}.size=1x.cycles=${cycles}.level=${L}"
	h5_filepath="${output_data_dir}/${arr}.L${L}.cycles=${cycles}.h5"
	if [ -f $h5_filepath ]; then
		h5dump --noindex -m ${precision} --width=400 -o "${arr_filepath}" -d p_${arr}_result_L${L} "${h5_filepath}" > /dev/null
		cat "${arr_filepath}" | tail -n+2 | tr -d ' ' | sed "s/,$//g" | tr -d "'" | sed "s/,/ /g" > "${arr_filepath}"2
		mv "${arr_filepath}"2 "${arr_filepath}"
	fi
	if [ ! -f "$arr_filepath" ]; then
		echo "ERROR: Can't find: ${arr_filepath}"
		exit 1
	fi
	mv "${arr_filepath}" "${output_data_dir}/${arr}.${suffix}.L$L"
}

execute() {
	cd "$test_dir"

	echo "${miniapp_plain_bin_dir}/euler3d_cpu_double_intel${OPTS_master}.b" -c "$config_plain"
	eval "${miniapp_plain_bin_dir}/euler3d_cpu_double_intel${OPTS_master}.b" -c "$config_plain"
	# for L in "${LEVELS[@]}"; do
		L=0
		for arr in ${arrays_to_compare_exact[@]}; do
			grab_output_dataset $L $arr master
		done
		for arr in ${arrays_to_compare_tolerable[@]}; do
			grab_output_dataset $L $arr master
		done
	# done

	cd "$test_dir"
	echo "${miniapp_op2_bin_dir}/mgcfd_seq" -l -c "$config_op2"
	eval "${miniapp_op2_bin_dir}/mgcfd_seq" -l -c "$config_op2"
	# for L in "${LEVELS[@]}"; do
		L=0
		for arr in ${arrays_to_compare_exact[@]}; do
			grab_output_dataset $L $arr op2
		done
		for arr in ${arrays_to_compare_tolerable[@]}; do
			grab_output_dataset $L $arr op2
		done
	# done
}

verify() {
	set -e

	cd "${output_data_dir}"
	for A in ${arrays_to_compare_exact[@]}; do
		# for ((l=0 ; l < 4; l++)); do
			l=0
			verify_level $A op2 0 $l 0.0
		# done
	done
	for A in ${arrays_to_compare_tolerable[@]}; do
		# for ((l=0 ; l < 4; l++)); do
			l=0
			verify_level $A op2 0 $l
		# done
	done
}

compile
# clean
execute
verify
