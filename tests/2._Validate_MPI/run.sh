#!/bin/bash

test_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# input_data_root_dir="$(cd ${test_dir}/../input_data && pwd)"
input_data_root_dir="../input_data"

####################
## Input settings ##
####################

input_data_dir="${input_data_root_dir}/m6wing/hdf5.original"
# input_data_dir="${input_data_root_dir}/rotor37/Rotor37_1M_OP2"
# input_data_dir="${input_data_root_dir}/rotor37/Rotor37_8M_OP2"

input_file=input.dat
LEVELS=(0 1 2 3)

bin_name=mgcfd_mpi_genseq

####################

###################
## Test settings ##
###################

miniapp_op2_dir=`cd "$test_dir"/../../ ; pwd`
miniapp_op2_bin_dir="${miniapp_op2_dir}/bin"
miniapp_op2_src_dir="${miniapp_op2_dir}/src"

output_data_dir="${test_dir}/data"
mkdir -p "${output_data_dir}"

config="${test_dir}/config"
master_config="${test_dir}/master.config"
test_config="${test_dir}/test.config"

cycles=25

generate_config() {
	set -e
	
	echo "input_file = $input_file" > "$config"
	echo "input_file_directory = ${input_data_dir}" >> "$config"
	
	# echo "output_file_prefix = ${output_data_dir}/" >> "$config"
	## NOTE: If I use the FULL filepath for 'output_file_prefix', it must 
	##       trigger a buffer overflow in HDF5 because it seg faults upon 
	##       final cleanup. So must use a shorter (relative) filepath:
	echo "output_file_prefix = ./data/" >> "$config"

	echo "cycles = $cycles" >> "$config"

	cp "$config" "$master_config"
	echo "output_variables = Y" >> "$master_config"

	cp "$config" "$test_config"
	echo "validate_result = Y" >> "$test_config"
}

compile() {
	set -e

	cd "${miniapp_op2_dir}"

	# make -j4 mgcfd_seq
	# make -j4 $bin_name
	make -j4 mgcfd_seq $bin_name
}

gen_solution() {
	set -e

	solution_needed=false
	# true_dd=`echo "$input_data_dir" | sed "s|${input_data_root_dir}|${HOME}/Datasets/Live|"`
	true_dd="$input_data_dir"
	for L in `seq 0 3`; do
		vf="variables.L${L}.cycles=${cycles}.h5"
		if [ ! -f "${true_dd}/solution.${vf}" ]; then
			solution_needed=true
			break
		fi
	done

	cd "$test_dir"
	if $solution_needed ; then
		eval "${miniapp_op2_bin_dir}"/mgcfd_seq OP_MAPS_BASE_INDEX=1 -c "$master_config"

		for L in `seq 0 3`; do
			vf="variables.L${L}.cycles=${cycles}.h5"
			mv "${output_data_dir}/${vf}" "${true_dd}/solution.${vf}"
		done
	fi

	for L in `seq 0 3`; do
		vf="variables.L${L}.cycles=${cycles}.h5"
		if [ ! -f "${input_data_dir}/solution.${vf}" ]; then
			ln -s "${true_dd}/solution.${vf}" "${input_data_dir}/solution.${vf}"
		fi
	done
}

execute() {
	set -e

	cd "$test_dir"
	mpirun -np 4 "${miniapp_op2_bin_dir}/${bin_name}" OP_MAPS_BASE_INDEX=1 -c "$test_config"
}

generate_config
compile
gen_solution
execute
