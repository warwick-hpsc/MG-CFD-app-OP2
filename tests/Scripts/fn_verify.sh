#!/bin/bash

function verify {
	dataset_master=$1
	dataset_test=$2
	header_line_length=$3

	for ((l=0 ; l < 4; l++)); do
		L="${LEVELS[$l]}"
		LC="${LEVELCODES[$l]}"

		if [ ! -f ${dataset_master}.L$L ]; then
			echo "Cannot find ${dataset_master}.L$L"
			exit 1
		fi
		if [ ! -f ${dataset_test}.L$L ]; then
			echo "Cannot find ${dataset_test}.L$L"
			exit 1
		fi

		# Verify the header, then if correct verify the contents:
		if [ "${header_line_length}" -gt "0" ] && [ "$(head -n 1 ${dataset_master}.L$L)" != "$(head -n 1 ${dataset_test}.L$L)" ]; then
			echo "Header of ${dataset_test} level $L incorrect"
			exit 1
		else
			if ! cmp -s ${dataset_master}.L$L ${dataset_test}.L$L
			then
				Rscript ../../Scripts/compare_points.R ${dataset_master}.L$L ${dataset_test}.L$L ${header_line_length} 0.0
				if [ -f "mismatch.out" ]; then
					rm mismatch.out
					echo "${dataset_test} mismatch on level $L"
					exit 1
				fi
			fi
		fi
	done
	echo "${dataset_test} correct on all levels"
}

function verify_level {
	dataset=$1
	opts=$2
	header_line_length=$3
	l=$4
	tolerance=$5

	fail="false"

	L="${LEVELS[$l]}"
	LC="${LEVELCODES[$l]}"

	if [ ! -f ${dataset}.master.L$L ]; then
		echo "Cannot find ${dataset}.master.L$L"
		exit 1
	fi
	if [ ! -f ${dataset}.${opts}.L$L ]; then
		echo "Cannot find ${dataset}.${opts}.L$L"
		exit 1
	fi

	# Verify the header, then if correct verify the contents:
	if [ "${header_line_length}" -gt "0" ] && [ "$(head -n 1 ${dataset}.master.L$L)" != "$(head -n 1 ${dataset}.${opts}.L$L)" ]; then
		echo "Header of ${dataset} level $L incorrect"
		fail="true"
		exit 1
	else
		if ! cmp -s ${dataset}.master.L$L ${dataset}.${opts}.L$L
		then
			Rscript ../../Scripts/compare_points.R ${dataset}.master.L$L ${dataset}.${opts}.L$L ${header_line_length} ${tolerance}
			if [ -f "mismatch.out" ]; then
				rm mismatch.out
				echo "${dataset} mismatch on level $L"
				fail="true"
				exit 1
			else
				echo "${dataset} correct on level $l within tolerance"
			fi
		else
			echo "${dataset} correct on level $l"
		fi
	fi
}
