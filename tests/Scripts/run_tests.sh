#!/bin/bash

BASE_DIR=`dirname $0`
cd $BASE_DIR
BASE_DIR=`pwd`

cd ../
TESTS_DIR=`pwd`

rm ../bin/*

for D in ./*; do
	cd $D
	if [ -f run.sh ]; then
		test_name=`basename $D`
		echo "Running test $test_name"
		./run.sh &> run.log
		if [ $? -eq 1 ]; then
			echo "Test fail detected"
			exit 1
		fi
	fi
	cd ../
done

echo "All tests passed"

exit 0