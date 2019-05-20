#!/bin/bash

cd `dirname $0`
projectDir=`pwd`

rm -rf "${projectDir}/cuda"
rm -rf "${projectDir}/openacc"
rm -rf "${projectDir}/openmp"
rm -rf "${projectDir}/openmp4"
rm -rf "${projectDir}/seq"
rm -rf "${projectDir}/src_op"
