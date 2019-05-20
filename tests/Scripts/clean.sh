#!/bin/bash

DIR=`dirname $0`

rm -fr ${DIR}/../*/data
rm -f ${DIR}/../*/run.log
rm -f ${DIR}/../*/*.config
rm -f ${DIR}/../input_data/*/*density*
rm -f ${DIR}/../input_data/*/*momentum*
rm -f ${DIR}/../input_data/*/*variables*
rm -f ${DIR}/../input_data/*/fluxes.*
rm -f ${DIR}/../input_data/*/edge_*
rm -f ${DIR}/../input_data/*/areas.*
rm -f ${DIR}/../input_data/*/step_factors.*
rm -f ${DIR}/../input_data/*/*.bin
