#!/bin/bash

patient="TCGA-CJ-4642"
paths_to_compare="data/patients/${patient}/matching_predictor_output/original_data"
echo "Comparing the results of the serial and of the parallel code for the patient \"${patient}\". WARNING: check that main.cpp is producing output in \"paths_to_compare\"."

## disable parallelization
cp omp_disabler_xcode.h omp_disabler.h
run

mkdir -p test_temporary_data
cp "${paths_to_compare}" test_temporary_data/serial_code

## enable parallelization
echo "" > omp_disabler.h
run
cp "${paths_to_compare}" test_temporary_data/parallel_code

output=$(diff -bur test_temporary_data/serial_code test_temporary_data/parallel_code)
if [[ $output ]]; then
    echo $output
else
    echo "the results coincide"
    rm -rf test_temporary_data
fi
