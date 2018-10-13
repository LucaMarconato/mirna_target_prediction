#!/bin/bash

patient="TCGA-CJ-4642"
echo "TODO: one path should refer to mirnas, another one to p_j_downregulated_given_c_bound, both of them must be binary data generated with this purpose, written in defined ordered way and without reference to pointers addresses"
exit
path_to_compare0="data/patients/${patient}/matchings_predictor_output/original_data"
echo "Comparing the results of the serial and of the parallel code for the patient \"${patient}\". WARNING: check that main.cpp is producing output in \"path_to_compare0\"."

## disable parallelization
cp omp_disabler_xcode.h omp_disabler.h
time ./expression_profiles

mkdir -p test_temporary_data
cp -r "${path_to_compare0}" test_temporary_data/serial_code

## enable parallelization
echo "" > omp_disabler.h
time ./expression_profiles

cp -r "${path_to_compare0}" test_temporary_data/parallel_code

output=$(diff -bur test_temporary_data/serial_code test_temporary_data/parallel_code)
if [[ $output ]]; then
    echo $output
else
    echo "the results coincide"
    rm -rf test_temporary_data
fi
