#!/bin/bash

patient="TCGA-CJ-4642"
path_to_compare0="data/patients/${patient}/matchings_predictor_output/original_data"
echo "Comparing the results of the serial and of the parallel code for the patient \"${patient}\". WARNING: check that main.cpp is producing output in \"path_to_compare0\"."

## makes the program output only data relevant for the testing
jq --tab '.debug_parallelization = true' global_parameters.json  > tmp.66015.json && mv tmp.66015.json global_parameters.json

## disables the parallelization
cp omp_disabler_xcode.h omp_disabler.h
time ./expression_profiles

mkdir -p test_temporary_data
cp -r "${path_to_compare0}" test_temporary_data/serial_code

## enables the parallelization
echo "" > omp_disabler.h
time ./expression_profiles

cp -r "${path_to_compare0}" test_temporary_data/parallel_code

## sort everyline in every file to compare, since the order of some output files is arbitrary (a std::unordered_map is showed)
to_sort = $(find test_temporary_data -name '*' -type f)
for f in ${to_sort}; do
    sort -o $f $f
done
output=$(diff -bur test_temporary_data/serial_code test_temporary_data/parallel_code)
if [[ $output ]]; then
    echo $output
else
    echo "the results coincide"
    rm -rf test_temporary_data
fi

## restore the normal output of the program
jq --tab '.debug_parallelization = false' global_parameters.json  > tmp.$$.json && mv tmp.$$.json global_parameters.json
