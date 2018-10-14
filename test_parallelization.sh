#!/bin/bash

patient="TCGA-CJ-4642"
path_to_compare="data/patients/${patient}/matchings_predictor_output/original_data"
echo "Comparing the results of the serial and of the parallel code for the patient \"${patient}\". WARNING: check that main.cpp is producing output in \"path_to_compare\"."

# echo "removing old output, if any"
# rm -rf test_temporary_data

# makes the program output only data relevant for the testing
my_uuid="dbe477e5-5959-4e45-b3ef-1ef8c5b73987"
jq --tab '.test_parallelization = true' global_parameters.json  > tmp.${my_uuid}.json && mv tmp.${my_uuid}.json global_parameters.json

# disables the parallelization
sh build.sh --serial
sh run.sh

mkdir -p test_temporary_data
cp -R "${path_to_compare}/" test_temporary_data/serial_code

# enables the parallelization
sh build.sh --parallel
sh run.sh

cp -R "${path_to_compare}/" test_temporary_data/parallel_code

# sort every line in every file to compare, since the order of some output files is arbitrary (a std::unordered_map is showed)
echo "sorting the files by line to make the files comparable"
to_sort=$(find test_temporary_data -name '*' -type f ! -name '*.png' ! -name '*.gif' ! -name '*_archived')
for f in ${to_sort}; do
    sort -o $f $f
    if [ "$?" -ne 0 ]; then
        echo "last file passed to sort: \"${f}\""
    fi
done
echo "finished sorting"

output=$(diff -bur test_temporary_data/serial_code test_temporary_data/parallel_code)
if [[ $output ]]; then
    echo $output
else
    echo "the results coincide"
    rm -rf test_temporary_data
fi

# restore the normal output of the program
jq --tab '.test_parallelization = false' global_parameters.json  > tmp.$$.json && mv tmp.$$.json global_parameters.json
