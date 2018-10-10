#!/bin/bash
# generate an Xcode project into the folder _build
cd xcode
# CC=/usr/local/Cellar/llvm/6.0.1/bin/clang
# CXX=/usr/local/Cellar/llvm/6.0.1/bin/clang++
# cmake -DCMAKE_CXX_COMPILER=/usr/local/Cellar/llvm/6.0.1/bin/clang++ -DCMAKE_C_COMPILER=/usr/local/Cellar/llvm/6.0.1/bin/clang -G Xcode -H.. -B_build
# CC=clang-omp
# CXX=clang-omp++
# CC=clang
# CXX=clang++
cp ../omp_disabler_xcode.h ../omp_disabler.h
cmake -G Xcode -H.. -B_build -DOMP_TEMPORARILY_DISABLED=ON
mkdir -p _build/Debug
mkdir -p _build/Debug/data
mkdir -p _build/Debug/data/processed
cd ../
rsync -avh global_parameters.json xcode/_build/Debug/global_parameters.json
rsync -avh data/processed/mirnas_with_scored_interactions.tsv xcode/_build/Debug/data/processed/mirnas_with_scored_interactions.tsv
rsync -avh data/processed/sites_with_scored_interactions.tsv xcode/_build/Debug/data/processed/sites_with_scored_interactions.tsv
rsync -avh data/processed/scored_interactions_processed.tsv xcode/_build/Debug/data/processed/scored_interactions_processed.tsv
rsync -avh data/patients xcode/_build/Debug/data/ --exclude=matchings_predictor_output
