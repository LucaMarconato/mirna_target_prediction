#!/bin/bash
# CC=clang
# CXX=clang++

# CC=/usr/local/Cellar/llvm/6.0.1/bin/clang
# CXX=/usr/local/Cellar/llvm/6.0.1/bin/clang++

CC=clang-omp
CXX=clang-omp++

# CC=/usr/local/opt/llvm/bin/clang
# CXX=/usr/local/opt/llvm/bin/clang++
# LDFLAGS=-L/usr/local/opt/llvm/lib/
# cmake -DCMAKE_BUILD_TYPE=Release .

echo "" > omp_disabler.h
cmake .
make
