#!/bin/bash
# CC=clang
# CXX=clang++

CC=/usr/local/Cellar/llvm/7.0.0/bin/clang
CXX=/usr/local/Cellar/llvm/7.0.0/bin/clang++

# CC=clang-omp
# CXX=clang-omp++

# CC=/usr/local/opt/llvm/bin/clang
# CXX=/usr/local/opt/llvm/bin/clang++
# LDFLAGS=-L/usr/local/opt/llvm/lib/
# cmake -DCMAKE_BUILD_TYPE=Release .

if [ "$1" = "--dont_change_omp_disabler_h" ]; then
    echo "not generating omp_disabler.h"
else
    if [ "$1" = "--serial" ]; then
        new_value="#define OMP_TEMPORARILY_DISABLED"
        echo "code compiled to avoid any openMP preprocessor directive, to comply with compiler not supporting openMP"
    else
        new_value="// #define OMP_TEMPORARILY_DISABLED"
        echo "code compiled to be parallel"
    fi
    echo ${new_value} > omp_disabler.h
fi

# CC=gcc-8
# CXX=g++-8
cmake .
# CC=gcc-8
# CXX=g++-8
make -j7
