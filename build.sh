#!/bin/bash
CC=clang
CXX=clang++
# CC=/usr/local/opt/llvm/bin/clang
# CXX=/usr/local/opt/llvm/bin/clang++
# LDFLAGS=-L/usr/local/opt/llvm/lib/
# cmake -DCMAKE_BUILD_TYPE=Release .
cmake .
make
