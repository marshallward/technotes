#!/bin/sh
FC=gfortran
FCFLAGS="-g -O3 -mavx -mfma -march=native"
FCFLAGS+=" -flto"
LD=gfortran
LDFLAGS="-flto"
export FC FCFLAGS LD LDFLAGS
./configure && make clean && make
