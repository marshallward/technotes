#!/bin/sh
FC=gfortran
#FCFLAGS="-g -O0"
#FCFLAGS="-g -O3"
FCFLAGS="-g -O3 -mavx -mfma -march=native"
FCFLAGS+=" -mavx512f"
FCFLAGS+=" -flto"
LD=gfortran
LDFLAGS="-flto"
export FC FCFLAGS LD LDFLAGS
./configure && make clean && make
