#!/bin/sh
FC=nvfortran
FCFLAGS="-g -O0"
#FCFLAGS=" -Mnofma"
FCFLAGS+=" -Minline"
FCFLAGS+=" -Minfo=inline"
LD=nvfortran
#LDFLAGS="-Mipa"
export FC FCFLAGS LD LDFLAGS
./configure && make clean && make
