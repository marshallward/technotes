#!/bin/sh
FC=ifx
FCFLAGS="-g -O3"
#FCFLAGS+=" -fp-model source"
FCFLAGS+=" -xHost"
#FCFLAGS+=" -fma"
FCFLAGS+=" -no-ftz"
FCFLAGS+=" -ipo"
LD=ifx
LDFLAGS="-ipo"
export FC FCFLAGS LD LDFLAGS
./configure && make clean && make
