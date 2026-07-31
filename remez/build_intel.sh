#!/bin/sh
FC=ifx
FCFLAGS="-g -O3 -ipo"
#FCFLAGS+=" -fp-model source"
FCFLAGS+=" -xHost"
LD=ifx
LDFLAGS="-ipo"
export FC FCFLAGS LD LDFLAGS
./configure && make clean && make
