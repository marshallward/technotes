## GCC
#FC = gfortran
#FCFLAGS = -g -O3 -mavx -mfma -march=native
#LD = gfortran
#LDFLAGS =

## ifx (Cray-wrapped)
#FC = ftn
#FCFLAGS = -g -O3 -qopenmp-simd
#LD = ftn
#LDFLAGS =

# NVHPC
FC = nvfortran
FCFLAGS = -g -O3 -stdpar=gpu -Mnofma -Minfo=all
LD = nvfortran
LDFLAGS = -stdpar=gpu

all: runexp

runexp: exp_repro.o profile.o
	$(LD) $(LDFLAGS) -o $@ $^

exp_repro.o: exp_repro.f90
	$(FC) $(FCFLAGS) -c $<

profile.o: profile.f90
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean
clean:
	rm -f exp_repro.o exp_repro.mod profile.o runexp
