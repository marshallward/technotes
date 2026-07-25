## GCC
#FC = gfortran
#FCFLAGS = -g -O3 \
#  -mavx -mfma -march=native \
#  -flto --param early-inlining-insns=64 \
#  -fopt-info-inline-all=inline.log \
#  -fno-range-check
#LD = gfortran
#LDFLAGS = -flto

# ifx (Cray-wrapped)
FC = ftn
FCFLAGS = -g -O3 -ipo
FCFLAGS += -assume protect_parens
#FCFLAGS += -fp-model source
LD = ftn
LDFLAGS = -ipo

## NVHPC
#FC = nvfortran
#FCFLAGS = -g -O3 -stdpar=gpu -Mnofma -Minfo=all
#LD = nvfortran
#LDFLAGS = -stdpar=gpu

all: runexp

runexp: exp_dev.o exp_repro.o profile.o
	$(LD) $(LDFLAGS) -o $@ $^

exp_dev.o: exp_dev.f90
	$(FC) $(FCFLAGS) -c $<

exp_repro.o: exp_repro.f90
	$(FC) $(FCFLAGS) -c $<

profile.o: profile.f90 exp_dev.o exp_repro.o
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean
clean:
	rm -f exp_dev.o exp_dev.mod exp_repro.o exp_repro_mod.mod profile.o runexp
