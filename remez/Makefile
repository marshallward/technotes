FC = gfortran
FCFLAGS = -g -O3
LD = gfortran

all: runexp

runexp: exp_repro.o profile.o
	$(LD) -o $@ $^

exp_repro.o: exp_repro.f90
	$(FC) $(FCFLAGS) -c $<

profile.o: profile.f90
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean
clean:
	rm -f exp_repro.o exp_repro.mod profile.o runexp
