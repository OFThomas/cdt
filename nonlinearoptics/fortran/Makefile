FC = gfortran 
FCFLAGS=-g -Wall -fbounds-check

FCFLAGS += -I/usr/include 
LDFLAGS=-L/usr/lib 
LOADLIBS= -llapack -lblas

PROGRAMS = olis_f90stdlib num_hom

all: $(PROGRAMS)

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F()
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
