
EXES= eigen eigen-extrap gral-matrix gral-explore

all: $(EXES)


eigen: eigen.o gral-simple.o
eigen-extrap: eigen-extrap.o gral-simple.o
gral-matrix: gral-matrix.o
gral-explore: gral-explore.o gral-simple.o

.c.o:
	cc -c -g -O2 $<

.o:
	cc -g -o $@ $^ -llapack -lf77blas -latlas -lg2c -lm


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png

