

EXES= diffy galois

all: $(EXES)


galois: galois.o
diffy: diffy.o


.o:
	cc -o $@ $^ -llapack -lf77blas -latlas -lg2c -lm

.c.o:
	cc -g -O2 -c $< 

clean:
	rm -f *.o junk a.out x

realclean: clean
	rm -f $(EXES)  *.png

