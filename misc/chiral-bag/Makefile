

all: bessel-test find_zero-test

bessel-test:	bessel-test.o bessel.o
	cc -o $@ $^ -lm


find_zero-test:	find_zero-test.o find_zero.o
	cc -o $@ $^ -lm


.c.o:
	cc -g -O2 -c -Wall $<

