

all: bessel-test find_zero-test find_eig sums asym

bessel-test:	bessel-test.o bessel.o
	cc -o $@ $^ -lm


find_zero-test:	find_zero-test.o find_zero.o
	cc -o $@ $^ -lm

find_eig:	find_eig.o bag_ener.o find_zero.o bessel.o
	cc -o $@ $^ -lm

sums:	sums.o bag_ener.o find_zero.o bessel.o
	cc -o $@ $^ -lm

asym: asym.o
	cc -o $@ $^ -lm

.c.o:
	cc -g -O2 -c -Wall $<

