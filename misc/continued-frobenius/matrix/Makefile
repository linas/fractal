
#
# Makefile that uses lapack for matrix operations
#

all: eigen comb

OBJS=../ache.o ../zetafn.o

eigen: eigen.o 

comb: comb.o

.o:
	cc -o $@ $< $(OBJS) -llapack -lf77blas -latlas -lg2c -lm

.c.o:
	cc -g -O2 -c $< -I..

