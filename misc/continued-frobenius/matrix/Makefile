
#
# Makefile that uses lapack for matrix operations
#

all: hyper eigen comb hyper

OBJS=../ache.o ../zetafn.o

comb: comb.o $(OBJS)
eigen: eigen.o $(OBJS)

hyper: hyper.o 

.o:
	cc -o $@ $^ -llapack -lf77blas -latlas -lg2c -lm

.c.o:
	cc -g -O2 -c $< -I..

