
#
# Makefile that uses lapack for matrix operations
#

TOP = ../../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib


all: hyper eigen comb hyper

OBJS=../ache.o ../zetafn.o $(LIB)/binomial.o $(LIB)/harmonic.o

comb: comb.o $(OBJS)
eigen: eigen.o $(OBJS)

hyper: hyper.o 

.o:
	cc -o $@ $^ -L $(LIB) -llapack -lf77blas -latlas -lg2c -lm

.c.o:
	cc -g -O2 -c $< -I.. -I$(INC)

