
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib


all: green

green:	green.o polynomial.o $(LIB)/bernoulli.o $(LIB)/binomial.o

green.o:	green.c
polynomial.o:	polynomial.c gaussian.h

.c.o:
	cc -c -O2 -I $(INC) $<

.o:
	cc -o $@ $^ -lm

