
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib


all: green polynomial

green:	green.o $(LIB)/bernoulli.o
polynomial:	polynomial.o

green.o:	green.c
polynomial.o:	polynomial.c

.c.o:
	cc -c -O2 -I $(INC) $<

.o:
	cc -o $@ $< -lm

